#include <limits>
#include <cstdlib>
#include <ctype.h>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "omxState.h"
#include "omxMatrix.h"
#include "glue.h"
#include "nloptcpp.h"
#include "finiteDifferences.h"

#include "nlopt.h"

struct fit_functional {
	GradientOptimizerContext &goc;

	fit_functional(GradientOptimizerContext &goc) : goc(goc) {};

	template <typename T1>
	double operator()(Eigen::MatrixBase<T1>& x) const {
		int mode = 0;
		return goc.solFun(x.derived().data(), &mode);
	}
};

static double nloptObjectiveFunction(unsigned n, const double *x, double *grad, void *f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	FitContext *fc = goc->fc;
	assert(n == fc->numParam);
	int mode = grad != 0;
	double fit = goc->solFun((double*) x, &mode);
	if (mode == -1) {
		if (!goc->feasible) {
			nlopt_opt opt = (nlopt_opt) goc->extraData;
			nlopt_force_stop(opt);
		}
		return nan("infeasible");
	}
	if (!grad) return fit;

	Eigen::Map< Eigen::VectorXd > Epoint((double*) x, n);
	Eigen::Map< Eigen::VectorXd > Egrad(grad, n);
	if (fc->wanted & FF_COMPUTE_GRADIENT) {
		Egrad = fc->grad;
	} else {
		if (goc->verbose >= 3) mxLog("fd_gradient start");
		fit_functional ff(*goc);
		fd_gradient(ff, Epoint, Egrad);
		fc->grad = Egrad;
	}
	if (goc->verbose >= 3) {
		mxPrintMat("gradient", Egrad);
	}
	return fit;
}

struct equality_functional {
	GradientOptimizerContext &goc;

	equality_functional(GradientOptimizerContext &goc) : goc(goc) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		memcpy(goc.fc->est, x.derived().data(), sizeof(double) * goc.fc->numParam);
		goc.fc->copyParamToModel();
		goc.solEqBFun();
		result = goc.equality;
	}
};

static void nloptEqualityFunction(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	assert(n == goc->fc->numParam);
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobian(grad, n, m);
	equality_functional ff(*goc);
	fd_jacobian(ff, Epoint, Eresult, grad==NULL, jacobian);
	//if (goc->verbose >= 3 && grad) std::cout << "equality:\n" << Eresult << "\n" << jacobian << "\n";
}

struct inequality_functional {
	GradientOptimizerContext &goc;

	inequality_functional(GradientOptimizerContext &goc) : goc(goc) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		memcpy(goc.fc->est, x.derived().data(), sizeof(double) * goc.fc->numParam);
		goc.fc->copyParamToModel();
		goc.myineqFun();
		result = goc.inequality;
	}
};

static void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	assert(n == goc->fc->numParam);
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobian(grad, n, m);
	inequality_functional ff(*goc);
	fd_jacobian(ff, Epoint, Eresult, grad==NULL, jacobian);
	if (grad && !std::isfinite(Eresult.sum())) {
		// infeasible at start of major iteration
		nlopt_opt opt = (nlopt_opt) goc->extraData;
		nlopt_force_stop(opt);
	}
	if (goc->verbose >= 3 && grad) {
		mxPrintMat("inequality jacobian", jacobian);
	}
}

void omxInvokeNLOPT(double *est, GradientOptimizerContext &goc)
{
	goc.optName = "SLSQP";
	goc.setupSimpleBounds();
	goc.useGradient = true;

	FitContext *fc = goc.fc;
	int oldWanted = fc->wanted;
	fc->wanted = 0;
	omxState *globalState = fc->state;
    
        nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, fc->numParam);
	goc.extraData = opt;
        //local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm
        
        //nlopt_set_local_optimizer(opt, local_opt);
        nlopt_set_lower_bounds(opt, goc.solLB.data());
        nlopt_set_upper_bounds(opt, goc.solUB.data());

	int eq, ieq;
	globalState->countNonlinearConstraints(eq, ieq);

	if (fc->CI) {
		nlopt_set_xtol_rel(opt, 1e-3);
		std::vector<double> tol(fc->numParam, 1e-3);
		//nlopt_set_xtol_abs(opt, tol.data());
	} else {
		// The 0.01 factor is a bit ridiculous. NPSOL doesn't
		// have the same definition of relative tolerance
		// compare to SLSQP.
		nlopt_set_ftol_rel(opt, Global->optimalityTolerance * 0.01);
		nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
	}
        
	nlopt_set_min_objective(opt, nloptObjectiveFunction, &goc);

        if (eq + ieq) {
		double feasibilityTolerance = Global->feasibilityTolerance;
                if (ieq > 0){
			goc.inequality.resize(ieq);
			std::vector<double> tol(ieq, feasibilityTolerance);
			nlopt_add_inequality_mconstraint(opt, ieq, nloptInequalityFunction, &goc, tol.data());
                }
                
                if (eq > 0){
			goc.equality.resize(eq);
			std::vector<double> tol(eq, feasibilityTolerance);
			nlopt_add_equality_mconstraint(opt, eq, nloptEqualityFunction, &goc, tol.data());
                }
	}
        
	int code = nlopt_optimize(opt, est, &fc->fit);
	if (goc.verbose >= 2) mxLog("nlopt_optimize returned %d", code);

        nlopt_destroy(opt);

	fc->wanted = oldWanted;

	// fatal errors
	if (code == NLOPT_INVALID_ARGS) {
		Rf_error("NLOPT invoked with invalid arguments");
	} else if (code == NLOPT_OUT_OF_MEMORY) {
		Rf_error("NLOPT ran out of memory");
	} else if (code == NLOPT_FORCED_STOP) {
		goc.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
	} else if (code < 0) {
		Rf_error("NLOPT fatal error %d", code);
	}

	// non-fatal errors
	if (eq + ieq == 0 && (fc->grad.array().abs() > 0.1).any()) {
		goc.informOut = INFORM_NOT_AT_OPTIMUM;
	} else if (code == NLOPT_MAXEVAL_REACHED) {
		goc.informOut = INFORM_ITERATION_LIMIT;
	} else if (code == NLOPT_ROUNDOFF_LIMITED) {
		goc.informOut = INFORM_UNCONVERGED_OPTIMUM;  // is this correct? TODO
	} else {
		goc.informOut = INFORM_CONVERGED_OPTIMUM;
	}
}
