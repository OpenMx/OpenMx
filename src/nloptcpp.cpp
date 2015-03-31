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

#include <iostream>

#ifdef HAS_NLOPT

#include <nlopt.h>

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
	assert(n == goc->fc->numParam);
	int mode = 1;
	double fit = goc->solFun((double*) x, &mode);
	//if (goc->verbose >= 3 && !grad) std::cout << "fit: " << fit << "\n";
	if (mode == -1) {
		if (!goc->feasible) {
			nlopt_opt opt = (nlopt_opt) goc->extraData;
			nlopt_force_stop(opt);
		}
		return std::numeric_limits<double>::max(); // this is wrong TODO
	}
	if (!grad) return fit;

	Eigen::Map< Eigen::VectorXd > Epoint((double*) x, n);
	Eigen::Map< Eigen::VectorXd > Egrad(grad, n);
	fit_functional ff(*goc);
	fd_gradient(ff, Epoint, Egrad);
	//if (goc->verbose >= 3) std::cout << "fit: " << fit << "\ngrad:\n" << Egrad << "\n";
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
	//if (goc->verbose >= 3 && grad) std::cout << "inequality:\n" << Eresult << "\n" << jacobian << "\n";
}

void omxInvokeNLOPT(double *est, GradientOptimizerContext &goc)
{
	goc.optName = "NLOPT";
	goc.setupSimpleBounds();
	goc.useGradient = false;  // not implemented yet, would be very helpful

	FitContext *fc = goc.fc;
	omxState *globalState = fc->state;
	int ncnln = globalState->ncnln;
    
        nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, fc->numParam);
	goc.extraData = opt;
        //local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm
        
        //nlopt_set_local_optimizer(opt, local_opt);
        nlopt_set_lower_bounds(opt, goc.solLB.data());
        nlopt_set_upper_bounds(opt, goc.solUB.data());
	nlopt_set_ftol_rel(opt, 1e-9);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
        
	nlopt_set_min_objective(opt, nloptObjectiveFunction, &goc);
        if (ncnln) {
		int neq = 0;
		for (int j = 0; j < globalState->numConstraints; j++) {
			if (globalState->conList[j].opCode == omxConstraint::EQUALITY)
				neq += globalState->conList[j].size;
		}
		int nineq = ncnln - neq;
            
                if(nineq > 0){
			goc.inequality.resize(nineq); // TODO remove
			std::vector<double> tol(nineq, sqrt(std::numeric_limits<double>::epsilon()));
			nlopt_add_inequality_mconstraint(opt, nineq, nloptInequalityFunction, &goc, tol.data());
                }
                
                if (neq > 0){
			goc.equality.resize(neq); // TODO remove
			std::vector<double> tol(neq, sqrt(std::numeric_limits<double>::epsilon()));
			nlopt_add_equality_mconstraint(opt, neq, nloptEqualityFunction, &goc, tol.data());
                }
	}
        
	int code = nlopt_optimize(opt, est, &fc->fit);

        nlopt_destroy(opt);

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
	if (code == NLOPT_MAXEVAL_REACHED) {
		goc.informOut = INFORM_ITERATION_LIMIT;
	} else if (code == NLOPT_ROUNDOFF_LIMITED) {
		goc.informOut = INFORM_UNCONVERGED_OPTIMUM;  // is this correct? TODO
	} else {
		goc.informOut = INFORM_CONVERGED_OPTIMUM;
	}
}

#endif
