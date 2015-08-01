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

namespace SLSQP {

	struct context {
		GradientOptimizerContext &goc;
		int origeq;
		int eqredundent;
		std::vector<bool> eqmask;
		context(GradientOptimizerContext &goc) : goc(goc) {
			eqredundent = 0;
		};
	};

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
	nlopt_opt opt = (nlopt_opt) goc->extraData;
	FitContext *fc = goc->fc;
	assert(n == fc->numParam);
	int mode = 0;
	double fit = goc->solFun((double*) x, &mode);
	if (grad) {
		fc->iterations += 1;
		if (goc->maxMajorIterations != -1 && fc->iterations >= goc->maxMajorIterations) {
			nlopt_force_stop(opt);
		}
	}
	if (grad && goc->verbose >= 2) {
		mxLog("major iteration fit=%.12f", fit);
	}
	if (mode == -1) {
		if (!goc->feasible) {
			nlopt_force_stop(opt);
		}
		return nan("infeasible");
	}
	if (!grad) return fit;

	Eigen::Map< Eigen::VectorXd > Epoint((double*) x, n);
	Eigen::Map< Eigen::VectorXd > Egrad(grad, n);
	if (fc->wanted & FF_COMPUTE_GRADIENT) {
		Egrad = fc->grad;
	} else if (fc->CI && fc->CI->varIndex >= 0) {
		Egrad.setZero();
		Egrad[fc->CI->varIndex] = fc->lowerBound? 1 : -1;
		fc->grad = Egrad;
	} else {
		if (goc->verbose >= 3) mxLog("fd_gradient start");
		fit_functional ff(*goc);
		gradient_with_ref(goc->gradientAlgo, goc->gradientIterations, goc->gradientStepSize,
				  ff, fit, Epoint, Egrad);
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
	context &ctx = *(context *)f_data;
	GradientOptimizerContext &goc = ctx.goc;
	assert(n == goc.fc->numParam);
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::VectorXd Eresult(ctx.origeq);
	Eigen::MatrixXd jacobian(n, ctx.origeq);
	equality_functional ff(goc);
	ff(Epoint, Eresult);
	if (grad) {
		fd_jacobian(goc.gradientAlgo, goc.gradientIterations, goc.gradientStepSize,
			    ff, Eresult, Epoint, jacobian);
		if (ctx.eqmask.size() == 0) {
			ctx.eqmask.assign(m, false);
			for (int c1=0; c1 < int(m-1); ++c1) {
				for (int c2=c1+1; c2 < int(m); ++c2) {
					bool match = (Eresult[c1] == Eresult[c2] &&
						      (jacobian.col(c1) == jacobian.col(c2)));
					if (match && !ctx.eqmask[c2]) {
						ctx.eqmask[c2] = match;
						++ctx.eqredundent;
						if (goc.verbose >= 2) {
							mxLog("nlopt: eq constraint %d is redundent with %d",
							      c1, c2);
						}
					}
				}
			}
			if (ctx.eqredundent) {
				if (goc.verbose >= 1) {
					mxLog("nlopt: detected %d redundent equality constraints; retrying",
					      ctx.eqredundent);
				}
				nlopt_opt opt = (nlopt_opt) goc.extraData;
				nlopt_force_stop(opt);
			}
		}
	}
	Eigen::Map< Eigen::VectorXd > Uresult(result, m);
	Eigen::Map< Eigen::MatrixXd > Ujacobian(grad, n, m);
	int dx=0;
	for (int cx=0; cx < int(m); ++cx) {
		if (ctx.eqmask[cx]) continue;
		Uresult[dx] = Eresult[cx];
		if (grad) {
			Ujacobian.col(dx) = jacobian.col(cx);
		}
		++dx;
	}
	if (goc.verbose >= 4 && grad) {
		mxPrintMat("eq result", Uresult);
		mxPrintMat("eq jacobian", Ujacobian);
	}
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
	if (grad && goc->verbose >= 2) {
		if (m == 1) {
			mxLog("major iteration ineq=%.12f", Eresult[0]);
		} else {
			mxPrintMat("major iteration ineq", Eresult);
		}
	}
	inequality_functional ff(*goc);
	ff(Epoint, Eresult);
	if (grad) {
		fd_jacobian(goc->gradientAlgo, goc->gradientIterations, goc->gradientStepSize,
			    ff, Eresult, Epoint, jacobian);
		if (!std::isfinite(Eresult.sum())) {
			// infeasible at start of major iteration
			nlopt_opt opt = (nlopt_opt) goc->extraData;
			nlopt_force_stop(opt);
		}
	}
	if (goc->verbose >= 3 && grad) {
		mxPrintMat("inequality jacobian", jacobian);
	}
}

};

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
		nlopt_set_xtol_rel(opt, 5e-3);
		std::vector<double> tol(fc->numParam, std::numeric_limits<double>::epsilon());
		nlopt_set_xtol_abs(opt, tol.data());
	} else {
		// The *2 is there to roughly equate accuracy with NPSOL.
		nlopt_set_ftol_rel(opt, goc.ControlTolerance * 2);
		nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
	}
        
	nlopt_set_min_objective(opt, SLSQP::nloptObjectiveFunction, &goc);

	double feasibilityTolerance = Global->feasibilityTolerance;
	SLSQP::context ctx(goc);
        if (eq + ieq) {
		ctx.origeq = eq;
                if (ieq > 0){
			goc.inequality.resize(ieq);
			std::vector<double> tol(ieq, feasibilityTolerance);
			nlopt_add_inequality_mconstraint(opt, ieq, SLSQP::nloptInequalityFunction, &goc, tol.data());
                }
                
                if (eq > 0){
			goc.equality.resize(eq);
			std::vector<double> tol(eq, feasibilityTolerance);
			nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());
                }
	}
        
	int priorIterations = fc->iterations;

	int code = nlopt_optimize(opt, est, &fc->fit);
	if (ctx.eqredundent) {
		nlopt_remove_equality_constraints(opt);
		eq -= ctx.eqredundent;
		std::vector<double> tol(eq, feasibilityTolerance);
		nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());

		code = nlopt_optimize(opt, est, &fc->fit);
	}

	if (goc.verbose >= 2) mxLog("nlopt_optimize returned %d", code);

        nlopt_destroy(opt);

	fc->wanted = oldWanted;

	if (code == NLOPT_INVALID_ARGS) {
		Rf_error("NLOPT invoked with invalid arguments");
	} else if (code == NLOPT_OUT_OF_MEMORY) {
		Rf_error("NLOPT ran out of memory");
	} else if (code == NLOPT_FORCED_STOP) {
		if (!goc.feasible) {
			goc.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
		} else {
			goc.informOut = INFORM_ITERATION_LIMIT;
		}
	} else if (code == NLOPT_ROUNDOFF_LIMITED) {
		if (fc->iterations - priorIterations <= 2) {
			Rf_error("%s: Failed due to singular matrix E or C in LSQ subproblem or "
				 "rank-deficient equality constraint subproblem or "
				 "positive directional derivative in line search", goc.optName);
		} else {
			goc.informOut = INFORM_NOT_AT_OPTIMUM;  // is this correct? TODO
		}
	} else if (code < 0) {
		Rf_error("NLOPT fatal error %d", code);
	} else if (code == NLOPT_MAXEVAL_REACHED) {
		goc.informOut = INFORM_ITERATION_LIMIT;
	} else {
		goc.informOut = INFORM_CONVERGED_OPTIMUM;
	}
}
