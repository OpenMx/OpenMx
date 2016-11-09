#include <limits>
#include <cstdlib>
#include <ctype.h>
#include "omxDefines.h"
#include "omxState.h"
#include "omxMatrix.h"
#include "glue.h"
#include "nloptcpp.h"
#include "finiteDifferences.h"

#include "nlopt.h"

#ifdef SHADOW_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#endif

namespace SLSQP {

	struct context {
		GradientOptimizerContext &goc;
		int origeq;
		int eqredundent;
		std::vector<bool> eqmask;
		context(GradientOptimizerContext &_goc) : goc(_goc) {
			eqredundent = 0;
		};
	};

struct fit_functional {
	GradientOptimizerContext &goc;

	fit_functional(GradientOptimizerContext &_goc) : goc(_goc) {};

	double operator()(double *x, int thrId) const {
		int mode = 0;
		return goc.evalFit(x, thrId, &mode);
	}
};

static double nloptObjectiveFunction(unsigned n, const double *x, double *grad, void *f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	nlopt_opt opt = (nlopt_opt) goc->extraData;
	int mode = grad != 0;
	double fit = goc->solFun((double*) x, &mode);
	if (grad) {
		if (goc->maxMajorIterations != -1 && goc->getIteration() >= goc->maxMajorIterations) {
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
	if (goc->getWanted() & FF_COMPUTE_GRADIENT) {
		Egrad = goc->grad;
	} else if (goc->hasKnownGradient()) {
		goc->setKnownGradient(Egrad);
		goc->grad = Egrad;
	} else {
		if (goc->verbose >= 3) mxLog("fd_gradient start");
		fit_functional ff(*goc);
		gradient_with_ref(goc->gradientAlgo, goc->numOptimizerThreads,
				  goc->gradientIterations, goc->gradientStepSize,
				  ff, fit, Epoint, Egrad);
		goc->grad = Egrad;
	}
	if (goc->verbose >= 3) {
		mxPrintMat("gradient", Egrad);
	}
	return fit;
}

struct equality_functional {
	GradientOptimizerContext &goc;

	equality_functional(GradientOptimizerContext &_goc) : goc(_goc) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.solEqBFun();
		result = goc.equality;
	}
};

static void nloptEqualityFunction(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
{
	context &ctx = *(context *)f_data;
	GradientOptimizerContext &goc = ctx.goc;
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::VectorXd Eresult(ctx.origeq);
	Eigen::MatrixXd jacobian(n, ctx.origeq);
	equality_functional ff(goc);
	ff(Epoint, Eresult);
	if (grad) {
		goc.eqNorm = Eresult.array().abs().sum();
		fd_jacobian(goc.gradientAlgo, goc.gradientIterations, goc.gradientStepSize,
			    ff, Eresult, Epoint, jacobian);
		if (ctx.eqmask.size() == 0) {
			ctx.eqmask.assign(m, false);
			for (int c1=0; c1 < int(m-1); ++c1) {
				for (int c2=c1+1; c2 < int(m); ++c2) {
					bool match = (Eresult[c1] == Eresult[c2] &&
						      (jacobian.col(c1) == jacobian.col(c2)));
					if (match && !ctx.eqmask[c2]) {
						ctx.eqmask[c2] = true;
						++ctx.eqredundent;
						if (goc.verbose >= 2) {
							mxLog("nlopt: eq constraint %d is redundent with %d",
							      c1, c2);
						}
					}
				}
			}
			for (int c1=0; c1 < int(m); ++c1) {
				if (ctx.eqmask[c1]) continue;
				if ((jacobian.col(c1).array() == 0.0).all()) {
					ctx.eqmask[c1] = true;
					++ctx.eqredundent;
					if (goc.verbose >= 2) {
						mxLog("nlopt: eq constraint %d is never active", c1);
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
	if (goc.verbose >= 3 && grad) {
		mxPrintMat("eq jacobian", Ujacobian);
	}
}

struct inequality_functional {
	GradientOptimizerContext &goc;

	inequality_functional(GradientOptimizerContext &_goc) : goc(_goc) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.myineqFun();
		result = goc.inequality;
	}
};

static void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobian(grad, n, m);
	inequality_functional ff(*goc);
	ff(Epoint, Eresult);
	if (grad && goc->verbose >= 2) {
		mxPrintMat("major iteration ineq", Eresult);
	}
	if (grad) {
		goc->ineqNorm = Eresult.array().abs().sum();
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

void omxInvokeNLOPT(GradientOptimizerContext &goc)
{
	double *est = goc.est.data();
	goc.optName = "SLSQP";
	goc.setupSimpleBounds();

	int oldWanted = goc.getWanted();
	goc.setWanted(0);
	omxState *globalState = goc.getState();
    
        nlopt_opt opt = nlopt_create(NLOPT_LD_SLSQP, goc.numFree);
	goc.extraData = opt;
        //local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm
        
        //nlopt_set_local_optimizer(opt, local_opt);
        nlopt_set_lower_bounds(opt, goc.solLB.data());
        nlopt_set_upper_bounds(opt, goc.solUB.data());

	int eq, ieq;
	globalState->countNonlinearConstraints(eq, ieq, false);

	// The *2 is there to roughly equate accuracy with NPSOL.
	nlopt_set_ftol_rel(opt, goc.ControlTolerance * 2);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
        
	nlopt_set_min_objective(opt, SLSQP::nloptObjectiveFunction, &goc);

	double feasibilityTolerance = Global->feasibilityTolerance;
	SLSQP::context ctx(goc);
        if (eq + ieq) {
		ctx.origeq = eq;
                if (ieq > 0){
			if (goc.verbose >= 2) mxLog("%d inequality constraints", ieq);
			goc.inequality.resize(ieq);
			std::vector<double> tol(ieq, feasibilityTolerance);
			nlopt_add_inequality_mconstraint(opt, ieq, SLSQP::nloptInequalityFunction, &goc, tol.data());
                }
                
                if (eq > 0){
			if (goc.verbose >= 2) mxLog("%d equality constraints", eq);
			goc.equality.resize(eq);
			std::vector<double> tol(eq, feasibilityTolerance);
			nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());
                }
	}
        
	int priorIterations = goc.getIteration();

	double fit = 0;
	int code = nlopt_optimize(opt, est, &fit);
	if (ctx.eqredundent) {
		nlopt_remove_equality_constraints(opt);
		eq -= ctx.eqredundent;
		std::vector<double> tol(eq, feasibilityTolerance);
		nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());

		code = nlopt_optimize(opt, est, &fit);
	}

	if (goc.verbose >= 2) mxLog("nlopt_optimize returned %d", code);

        nlopt_destroy(opt);

	goc.setWanted(oldWanted);

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
		if (goc.eqNorm > feasibilityTolerance || goc.ineqNorm > feasibilityTolerance) {
			goc.informOut = INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE;
		} else if (goc.getIteration() - priorIterations <= 2) {
			Rf_error("%s: Failed due to singular matrix E or C in LSQ subproblem or "
				 "rank-deficient equality constraint subproblem or "
				 "positive directional derivative in line search "
				 "(eq %.4g ineq %.4g)", goc.optName, goc.eqNorm, goc.ineqNorm);
		} else {
			goc.informOut = INFORM_NOT_AT_OPTIMUM;  // is this correct? TODO
		}
	} else if (code < 0) {
		Rf_error("NLOPT unrecognized error %d; please report to developers", code);
	} else if (code == NLOPT_MAXEVAL_REACHED) {
		goc.informOut = INFORM_ITERATION_LIMIT;
	} else {
		goc.informOut = INFORM_CONVERGED_OPTIMUM;
	}
}
