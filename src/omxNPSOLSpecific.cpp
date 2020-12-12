/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <ctype.h>

#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "omxImportFrontendState.h"
#include "Compute.h"
#include "ComputeGD.h"
#include "npsolswitch.h"
#include "EnableWarnings.h"

#ifdef  __cplusplus
extern "C" {
#endif

/* NPSOL-related functions */
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,
                            double *bl, double *bu, void* funcon, void* funobj, int *inform, int *iter,
                            int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,
                            double *x, int *iw, int *leniw, double *w, int *lenw);
extern void F77_SUB(npoptn)(char* string, int Rf_length);

#ifdef  __cplusplus
}
#endif

#if HAS_NPSOL

static class GradientOptimizerContext *NPSOL_GOpt;

void F77_SUB(npsolObjectiveFunction)(int* mode, int* n, double* x,
				     double* f, double* g, int* nstate)
{
	if (x != NPSOL_GOpt->est.data()) return;  // this is strange but necessary
	double fit = NPSOL_GOpt->recordFit(x, mode);
	*f = fit;

  if (*mode > 0) {
    for (int vx = 0; vx < *n; ++vx) {
      if (!std::isfinite(NPSOL_GOpt->grad[vx])) continue;;
      g[vx] = NPSOL_GOpt->grad[vx];
    }
  }
}

/* (Non)Linear Constraint Functions */
void F77_SUB(npsolConstraintFunction)
	(	int *mode, int *ncnln, int *n,
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{
	// "Note that if there are any nonlinear constraints then the
	// first call to CONFUN will precede the first call to
	// OBJFUN." Hence, we must copyParamToModel here.

	NPSOL_GOpt->copyFromOptimizer(x);
	Eigen::Map< Eigen::VectorXd > cE(c, *ncnln);

  // needc could improve performance on some models
	//Eigen::Map< Eigen::VectorXi > needcE(needc, *ncnln);

  int verbose = NPSOL_GOpt->verbose;
	if (mode == 0) {
    NPSOL_GOpt->evalAllC(c);
    if (verbose >= 3) mxPrintMat("constraints", cE);
  } else {
    Eigen::Map< Eigen::MatrixXd > cJacE(cJac, *ldJ, *n);
    if (NPSOL_GOpt->AllC.anyAnalyticJac()) {
      NPSOL_GOpt->evalAllC(c, cJac);
    } else {
      // Let NPSOL figure out the Jacobian, for backward compatibility
      NPSOL_GOpt->evalAllC(c);
    }
    if (verbose >= 3) {
      mxPrintMat("constraints", cE);
			mxPrintMat("jacobian", cJacE);
		}
	}
}

static double getNPSOLFeasibilityTolerance()
{
	// convert units (who knows what the units are)
	return Global->feasibilityTolerance * 2. / 5.;
}

struct NPSOLModeSwitch {
  NPSOLModeSwitch() { Global->NPSOL_HACK += 1; }
  ~NPSOLModeSwitch() { Global->NPSOL_HACK -= 1; }
};

static void omxNPSOL1(double *est, GradientOptimizerContext &rf, int nclin, int ncnln)
{
  NPSOLModeSwitch modeSwitch;
	rf.setEngineName("NPSOL");
	rf.setupAllBounds();
	{
		double ft = (nclin + ncnln)? getNPSOLFeasibilityTolerance() : 1e-5;
		std::string opt = string_snprintf("Feasibility tolerance %.8g", ft);
		F77_CALL(npoptn)((char*) opt.c_str(), opt.size());
	}
	if (std::isfinite(rf.ControlTolerance)) {
		std::string opt = string_snprintf("Optimality tolerance %.8g", rf.ControlTolerance);
		F77_CALL(npoptn)((char*) opt.c_str(), opt.size());
	}
	if (rf.maxMajorIterations != -1) {
		std::string opt = string_snprintf("Major iterations %d",
						  rf.maxMajorIterations);
		F77_CALL(npoptn)((char*) opt.c_str(), opt.size());
	}
	if (rf.warmStart) {
		std::string opt = string_snprintf("Warm start");
		F77_CALL(npoptn)((char*) opt.c_str(), opt.size());
	} else {
		std::string opt = string_snprintf("Cold start");
		F77_CALL(npoptn)((char*) opt.c_str(), opt.size());
	}

	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) mxThrow("NPSOL is not reentrant");
	NPSOL_GOpt = &rf;

	if (ncnln + nclin == 0) { //<--We might have to move to a worse fit to satisfy even linear constraints.
		// ensure we never move to a worse point
		int mode = 0;
		double fit = rf.recordFit(rf.est.data(), &mode);
		if (!std::isfinite(fit)) {
			rf.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
			NPSOL_GOpt = NULL;
			return;
		}
		if (mode == -1) {
			NPSOL_GOpt = NULL;
			return;
		}
	}

	int n = rf.numFree;
    int nlinwid = std::max(1, nclin);
    int nlnwid = std::max(1, ncnln);
        int nctotl = n + nlinwid + nlnwid;

        int leniw = 3 * n + nclin + 2 * ncnln;
        int lenw = 2 * n * n + n * nclin + 2 * n * ncnln + 20 * n + 11 * nclin + 21 * ncnln;

        int ldA = nlinwid;          // NPSOL specifies this should be either 1 or nclin, whichever is greater
        int ldJ = nlnwid;           // NPSOL specifies this should be either 1 or nclin, whichever is greater
        int ldR = n;                // TODO: Test alternative versions of the size of R to see what's best.

    /* Allocate arrays */
	Eigen::MatrixXd A(ldA, n);  // maybe transposed?
	if(nclin){
    // Not implemented yet
		// A.setZero();
		// rf.linearConstraintCoefficients(A);
	}
	rf.constraintFunValsOut.resize(nlnwid);//Eigen::VectorXd c(nlnwid);
	rf.constraintJacobianOut.resize(ldJ, n);//Eigen::MatrixXd cJac(ldJ, n);
	rf.LagrMultipliersOut.resize(nctotl);//Eigen::VectorXd clambda(nctotl);
	Eigen::VectorXd w(lenw);
	rf.constraintStatesOut.resize(nctotl);//Eigen::VectorXi istate(nctotl);
	Eigen::VectorXi iw(leniw);

	if (rf.warmStart) {
		rf.constraintStatesOut.setZero();
		rf.LagrMultipliersOut.setZero();
	}

    /*  F77_CALL(npsol)
        (   int *n,                 -- Number of variables
            int *nclin,             -- Number of linear constraints
            int *ncnln,             -- Number of nonlinear constraints
            int *ldA,               -- Row dimension of A (Linear Constraints)
            int *ldJ,               -- Row dimension of cJac (Jacobian)
            int *ldR,               -- Row dimension of R (Hessian)
            double *A,              -- Linear Constraints Array A (in Column-major order)
            double *bl,             -- Lower Bounds Array (at least n + nclin + ncnln long)
            double *bu,             -- Upper Bounds Array (at least n + nclin + ncnln long)
            function funcon,        -- Nonlinear constraint function
            function funobj,        -- Objective function
            int *inform,            -- Used to report state.  Need not be initialized.
            int *iter,              -- Used to report number of major iterations performed.  Need not be initialized.
            int *istate,            -- Initial State.  Need not be initialized unless using Warm Start.
            double *c,              -- Array of length ncnln.  Need not be initialized.  Reports nonlinear constraints at final iteration.
            double *cJac,           -- Array of Row-length ldJ.  Unused if ncnln = 0. Generally need not be initialized.
            double *clambda,        -- Array of length n+nclin+ncnln.  Need not be initialized unless using Warm Start. Reports final QP multipliers.
            double *f,              -- Used to report final objective value.  Need not be initialized.
            double *g,              -- Array of length n. Used to report final objective gradient.  Need not be initialized.
            double *R,              -- Array of length ldR.  Need not be initialized unless using Warm Start.
            double *x,              -- Array of length n.  Contains initial solution estimate.
            int *iw,                -- Array of length leniw. Need not be initialized.  Provides workspace.
            int *leniw,             -- Length of iw.  Must be at least 3n + nclin + ncnln.
            double *w,              -- Array of length lenw. Need not be initialized.  Provides workspace.
            int *lenw               -- Length of w.  Must be at least 2n^2 + n*nclin + 2*n*ncnln + 20*n + 11*nclin +21*ncnln
        )

        bl, bu, istate, and clambda are all length n+nclin+ncnln.
            First n elements refer to the vars, in order.
            Next nclin elements refer to bounds on Ax
            Last ncnln elements refer to bounds on c(x)

        All arrays must be in column-major order.
        */

  Eigen::ArrayXd workingGrad(n);
  workingGrad.setZero();
	rf.hessOut.resize(n, n);
	double fit; // do not pass in &fc->fit
	int iter_out; // ignored
	F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR,
          A.data(), rf.solLB.data(), rf.solUB.data(),
          (void*)F77_SUB(npsolConstraintFunction), (void*) F77_SUB(npsolObjectiveFunction),
          &rf.informOut, &iter_out, rf.constraintStatesOut.data(),
          rf.constraintFunValsOut.data(), rf.constraintJacobianOut.data(), rf.LagrMultipliersOut.data(), &fit, workingGrad.data(), rf.hessOut.data(), rf.est.data(),
          iw.data(), &leniw, w.data(), &lenw);

	/*We tell NPSOL that all of the MxConstraints are nonlinear, because NPSOL's special handling of linear constraints doesn't seem to
	work as documented.  However, nlinwid = max(1,nclin), and nctotl = n + nlinwid + nlnwid, so arrays that are sized to nctotl will always
	have at least one element that NPSOL never uses.  Obviously, if ncnln==0, the only constraints will be box constraints, and there will be a
	lot of unused elements in the constraint-related arrays:*/
	if(ncnln==0){
		rf.constraintFunValsOut.resize(0);
		rf.constraintJacobianOut.resize(0,0);
		rf.LagrMultipliersOut.conservativeResize(n);
		rf.constraintStatesOut.conservativeResize(n);
	}
	else{
		rf.LagrMultipliersOut.conservativeResize(n+ncnln);
		rf.constraintStatesOut.conservativeResize(n+ncnln);
	}

	// NPSOL can return the wrong fit and estimates, but hard to
	// know what to do if there are constraints.
	if (rf.bestEst.size() && ncnln+nclin == 0) {
		rf.useBestFit();
	}

	// With constraints, it's probably garbage,
	//if (ncnln) rf.grad.setConstant(NA_REAL);

    NPSOL_GOpt = NULL;
}

void omxNPSOL(GradientOptimizerContext &rf)
{
	double *est = rf.est.data();
	Eigen::Map< Eigen::ArrayXd > Est(est, rf.numFree);
	Eigen::ArrayXd startingPoint = Est;

	int nclin=0, ncnln = rf.AllC.getCount();

	omxNPSOL1(est, rf, nclin, ncnln);

	if (nclin + ncnln == 0) return;

	const int maxRetries = 10;
	int retry = 0;
	double best = std::numeric_limits<double>::max();
	while (++retry < maxRetries) {
		Eigen::VectorXd cE(ncnln);
    rf.evalAllC(cE.data());

		double norm = cE.norm();
		if (rf.verbose >= 1) {
			mxLog("NPSOL[%d]: fit %f feasibility constraints at %g (best %g)",
			      retry, rf.getFit(), norm, best);
		}
		if (norm >= best) {
			// NPSOL can jump off a cliff and get lost.
			// Try to move back toward a more feasible region.
			Est = (Est + startingPoint) / 2;
		} else {
			best = norm;
		}
		if (!(cE.array().abs() < getNPSOLFeasibilityTolerance()).all()) {
			omxNPSOL1(est, rf, nclin, ncnln);
		} else {
			break;
		}
	}
}

void omxSetNPSOLOpts(SEXP options)
{
    ProtectAutoBalanceDoodad mpi;

    static const char *whitelist[] = {
        "Central Difference Interval",
        "Crash Tolerance",
        "Derivative level",
        "Difference interval",
        "Feasibility tolerance",
        "Function precision",
        "Hessian",
        "Infinite bound size",
        "Infinite step size",
        "Iteration limit",
        "Iters",
        "Line search tolerance",
        "Major iteration limit",
        "Major iterations",
        "Print level",
        "Print file",
        "Minor iteration limit",
        "Minor print level",
        "Nolist",
        "Optimality tolerance",
        "Step limit",
        "Summary file",
        "Verify level",
        0
    };

    const int opBufLen = 250;
    char optionCharArray[opBufLen];
    int numOptions = Rf_length(options);
    SEXP optionNames;
    Rf_protect(optionNames = Rf_getAttrib(options, R_NamesSymbol));
    for(int i = 0; i < numOptions; i++) {
        const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
        const char *nextOptionValue = CHAR(Rf_asChar(VECTOR_ELT(options, i)));
        bool ok=false;
        for (int wx=0; whitelist[wx]; ++wx) {
            if (matchCaseInsensitive(nextOptionName, whitelist[wx])) {
                ok=true;
                break;
            }
        }
        if (!ok) continue;
        snprintf(optionCharArray, opBufLen, "%s %s", nextOptionName, nextOptionValue);
        F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
        if(OMX_DEBUG) { mxLog("Option %s ", optionCharArray);}
    }
}
#endif
