/*
 *  Copyright 2007-2015 The OpenMx Project
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
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "glue.h"
#include "omxImportFrontendState.h"
#include "Compute.h"
#include "npsolswitch.h"
#include "omxBuffer.h"

template <typename T1>
void GradientOptimizerContext::allConstraintsFun(Eigen::MatrixBase<T1> &constraintOut)
{
	omxState *globalState = fc->state;
	int l=0;
	for(int j = 0; j < (int) globalState->conList.size(); j++) {
		omxConstraint &cs = *globalState->conList[j];
		cs.refreshAndGrab(fc, omxConstraint::LESS_THAN, &constraintOut(l));
		l += cs.size;
	}
}

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

static struct GradientOptimizerContext *NPSOL_GOpt;

void F77_SUB(npsolObjectiveFunction)(int* mode, int* n, double* x,
				     double* f, double* g, int* nstate)
{
	if (x != NPSOL_GOpt->fc->est) return;  // this is strange but necessary
	double fit = NPSOL_GOpt->recordFit(x, mode);
	*f = fit;
}

/* (Non)Linear Constraint Functions */
void F77_SUB(npsolConstraintFunction)
	(	int *mode, int *ncnln, int *n,
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{
	if(*mode==1) return;

	// "Note that if there are any nonlinear constraints then the
	// first call to CONFUN will precede the first call to
	// OBJFUN." Hence, we must copyParamToModel here.

	NPSOL_GOpt->fc->copyParamToModel();

	Eigen::Map< Eigen::VectorXd > cE(c, *ncnln);
	NPSOL_GOpt->allConstraintsFun(cE);
}

void omxNPSOL(double *est, GradientOptimizerContext &rf)
{
	rf.optName = "NPSOL";
	rf.setupAllBounds();
	if (std::isfinite(rf.ControlTolerance)) {
		std::string opt = string_snprintf("Optimality tolerance %.8g", rf.ControlTolerance);
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
	//if (NPSOL_fitMatrix) Rf_error("NPSOL is not reentrant");
	FitContext *fc = rf.fc;
	NPSOL_GOpt = &rf;
	{
		// ensure we never move to a worse point
		int mode = 0;
		double fit = rf.recordFit(fc->est, &mode);
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
	fc->grad.resize(fc->numParam); // ensure memory is allocated

    omxState *globalState = fc->state;
    int nclin = 0;
    int nlinwid = std::max(1, nclin);
    int equality, inequality;
    globalState->countNonlinearConstraints(equality, inequality);
    int ncnln = equality + inequality;
    int nlnwid = std::max(1, ncnln);
 
	int n = int(fc->numParam);
 
        int nctotl = n + nlinwid + nlnwid;
 
        int leniw = 3 * n + nclin + 2 * ncnln;
        int lenw = 2 * n * n + n * nclin + 2 * n * ncnln + 20 * n + 11 * nclin + 21 * ncnln;
 
        int ldA = nlinwid;          // NPSOL specifies this should be either 1 or nclin, whichever is greater
        int ldJ = nlnwid;           // NPSOL specifies this should be either 1 or nclin, whichever is greater
        int ldR = n;                // TODO: Test alternative versions of the size of R to see what's best.
 
    /* Allocate arrays */
	Eigen::ArrayXXd A(ldA, n);  // maybe transposed?
	Eigen::VectorXd c(nlnwid);
	Eigen::MatrixXd cJac(ldJ, n); // maybe transposed?
	Eigen::VectorXd clambda(nctotl);
	Eigen::VectorXd w(lenw);
	Eigen::VectorXi istate(nctotl);
	Eigen::VectorXi iw(leniw);
 
	if (rf.warmStart) {
		istate.setZero();
		clambda.setZero();
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
            double *R,              -- Array of length ldR.  Need not be intialized unless using Warm Start.
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
 
	fc->grad.resize(n);
	rf.hessOut.resize(n, n);
	double fit; // do not pass in &fc->fit
	int iter_out; // ignored
	F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A.data(),
			rf.solLB.data(), rf.solUB.data(), (void*)F77_SUB(npsolConstraintFunction),
			(void*) F77_SUB(npsolObjectiveFunction), &rf.informOut, &iter_out,
			istate.data(), c.data(), cJac.data(),
			clambda.data(), &fit, fc->grad.data(), rf.hessOut.data(), fc->est, iw.data(), &leniw, w.data(), &lenw);

	// NPSOL can return the wrong fit and estimates, but hard to
	// know what to do if there are constraints.
	if (rf.bestEst.size() && ncnln == 0) {
		fc->fit = rf.bestFit;
		memcpy(fc->est, rf.bestEst.data(), sizeof(double) * fc->numParam);
	}

    NPSOL_GOpt = NULL;
}

void omxSetNPSOLOpts(SEXP options)
{
    omxManageProtectInsanity mpi;
    
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
