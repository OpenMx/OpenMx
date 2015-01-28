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
void CSOLNPFit::allConstraintsFun(Eigen::MatrixBase<T1> &constraintOut)
{
	omxState *globalState = fc->state;
	int ncnln = globalState->ncnln;
	int l=0;
	for(int j = 0; j < globalState->numConstraints; j++) {
		omxRecompute(globalState->conList[j].result, fc);
		for(int k = 0; k < globalState->conList[j].size; k++){
			constraintOut[l++] = globalState->conList[j].result->data[k];
		}
	}
}

#if HAS_NPSOL
static const char* anonMatrix = "anonymous matrix";
static omxMatrix *NPSOL_fitMatrix = NULL;
static int NPSOL_currentInterval = -1;
static FitContext *NPSOL_fc = NULL;
static bool NPSOL_useGradient;
static int NPSOL_verbose;
static struct CSOLNPFit *NPSOL_GOpt;
#endif

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

// NOTE: All non-linear constraints are applied regardless of free TODO
// variable group.  This is probably wrong. TODO
static void omxSetupBoundsAndConstraints(FitContext *fc, double * bl, double * bu)
{
	FreeVarGroup *freeVarGroup = fc->varGroup;
	omxState *globalState = fc->state;
	size_t n = freeVarGroup->vars.size();

	/* Set min and max limits */
	for(size_t index = 0; index < n; index++) {
		bl[index] = freeVarGroup->vars[index]->lbound;
		bu[index] = freeVarGroup->vars[index]->ubound;
	}

	int index = n;
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {		// Nonlinear constraints:
		if(OMX_DEBUG) { mxLog("Constraint %d: ", constraintIndex);}
		omxConstraint::Type type = globalState->conList[constraintIndex].opCode;
		switch(type) {
		case omxConstraint::LESS_THAN:
			if(OMX_DEBUG) { mxLog("Bounded at (-INF, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = -0.0;
				index++;
			}
			break;
		case omxConstraint::EQUALITY:
			if(OMX_DEBUG) { mxLog("Bounded at (-0, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = -0.0;
				bu[index] = 0.0;
				index++;
			}
			break;
		case omxConstraint::GREATER_THAN:
			if(OMX_DEBUG) { mxLog("Bounded at (0, INF).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				if(OMX_DEBUG) { mxLog("\tBounds set for constraint %d.%d.", constraintIndex, offset);}
				bl[index] = 0.0;
				bu[index] = INF;
				index++;
			}
			break;
		default:
			Rf_error("Unknown constraint type %d", type);
		}
	}
}

void F77_SUB(npsolObjectiveFunction)(int* mode, int* n, double* x,
				     double* f, double* g, int* nstate)
{
	if (x != NPSOL_GOpt->fc->est) return;  // this is strange but necessary
	double fit = NPSOL_GOpt->solFun(x, mode);
	*f = fit;
}

/****** Objective Function *********/
static void
npsolObjectiveFunction1(int* mode, int* n, double* x,
                        double* f, double* g, double *hessian, int* nstate )
{
        if (*mode == 1) NPSOL_fc->iterations += 1;  //major iteration

        omxMatrix* fitMatrix = NPSOL_fitMatrix;

        // x == fc->est
        NPSOL_fc->copyParamToModel();

        int want = FF_COMPUTE_FIT;

        if (*mode > 0 && NPSOL_useGradient &&
            fitMatrix->fitFunction->gradientAvailable && NPSOL_currentInterval < 0) {
                NPSOL_fc->grad = Eigen::VectorXd::Zero(NPSOL_fc->numParam);
                want |= FF_COMPUTE_GRADIENT;
        }

        ComputeFit("NPSOL", fitMatrix, want, NPSOL_fc);

        *f = NPSOL_fc->fit;

        if (!std::isfinite(*f) || isErrorRaised()) {
                *mode = -1;
        }
}

void F77_SUB(npsolObjectiveFunctionOld)
        (       int* mode, int* n, double* x,
                double* f, double* g, int* nstate )
{
        npsolObjectiveFunction1(mode, n, x, f, g, NULL, nstate);
}

/* Objective function for confidence interval limit finding. 
 * Replaces the standard objective function when finding confidence intervals. */
void F77_SUB(npsolLimitObjectiveFunction)
	(	int* mode, int* n, double* x, double* f, double* g, int* nstate ) {
		
		F77_CALL(npsolObjectiveFunctionOld)(mode, n, x, f, g, nstate);	// Standard objective function call

		omxConfidenceInterval *oCI = Global->intervalList[NPSOL_currentInterval];
		
		omxRecompute(oCI->matrix, NPSOL_fc);
		
		double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);

		/* Catch boundary-passing condition */
		if(std::isnan(CIElement) || std::isinf(CIElement)) {
			NPSOL_fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			*mode = -1;
			return;
		}

		if(oCI->calcLower) {
			double diff = oCI->lbound - *f;		// Offset - likelihood
			*f = diff * diff + CIElement;
				// Minimize element for lower bound.
		} else {
			double diff = oCI->ubound - *f;			// Offset - likelihood
			*f = diff * diff - CIElement;
				// Maximize element for upper bound.
		}
}

/* (Non)Linear Constraint Functions */
void F77_SUB(npsolConstraintFunction)
	(	int *mode, int *ncnln, int *n,
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	if(OMX_DEBUG) { mxLog("Constraint function called.");}

	if(*mode==1) {
		if(OMX_DEBUG) {
			mxLog("But only gradients requested.  Returning.");
			mxLog("-=====================================================-");
		}

		return;
	}

	int j, k, l = 0;

	NPSOL_fc->copyParamToModel();

	omxState *globalState = NPSOL_fc->state;
	for(j = 0; j < globalState->numConstraints; j++) {
		omxRecompute(globalState->conList[j].result, NPSOL_fc);
		if(OMX_DEBUG) { omxPrint(globalState->conList[j].result, "Constraint evaluates as:"); }
		for(k = 0; k < globalState->conList[j].size; k++){
			c[l++] = globalState->conList[j].result->data[k];
		}
	}
	if(OMX_DEBUG) { mxLog("-=======================================================-"); }
}

void omxNPSOL(double *est, RegularFit &rf)
{
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
	NPSOL_fitMatrix = rf.fitMatrix;
	NPSOL_verbose = rf.verbose;

	NPSOL_useGradient = rf.useGradient;
	NPSOL_fc = rf.fc;
	NPSOL_GOpt = &rf;
	FitContext *fc = rf.fc;
	double *x = fc->est;
	fc->grad.resize(fc->numParam); // ensure memory is allocated
	double *g = fc->grad.data();

    omxState *globalState = fc->state;
    int nclin = 0;
    int nlinwid = std::max(1, nclin);
    int ncnln = globalState->ncnln;
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
 
	rf.hessOut.resize(n, n);
	double fit; // do not pass in &fc->fit
	int iter_out; // ignored
	F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A.data(),
			rf.solLB.data(), rf.solUB.data(), (void*)F77_SUB(npsolConstraintFunction),
			(void*) F77_SUB(npsolObjectiveFunction), &rf.informOut, &iter_out,
			istate.data(), c.data(), cJac.data(),
			clambda.data(), &fit, g, rf.hessOut.data(), x, iw.data(), &leniw, w.data(), &lenw);

    NPSOL_fitMatrix = NULL;
    NPSOL_fc = NULL;
    NPSOL_GOpt = NULL;
}
 
 
// Mostly duplicated code in omxCSOLNPConfidenceIntervals
// needs to be refactored so there is only 1 copy of CI
// code that can use whatever optimizer is provided.
void omxNPSOLConfidenceIntervals(omxMatrix *fitMatrix, FitContext *opt, double tolerance)
{
	if (std::isfinite(tolerance)) {
		std::string option = string_snprintf("Optimality tolerance %.8g", tolerance);
		F77_CALL(npoptn)((char*) option.c_str(), option.size());
	}
	{
		std::string option = string_snprintf("Cold start");
		F77_CALL(npoptn)((char*) option.c_str(), option.size());
	}

    FitContext fc(opt, opt->varGroup);
    fc.createChildren();

	int ciMaxIterations = Global->ciMaxIterations;
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) Rf_error("NPSOL is not reentrant");
	NPSOL_fitMatrix = fitMatrix;
	NPSOL_fc = &fc;
	FreeVarGroup *freeVarGroup = opt->varGroup;

    double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *w=NULL; //  *g, *R, *cJac,
 
    int ldA, ldJ, ldR, inform, iter, leniw, lenw; 
 
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
 
    int *iw = NULL;
 
    int *istate = NULL;                 // Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))
 
    int nctotl, nlinwid, nlnwid;    // Helpful side variables.
 
    int n = int(freeVarGroup->vars.size());
    double f = opt->fit;
    Eigen::VectorXd gradient(n);
    Eigen::MatrixXd hessian(n, n);

    /* NPSOL Arguments */
    void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);
 
    funcon = F77_SUB(npsolConstraintFunction);
 
    int nclin = 0;
    omxState *globalState = NPSOL_fc->state;
    int ncnln = globalState->ncnln;
 
    /* Set boundaries and widths. */
    if(nclin <= 0) {
        nclin = 0;                  // This should never matter--nclin should always be non-negative.
        nlinwid = 1;                // For memory allocation purposes, nlinwid > 0
    } else {                        // nlinwid is  used to calculate ldA, and eventually the size of A.
        nlinwid = nclin;
    }
 
    if(ncnln <= 0) {
        ncnln = 0;                  // This should never matter--ncnln should always be non-negative.
        nlnwid = 1;                 // For memory allocation purposes nlnwid > 0
    } else {                        // nlnwid is used to calculate ldJ, and eventually the size of J.
        nlnwid = ncnln;
    }
 
    nctotl = n + nlinwid + nlnwid;
 
    leniw = 3 * n + nclin + 2 * ncnln;
    lenw = 2 * n * n + n * nclin + 2 * n * ncnln + 20 * n + 11 * nclin + 21 * ncnln;
 
    ldA = nlinwid;          // NPSOL specifies this should be either 1 or nclin, whichever is greater
    ldJ = nlnwid;           // NPSOL specifies this should be either 1 or nclin, whichever is greater
    ldR = n;                // TODO: Test alternative versions of the size of R to see what's best.
 
    /* Allocate arrays */
    A       = (double*) R_alloc (ldA * n, sizeof ( double )  );
    bl      = (double*) R_alloc ( nctotl, sizeof ( double ) );
    bu      = (double*) R_alloc (nctotl, sizeof ( double ) );
    c       = (double*) R_alloc (nlnwid, sizeof ( double ));
    cJac    = (double*) R_alloc (ldJ * n, sizeof ( double ) );
    clambda = (double*) R_alloc (nctotl, sizeof ( double )  );
    w       = (double*) R_alloc (lenw, sizeof ( double ));
    istate  = (int*) R_alloc (nctotl, sizeof ( int ) );
    iw      = (int*) R_alloc (leniw, sizeof ( int ));
 
 
    omxSetupBoundsAndConstraints(opt, bl, bu);
 
        if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }

	const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?

        for(int i = 0; i < (int) Global->intervalList.size(); i++) {

		omxConfidenceInterval *currentCI = Global->intervalList[i];

		const char *matName = anonMatrix;
		if (currentCI->matrix->name) {
			matName = currentCI->matrix->name;
		}
		Global->checkpointMessage(opt, opt->est, "%s[%d, %d] begin lower interval",
					  matName, currentCI->row + 1, currentCI->col + 1);
 
 
		memcpy(fc.est, opt->est, n * sizeof(double)); // Reset to previous optimum
			NPSOL_currentInterval = i;

            currentCI->lbound += opt->fit;          // Convert from offsets to targets
            currentCI->ubound += opt->fit;          // Convert from offsets to targets
 
	    if (std::isfinite(currentCI->lbound)) {
            /* Set up for the lower bound */
            inform = -1;
            // Number of times to keep trying.
            int cycles = ciMaxIterations;
            double value = INF;
            while(inform != 0 && cycles > 0) {
                /* Find lower limit */
                currentCI->calcLower = TRUE;
                F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
                    (void*) F77_SUB(npsolLimitObjectiveFunction), &inform, &iter, istate, c, cJac,
				clambda, &f, gradient.data(), hessian.data(), fc.est, iw, &leniw, w, &lenw);
 
                currentCI->lCode = inform;
                if(f < value) {
                    currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
                }
 
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of lower interval %d failed: Bad inform value of %d",
                            i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opt->est[j]) > objDiff) {
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opt->est[j] + objDiff;
                        }
                    }
                }
            }
	    }
 
	    if (std::isfinite(currentCI->ubound)) {
		Global->checkpointMessage(opt, opt->est, "%s[%d, %d] begin upper interval",
					  matName, currentCI->row + 1, currentCI->col + 1);

		memcpy(fc.est, opt->est, n * sizeof(double)); // Reset to previous optimum
 
            /* Reset for the upper bound */
		double value = INF;
            inform = -1;
            double cycles = ciMaxIterations;
 
            while(inform != 0 && cycles > 0) {
                /* Find upper limit */
                currentCI->calcLower = FALSE;
                F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
                                    (void*) F77_SUB(npsolLimitObjectiveFunction), &inform, &iter, istate, c, cJac,
				clambda, &f, gradient.data(), hessian.data(), fc.est, iw, &leniw, w, &lenw);
 
                currentCI->uCode = inform;
                if(f < value) {
                    currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
                }
 
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of upper interval %d failed: Bad inform value of %d",
                            i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opt->est[j]) > objDiff){
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opt->est[j] + objDiff;
                        }
                    }
                }
            }
            if(OMX_DEBUG) {mxLog("Found Upper bound %d.", i);}
        }
	}

	NPSOL_fc = NULL;
	NPSOL_fitMatrix = NULL;
	NPSOL_currentInterval = -1;
}
 
void omxInvokeNPSOL(omxMatrix *fitMatrix, FitContext *fc,
                   int *inform_out, bool useGradient, FreeVarGroup *freeVarGroup,
                   int verbose, double *hessOut, double tolerance, bool warmStart)
{
       RegularFit rf("NPSOL", fc, fitMatrix, verbose);
       //rf.ControlMajorLimit = majIter;
       //rf.ControlMinorLimit = minIter;
       //rf.ControlFuncPrecision = funcPrecision;
       rf.ControlTolerance = tolerance;
       rf.warmStart = warmStart;
       Eigen::Map< Eigen::MatrixXd > hessWrap(hessOut, fc->numParam, fc->numParam);
       if (rf.warmStart) {
               rf.hessOut = hessWrap;
       }
       rf.useGradient = useGradient;
       rf.verbose = verbose;
       rf.setupAllBounds();
       omxNPSOL(fc->est, rf);
       *inform_out = rf.informOut;
       hessWrap = rf.hessOut;
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
