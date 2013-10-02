/*
 *  Copyright 2007-2013 The OpenMx Project
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
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "glue.h"
#include "omxImportFrontendState.h"
#include "Compute.h"
#include "npsolswitch.h"

/* NPSOL-specific globals */
const double NPSOL_BIGBND = 1e20;

const char* anonMatrix = "anonymous matrix";
static omxMatrix *NPSOL_fitMatrix = NULL;
static int NPSOL_currentInterval = -1;
static FitContext *NPSOL_fc = NULL;
static bool NPSOL_useGradient;
static int NPSOL_verbose;

#ifdef  __cplusplus
extern "C" {
#endif

/* NPSOL-related functions */
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,
                            double *bl, double *bu, void* funcon, void* funobj, int *inform, int *iter, 
                            int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,
                            double *x, int *iw, int *leniw, double *w, int *lenw);
extern void F77_SUB(npoptn)(char* string, int length);

#ifdef  __cplusplus
}
#endif

// NOTE: All non-linear constraints are applied regardless of free TODO
// variable group.  This is probably wrong. TODO
void omxSetupBoundsAndConstraints(FreeVarGroup *freeVarGroup, double * bl, double * bu)
{
	size_t n = freeVarGroup->vars.size();

	/* Set min and max limits */
	for(size_t index = 0; index < n; index++) {
		bl[index] = freeVarGroup->vars[index]->lbound;
		bu[index] = freeVarGroup->vars[index]->ubound;
	}

	int index = n;
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {		// Nonlinear constraints:
		if(OMX_DEBUG) { mxLog("Constraint %d: ", constraintIndex);}
		switch(globalState->conList[constraintIndex].opCode) {
		case 0:									// Less than: Must be strictly less than 0.
			if(OMX_DEBUG) { mxLog("Bounded at (-INF, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = -0.0;
				index++;
			}
			break;
		case 1:									// Equal: Must be roughly equal to 0.
			if(OMX_DEBUG) { mxLog("Bounded at (-0, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = -0.0;
				bu[index] = 0.0;
				index++;
			}
			break;
		case 2:									// Greater than: Must be strictly greater than 0.
			if(OMX_DEBUG) { mxLog("Bounded at (0, INF).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				if(OMX_DEBUG) { mxLog("\tBounds set for constraint %d.%d.", constraintIndex, offset);}
				bl[index] = 0.0;
				bu[index] = INF;
				index++;
			}
			break;
		default:
			if(OMX_DEBUG) { mxLog("Bounded at (-INF, INF).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = INF;
				index++;
			}
			break;
		}
	}
}

#if HAS_NPSOL

/****** Objective Function *********/
static void
npsolObjectiveFunction1(int* mode, int* n, double* x,
			double* f, double* g, double *hessian, int* nstate )
{
	unsigned short int checkpointNow = FALSE;

	if(OMX_DEBUG) {mxLog("Starting Objective Run.");}

	if(*mode == 1) {
		omxSetMajorIteration(globalState, globalState->majorIteration + 1);
		omxSetMinorIteration(globalState, 0);
		checkpointNow = TRUE;					// Only checkpoint at major iterations.
	} else omxSetMinorIteration(globalState, globalState->minorIteration + 1);

	omxMatrix* fitMatrix = NPSOL_fitMatrix;
	omxResetStatus(globalState);						// Clear Error State recursively
	/* Interruptible? */
	R_CheckUserInterrupt();

	NPSOL_fc->copyParamToModel(globalState, x);

	if (*mode > 0 && NPSOL_useGradient &&
	    fitMatrix->fitFunction->gradientAvailable && NPSOL_currentInterval < 0) {
		size_t numParams = NPSOL_fc->varGroup->vars.size();
		OMXZERO(NPSOL_fc->grad, numParams);

		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT|FF_COMPUTE_GRADIENT, NPSOL_fc);
		if (NPSOL_verbose) {
			NPSOL_fc->log("NPSOL", FF_COMPUTE_FIT|FF_COMPUTE_ESTIMATE|FF_COMPUTE_GRADIENT);
		}
	} else {
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, NPSOL_fc);
		if (NPSOL_verbose) {
			NPSOL_fc->log("NPSOL", FF_COMPUTE_FIT|FF_COMPUTE_ESTIMATE);
		}
	}

	omxExamineFitOutput(globalState, fitMatrix, mode);

	if (isErrorRaised(globalState)) {
		if(OMX_DEBUG) {
			mxLog("Error status reported.");
		}
		*mode = -1;
	}

	*f = fitMatrix->data[0];
	if(OMX_VERBOSE) {
		mxLog("Fit function value is: %f, Mode is %d.", fitMatrix->data[0], *mode);
	}

	if(OMX_DEBUG) { mxLog("-======================================================-"); }

	if(checkpointNow && globalState->numCheckpoints != 0) {	// If it's a new major iteration
		omxSaveCheckpoint(x, *f, FALSE);		// Check about saving a checkpoint
	}

}

void F77_SUB(npsolObjectiveFunction)
	(	int* mode, int* n, double* x,
		double* f, double* g, int* nstate )
{
	npsolObjectiveFunction1(mode, n, x, f, g, NULL, nstate);
}

/* Objective function for confidence interval limit finding. 
 * Replaces the standard objective function when finding confidence intervals. */
void F77_SUB(npsolLimitObjectiveFunction)
	(	int* mode, int* n, double* x, double* f, double* g, int* nstate ) {
		
		if(OMX_VERBOSE) mxLog("Calculating interval %d, %s boundary:", NPSOL_currentInterval, (Global->intervalList[NPSOL_currentInterval].calcLower?"lower":"upper"));

		F77_CALL(npsolObjectiveFunction)(mode, n, x, f, g, nstate);	// Standard objective function call

		omxConfidenceInterval *oCI = &(Global->intervalList[NPSOL_currentInterval]);
		
		omxRecompute(oCI->matrix);
		
		double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);

		if(OMX_DEBUG) {
			mxLog("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.",
				oCI->lbound, oCI->ubound, *f, CIElement);
		}

		/* Catch boundary-passing condition */
		if(isnan(CIElement) || isinf(CIElement)) {
			omxRaiseError(globalState, -1, 
				"Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
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

		if(OMX_DEBUG) {
			mxLog("Interval fit function in previous iteration was calculated to be %f.", *f);
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

	NPSOL_fc->copyParamToModel(globalState, x);

	for(j = 0; j < globalState->numConstraints; j++) {
		omxRecompute(globalState->conList[j].result);
		if(OMX_VERBOSE) { omxPrint(globalState->conList[j].result, "Constraint evaluates as:"); }
		for(k = 0; k < globalState->conList[j].size; k++){
			c[l++] = globalState->conList[j].result->data[k];
		}
	}
	if(OMX_DEBUG) { mxLog("-=======================================================-"); }
}

void omxInvokeNPSOL(omxMatrix *fitMatrix, FitContext *fc,
		    int *inform_out, int *iter_out, bool useGradient, FreeVarGroup *freeVarGroup,
		    int verbose)
{
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) error("NPSOL is not reentrant");
	NPSOL_fitMatrix = fitMatrix;
	NPSOL_verbose = verbose;

	NPSOL_useGradient = useGradient;
	NPSOL_fc = fc;
	double *x = fc->est;
	double *g = fc->grad;

    double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *w=NULL; //  *g, *R, *cJac,
 
    int k, ldA, ldJ, ldR, inform, iter, leniw, lenw; 
 
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
 
    int *iw = NULL;
 
    int *istate = NULL;                 // Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))
 
    int nctotl, nlinwid, nlnwid;    // Helpful side variables.
 
    int nclin = 0;
    int ncnln = globalState->ncnln;
 
    /* NPSOL Arguments */
    void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);
 
    funcon = F77_SUB(npsolConstraintFunction);
 
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
 
	int n = int(freeVarGroup->vars.size());
 
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
 
        /* Set up actual run */
 
        omxSetupBoundsAndConstraints(freeVarGroup, bl, bu);
 
        /* Initialize Starting Values */
        if(OMX_VERBOSE) {
            mxLog("--------------------------");
            mxLog("Starting Values (%d) are:", n);
        }
        for(k = 0; k < n; k++) {
            if((x[k] == 0.0)) {
                x[k] += 0.1;
            }
            if(OMX_VERBOSE) { mxLog("%d: %f", k, x[k]); }
        }
        if(OMX_DEBUG) {
            mxLog("--------------------------");
            mxLog("Setting up optimizer...");
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
 
        if(OMX_DEBUG) {
            mxLog("Set.");
        }
 
	std::vector<double> hessOut(n * n);

	F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
			(void*) F77_SUB(npsolObjectiveFunction), &inform, &iter, istate, c, cJac,
			clambda, &fc->fit, g, hessOut.data(), x, iw, &leniw, w, &lenw);

	memcpy(fc->hess, hessOut.data(), sizeof(double) * n * n);

        if(OMX_DEBUG) { mxLog("Final Objective Value is: %f", fc->fit); }
 
        omxSaveCheckpoint(x, fc->fit, TRUE);
 
	NPSOL_fc->copyParamToModel(globalState);
 
    *inform_out = inform;
    *iter_out   = iter;
 
    NPSOL_fitMatrix = NULL;
    NPSOL_fc = NULL;
}
 
 
void omxNPSOLConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc)
{
	int ciMaxIterations = Global->ciMaxIterations;
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) error("NPSOL is not reentrant");
	NPSOL_fitMatrix = fitMatrix;
	NPSOL_fc = fc;
	FreeVarGroup *freeVarGroup = fitMatrix->fitFunction->freeVarGroup;

    double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *w=NULL; //  *g, *R, *cJac,
 
    int ldA, ldJ, ldR, inform, iter, leniw, lenw; 
 
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
 
    int *iw = NULL;
 
    int *istate = NULL;                 // Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))
 
    int nctotl, nlinwid, nlnwid;    // Helpful side variables.
 
    int n = int(freeVarGroup->vars.size());
    double optimum = fc->fit;
    double *optimalValues = fc->est;
    double f = optimum;
    std::vector< double > x(n, *optimalValues);
    std::vector< double > gradient(n);
    std::vector< double > hessian(n * n);

    /* NPSOL Arguments */
    void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);
 
    funcon = F77_SUB(npsolConstraintFunction);
 
    int nclin = 0;
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
 
 
    omxSetupBoundsAndConstraints(freeVarGroup, bl, bu);
 
        if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }

        for(int i = 0; i < Global->numIntervals; i++) {

			omxConfidenceInterval *currentCI = &(Global->intervalList[i]);

			int msgLength = 45;
 
			if (currentCI->matrix->name == NULL) {
				msgLength += strlen(anonMatrix);
			} else {
				msgLength += strlen(currentCI->matrix->name);
			}
            
			char *message = Calloc(msgLength, char);
 
			if (currentCI->matrix->name == NULL) {
				snprintf(message, msgLength, "%s[%d, %d] begin lower interval",
					anonMatrix, currentCI->row + 1, currentCI->col + 1);
			} else {
				snprintf(message, msgLength, "%s[%d, %d] begin lower interval",
					currentCI->matrix->name, currentCI->row + 1, currentCI->col + 1);
			}
 
			omxWriteCheckpointMessage(message);
 
			memcpy(x.data(), optimalValues, n * sizeof(double)); // Reset to previous optimum
			NPSOL_currentInterval = i;

            currentCI->lbound += optimum;          // Convert from offsets to targets
            currentCI->ubound += optimum;          // Convert from offsets to targets
 
            /* Set up for the lower bound */
            inform = -1;
            // Number of times to keep trying.
            int cycles = ciMaxIterations;
            double value = INF;
            double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
            while(inform != 0 && cycles > 0) {
                /* Find lower limit */
                currentCI->calcLower = TRUE;
                F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
                    (void*) F77_SUB(npsolLimitObjectiveFunction), &inform, &iter, istate, c, cJac,
				clambda, &f, gradient.data(), hessian.data(), x.data(), iw, &leniw, w, &lenw);
 
                currentCI->lCode = inform;
                if(f < value) {
                    currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
		    omxSaveCheckpoint(x.data(), f, TRUE);
                }
 
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of lower interval %d failed: Bad inform value of %d",
                            i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(x[j] - optimalValues[j]) > objDiff) {
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            x[j] = optimalValues[j] + objDiff;
                        }
                    }
                }
            }
 
            if(OMX_DEBUG) { mxLog("Found lower bound %d.  Seeking upper.", i); }
            // TODO: Repopulate original optimizer state in between CI calculations

			if (currentCI->matrix->name == NULL) {
				snprintf(message, msgLength, "%s[%d, %d] begin upper interval", 
					anonMatrix, currentCI->row + 1, currentCI->col + 1);
			} else {
				snprintf(message, msgLength, "%s[%d, %d] begin upper interval",
					currentCI->matrix->name, currentCI->row + 1, currentCI->col + 1);
			}
 
			omxWriteCheckpointMessage(message);
 
			Free(message);
 
			memcpy(x.data(), optimalValues, n * sizeof(double));
 
            /* Reset for the upper bound */
            value = INF;
            inform = -1;
            cycles = ciMaxIterations;
 
            while(inform != 0 && cycles >= 0) {
                /* Find upper limit */
                currentCI->calcLower = FALSE;
                F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
                                    (void*) F77_SUB(npsolLimitObjectiveFunction), &inform, &iter, istate, c, cJac,
				clambda, &f, gradient.data(), hessian.data(), x.data(), iw, &leniw, w, &lenw);
 
                currentCI->uCode = inform;
                if(f < value) {
                    currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
		    omxSaveCheckpoint(x.data(), f, TRUE);
                }
 
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of upper interval %d failed: Bad inform value of %d",
                            i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(x[j] - optimalValues[j]) > objDiff){
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            x[j] = optimalValues[j] + objDiff;
                        }
                    }
                }
            }
            if(OMX_DEBUG) {mxLog("Found Upper bound %d.", i);}
        }

	NPSOL_fc = NULL;
	NPSOL_fitMatrix = NULL;
	NPSOL_currentInterval = -1;
}
 
void omxSetNPSOLOpts(SEXP options)
{
		char optionCharArray[250] = "";			// For setting options
		int numOptions = length(options);
		SEXP optionNames;
		PROTECT(optionNames = GET_NAMES(options));
		for(int i = 0; i < numOptions; i++) {
			const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
			const char *nextOptionValue = STRING_VALUE(VECTOR_ELT(options, i));
			if (matchCaseInsensitive(nextOptionName, "CI Max Iterations")) {
			} else if(matchCaseInsensitive(nextOptionName, "Analytic Gradients")) {
			} else if(matchCaseInsensitive(nextOptionName, "Number of Threads")) {
			} else {
				sprintf(optionCharArray, "%s %s", nextOptionName, nextOptionValue);
				F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
				if(OMX_DEBUG) { mxLog("Option %s ", optionCharArray); }
			}
		}
		UNPROTECT(1); // optionNames
}

#endif
