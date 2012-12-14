/*
 *  Copyright 2007-2012 The OpenMx Project
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

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include "npsolWrap.h"
#include "omxOpenmpWrap.h"

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>
#include "omxState.h"
#include "omxGlobalState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxExpectation.h"
#include "omxNPSOLSpecific.h"
#include "omxImportFrontendState.h"
#include "omxExportBackendState.h"
#include "omxHessianCalculation.h"
#include "omxOptimizer.h"

//#include "omxSymbolTable.h"

/* Set up R .Call info */
R_CallMethodDef callMethods[] = {
{"omxBackend", (void*(*)())&omxBackend, 12},
{"omxCallAlgebra", (void*(*)())&omxCallAlgebra, 3},
{"findIdenticalRowsData", (void*(*)())&findIdenticalRowsData, 5},
{NULL, NULL, 0}
};

void R_init_mylib(DllInfo *info) {
/* Register routines, allocate resources. */
R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_mylib(DllInfo *info) {
/* Release resources. */
}

/* Globals for function evaluation */
SEXP RObjFun, RConFun;			// Pointers to the functions NPSOL calls
SEXP env;						// Environment for evaluation and object hunting

/* Main functions */
SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options) {

	if(OMX_DEBUG) { Rprintf("-----------------------------------------------------------------------\n");}
	if(OMX_DEBUG) { Rprintf("Explicit call to algebra %d.\n", INTEGER(algNum));}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;
	char output[250];
	int errOut = 0;

	/* Create new omxState for current state storage and initialize it. */
	
	globalState = (omxState*) R_alloc(1, sizeof(omxState));
	omxInitState(globalState, NULL, 1);
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", globalState);}

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}
	globalState->numMats = length(matList);
	globalState->matrixList = (omxMatrix**) R_alloc(length(matList), sizeof(omxMatrix*));

	for(k = 0; k < length(matList); k++) {
		PROTECT(nextMat = VECTOR_ELT(matList, k));	// This is the matrix + populations
		globalState->matrixList[k] = omxNewMatrixFromRPrimitive(nextMat, globalState, 1, - k - 1);
		if(OMX_DEBUG) {
			Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n",
				globalState->matrixList[k], globalState->matrixList[k]->rows, globalState->matrixList[k]->cols);
		}
		UNPROTECT(1); // nextMat
	}

	algebra = omxNewAlgebraFromOperatorAndArgs(algebraNum, globalState->matrixList, globalState->numMats, globalState);

	if(algebra==NULL) {
		error(globalState->statusMsg);
	}

	if(OMX_DEBUG) {Rprintf("Completed Algebras and Matrices.  Beginning Initial Compute.\n");}
	omxStateNextEvaluation(globalState);

	omxRecompute(algebra);

	PROTECT(ans = allocMatrix(REALSXP, algebra->rows, algebra->cols));
	for(l = 0; l < algebra->rows; l++)
		for(j = 0; j < algebra->cols; j++)
			REAL(ans)[j * algebra->rows + l] =
				omxMatrixElement(algebra, l, j);

	UNPROTECT(1);	/* algebra */

	if(OMX_DEBUG) { Rprintf("All Algebras complete.\n"); }

	if(globalState->statusCode != 0) {
		errOut = globalState->statusCode;
		strncpy(output, globalState->statusMsg, 250);
	}

	omxFreeAllMatrixData(algebra);
	omxFreeState(globalState);

	if(errOut != 0) {
		error(output);
	}

	return ans;
}

SEXP omxBackend(SEXP fitfunction, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList, SEXP expectList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state) {

	double *x, *g, *R;
	double f;
	double *est, *grad, *hess;

	/* Helpful variables */

	int k, l;					// Index Vars
	
	int errOut = 0;                 // Error state: Clear

	SEXP nextLoc;

	int n;
	int calculateStdErrors = FALSE;
	int numHessians = 0;
	int ciMaxIterations = 5;
	int disableOptimizer = 0;
	int numThreads = 1;
	int analyticGradients = 0;

	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these.  For now, all checking is done on the front-end. */
//	if(!isVector(startVals)) error ("startVals must be a vector");
//	if(!isVector(matList)) error ("matList must be a list");
//	if(!isVector(algList)) error ("algList must be a list");

	omx_omp_init();

	/* 	Set NPSOL options */
	omxSetNPSOLOpts(options, &numHessians, &calculateStdErrors, 
		&ciMaxIterations, &disableOptimizer, &numThreads, 
		&analyticGradients, length(startVals));

	/* Create new omxState for current state storage and initialize it. */
	globalState = (omxState*) R_alloc(1, sizeof(omxState));
	omxInitState(globalState, NULL, numThreads);
	globalState->numFreeParams = length(startVals);
	globalState->analyticGradients = analyticGradients;
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", globalState);}

	/* Retrieve Data Objects */
	if(!errOut) errOut = omxProcessMxDataEntities(data);
    
	/* Retrieve All Matrices From the MatList */
	if(!errOut) errOut = omxProcessMxMatrixEntities(matList);
	
	/* Initialize all Expectations Here */
	if(!errOut) errOut = omxProcessMxExpectationEntities(expectList);

	/* Process Algebras Here */
	if(!errOut) errOut = omxProcessMxAlgebraEntities(algList);

	/* Complete Expectations */
	if(!errOut) errOut = omxCompleteMxExpectationEntities();

	/* Initial Matrix and Algebra Calculations */
    if(!errOut) {
		// disable parallelism until omxDuplicateState() can be invoked
    	globalState->numChildren = 0;
		errOut = omxInitialMatrixAlgebraCompute();
		globalState->numChildren = (numThreads > 1) ? numThreads : 0;
	}

	/* Process Fit Function */
	if(!errOut) errOut = omxProcessFitFunction(fitfunction);
	
	// TODO: Make calculateHessians an option instead.

	if(!errOut) {	// In the event of an initialization error, skip all this.

		int numChildren = globalState->numChildren;

		/* Process Matrix and Algebra Population Function */
		/*
		 Each matrix is a list containing a matrix and the other matrices/algebras that are
		 populated into it at each iteration.  The first element is already processed, above.
		 The rest of the list will be processed here.
		*/
		for(int j = 0; j < globalState->numMats; j++) {
			PROTECT(nextLoc = VECTOR_ELT(matList, j));		// This is the matrix + populations
			omxProcessMatrixPopulationList(globalState->matrixList[j], nextLoc);
			UNPROTECT(1);
		}

		/* Process Free Var List */
		omxProcessFreeVarList(varList, globalState->numFreeParams);

		/* Processing Constraints */
		omxProcessConstraints(constraints);

		/* Process Confidence Interval List */
		omxProcessConfidenceIntervals(intervalList);

		/* Process Checkpoint List */
		omxProcessCheckpointOptions(checkpointList);

		for(int i = 0; i < numChildren; i++) {
			omxDuplicateState(globalState->childList[i], globalState);
		}
		
	} else { // End if(errOut)
		error(globalState->statusMsg);
	}

	n = globalState->numFreeParams;

	if(n == 0) {			// Special Case for the evaluation-only condition
		g = NULL;
		x = NULL;
		R = NULL;
	} else {
		g		= (double*) R_alloc (n, sizeof ( double ) );
		x		= (double*) R_alloc ((n+1), sizeof ( double ));
		R		= (double*) R_alloc (n * n, sizeof ( double ));

		for(k = 0; k < n; k++) {
			x[k] = REAL(startVals)[k];
		}

	}
	
	omxInvokeNPSOL(&f, x, g, R, disableOptimizer);

	SEXP minimum, estimate, gradient, hessian, code, status, statusMsg, iterations;
	SEXP evaluations, ans=NULL, names=NULL, algebras, matrices, expectations, optimizer;
	SEXP intervals, NAmat, intervalCodes, calculatedHessian, stdErrors;

	int numReturns = 14;

	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(status = allocVector(VECSXP, 3));
	PROTECT(iterations = NEW_NUMERIC(1));
	PROTECT(evaluations = NEW_NUMERIC(2));
	PROTECT(matrices = NEW_LIST(globalState->numMats));
	PROTECT(algebras = NEW_LIST(globalState->numAlgs));
	PROTECT(expectations = NEW_LIST(globalState->numExpects));

	/* N-dependent SEXPs */
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(optimizer = allocVector(VECSXP, 2));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));
	PROTECT(calculatedHessian = allocMatrix(REALSXP, n, n));
	PROTECT(stdErrors = allocMatrix(REALSXP, n, 1)); // for optimizer
	PROTECT(names = allocVector(STRSXP, 2)); // for optimizer
	PROTECT(intervals = allocMatrix(REALSXP, globalState->numIntervals, 2)); // for optimizer
	PROTECT(intervalCodes = allocMatrix(INTSXP, globalState->numIntervals, 2)); // for optimizer
	PROTECT(NAmat = allocMatrix(REALSXP, 1, 1)); // In case of missingness
	REAL(NAmat)[0] = R_NaReal;

	/* Store outputs for return */
	if(fitfunction != NULL) {
		REAL(minimum)[0] = f;
		globalState->optimum = f;
	} else {
		REAL(minimum)[0] = R_NaReal;
	}

	est = REAL(estimate);  // Aliases to avoid repeated function calls.
	grad = REAL(gradient);
	hess = REAL(hessian);

	for(k = 0; k < n; k++) {		// Do these with memcpy instead, probably.
		est[k] = x[k];
		grad[k] = g[k];
		for(l = 0; l < n; l++) {			// Save the NPSOL hessian, in case somebody wants it
			hess[k*n + l] = R[k*n + l];		// This is ok, because they're both in Col-Major already.
		}
	}

	omxSaveState(globalState, x, f);		// Keep the current values for the globalState.

	/* Fill in details from the optimizer */
	SET_VECTOR_ELT(optimizer, 0, gradient);
	SET_VECTOR_ELT(optimizer, 1, hessian);

	SET_STRING_ELT(names, 0, mkChar("minimum"));
	SET_STRING_ELT(names, 1, mkChar("estimate"));
	namesgets(optimizer, names);

	REAL(code)[0] = globalState->inform;
	REAL(iterations)[0] = globalState->iter;
	REAL(evaluations)[0] = globalState->computeCount;

	/* Fill Status code. */
	SET_VECTOR_ELT(status, 0, code);
	PROTECT(code = NEW_NUMERIC(1));
	REAL(code)[0] = globalState->statusCode;
	SET_VECTOR_ELT(status, 1, code);
	PROTECT(statusMsg = allocVector(STRSXP, 1));
	SET_STRING_ELT(statusMsg, 0, mkChar(globalState->statusMsg));
	SET_VECTOR_ELT(status, 2, statusMsg);

	if(numHessians && globalState->fitMatrix != NULL && globalState->optimumStatus >= 0) {		// No hessians or standard errors if the optimum is invalid
		if(globalState->numConstraints == 0) {
			if(OMX_DEBUG) { Rprintf("Calculating Hessian for Fit Function.\n");}
			int gotHessians = omxEstimateHessian(numHessians, .0001, 4, globalState);
			if(gotHessians) {
				if(calculateStdErrors) {
					for(int j = 0; j < numHessians; j++) {		//TODO: Fix Hessian calculation to allow more if requested
						if(OMX_DEBUG) { Rprintf("Calculating Standard Errors for Fit Function.\n");}
						omxFitFunction* oo = globalState->fitMatrix->fitFunction;
						if(oo->getStandardErrorFun != NULL) {
							oo->getStandardErrorFun(oo);
						} else {
							omxCalculateStdErrorFromHessian(2.0, oo);
						}
					}
				}
			} else {
				numHessians = 0;
			}
		} else {
			numHessians = 0;
		}
	}

	/* Likelihood-based Confidence Interval Calculation */
	if(globalState->numIntervals) {
		omxNPSOLConfidenceIntervals(&f, x, g, R, ciMaxIterations);
	}  

	handleFreeVarList(globalState, globalState->optimalValues, n);  // Restore to optima for final compute
	if(!errOut) omxFinalAlgebraCalculation(globalState, matrices, algebras, expectations); 

	omxPopulateFitFunction(globalState, numReturns, &ans, &names);

	if(numHessians) {
		omxPopulateHessians(numHessians, globalState->fitMatrix, 
			calculatedHessian, stdErrors, calculateStdErrors, n);
	}

	if(globalState->numIntervals) {	// Populate CIs
		omxPopulateConfidenceIntervals(globalState, intervals, intervalCodes);
	}
	
	REAL(evaluations)[1] = globalState->computeCount;

	int nextEl = 0;

	SET_STRING_ELT(names, nextEl++, mkChar("minimum"));
	SET_STRING_ELT(names, nextEl++, mkChar("estimate"));
	SET_STRING_ELT(names, nextEl++, mkChar("gradient"));
	SET_STRING_ELT(names, nextEl++, mkChar("hessianCholesky"));
	SET_STRING_ELT(names, nextEl++, mkChar("status"));
	SET_STRING_ELT(names, nextEl++, mkChar("iterations"));
	SET_STRING_ELT(names, nextEl++, mkChar("evaluations"));
	SET_STRING_ELT(names, nextEl++, mkChar("matrices"));
	SET_STRING_ELT(names, nextEl++, mkChar("algebras"));
	SET_STRING_ELT(names, nextEl++, mkChar("expectations"));
	SET_STRING_ELT(names, nextEl++, mkChar("confidenceIntervals"));
	SET_STRING_ELT(names, nextEl++, mkChar("confidenceIntervalCodes"));
	SET_STRING_ELT(names, nextEl++, mkChar("calculatedHessian"));
	SET_STRING_ELT(names, nextEl++, mkChar("standardErrors"));

	nextEl = 0;

	SET_VECTOR_ELT(ans, nextEl++, minimum);
	SET_VECTOR_ELT(ans, nextEl++, estimate);
	SET_VECTOR_ELT(ans, nextEl++, gradient);
	SET_VECTOR_ELT(ans, nextEl++, hessian);
	SET_VECTOR_ELT(ans, nextEl++, status);
	SET_VECTOR_ELT(ans, nextEl++, iterations);
	SET_VECTOR_ELT(ans, nextEl++, evaluations);
	SET_VECTOR_ELT(ans, nextEl++, matrices);
	SET_VECTOR_ELT(ans, nextEl++, algebras);
	SET_VECTOR_ELT(ans, nextEl++, expectations);
	SET_VECTOR_ELT(ans, nextEl++, intervals);
	SET_VECTOR_ELT(ans, nextEl++, intervalCodes);
	if(numHessians == 0) {
		SET_VECTOR_ELT(ans, nextEl++, NAmat);
	} else {
		SET_VECTOR_ELT(ans, nextEl++, calculatedHessian);
	}
	if(!calculateStdErrors) {
		SET_VECTOR_ELT(ans, nextEl++, NAmat);
	} else {
		SET_VECTOR_ELT(ans, nextEl++, stdErrors);
	}
	namesgets(ans, names);

	if(OMX_VERBOSE) {
		Rprintf("Inform Value: %d\n", globalState->optimumStatus);
		Rprintf("--------------------------\n");
	}

	/* Free data memory */
	omxFreeState(globalState);

	UNPROTECT(numReturns);						// Unprotect Output Parameters
	UNPROTECT(8);								// Unprotect internals

	if(OMX_DEBUG) {Rprintf("All vectors freed.\n");}

	return(ans);

}

