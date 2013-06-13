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

#include <R.h>
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

omp_lock_t GlobalRLock;

static R_CallMethodDef callMethods[] = {
	{"omxBackend", (DL_FUNC) omxBackend, 12},
	{"omxCallAlgebra", (DL_FUNC) omxCallAlgebra, 3},
	{"findIdenticalRowsData", (DL_FUNC) findIdenticalRowsData, 5},
	{NULL, NULL, 0}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_OpenMx(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	omx_omp_init_lock(&GlobalRLock);

	// There is no code that will change behavior whether openmp
	// is set for nested or not. I'm just keeping this in case it
	// makes a difference with older versions of openmp. 2012-12-24 JNP
#if defined(_OPENMP) && _OPENMP <= 200505
	omp_set_nested(0);
#endif
}

void R_unload_OpenMx(DllInfo *info) {
	omx_omp_destroy_lock(&GlobalRLock);
}

#ifdef  __cplusplus
}
#endif

void string_to_try_error( const std::string& str )
{
	error("%s", str.c_str());
}

void exception_to_try_error( const std::exception& ex )
{
	string_to_try_error(ex.what());
}

/* Main functions */
SEXP omxCallAlgebra2(SEXP matList, SEXP algNum, SEXP options) {

	omxManageProtectInsanity protectManager;

	if(OMX_DEBUG) { Rprintf("-----------------------------------------------------------------------\n");}
	if(OMX_DEBUG) { Rprintf("Explicit call to algebra %d.\n", INTEGER(algNum));}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;
	char output[250];
	int errOut = 0;

	/* Create new omxState for current state storage and initialize it. */
	
	globalState = new omxState;
	omxInitState(globalState, NULL);
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", globalState);}

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}

	omxMatrix *args[length(matList)];
	for(k = 0; k < length(matList); k++) {
		PROTECT(nextMat = VECTOR_ELT(matList, k));	// This is the matrix + populations
		args[k] = omxNewMatrixFromRPrimitive(nextMat, globalState, 1, - k - 1);
		globalState->matrixList.push_back(args[k]);
		if(OMX_DEBUG) {
			Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n",
				globalState->matrixList[k], globalState->matrixList[k]->rows, globalState->matrixList[k]->cols);
		}
	}

	algebra = omxNewAlgebraFromOperatorAndArgs(algebraNum, args, length(matList), globalState);

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

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options)
{
	try {
		return omxCallAlgebra2(matList, algNum, options);
	} catch( std::exception& __ex__ ) {
		exception_to_try_error( __ex__ );
	} catch(...) {
		string_to_try_error( "c++ exception (unknown reason)" );
	}
}

SEXP omxBackend2(SEXP fitfunction, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList, SEXP expectList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state) {

	/* Helpful variables */

	SEXP nextLoc;

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

	omxManageProtectInsanity protectManager;

	/* 	Set NPSOL options */
	omxSetNPSOLOpts(options, &numHessians, &calculateStdErrors, 
		&ciMaxIterations, &disableOptimizer, &numThreads, 
		&analyticGradients, length(startVals));

	/* Create new omxState for current state storage and initialize it. */
	globalState = new omxState;
	omxInitState(globalState, NULL);
	globalState->numThreads = numThreads;
	globalState->numFreeParams = length(startVals);
	globalState->analyticGradients = analyticGradients;
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", globalState);}

	/* Retrieve Data Objects */
	omxProcessMxDataEntities(data);
	if (globalState->statusMsg[0]) error(globalState->statusMsg);
    
	/* Retrieve All Matrices From the MatList */
	omxProcessMxMatrixEntities(matList);
	if (globalState->statusMsg[0]) error(globalState->statusMsg);

	if (length(startVals) != length(varList)) error("varList and startVals must be the same length");

	/* Process Free Var List */
	omxProcessFreeVarList(varList);
	if (globalState->statusMsg[0]) error(globalState->statusMsg);

	omxProcessMxExpectationEntities(expectList);
	if (globalState->statusMsg[0]) error(globalState->statusMsg);

	omxProcessMxAlgebraEntities(algList);
	if (globalState->statusMsg[0]) error(globalState->statusMsg);

	omxCompleteMxExpectationEntities();
	if (globalState->statusMsg[0]) error(globalState->statusMsg);

	for (size_t ax=0; ax < globalState->algebraList.size(); ax++) {
		omxMatrix *alg = globalState->algebraList[ax];
		if (!alg->fitFunction) continue;
		omxInitializeFitFunction(alg);
		if (globalState->statusMsg[0]) error(globalState->statusMsg);
	}

	// This is the chance to check for matrix
	// conformability, etc.  Any errors encountered should
	// be reported using R's error() function, not
	// omxRaiseErrorf.

	omxInitialMatrixAlgebraCompute();
	omxResetStatus(globalState);

	if(!isNull(fitfunction)) {
		if(OMX_DEBUG) { Rprintf("Processing fit function.\n"); }
		globalState->fitMatrix = omxMatrixLookupFromState1(fitfunction, globalState);
	}
	if (globalState->statusMsg[0]) error(globalState->statusMsg);
	
	// TODO: Make calculateHessians an option instead.

	/* Process Matrix and Algebra Population Function */
	/*
	  Each matrix is a list containing a matrix and the other matrices/algebras that are
	  populated into it at each iteration.  The first element is already processed, above.
	  The rest of the list will be processed here.
	*/
	for(int j = 0; j < length(matList); j++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, j));		// This is the matrix + populations
		omxProcessMatrixPopulationList(globalState->matrixList[j], nextLoc);
	}

	/* Processing Constraints */
	omxProcessConstraints(constraints);

	/* Process Confidence Interval List */
	omxProcessConfidenceIntervals(intervalList);

	/* Process Checkpoint List */
	omxProcessCheckpointOptions(checkpointList);

	// Probably, this is always the same for all children and
	// doesn't need to be copied to child states.
	cacheFreeVarDependencies(globalState);

	omxFitFunctionCreateChildren(globalState, numThreads);

	int n = globalState->numFreeParams;

	SEXP minimum, estimate, gradient, hessian;
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));

	if (n>0) { memcpy(REAL(estimate), REAL(startVals), sizeof(double)*n); }
	
	omxInvokeNPSOL(REAL(minimum), REAL(estimate), REAL(gradient), REAL(hessian), disableOptimizer);

	SEXP code, status, statusMsg, iterations;
	SEXP evaluations, ans=NULL, names=NULL, algebras, matrices, expectations, optimizer;
	SEXP intervals, NAmat, intervalCodes, calculatedHessian, stdErrors;

	int numReturns = 14;

	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(status = allocVector(VECSXP, 3));
	PROTECT(iterations = NEW_NUMERIC(1));
	PROTECT(evaluations = NEW_NUMERIC(2));
	PROTECT(matrices = NEW_LIST(globalState->matrixList.size()));
	PROTECT(algebras = NEW_LIST(globalState->algebraList.size()));
	PROTECT(expectations = NEW_LIST(globalState->numExpects));

	PROTECT(optimizer = allocVector(VECSXP, 2));
	PROTECT(calculatedHessian = allocMatrix(REALSXP, n, n));
	PROTECT(stdErrors = allocMatrix(REALSXP, n, 1)); // for optimizer
	PROTECT(names = allocVector(STRSXP, 2)); // for optimizer
	PROTECT(intervals = allocMatrix(REALSXP, globalState->numIntervals, 2)); // for optimizer
	PROTECT(intervalCodes = allocMatrix(INTSXP, globalState->numIntervals, 2)); // for optimizer
	PROTECT(NAmat = allocMatrix(REALSXP, 1, 1)); // In case of missingness
	REAL(NAmat)[0] = R_NaReal;

	omxSaveState(globalState, REAL(estimate), REAL(minimum)[0]);

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
						omxCalculateStdErrorFromHessian(2.0, oo);
					}
				}
			} else {
				numHessians = 0;
			}
		} else {
			numHessians = 0;
		}
	} else {
		numHessians = 0;
	}

	/* Likelihood-based Confidence Interval Calculation */
	if(globalState->numIntervals) {
		omxNPSOLConfidenceIntervals(REAL(minimum), REAL(estimate), REAL(gradient), REAL(hessian), ciMaxIterations);
	}  

	// What if fitfunction has its own repopulateFun? TODO
	handleFreeVarListHelper(globalState, globalState->optimalValues, n);

	omxFinalAlgebraCalculation(globalState, matrices, algebras, expectations); 

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

	if(OMX_DEBUG) {Rprintf("All vectors freed.\n");}

	return(ans);

}

SEXP omxBackend(SEXP fitfunction, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList, SEXP expectList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state)
{
	try {
		return omxBackend2(fitfunction, startVals, constraints,
				   matList, varList, algList, expectList,
				   data, intervalList, checkpointList, options, state);
	} catch( std::exception& __ex__ ) {
		exception_to_try_error( __ex__ );
	} catch(...) {
		string_to_try_error( "c++ exception (unknown reason)" );
	}
}

