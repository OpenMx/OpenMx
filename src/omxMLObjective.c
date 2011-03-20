/*
 *  Copyright 2007-2009 The OpenMx Project
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
#include "omxAlgebraFunctions.h"

extern void omxInitFIMLObjective(omxObjective* oo, SEXP rObj); // Included in case ML gets a raw data object.

#ifndef _OMX_ML_OBJECTIVE_
#define _OMX_ML_OBJECTIVE_ TRUE

typedef struct omxMLObjective {

	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* localCov;
	omxMatrix* localProd;
	omxMatrix* P;
	omxMatrix* C;
	omxMatrix* I;
	double n;
	double logDetObserved;

	double* work;
	int lwork;

} omxMLObjective;

void omxDestroyMLObjective(omxObjective *oo) {

	if(OMX_DEBUG) {Rprintf("Freeing ML Objective.");}
	omxMLObjective* omlo = ((omxMLObjective*)oo->argStruct);

	if(omlo->localCov != NULL)	omxFreeMatrixData(omlo->localCov);
	if(omlo->localProd != NULL)	omxFreeMatrixData(omlo->localProd);
	if(omlo->P != NULL)			omxFreeMatrixData(omlo->P);
	if(omlo->C != NULL)			omxFreeMatrixData(omlo->C);
	if(omlo->I != NULL)			omxFreeMatrixData(omlo->I);
}

omxRListElement* omxSetFinalReturnsMLObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 2;
	omxRListElement* retVal = (omxRListElement*) R_alloc(2, sizeof(omxRListElement));

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);

	retVal[1].numValues = 1;
	retVal[1].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[1].label, "SaturatedLikelihood", 20);
	retVal[1].values[0] = (((omxMLObjective*)oo->argStruct)->logDetObserved + ((omxMLObjective*)oo->argStruct)->observedCov->cols) * (((omxMLObjective*)oo->argStruct)->n - 1);

	return retVal;
}

void omxCallMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning ML Evaluation.\n");}
	// Requires: Data, means, covariances.

	double sum = 0.0, det = 1.0;
	char u = 'U';
	char r = 'R';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	double minusoned = -1.0;
	int onei = 1;
	double fmean = 0.0;

	omxMatrix *scov, *smeans, *cov, *means, *localCov, *localProd, *P, *C, *I;


    /* Locals for readability.  Compiler should cut through this. */
	scov 		= ((omxMLObjective*)oo->argStruct)->observedCov;
	smeans		= ((omxMLObjective*)oo->argStruct)->observedMeans;
	cov			= ((omxMLObjective*)oo->argStruct)->expectedCov;
	means 		= ((omxMLObjective*)oo->argStruct)->expectedMeans;
	localCov 	= ((omxMLObjective*)oo->argStruct)->localCov;
	localProd 	= ((omxMLObjective*)oo->argStruct)->localProd;
	P		 	= ((omxMLObjective*)oo->argStruct)->P;
	C		 	= ((omxMLObjective*)oo->argStruct)->C;
	I		 	= ((omxMLObjective*)oo->argStruct)->I;
	double n 	= ((omxMLObjective*)oo->argStruct)->n;
	double Q	= ((omxMLObjective*)oo->argStruct)->logDetObserved;

    /* Recompute and recopy */
	omxRecompute(cov);							// We assume data won't need to be recomputed
	omxCopyMatrix(localCov, cov);				// But expected cov is destroyed in inversion

	if(OMX_DEBUG_ALGEBRA) {
		omxPrint(scov, "Observed Covariance is");
		omxPrint(localCov, "Implied Covariance Is");
	}

	/* Calculate |expected| */

//	F77_CALL(dgetrf)(&(localCov->cols), &(localCov->rows), localCov->data, &(localCov->cols), ipiv, &info);
	F77_CALL(dpotrf)(&u, &(localCov->cols), localCov->data, &(localCov->cols), &info);

	if(OMX_DEBUG_ALGEBRA) { Rprintf("Info on LU Decomp: %d\n", info);
	omxPrint(localCov, "After Decomp:");}
	if(info > 0) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Expected covariance matrix is non-positive-definite");
		if(oo->matrix->currentState->computeCount <= 0) {
			strncat(errstr, " at starting values", 20);
		}
		strncat(errstr, ".\n", 3);
		omxRaiseError(oo->matrix->currentState, -1, errstr);						// Raise error
		free(errstr);
		return;																		// Leave output untouched
	}

	for(info = 0; info < localCov->cols; info++) { 	    	// |cov| is the square of the product of the diagonal elements of U from the LU factorization.
		det *= localCov->data[info+localCov->rows*info];
	}
	det *= det;

	if(OMX_DEBUG_ALGEBRA) { Rprintf("Determinant of Expected Cov: %f\n", det); }
	det = log(fabs(det));
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Log of Determinant of Expected Cov: %f\n", det); }

	/* Calculate Expected^(-1) */
//	F77_CALL(dgetri)(&(localCov->rows), localCov->data, &(localCov->cols), ipiv, work, lwork, &info);
	F77_CALL(dpotri)(&u, &(localCov->rows), localCov->data, &(localCov->cols), &info);
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Info on Invert: %d\n", info); }

	if(OMX_DEBUG_ALGEBRA) {omxPrint(cov, "Expected Covariance Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localCov, "Inverted Matrix:");}

	/* Calculate C = Observed * expected^(-1) */

	if(OMX_DEBUG_ALGEBRA) {Rprintf("Call is: DSYMM(%d, %d, %f, %0x, %d, %0x, %d, %f, %0x, %d)",
					(scov->rows), (localCov->cols), oned, scov->data, (localCov->leading),
					localCov->data, (localCov->leading), zerod, localProd->data, (localProd->leading));}


	// Stop gcc from issuing a warning
	int majority = *(scov->majority) == 'n' ? scov->rows : scov->cols;

	/*  TODO:  Make sure leading edges are being appropriately calculated, and sub them back into this */
	F77_CALL(dsymm)(&r, &u, &(localCov->rows), &(scov->cols),
					&oned, localCov->data, &(majority),
 					scov->data, &(majority),
					&zerod, localProd->data, &(localProd->leading));

    /* And get the trace of the result */

	for(info = 0; info < localCov->cols; info++) {
		sum += localProd->data[info*localCov->cols + info];
	}

//	for(info = 0; info < (localCov->cols * localCov->rows); info++) {
//		sum += localCov->data[info] * scov->data[info];
//	}

	if(OMX_DEBUG_ALGEBRA) {omxPrint(scov, "Observed Covariance Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localCov, "Inverse Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localProd, "Product Matrix:");}

	if(means != NULL) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Means Likelihood Calculation"); }
		omxRecompute(means);
		omxCopyMatrix(P, means);
		// P = means - smeans
		if(OMX_DEBUG_ALGEBRA) {omxPrint(means, "means");}
		if(OMX_DEBUG_ALGEBRA) {omxPrint(smeans, "smeans");}
		F77_CALL(daxpy)(&(smeans->cols), &minusoned, smeans->data, &onei, P->data, &onei);
		if(OMX_DEBUG_ALGEBRA) {omxPrint(P, "means - smeans");}
		// C = P * Cov
		F77_CALL(dsymv)(&u, &(localCov->rows), &oned, localCov->data, &(localCov->leading), P->data, &onei, &zerod, C->data, &onei);
		// P = C * P'
		fmean = F77_CALL(ddot)(&(C->cols), P->data, &onei, C->data, &onei);

		if(OMX_DEBUG_ALGEBRA) { Rprintf("Mean contribution to likelihood is %f per row.\n", fmean); }
		if(fmean < 0.0) fmean = 0.0;
	}

	oo->matrix->data[0] = (sum + det) * (n - 1) + fmean * (n);

	if(OMX_DEBUG) { Rprintf("MLObjective value comes to: %f (Chisq: %f).\n", oo->matrix->data[0], (sum + det) - Q - cov->cols); }

}

unsigned short int omxNeedsUpdateMLObjective(omxObjective* oo) {
 	if(omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedCov))
		return 1;
	if(((omxMLObjective*)oo->argStruct)->expectedMeans != NULL)
		return omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedMeans);
	return 0;
}

void omxInitMLObjective(omxObjective* oo, SEXP rObj) {

	if(OMX_DEBUG) { Rprintf("Initializing ML objective function.\n"); }

	SEXP nextMatrix;
	int info=0;
	double det=1.0;
	char u = 'U';
	omxMLObjective *newObj = (omxMLObjective*) R_alloc(1, sizeof(omxMLObjective));
	
	/* Set Objective Calls to ML Objective Calls */
	oo->objectiveFun = omxCallMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateMLObjective;
	oo->destructFun = omxDestroyMLObjective;
	oo->setFinalReturns = omxSetFinalReturnsMLObjective;
	oo->repopulateFun = NULL;
	oo->argStruct = (void*) newObj;


	if(OMX_DEBUG) { Rprintf("Retrieving data.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));
	omxData* dataMat = omxNewDataFromMxDataPtr(nextMatrix, oo->matrix->currentState);
	if(strncmp(omxDataType(dataMat), "cov", 3) != 0 && strncmp(omxDataType(dataMat), "cor", 3) != 0) {
		if(strncmp(omxDataType(dataMat), "raw", 3) == 0) {
			if(OMX_DEBUG) { Rprintf("Raw Data: Converting to FIML.\n"); }
			omxInitFIMLObjective(oo, rObj);
			return;
		}
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "ML Objective unable to handle data type %s.\n", omxDataType(dataMat));
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
		if(OMX_DEBUG) { Rprintf("ML Objective unable to handle data type %s.  Aborting.\n", omxDataType(dataMat)); }
		return;
	}
	if(OMX_DEBUG) { Rprintf("Processing Observed Covariance.\n"); }
	newObj->observedCov = omxDataMatrix(dataMat, NULL);
	if(OMX_DEBUG) { Rprintf("Processing Observed Means.\n"); }
	newObj->observedMeans = omxDataMeans(dataMat, NULL, NULL);
	if(OMX_DEBUG && newObj->observedMeans == NULL) { Rprintf("ML: No Observed Means.\n"); }
	if(OMX_DEBUG) { Rprintf("Processing n.\n"); }
	newObj->n = omxDataNumObs(dataMat);
	UNPROTECT(1); // nextMatrix

	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	if(OMX_DEBUG) { Rprintf("Processing Expected Means.\n"); }
	if(!R_FINITE(INTEGER(nextMatrix)[0])) {
		if(OMX_DEBUG) {
			Rprintf("ML: No Expected Means.\n");
		}
		newObj->expectedMeans = NULL;
	} else {
		newObj->expectedMeans = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	}
	UNPROTECT(1);

	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	if(OMX_DEBUG) { Rprintf("Processing Expected Covariance.\n"); }
	newObj->expectedCov = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);

	/* Temporary storage for calculation */
	int rows = newObj->observedCov->rows;
	int cols = newObj->observedCov->cols;
	newObj->localCov = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->localProd = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->P = omxInitMatrix(NULL, 1, cols, TRUE, oo->matrix->currentState);
	newObj->C = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->I = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);

	for(int i = 0; i < rows; i++) omxSetMatrixElement(newObj->I, i, i, 1.0);

	omxCopyMatrix(newObj->localCov, newObj->observedCov);

	newObj->lwork = newObj->expectedCov->rows;
	newObj->work = (double*)R_alloc(newObj->lwork, sizeof(double));


	F77_CALL(dpotrf)(&u, &(newObj->localCov->cols), newObj->localCov->data, &(newObj->localCov->cols), &info);

	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", info); }
	if(info != 0) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Observed Covariance Matrix is non-positive-definite.\n");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
		return;
	}
	for(info = 0; info < newObj->localCov->cols; info++) {
		det *= omxMatrixElement(newObj->localCov, info, info);
	}
	det *= det;					// Product of squares.

	if(OMX_DEBUG) { Rprintf("Determinant of Observed Cov: %f\n", det); }
	newObj->logDetObserved = log(det);
	if(OMX_DEBUG) { Rprintf("Log Determinant of Observed Cov: %f\n", newObj->logDetObserved); }

	omxCopyMatrix(newObj->localCov, newObj->expectedCov);

}

#endif /* _OMX_ML_OBJECTIVE_ */
