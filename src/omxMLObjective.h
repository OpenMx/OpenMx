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

#ifndef _OMX_ML_OBJECTIVE_
#define _OMX_ML_OBJECTIVE_ TRUE

extern omxMatrix** matrixList;

typedef struct omxMLObjective {

	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* localCov;
	omxMatrix* localProd;
	double logDetObserved;
	
	double* work;
	int lwork;

} omxMLObjective;

void omxDestroyMLObjective(omxObjective *oo) {
	
	omxMLObjective* omlo = (omxMLObjective*)oo;
/*
	omxFreeMatrixData(omlo->observedCov);
	omxFreeMatrixData(omlo->observedMeans);
	omxFreeMatrixData(omlo->expectedCov);
	omxFreeMatrixData(omlo->expectedMeans);
	omxFreeMatrixData(omlo->localCov);
	*/
}

void omxCallMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning ML Evaluation.\n");}
	// Requires: Data, means, covariances.

	SEXP matrixDims;
	int *dimList;
	double sum = 0.0, det = 1.0;
	char u = 'U';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	int onei = 1;
	int mainDist = 0;
//	double Q = 0.0;
	double logDet = 0;
	int nextRow, nextCol, numCols, numRemoves;

	omxMatrix *scov, *smeans, *cov, *means, *localCov, *localProd;


    /* Locals for readability.  Compiler should cut through this. */
	scov 		= ((omxMLObjective*)oo->argStruct)->observedCov;
	smeans		= ((omxMLObjective*)oo->argStruct)->observedMeans;
	cov			= ((omxMLObjective*)oo->argStruct)->expectedCov;
	means 		= ((omxMLObjective*)oo->argStruct)->expectedMeans;
	localCov 	= ((omxMLObjective*)oo->argStruct)->localCov;
	localProd 	= ((omxMLObjective*)oo->argStruct)->localProd;
	double Q	= ((omxMLObjective*)oo->argStruct)->logDetObserved;
	double* work = ((omxMLObjective*)oo->argStruct)->work;
	int* lwork = &(((omxMLObjective*)oo->argStruct)->lwork);
	int ipiv[scov->rows];

    /* Recompute and recopy */
	omxRecomputeMatrix(cov);							// We assume data won't need to be recomputed
	omxCopyMatrix(localCov, cov);						// But expected cov is destroyed in inversion
	
	if(OMX_DEBUG) {
		omxPrintMatrix(scov, "Observed Covariance is");
		omxPrintMatrix(localCov, "Implied Covariance Is");
	}
	
	/* Calculate |expected| */
	
	F77_CALL(dgetrf)(&(localCov->cols), &(localCov->rows), localCov->data, &(localCov->cols), ipiv, &info);

	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", info); }
	if(info > 0) {
	//	error("Expected Covariance Matrix is exactly singular.  Maybe a variance estimate has dropped to zero?\n");
	}

	for(info = 0; info < localCov->cols; info++) { 	    	// |cov| is the product of the diagonal elements of U from the LU factorization.
		det *= localCov->data[info+localCov->rows*info];	// Normally, we'd need to worry about transformations made during LU, but
	}														// we're safe here because the determinant of a covariance matrix > 0.	
															// TODO: Prove this for negative estimated variances.
	if(det <= 0) error("Non-positive-definite.\n");
	if(OMX_DEBUG) { Rprintf("Determinant of Expected Cov: %f\n", det); }
	det = log(fabs(det));
	if(OMX_DEBUG) { Rprintf("Log of Determinant of Expected Cov: %f\n", det); }
	
	/* Calculate Expected^(-1) */
	F77_CALL(dgetri)(&(localCov->rows), localCov->data, &(localCov->cols), ipiv, work, lwork, &info);
	if(OMX_DEBUG) { Rprintf("Info on Invert: %d\n", info); }
	
	if(OMX_DEBUG) {omxPrintMatrix(cov, "Expected Covariance Matrix:");}
	if(OMX_DEBUG) {omxPrintMatrix(localCov, "Inverted Matrix:");}
	
	/* Calculate Observed * expected */
	
	if(OMX_DEBUG) {Rprintf("Call is: DGEMM(%c, %c, %d, %d, %d, %f, %0x, %d, %0x, %d, %f, %0x, %d)", *(scov->majority), *(localCov->majority), (scov->rows), (localCov->cols), (scov->cols), oned, scov->data, (localCov->leading), localCov->data, (localCov->leading), zerod, localCov->data, (localCov->leading));}
	
	F77_CALL(dgemm)((scov->majority), (localCov->majority), &(scov->rows), &(localCov->cols), &(scov->cols), &oned, scov->data, &(localCov->leading), localCov->data, &(localCov->leading), &zerod, localProd->data, &(localCov->leading));

    /* And get the trace of the result */

	for(info = 0; info < localCov->cols; info++) {
		sum += localProd->data[info*localCov->cols + info];
	}
	
//	for(info = 0; info < (localCov->cols * localCov->rows); info++) {
//		sum += localCov->data[info] * scov->data[info];
//	}
	
	if(OMX_DEBUG) {omxPrintMatrix(scov, "Observed Covariance Matrix:");}
	if(OMX_DEBUG) {omxPrintMatrix(localCov, "Inverse Matrix:");}
	if(OMX_DEBUG) {omxPrintMatrix(localProd, "Product Matrix:");}
	if(OMX_DEBUG) {Rprintf("k is %d and Q is %f.\n", scov->rows, Q);}
	
    oo->myMatrix->data[0] = fabs(sum + det - scov->rows - Q);

	if(OMX_DEBUG) { Rprintf("MLObjective value comes to: %f (was: %f).\n", oo->myMatrix->data[0], (sum + det)); }

}

unsigned short int omxNeedsUpdateMLObjective(omxObjective* oo) {
 	if(omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedCov))
		return 1;
	if(((omxMLObjective*)oo->argStruct)->expectedMeans != NULL) 
		return omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedMeans);
	return 0;
}

void omxInitMLObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	SEXP nextMatrix, itemList, nextItem, dataElt;
	int nextDef, index, data, column;
	int *items, info=0;
	double det=1.0;
	omxMLObjective *newObj = (omxMLObjective*) R_alloc(1, sizeof(omxMLObjective));
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	if(!R_FINITE(REAL(nextMatrix)[0])) {
		if(OMX_DEBUG) {
			Rprintf("ML: No Expected Means.\n");
		}
		newObj->expectedMeans = NULL;
	} else {
		newObj->expectedMeans = omxNewMatrixFromMxMatrixPtr(nextMatrix);
	}
	UNPROTECT(1);
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	newObj->expectedCov = omxNewMatrixFromMxMatrixPtr(nextMatrix);
	UNPROTECT(1);
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));   // TODO: Need better way to process data elements.
	index = (int) REAL(nextMatrix)[0];
	PROTECT(nextMatrix = VECTOR_ELT(dataList, index));
	PROTECT(dataElt = GET_SLOT(nextMatrix, install("matrix")));
	newObj->observedCov = omxNewMatrixFromMxMatrix(dataElt);
	PROTECT(dataElt = GET_SLOT(nextMatrix, install("vector")));
	if(!R_FINITE(REAL(dataElt)[0])) {
		if(OMX_DEBUG) {
			Rprintf("ML: No Observed Means.\n");
		}
		newObj->observedMeans = NULL;
	} else {
		newObj->observedMeans = omxNewMatrixFromMxMatrixPtr(nextMatrix);
	}

	UNPROTECT(4);
	
	/* Temporary storage for calculation */
	newObj->localCov = omxInitMatrix(NULL, newObj->observedCov->rows, newObj->observedCov->cols, TRUE);
	newObj->localProd = omxInitMatrix(NULL, newObj->observedCov->rows, newObj->observedCov->cols, TRUE);
	
	int ipiv[newObj->observedCov->rows];
	
	omxCopyMatrix(newObj->localCov, newObj->observedCov);
	
	newObj->lwork = newObj->expectedCov->rows;
	newObj->work = (double*)R_alloc(newObj->lwork, sizeof(double));
	
	F77_CALL(dgetrf)(&(newObj->localCov->cols), &(newObj->localCov->rows), newObj->localCov->data, &(newObj->localCov->cols), ipiv, &info);

	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", info); }
	if(info > 0) {
		error("Observed Covariance Matrix is non-positive-definite. Collinearity may be an issue.\n");
	}

	for(info = 0; info < newObj->localCov->cols; info++) { 
		det *= newObj->localCov->data[info+newObj->localCov->rows*info];
	}

	if(OMX_DEBUG) { Rprintf("Determinant of Observed Cov: %f\n", det); }
	newObj->logDetObserved = log(fabs(det));
	if(OMX_DEBUG) { Rprintf("Log Determinant of Observed Cov: %f\n", newObj->logDetObserved); }
	
	omxCopyMatrix(newObj->localCov, newObj->expectedCov);
	
	oo->needsUpdateFun = omxNeedsUpdateMLObjective;
	oo->argStruct = (void*) newObj;
	
}

#endif /* _OMX_ML_OBJECTIVE_ */
