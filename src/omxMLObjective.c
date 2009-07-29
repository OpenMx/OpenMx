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

	omxMLObjective* omlo = ((omxMLObjective*)oo->argStruct);

	omxFreeMatrixData(omlo->localCov);
	omxFreeMatrixData(omlo->localProd);
	omxFreeMatrixData(omlo->P);
	omxFreeMatrixData(omlo->C);
	omxFreeMatrixData(omlo->I);
}

omxRListElement* omxSetFinalReturnsMLObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 1;
	omxRListElement* retVal = (omxRListElement*) Calloc(1, omxRListElement);
	retVal->numValues = 1;
	retVal->values = (double*) Calloc(1, double);
	strncpy(retVal->label, "Saturated", 10);
	*(retVal->values) = ((omxMLObjective*)oo->argStruct)->logDetObserved + ((omxMLObjective*)oo->argStruct)->observedCov->cols;
	return retVal;
}

void omxCallMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning ML Evaluation.\n");}
	// Requires: Data, means, covariances.

	SEXP matrixDims;
	int *dimList;
	double sum = 0.0, det = 1.0;
	char u = 'U';
	char r = 'R';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	double minusoned = -1.0;
	int onei = 1;
	int mainDist = 0;
	double fmean = 0.0;
	double logDet = 0;
	int nextRow, nextCol, numCols, numRemoves;

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
	double* work = ((omxMLObjective*)oo->argStruct)->work;
	int* lwork = &(((omxMLObjective*)oo->argStruct)->lwork);
	int ipiv[scov->rows];

    /* Recompute and recopy */
	omxRecompute(cov);							// We assume data won't need to be recomputed
	omxCopyMatrix(localCov, cov);						// But expected cov is destroyed in inversion
	
	if(OMX_DEBUG) {
		omxPrint(scov, "Observed Covariance is");
		omxPrint(localCov, "Implied Covariance Is");
	}
	
	/* Calculate |expected| */
	
//	F77_CALL(dgetrf)(&(localCov->cols), &(localCov->rows), localCov->data, &(localCov->cols), ipiv, &info);
	F77_CALL(dpotrf)(&u, &(localCov->cols), localCov->data, &(localCov->cols), &info);
	
	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", info); 
	omxPrint(localCov, "After Decomp:");}
	if(info > 0) {
		int n = oo->matrix->currentState->numFreeParams;
			/* This section needs to be replaced.  Once we have a back-end-to-front-end error protocol. */
		int errlen = 50+10*n;
		char errstr[50+10*n];
		char shortstr[30];
		sprintf(errstr, "Non-positive-definite at free parameters:");
		for(int i = 0; i < n; i++) {
			sprintf(shortstr, "\t%3.3f", *(oo->matrix->currentState->freeVarList[i].location[0]));
			strncat(errstr, shortstr, 10);
		}
		strncat(errstr, "\n", 1);
		if(errlen > 250) errlen = 250;
		omxRaiseError(oo->matrix->currentState, -1, errstr);						// Raise error
		return;																		// Leave output untouched
	}

	for(info = 0; info < localCov->cols; info++) { 	    	// |cov| is the square of the product of the diagonal elements of U from the LU factorization.
		det *= localCov->data[info+localCov->rows*info];
	}	
	det *= det;
	// TODO: Prove this will still work for negative estimated variances (Heywood cases).
	
	if(OMX_DEBUG) { Rprintf("Determinant of Expected Cov: %f\n", det); }
	det = log(fabs(det));
	if(OMX_DEBUG) { Rprintf("Log of Determinant of Expected Cov: %f\n", det); }
	
	/* Calculate Expected^(-1) */
//	F77_CALL(dgetri)(&(localCov->rows), localCov->data, &(localCov->cols), ipiv, work, lwork, &info);
	F77_CALL(dpotri)(&u, &(localCov->rows), localCov->data, &(localCov->cols), &info);
	if(OMX_DEBUG) { Rprintf("Info on Invert: %d\n", info); }
	
	if(OMX_DEBUG) {omxPrint(cov, "Expected Covariance Matrix:");}
	if(OMX_DEBUG) {omxPrint(localCov, "Inverted Matrix:");}
	
	/* Calculate C = Observed * expected^(-1) */
	
	if(OMX_DEBUG) {Rprintf("Call is: DSYMM(%d, %d, %f, %0x, %d, %0x, %d, %f, %0x, %d)",  
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
	
	if(OMX_DEBUG) {omxPrint(scov, "Observed Covariance Matrix:");}
	if(OMX_DEBUG) {omxPrint(localCov, "Inverse Matrix:");}
	if(OMX_DEBUG) {omxPrint(localProd, "Product Matrix:");}
		
	if(means != NULL) {
		if(OMX_DEBUG) { Rprintf("Means Likelihood Calculation"); }
		omxRecompute(smeans);
		omxCopyMatrix(P, means);
		// P = means - smeans
		if(OMX_DEBUG) {omxPrint(means, "means");}
		if(OMX_DEBUG) {omxPrint(smeans, "smeans");}
		F77_CALL(dgemm)(I->majority, smeans->majority, &(smeans->rows), &(smeans->cols), &(smeans->rows), &minusoned, I->data, &(I->leading), smeans->data, &(smeans->leading), &oned, P->data, &(P->leading));
		// C = P * Cov
		if(OMX_DEBUG) {omxPrint(P, "means - smeans");}
		F77_CALL(dgemm)(P->majority, localCov->majority, &(P->rows), &(P->cols), &(localCov->cols), &oned, P->data, &(P->leading), localCov->data, &(localCov->leading), &zerod, C->data, &(C->leading));
		if(OMX_DEBUG) {omxPrint(C, "P*Cov");}
		// P = C * P'
		F77_CALL(dgemm)(C->majority, P->minority, &onei, &(C->cols), &onei, &oned, C->data, &(C->leading), P->data, &(P->lagging), &zerod, &fmean, &onei);
		if(OMX_DEBUG) { Rprintf("Mean contribution to likelihood is %f per row.\n", fmean); } 
		fmean = fmean * n;
		if(fmean < 0.0) fmean = 0.0;
	}
	
	oo->matrix->data[0] = sum + det + fmean;

	if(OMX_DEBUG) { Rprintf("MLObjective value comes to: %f (Chisq: %f).\n", oo->matrix->data[0], (sum + det) - Q - cov->cols); }

}

unsigned short int omxNeedsUpdateMLObjective(omxObjective* oo) {
 	if(omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedCov))
		return 1;
	if(((omxMLObjective*)oo->argStruct)->expectedMeans != NULL) 
		return omxNeedsUpdate(((omxMLObjective*)oo->argStruct)->expectedMeans);
	return 0;
}

void omxInitMLObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	if(OMX_DEBUG) { Rprintf("Initializing ML objective function.\n"); }
	
	SEXP nextMatrix, itemList, nextItem, dataElt;
	int nextDef, index, data, column;
	int *items, info=0;
	double det=1.0;
	char u = 'U';
	omxMLObjective *newObj = (omxMLObjective*) R_alloc(1, sizeof(omxMLObjective));
	
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
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));   // TODO: Need better way to process data elements.
	if(OMX_DEBUG) { Rprintf("Processing Observed Covariance.\n"); }
	index = (int) REAL(nextMatrix)[0];
	PROTECT(nextMatrix = VECTOR_ELT(dataList, index));
	PROTECT(dataElt = GET_SLOT(nextMatrix, install("data")));
	newObj->observedCov = omxNewMatrixFromMxMatrix(dataElt, oo->matrix->currentState);
	if(OMX_DEBUG) { Rprintf("Processing Observed Means.\n"); }
	PROTECT(dataElt = GET_SLOT(nextMatrix, install("means")));
	if(!R_FINITE(REAL(dataElt)[0])) {
		if(OMX_DEBUG) {
			Rprintf("ML: No Observed Means.\n");
		}
		newObj->observedMeans = NULL;
	} else {
		newObj->observedMeans = omxNewMatrixFromMxMatrix(dataElt, oo->matrix->currentState);
	}
	UNPROTECT(1); // dataElt
	if(OMX_DEBUG && newObj->observedMeans == NULL) { Rprintf("ML: No Observed Means.\n"); }
	if(OMX_DEBUG) { Rprintf("Processing n.\n"); }
	PROTECT(dataElt = GET_SLOT(nextMatrix, install("numObs")));
	newObj->n = REAL(dataElt)[0];
	UNPROTECT(4);
	
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
	if(info > 0) {
		error("Observed Covariance Matrix is non-positive-definite. Collinearity may be an issue.\n");
	}

	for(info = 0; info < newObj->localCov->cols; info++) { 
		det *= omxMatrixElement(newObj->localCov, info, info);
	}
	det *= det;					// Product of squares.

	if(OMX_DEBUG) { Rprintf("Determinant of Observed Cov: %f\n", det); }
	newObj->logDetObserved = log(det);
	if(OMX_DEBUG) { Rprintf("Log Determinant of Observed Cov: %f\n", newObj->logDetObserved); }
	
	omxCopyMatrix(newObj->localCov, newObj->expectedCov);
	
	oo->objectiveFun = omxCallMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateMLObjective;
	oo->destructFun = omxDestroyMLObjective;
	oo->setFinalReturns = omxSetFinalReturnsMLObjective;
	oo->repopulateFun = NULL;
	oo->argStruct = (void*) newObj;
	
}

#endif /* _OMX_ML_OBJECTIVE_ */
