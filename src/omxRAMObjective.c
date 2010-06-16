/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxObjective.h"

#ifndef _OMX_RAM_OBJECTIVE_
#define _OMX_RAM_OBJECTIVE_ TRUE

typedef struct {

	omxMatrix *cov, *means, *I;
	omxMatrix *A, *S, *F, *M;
	omxMatrix *C, *X, *Y, *Z, *Ax, *P, *V, *mCov;

	double logDetObserved;
	double n;
	double* work;
	int lwork;

} omxRAMObjective;

omxRListElement* omxSetFinalReturnsRAMObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 2;
	omxRListElement* retVal = (omxRListElement*) R_alloc(2, sizeof(omxRListElement));
	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);

	retVal[1].numValues = 1;
	retVal[1].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[1].label, "SaturatedLikelihood", 20);
	retVal[1].values[0] = ((omxRAMObjective*)oo->argStruct)->logDetObserved;
	return retVal;
}

void omxDestroyRAMObjective(omxObjective *oo) {
	omxRAMObjective *argStruct = ((omxRAMObjective*)oo->argStruct);

	oo->matrix->currentState->saturatedModel = argStruct->logDetObserved;

	omxFreeMatrixData(argStruct->I);
	omxFreeMatrixData(argStruct->C);
	omxFreeMatrixData(argStruct->X);
	omxFreeMatrixData(argStruct->Y);
	omxFreeMatrixData(argStruct->Z);
	omxFreeMatrixData(argStruct->Ax);
	omxFreeMatrixData(argStruct->P);
	omxFreeMatrixData(argStruct->V);
	omxFreeMatrixData(argStruct->mCov);
}

void omxCallRAMObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG_ALGEBRA) { Rprintf("RAM Objective Running.\n"); }
	omxMatrix *cov = ((omxRAMObjective*)oo->argStruct)->cov;	// The compiler should cut through this.
	omxMatrix *I = ((omxRAMObjective*)oo->argStruct)->I;
	omxMatrix *A = ((omxRAMObjective*)oo->argStruct)->A;
	omxMatrix *S = ((omxRAMObjective*)oo->argStruct)->S;
	omxMatrix *F = ((omxRAMObjective*)oo->argStruct)->F;
	omxMatrix *C = ((omxRAMObjective*)oo->argStruct)->C;
	omxMatrix *X = ((omxRAMObjective*)oo->argStruct)->X;
	omxMatrix *Y = ((omxRAMObjective*)oo->argStruct)->Y;
	omxMatrix *Z = ((omxRAMObjective*)oo->argStruct)->Z;
	omxMatrix *Ax= ((omxRAMObjective*)oo->argStruct)->Ax;
	omxMatrix *M = ((omxRAMObjective*)oo->argStruct)->M;				// Expected means
	omxMatrix *P = ((omxRAMObjective*)oo->argStruct)->P;				// Calculation space
	omxMatrix *V = ((omxRAMObjective*)oo->argStruct)->V;				// Calculation space
	omxMatrix *mCov = ((omxRAMObjective*)oo->argStruct)->mCov;			// Calculation space
	omxMatrix *means = ((omxRAMObjective*)oo->argStruct)->means;		// Observed means
	double Q = ((omxRAMObjective*)oo->argStruct)->logDetObserved;
	double* work = ((omxRAMObjective*)oo->argStruct)->work;
	int* lwork = &(((omxRAMObjective*)oo->argStruct)->lwork);
	double n = (((omxRAMObjective*)oo->argStruct)->n);

	/* Since we're not using AlgebraFunctions, we need to recompute these by hand. */
	omxRecompute(A);
	omxRecompute(S);
	omxRecompute(F);

	if(OMX_DEBUG_ALGEBRA) {
		omxPrint(A, "A");
		omxPrint(S, "S");
		omxPrint(F, "F");
	}

	const char NoTrans = 'n';
	double MinusOne = -1.0;
	double fmean = 0.0;
	double Zero = 0.0;
	double One = 1.0;
	int OneI = 1;
	int ipiv[I->rows], j, k;
	char u = 'U';

	/* Z = (I-A)^-1 */
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Beginning RAM Objective Calculation.\n"); }

	/* Taylor Expansion optimized I-A calculation */
	if(I->colMajor != A->colMajor) {
		omxTransposeMatrix(I);
	}
	omxCopyMatrix(Z, A);
	/* Optimized I-A inversion */
	F77_CALL(dgemm)(&NoTrans, A->majority, &(I->cols), &(I->rows), &(A->rows), &One, I->data, &(I->cols), I->data, &(I->cols), &One, Z->data, &(Z->cols));  // Z = I + A = A^0 + A^1

	for(int i = 1; i < (A->rows - 1); i++) { // TODO: Efficiently determne how many times to do this
		omxCopyMatrix(Ax, I);
		F77_CALL(dgemm)(A->majority, A->majority, &(Z->cols), &(Z->rows), &(A->rows), &One, Z->data, &(Z->cols), A->data, &(A->cols), &One, Ax->data, &(Ax->cols));  // Ax = Z %*% A + I
		omxMatrix* m = Z; Z = Ax; Ax = m;	// Juggle to make Z equal to Ax
	}

	if(OMX_DEBUG_ALGEBRA) {omxPrint(Z, "Z");}
	
	/* C = FZSZ'F' */
	if(OMX_DEBUG_ALGEBRA) { Rprintf("DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(F->majority), *(Z->majority), (F->rows), (Z->cols), (Z->rows), One, F->data, (F->leading), Z->data, (Z->leading), Zero, Y->data, (Y->leading));}
	F77_CALL(dgemm)(F->majority, Z->majority, &(F->rows), &(Z->cols), &(Z->rows), &One, F->data, &(F->leading), Z->data, &(Z->leading), &Zero, Y->data, &(Y->leading)); 	// Y = FZ

	if(OMX_DEBUG_ALGEBRA) {omxPrint(Y, "Y=FZ");}

	if(OMX_DEBUG_ALGEBRA) { Rprintf("DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(Y->majority), *(S->majority), (Y->rows), (S->cols), (S->rows), One, Y->data, (Y->leading), S->data, (S->leading), Zero, X->data, (X->leading));}
	F77_CALL(dgemm)(Y->majority, S->majority, &(Y->rows), &(S->cols), &(S->rows), &One, Y->data, &(Y->leading), S->data, &(S->leading), &Zero, X->data, &(X->leading)); 	// X = FZS

	if(OMX_DEBUG_ALGEBRA) {omxPrint(X, "X = FZS");}

	if(OMX_DEBUG_ALGEBRA) { Rprintf("DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(X->majority), *(Y->minority), (X->rows), (Y->rows), (Y->cols), One, X->data, (X->leading), Y->data, (Y->lagging), Zero, C->data, (C->leading));}
	F77_CALL(dgemm)(X->majority, Y->minority, &(X->rows), &(Y->rows), &(Y->cols), &One, X->data, &(X->leading), Y->data, &(Y->leading), &Zero, C->data, &(C->leading)); 	// C = FZSZ'F' (Because (FZ)' = Z'F')

	/* Val = sum(diag(tempCov %*% solve(PredictedCov))) + log(det(PredictedCov)) */
	/* Alternately, Val = sum (tempCov .* PredictedCov^-1) + log(det(PredictedCov)) */

	omxCopyMatrix(mCov, C);

	if(OMX_DEBUG_ALGEBRA) {omxPrint(C, "Model-Implied Covariance Matrix");}

	F77_CALL(dpotrf)(&u, &(C->cols), C->data, &(C->cols), &k);

	if(OMX_DEBUG_ALGEBRA) { Rprintf("Info on LU Decomp: %d\n", k); }
	if(k > 0) {
		char errStr[250];
		strncpy(errStr, "Backing out of parameter space region where expected covariance is non-positive-definite.", 100);
		omxRaiseError(oo->matrix->currentState, -1, errStr);			// Raise Error
		return;
	}
	double det = 1.0;
	double sum = 0;

	for(k = 0; k < C->cols; k++) { 		// |A| is the sum of the diagonal elements of U from the LU factorization.
		det *= C->data[k+C->rows*k];	// Normally, we'd need to worry about transformations made during LU, but
	}									// we're safe here because the determinant of a covariance matrix > 0.
	det *= det;
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Determinant of F(I-A)^-1*S*(I-A)^1'*F': %f\n", det); }
	det = log(det);
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Log of Determinant of F(I-A)^-1*S*(I-A)^1'*F': %f (lwork: %d/%d)\n", det, *lwork, C->rows); }

	F77_CALL(dpotri)(&u, &(C->rows), C->data, &(C->leading), &k);
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Info on Invert: %d\n", k); }

	for(k = 0; k < C->cols; k++) {
		for(j = 0; j <=k ; j++) {
			if(j == k) {
				sum += omxMatrixElement(C, j, k) * omxMatrixElement(cov, j, k);
			} else {
				sum += 2 * omxMatrixElement(C, j, k) * omxMatrixElement(cov, j, k);
			}

		}
	}
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Trace of (observed .* expected^-1) = %f\n", sum); }

	if(means != NULL && M != NULL) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Means Likelihood Calculation"); }
		// For now, assume M has only one column.  If that changes, we need to multiply M by U, where U is an nx1 unit vector.
		if(M->rows > 1) { error("NYI: Back-end currently only supports one row of means.");}
		omxRecompute(Y);
		// omxRecompute(C);
		// omxRecompute(P);
		omxRecompute(M);
		omxCopyMatrix(P, means);

		F77_CALL(dgemv)(Y->majority, &(Y->rows), &(Y->cols), &MinusOne, Y->data, &(Y->leading), M->data, &OneI, &One, P->data, &OneI);
		// C = P * Cov^-1
		F77_CALL(dsymv)(&u, &(C->rows), &One, C->data, &(C->leading), P->data, &OneI, &Zero, V->data, &OneI);
		// P = C * P'
		fmean = F77_CALL(ddot)(&(C->cols), P->data, &OneI, V->data, &OneI);

		if(OMX_DEBUG_ALGEBRA) { Rprintf("Mean contribution to likelihood is %f per row, total %f.\n", fmean, fmean/n); }
		if(fmean < 0.0) fmean = 0.0;
	}

	if(OMX_DEBUG_ALGEBRA) { Rprintf("RAMObjective value: %f + %f - %f= %f.\n", (sum + det), fmean, Q, (sum + det) + fmean - Q); }

	oo->matrix->data[0] = ((sum + det) * (n-1) + fmean * n);

	if(OMX_DEBUG) { Rprintf("RAMObjective value comes to: %f (was: %f).\n", oo->matrix->data[0], (sum + det)); }

}

unsigned short int omxNeedsUpdateRAMObjective(omxObjective* oo) {
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Checking if RAM needs update.  RAM uses A:0x%x, S:0x%x, F:0x%x, M:0x%x.\n",((omxRAMObjective*)oo->argStruct)->A, ((omxRAMObjective*)oo->argStruct)->S, ((omxRAMObjective*)oo->argStruct)->F, ((omxRAMObjective*)oo->argStruct)->M); }
	return(omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->A)
	 	|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->S)
	 	|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->F)
		|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->M));

	// Note: cov is data, and should never need updating.
}

void omxInitRAMObjective(omxObjective* oo, SEXP rObj) {

	int l, k;

	omxRAMObjective *newObj = (omxRAMObjective*) R_alloc(1, sizeof(omxRAMObjective));
	
	omxState* currentState = oo->matrix->currentState;

	if(OMX_DEBUG) { Rprintf("Initializing RAM objective function.\n"); }

	// Read the observed covariance matrix from the data argument.
	omxData* dataElt = omxNewDataFromDataSlot(rObj, currentState, "data");
	newObj->cov = omxDataMatrix(dataElt, NULL);
	
	if(newObj->cov->rows != newObj->cov->cols) {
		error("Covariance/Correlation matrix is not square.  Perhaps this should be a FIML evaluation.");
	}
	
	if(OMX_DEBUG) { Rprintf("Processing observed means.\n"); }
	newObj->means = omxDataMeans(dataElt, 0, NULL);
	if(OMX_DEBUG && newObj->means == NULL) { Rprintf("RAM: No Observed Means.\n"); }
	
	if(OMX_DEBUG) { Rprintf("Processing n.\n"); }
	newObj->n = omxDataNumObs(dataElt);

	if(OMX_DEBUG) { Rprintf("Processing M.\n"); }
	newObj->M = omxNewMatrixFromIndexSlot(rObj, currentState, "M");
	if(newObj->M != NULL) omxRecompute(newObj->M);
	else if(OMX_DEBUG) Rprintf("No M found.\n");

	if(OMX_DEBUG) { Rprintf("Processing A.\n"); }
	newObj->A = omxNewMatrixFromIndexSlot(rObj, currentState, "A");
	omxRecompute(newObj->A);

	if(OMX_DEBUG) { Rprintf("Processing S.\n"); }
	newObj->S = omxNewMatrixFromIndexSlot(rObj, currentState, "S");
	omxRecompute(newObj->S);

	if(OMX_DEBUG) { Rprintf("Processing F.\n"); }
	newObj->F = omxNewMatrixFromIndexSlot(rObj, currentState, "F");
	omxRecompute(newObj->F);

	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { Rprintf("Generating I.\n"); }
	newObj->I = omxNewIdentityMatrix(newObj->A->rows, currentState);
	omxRecompute(newObj->I);

	l = newObj->cov->rows;
	k = newObj->A->cols;


	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }
	newObj->Z = 	omxInitMatrix(NULL, k, k, TRUE, oo->matrix->currentState);
	newObj->Ax = 	omxInitMatrix(NULL, k, k, TRUE, oo->matrix->currentState);
	newObj->Y = 	omxInitMatrix(NULL, l, k, TRUE, oo->matrix->currentState);
	newObj->X = 	omxInitMatrix(NULL, l, k, TRUE, oo->matrix->currentState);
	newObj->C = 	omxInitMatrix(NULL, l, l, TRUE, oo->matrix->currentState);
	newObj->mCov = 	omxInitMatrix(NULL, l, l, TRUE, oo->matrix->currentState);
	newObj->P = 	omxInitMatrix(NULL, 1, l, TRUE, oo->matrix->currentState);
	newObj->V = 	omxInitMatrix(NULL, 1, l, TRUE, oo->matrix->currentState);
	newObj->lwork = k;
	newObj->work = 	(double*)R_alloc(newObj->lwork, sizeof(double));

	int info;
	char u = 'U';
	double det = 1.0, sum = 0.0;
	
	if(OMX_DEBUG) { Rprintf("Precomputing necessary structures.\n"); }
	omxCopyMatrix(newObj->C, newObj->cov);

	/* Call BLAS routines for inverting observed covariance matrix */
	F77_CALL(dpotrf)(&u, &(newObj->C->cols), newObj->C->data, &(newObj->C->cols), &info);
	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", info); }
	if(info > 0) {
		error("Observed Covariance Matrix is non-positive-definite. Collinearity may be an issue.\n");
	}
	for(info = 0; info < newObj->C->cols; info++) {
		det *= omxMatrixElement(newObj->C, info, info);			// Determinant
	}

	/* Precalculate offset */
	if(OMX_DEBUG) { Rprintf("Determinant of Observed Cov: %f\n", det); }
	newObj->logDetObserved = (log(det * det) + newObj->cov->rows) * (newObj->n - 1);
	if(OMX_DEBUG) { Rprintf("Log Determinant %f + %f = : %f\n", log(fabs(det)), sum, newObj->logDetObserved); }

	/* Register functions */
	oo->objectiveFun = omxCallRAMObjective;
	oo->destructFun = omxDestroyRAMObjective;
	oo->setFinalReturns = omxSetFinalReturnsRAMObjective;
	oo->needsUpdateFun = omxNeedsUpdateRAMObjective;
	oo->repopulateFun = NULL;

	oo->argStruct = (void*) newObj;

}


#endif /* _OMX_RAM_OBJECTIVE_ */
