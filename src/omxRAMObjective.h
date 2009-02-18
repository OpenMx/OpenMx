#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_RAM_OBJECTIVE_
#define _OMX_RAM_OBJECTIVE_ TRUE

typedef struct {

	omxMatrix *cov, *I;
	omxMatrix *A, *S, *F;
	omxMatrix *C, *Y, *Z;
	
	double* work;
	int lwork;

} omxRAMObjective;

void omxDestroyRAMObjective(omxObjective *oo) {


}

void omxCallRAMObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
	
	omxMatrix *cov = ((omxRAMObjective*)oo->argStruct)->cov;	// The compiler should cut through this.
	omxMatrix *I = ((omxRAMObjective*)oo->argStruct)->I;
	omxMatrix *A = ((omxRAMObjective*)oo->argStruct)->A;
	omxMatrix *S = ((omxRAMObjective*)oo->argStruct)->S;
	omxMatrix *F = ((omxRAMObjective*)oo->argStruct)->F;
	omxMatrix *C = ((omxRAMObjective*)oo->argStruct)->C;
	omxMatrix *Y = ((omxRAMObjective*)oo->argStruct)->Y;
	omxMatrix *Z = ((omxRAMObjective*)oo->argStruct)->Z;
	double* work = ((omxRAMObjective*)oo->argStruct)->work;
	int* lwork = &(((omxRAMObjective*)oo->argStruct)->lwork);

	/* Since we're not using AlgebraFunctions, we need to recompute these by hand. */
	omxRecomputeMatrix(A);
	omxRecomputeMatrix(S);
	omxRecomputeMatrix(F);

	omxCopyMatrix(Z, A);
	omxCopyMatrix(Y, I);
	omxCopyMatrix(C, I);

	const char NoTrans = 'n';
	const char Trans = 'T';
	double MinusOne = -1.0;
	double Zero = 0.0;
	double One = 1.0;
	double Two = 2.0;
	int ipiv[I->rows], k;
	
	/* Z = (I-A)^-1 */
	if(OMX_DEBUG) { Rprintf("Beginning Objective Calculation.\n"); }
	
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(I->cols), &(I->rows), &(Z->rows), &One, I->data, &(I->cols), I->data, &(I->cols), &MinusOne, Z->data, &(Z->cols));
	F77_CALL(dgetrf)(&(Z->cols), &(Z->rows), Z->data, &(Z->cols), ipiv, &k);
	F77_CALL(dgetri)(&(Z->rows), Z->data, &(Z->cols), ipiv, work, lwork, &k);

	/* C = FZSZ'F' */ // There MUST be an easier way to do this.  I'm thinking matrix class->
	F77_CALL(dgemm)(&Trans, &Trans, &(Z->cols), &(Z->rows), &(F->cols), &One, Z->data, &(Z->cols), F->data, &(F->cols), &Zero, Y->data, &(Y->cols)); 		// C = ...Z'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(S->cols), &(S->rows), &(Y->rows), &One, S->data, &(S->cols), Y->data, &(Y->cols), &Zero, C->data, &(C->cols)); // C = ..SZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(Z->cols), &(Z->rows), &(C->rows), &One, Z->data, &(Z->cols), C->data, &(C->cols), &Zero, Y->data, &(Y->cols)); 	// C = .ZSZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(F->cols), &(F->rows), &(Y->rows), &One, F->data, &(F->cols), Y->data, &(Y->cols), &Zero, C->data, &(C->cols));	// C = FZSZ'F'
	
	/* Val = sum(diag(tempCov %*% solve(PredictedCov))) + log(det(PredictedCov)) */
	/* Alternately, Val = sum (tempCov .* PredictedCov^-1) + log(det(PredictedCov)) */	
	
	F77_CALL(dgetrf)(&(C->cols), &(C->rows), C->data, &(C->cols), ipiv, &k);

	if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", k); }
	double det = 1.0;
	double sum = 0;

	for(k = 0; k < C->cols; k++) { 		// |A| is the sum of the diagonal elements of U from the LU factorization.
		det *= C->data[k+C->rows*k];  	// Normally, we'd need to worry about transformations made during LU, but
	}									// we're safe here because the determinant of a covariance matrix > 0.

	if(OMX_DEBUG) { Rprintf("Determinant of F(I-A)^-1*S*(I-A)^1'*F': %f\n", det); }
	det = log(det);
	
	F77_CALL(dgetri)(&(C->rows), C->data, &(C->cols), ipiv, work, lwork, &k);
	if(OMX_DEBUG) { Rprintf("Info on Invert: %d\n", k); }
	
	for(k = 0; k < (C->cols * C->rows); k++) {
		sum += C->data[k] * cov->data[k];
	}

	if(OMX_DEBUG) {omxPrintMatrix(C, "Inverted Matrix:");}
	if(OMX_DEBUG) {omxPrintMatrix(cov, "Covariance Matrix:");}

	oo->myMatrix->data[0] = (sum + det);

	if(OMX_DEBUG) { Rprintf("RAMObjective value comes to: %f (was: %f).\n", oo->myMatrix->data[0], (sum + det)); }

}

void omxInitRAMObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	int l, k;
	
	omxRAMObjective *newObj = (omxRAMObjective*) R_alloc(sizeof(omxRAMObjective), 1);
	
	if(OMX_DEBUG) { Rprintf("Using RAM objective function.\n"); }

	// Read the observed covariance matrix from the data argument.
	
	SEXP newMatrix;
	int index;

	PROTECT(newMatrix = GET_SLOT(rObj, install("data")));
	index = round(REAL(newMatrix)[0]);
	PROTECT(newMatrix = VECTOR_ELT(dataList, index));
	Rprintf("%d.\n", index);
	newObj->cov = omxNewMatrixFromMxMatrix(newMatrix); 	// Covariance matrix is the data arg.
	UNPROTECT(2);

	PROTECT(newMatrix = GET_SLOT(rObj, install("A")));
	newObj->A = omxNewMatrixFromMxMatrixPtr(newMatrix);
	UNPROTECT(1);
	
	PROTECT(newMatrix = GET_SLOT(rObj, install("S")));
	newObj->S = omxNewMatrixFromMxMatrixPtr(newMatrix);
	UNPROTECT(1);
	
	PROTECT(newMatrix = GET_SLOT(rObj, install("F")));
	newObj->F = omxNewMatrixFromMxMatrixPtr(newMatrix);
	UNPROTECT(1);
	
	/* Identity Matrix, Size Of A */
	newObj->I = (omxMatrix*) R_alloc(1, sizeof(omxMatrix));
	omxInitMatrix(newObj->I, newObj->A->rows, newObj->A->cols, FALSE);
	for(k =0; k < newObj->I->rows; k++) {
		for(l = 0; l < newObj->I->cols; l++) {
			if(l == k) {
				omxSetMatrixElement(newObj->I, k, l, 1);
			} else {
				omxSetMatrixElement(newObj->I, k, l, 0);
			}
		}
	}
	
	k = newObj->A->rows;

	newObj->Z = omxInitMatrix(NULL, k, k, TRUE);
	newObj->Y = omxInitMatrix(NULL, k, k, TRUE);
	newObj->C = omxInitMatrix(NULL, k, k, TRUE);
	newObj->lwork = k;
	newObj->work = (double*)R_alloc(sizeof(double), newObj->lwork);
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_RAM_OBJECTIVE_ */