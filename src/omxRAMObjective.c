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
#include "omxBLAS.h"
#include "omxObjectiveMetadata.h"
#include "omxFIMLObjective.h"
#include "omxRAMObjective.h"

// Forward declarations
void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);

void omxInitRAMObjectiveWithSummaryData(omxObjective* oo, SEXP rObj);
void omxInitRAMObjectiveWithRawData(omxObjective* oo, SEXP rObj);


omxRListElement* omxSetFinalReturnsRAMObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 3;
	omxRListElement* retVal = (omxRListElement*) R_alloc(*numReturns, sizeof(omxRListElement));

	char u = 'U';
	char r = 'R';
	double sum = 0;
	double det = 1;
    double oned = 1.0, zerod = 0.0;
	omxMatrix* cov = ((omxRAMObjective*)oo->argStruct)->cov;
	int ncols = ((omxRAMObjective*)oo->argStruct)->cov->cols;
	omxMatrix* diag = omxInitTemporaryMatrix(NULL, ncols, ncols, TRUE, oo->matrix->currentState);
	
	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);

	retVal[1].numValues = 1;
	retVal[1].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[1].label, "SaturatedLikelihood", 20);
	retVal[1].values[0] = ((omxRAMObjective*)oo->argStruct)->logDetObserved;
	
	retVal[2].numValues = 1;
	retVal[2].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[2].label, "IndependenceLikelihood", 23);
	// Independence model assumes all-zero manifest covariances.
	for(int i = 0; i < ncols; i++) {
		double value = omxMatrixElement(cov, i, i);
		omxSetMatrixElement(diag, i, i, 1./value);
		det *= value;
	}
    det = log(det);
    if(OMX_DEBUG) { omxPrint(diag, "Diag:"); }
	// (det(expected) + tr(observed * expected^-1)) * (n - 1);
	F77_CALL(dsymm)(&r, &u, &(diag->rows), &(diag->cols),
					&oned, diag->data, &(diag->leading),
 					cov->data, &(cov->leading),
					&zerod, diag->data, &(diag->leading));
	for(int i = 0; i < ncols; i++) {
		sum += omxMatrixElement(diag, i, i);
	}
	if(OMX_DEBUG) { omxPrint(cov, "Observed:"); }
	retVal[2].values[0] = (sum + det) * (((omxRAMObjective*)oo->argStruct)->n - 1);
    omxFreeMatrixData(diag);
	
	return retVal;
}

void omxPopulateRAMAttributesRawData(omxObjective *oo, SEXP algebra) {
	omxFIMLObjective *argStruct = ((omxFIMLObjective*)oo->argStruct);
	SEXP expCovExt, expMeanExt, rowLikelihoodsExt;
	omxMatrix *expCovInt, *expMeanInt, *rowLikelihoodsInt;
	expCovInt = argStruct->cov;
	expMeanInt = argStruct->means;
	rowLikelihoodsInt = argStruct->rowLikelihoods;

	
	PROTECT(expCovExt = allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL) {
		PROTECT(expMeanExt = allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		PROTECT(expMeanExt = allocMatrix(REALSXP, 0, 0));		
	}
	PROTECT(rowLikelihoodsExt = allocVector(REALSXP, rowLikelihoodsInt->rows));
	for(int row = 0; row < rowLikelihoodsInt->rows; row++)
		REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);

	setAttrib(algebra, install("expCov"), expCovExt);
	setAttrib(algebra, install("expMean"), expMeanExt);
	setAttrib(algebra, install("likelihoods"), rowLikelihoodsExt);

	UNPROTECT(3);
}

void omxPopulateRAMAttributesSummaryData(omxObjective *oo, SEXP algebra) {
	omxRAMObjective *argStruct = ((omxRAMObjective*)oo->argStruct);
    omxMatrix *I = argStruct->I;
    omxMatrix *A = argStruct->A;
    omxMatrix *S = argStruct->S;
    omxMatrix *F = argStruct->F;
    omxMatrix *Mns = argStruct->Mns;
    int numIters = argStruct->numIters;
    omxMatrix *X = argStruct->X;
    omxMatrix *Y = argStruct->Y;
    omxMatrix *Z = argStruct->Z;
    omxMatrix *Ax= argStruct->Ax;
	omxMatrix *expCovInt = argStruct->C;			// Expected covariance
	omxMatrix *expMeanInt = argStruct->M;			// Expected means

    omxRecompute(A);
    omxRecompute(S);
    omxRecompute(F);
    if(expMeanInt != NULL) omxRecompute(expMeanInt);

    omxCalculateRAMCovarianceAndMeans(A, S, F, expMeanInt, expCovInt, Mns, numIters, I, Z, Y, X, Ax);

	SEXP expCovExt, expMeanExt;
	PROTECT(expCovExt = allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL) {
		PROTECT(expMeanExt = allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		PROTECT(expMeanExt = allocMatrix(REALSXP, 0, 0));		
	}   
	setAttrib(algebra, install("expCov"), expCovExt);
	setAttrib(algebra, install("expMean"), expMeanExt);
    setAttrib(algebra, install("likelihoods"), PROTECT(allocVector(REALSXP, 0)));
	UNPROTECT(3);
}

void omxDestroyRAMObjective(omxObjective *oo) {
	omxRAMObjective *argStruct = ((omxRAMObjective*)oo->argStruct);

	oo->matrix->currentState->saturatedModel = argStruct->logDetObserved;

	if(argStruct->I   != NULL) omxFreeMatrixData(argStruct->I);
	if(argStruct->C   != NULL) omxFreeMatrixData(argStruct->C);
	if(argStruct->X   != NULL) omxFreeMatrixData(argStruct->X);
	if(argStruct->Y   != NULL) omxFreeMatrixData(argStruct->Y);
	if(argStruct->Z   != NULL) omxFreeMatrixData(argStruct->Z);
	if(argStruct->Ax  != NULL) omxFreeMatrixData(argStruct->Ax);
	if(argStruct->P   != NULL) omxFreeMatrixData(argStruct->P);
	if(argStruct->V   != NULL) omxFreeMatrixData(argStruct->V);
	if(argStruct->Mns != NULL) omxFreeMatrixData(argStruct->Mns);
	// omxFreeMatrixData(argStruct->mCov);
}

/*
 * omxCalculateRAMCovarianceAndMeans
 * 			Just like it says on the tin.  Calculates the mean and covariance matrices
 * for a RAM model.  M is the number of total variables, latent and manifest. N is
 * the number of manifest variables.
 *
 * params:
 * omxMatrix *A, *S, *F 	: matrices as specified in the RAM model.  MxM, MxM, and NxM
 * omxMatrix *M				: vector containing model implied means. 1xM
 * omxMatrix *Cov			: On output: model-implied manifest covariance.  NxN.
 * omxMatrix *Means			: On output: model-implied manifest means.  1xN.
 * int numIterations		: Precomputed number of iterations of taylor series expansion.
 * omxMatrix *I				: Identity matrix.  If left NULL, will be populated.  MxM.
 * omxMatrix *Z				: On output: Computed (I-A)^-1. MxM.
 * omxMatrix *Y, *X, *Ax	: Space for computation. NxM, NxM, MxM.  On exit, populated.
 */

void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax) {
	
	if(OMX_DEBUG) { Rprintf("Running RAM computation."); }
	
	double oned = 1.0, zerod=0.0, minusOned = -1.0;
	int onei = 1;
	
	if(Ax == NULL || I == NULL || Z == NULL || Y == NULL || X == NULL) {
		error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
	}
	
	if(Cov == NULL || Means == NULL) {
		return; // We're not populating results, so why are we bothering, here?
	}
	
	// if(   (Cov->rows != Cov->cols)  || (A->rows  != A->cols)  // Conformance check
	// 	|| (X->rows  != Cov->cols)  || (X->cols  != A->rows)
	// 	|| (Y->rows  != Cov->cols)  || (Y->cols  != A->rows)
	// 	|| (Ax->rows != Cov->cols)  || (Ax->cols != A->rows)
	// 	|| (I->rows  != Cov->cols)  || (I->cols  != Cov->rows)
	// 	|| (Y->rows  != Cov->cols)  || (Y->cols  != A->rows)
	// 	|| (M->cols  != Cov->cols)  || (M->rows  != 1)
	// 	|| (Means->rows != 1)       || (Means->cols != Cov->cols) ) {
	// 		error("INTERNAL ERROR: Incorrectly sized matrices being passed to omxRAMObjective Calculation.\n Please report this to the OpenMx development team.");
	// }
	
	if(numIters == NA_INTEGER) {
		int ipiv[A->rows], lwork = 4 * A->rows * A->cols, k;		// TODO: Speedups can be made by preallocating this.
		double work[lwork];
		if(OMX_DEBUG_ALGEBRA) { Rprintf("RAM Algebra (I-A) inversion using standard (general) inversion.\n"); }
		

		/* Z = (I-A)^-1 */
		if(I->colMajor != A->colMajor) {
			omxTransposeMatrix(I);
		}
		omxCopyMatrix(Z, A);

		/* Z = (I-A)^-1 */
		F77_CALL(omxunsafedgemm)(I->majority, Z->majority, &(I->cols), &(I->rows), &(Z->rows), &oned, I->data, &(I->cols), I->data, &(I->cols), &minusOned, Z->data, &(Z->cols));

		F77_CALL(dgetrf)(&(Z->rows), &(Z->cols), Z->data, &(Z->leading), ipiv, &k);
		if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", k); }
		if(k > 0) {
		        char errStr[250];
		        strncpy(errStr, "(I-A) is exactly singular.", 100);
		        omxRaiseError(A->currentState, -1, errStr);                    // Raise Error
		        return;
		}
		F77_CALL(dgetri)(&(Z->rows), Z->data, &(Z->leading), ipiv, work, &lwork, &k);
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Info on Invert: %d\n", k); }

		if(OMX_DEBUG_ALGEBRA) {omxPrint(Z, "Z");}
		
	} else {
	
		if(OMX_DEBUG_ALGEBRA) { Rprintf("RAM Algebra (I-A) inversion using optimized expansion.\n"); }
		
		/* Taylor Expansion optimized I-A calculation */
		if(I->colMajor != A->colMajor) {
			omxTransposeMatrix(I);
		}

		if(I->colMajor != Ax->colMajor) {
			omxTransposeMatrix(Ax);
		}
	
		omxCopyMatrix(Z, A);
	
		/* Optimized I-A inversion: Z = (I-A)^-1 */
		F77_CALL(omxunsafedgemm)(I->majority, A->majority, &(I->cols), &(I->rows), &(A->rows), &oned, I->data, &(I->cols), I->data, &(I->cols), &oned, Z->data, &(Z->cols));  // Z = I + A = A^0 + A^1

		for(int i = 1; i <= numIters; i++) { // TODO: Efficiently determine how many times to do this
			// The sequence goes like this: (I + A), I + (I + A) * A, I + (I + (I + A) * A) * A, ...
			// Which means only one DGEMM per iteration.
			if(OMX_DEBUG_ALGEBRA) { Rprintf("....RAM: Iteration #%d/%d\n", i, numIters);}
			omxCopyMatrix(Ax, I);
			F77_CALL(omxunsafedgemm)(A->majority, A->majority, &(Z->cols), &(Z->rows), &(A->rows), &oned, Z->data, &(Z->cols), A->data, &(A->cols), &oned, Ax->data, &(Ax->cols));  // Ax = Z %*% A + I
			omxMatrix* m = Z; Z = Ax; Ax = m;	// Juggle to make Z equal to Ax
		}
	}
		
	/* Cov = FZSZ'F' */
	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(F->majority), *(Z->majority), (F->rows), (Z->cols), (Z->rows), onei, F->data, (F->leading), Z->data, (Z->leading), zerod, Y->data, (Y->leading));}
	F77_CALL(omxunsafedgemm)(F->majority, Z->majority, &(F->rows), &(Z->cols), &(Z->rows), &oned, F->data, &(F->leading), Z->data, &(Z->leading), &zerod, Y->data, &(Y->leading)); 	// Y = FZ

	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(Y->majority), *(S->majority), (Y->rows), (S->cols), (S->rows), onei, Y->data, (Y->leading), S->data, (S->leading), zerod, X->data, (X->leading));}
	F77_CALL(omxunsafedgemm)(Y->majority, S->majority, &(Y->rows), &(S->cols), &(S->rows), &oned, Y->data, &(Y->leading), S->data, &(S->leading), &zerod, X->data, &(X->leading)); 	// X = FZS

	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(X->majority), *(Y->minority), (X->rows), (Y->rows), (Y->cols), onei, X->data, (X->leading), Y->data, (Y->lagging), zerod, Cov->data, (Cov->leading));}
	F77_CALL(omxunsafedgemm)(X->majority, Y->minority, &(X->rows), &(Y->rows), &(Y->cols), &oned, X->data, &(X->leading), Y->data, &(Y->leading), &zerod, Cov->data, &(Cov->leading));
	 // Cov = FZSZ'F' (Because (FZ)' = Z'F')
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Cov, "....RAM: Model-implied Covariance Matrix:");}
	
	if(M != NULL && Means != NULL) {
		F77_CALL(omxunsafedgemv)(Y->majority, &(Y->rows), &(Y->cols), &oned, Y->data, &(Y->leading), M->data, &onei, &zerod, Means->data, &onei);
	}
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....RAM: Model-implied Means Vector:");}
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
	omxMatrix *Mns = ((omxRAMObjective*)oo->argStruct)->Mns;			// Calculation space
	omxMatrix *means = ((omxRAMObjective*)oo->argStruct)->means;		// Observed means
	int numIters = ((omxRAMObjective*)oo->argStruct)->numIters;			// Precalculated # of iterations
	double Q = ((omxRAMObjective*)oo->argStruct)->logDetObserved;
	int* lwork = &(((omxRAMObjective*)oo->argStruct)->lwork);
	double n = (((omxRAMObjective*)oo->argStruct)->n);

	/* Since we're not using AlgebraFunctions, we need to recompute these by hand. */
	omxRecompute(A);
	omxRecompute(S);
	omxRecompute(F);
	if(M != NULL) 
		omxRecompute(M);

	if(OMX_DEBUG_ALGEBRA) {
		omxPrint(A, "A");
		omxPrint(S, "S");
		omxPrint(F, "F");
	}

	double fmean = 0.0;
	double MinusOne = -1.0;
	double Zero = 0.0;
	double One = 1.0;
	int OneI = 1;
	int j, k;
	char u = 'U';	
	
	/* C = F((I-A)^1)S((I-A)^1)'F'   Also, if there's an M, means = M((I-A)^-1)'F' */
	omxCalculateRAMCovarianceAndMeans(A, S, F, M, C, Mns, numIters, I, Z, Y, X, Ax);

	/* Val = sum(diag(tempCov %*% solve(PredictedCov))) + log(det(PredictedCov)) */
	/* Alternately, Val = sum (tempCov .* PredictedCov^-1) + log(det(PredictedCov)) */

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
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "Model-implied Means Vector:");}
		omxCopyMatrix(P, means);

		// P = Expected - Observed
		F77_CALL(daxpy)(&(P->cols), &MinusOne, Mns->data, &OneI, P->data, &OneI);

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
	
	omxState* currentState = oo->matrix->currentState;	

	omxData* dataElt = omxNewDataFromDataSlot(rObj, currentState, "data");	
	if(strncmp(dataElt->type, "raw", 3) == 0) {
		omxInitRAMObjectiveWithRawData(oo, rObj);
	} else {
		omxInitRAMObjectiveWithSummaryData(oo, rObj);
	}	
}

void omxInitRAMObjectiveWithSummaryData(omxObjective* oo, SEXP rObj) {
	
	int l, k;
	
	SEXP slotValue;

	omxRAMObjective *newObj = (omxRAMObjective*) R_alloc(1, sizeof(omxRAMObjective));
	
	omxState* currentState = oo->matrix->currentState;

	if(OMX_DEBUG) { Rprintf("Initializing RAM objective function with covariance/correlation data.\n"); }

	/* Set Objective Functions to RAM Objective Functions*/
	oo->objectiveFun = omxCallRAMObjective;
	oo->destructFun = omxDestroyRAMObjective;
	oo->setFinalReturns = omxSetFinalReturnsRAMObjective;
	oo->needsUpdateFun = omxNeedsUpdateRAMObjective;
	oo->populateAttrFun = omxPopulateRAMAttributesSummaryData;
	oo->repopulateFun = NULL;

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
	
	if(OMX_DEBUG) { Rprintf("Processing expansion iteration depth.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("depth")));
	newObj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { Rprintf("Using %d iterations.", newObj->numIters); }
	UNPROTECT(1);

	l = newObj->cov->rows;
	k = newObj->A->cols;


	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }
	newObj->Z = 	omxInitMatrix(NULL, k, k, TRUE, oo->matrix->currentState);
	newObj->Ax = 	omxInitMatrix(NULL, k, k, TRUE, oo->matrix->currentState);
	newObj->Y = 	omxInitMatrix(NULL, l, k, TRUE, oo->matrix->currentState);
	newObj->X = 	omxInitMatrix(NULL, l, k, TRUE, oo->matrix->currentState);
	newObj->C = 	omxInitMatrix(NULL, l, l, TRUE, oo->matrix->currentState);
	newObj->Mns = 	omxInitMatrix(NULL, 1, l, TRUE, oo->matrix->currentState);
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

	oo->argStruct = (void*) newObj;	
}

void omxInitRAMObjectiveWithRawData(omxObjective* oo, SEXP rObj) {

	if(OMX_DEBUG) { Rprintf("Initializing RAM objective function with raw data.\n"); }

	omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(1, sizeof(omxFIMLObjective));
		
	/* Set default Objective calls to FIML Objective Calls */
	oo->objectiveFun = omxCallFIMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateFIMLObjective;
	oo->setFinalReturns = omxSetFinalReturnsFIMLObjective;
	oo->destructFun = omxDestroyFIMLObjective;
	oo->populateAttrFun = omxPopulateRAMAttributesRawData;	
	oo->repopulateFun = NULL;
	
	omxObjectiveMetadataContainer oomc = {NULL, NULL, NULL, NULL, NULL};
	omxObjectiveMetadataContainer *poomc = &oomc;
	
	omxInitRAMMetaData(rObj, poomc, oo->matrix->currentState);
	
	newObj->cov = oomc.cov;
	newObj->means = oomc.means;
	newObj->subObjective = oomc.subObjective;
	newObj->covarianceMeansFunction = oomc.covarianceMeansFunction;
	newObj->destroySubObjective = oomc.destroySubObjective;

	initFIMLObjectiveHelper(oo, rObj, newObj);
	
	oo->argStruct = (void*) newObj;	
}
