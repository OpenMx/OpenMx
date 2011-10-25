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
#include "omxFIMLObjective.h"
#include "omxRAMObjective.h"

extern void omxCreateMLObjective(omxObjective* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means);

void omxCallRAMObjective(omxObjective* oo) {
    if(OMX_DEBUG) { Rprintf("RAM Subobjective called.\n"); }
	omxRAMObjective* oro = (omxRAMObjective*)(oo->argStruct);
	
	omxRecompute(oro->A);
	omxRecompute(oro->S);
	omxRecompute(oro->F);
	if(oro->M != NULL)
	    omxRecompute(oro->M);
	    
	omxCalculateRAMCovarianceAndMeans(oro->A, oro->S, oro->F, oro->M, oro->cov, oro->means, oro->numIters, oro->I, oro->Z, oro->Y, oro->X, oro->Ax);
}

void omxDestroyRAMObjective(omxObjective* oo) {

	if(OMX_DEBUG) { Rprintf("Destroying RAM Objective.\n"); }
	
	omxRAMObjective* argStruct = (omxRAMObjective*)(oo->argStruct);

	/* We allocated 'em, so we destroy 'em. */
	if(argStruct->cov != NULL)
		omxFreeMatrixData(argStruct->cov);

	if(argStruct->means != NULL)
		omxFreeMatrixData(argStruct->means);
	
	omxFreeMatrixData(argStruct->I);
	omxFreeMatrixData(argStruct->X);
	omxFreeMatrixData(argStruct->Y);
	omxFreeMatrixData(argStruct->Z);
	omxFreeMatrixData(argStruct->Ax);

	if(argStruct->ppmlData != NULL) 
		omxFreeData(argStruct->ppmlData);
	
}

void omxPopulateRAMAttributes(omxObjective *oo, SEXP algebra) {
    if(OMX_DEBUG) { Rprintf("Populating RAM Attributes.\n"); }

	omxRAMObjective* oro = (omxRAMObjective*) (oo->argStruct);
	omxMatrix* A = oro->A;
	omxMatrix* S = oro->S;
	omxMatrix* X = oro->X;
	omxMatrix* Ax= oro->Ax;
	omxMatrix* Z = oro->Z;
	omxMatrix* I = oro->I;
    int numIters = oro->numIters;
    double oned = 1.0, zerod = 0.0;
    
    omxRecompute(A);
    omxRecompute(S);
	
	omxFastRAMInverse(numIters, A, Z, Ax, I ); // Z = (I-A)^-1
	
	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c \n %d %d %d \n %f \n %x %d %x %d \n %f %x %d.\n", *(Z->majority), *(S->majority), (Z->rows), (S->cols), (S->rows), oned, Z->data, (Z->leading), S->data, (S->leading), zerod, Ax->data, (Ax->leading));}
	// F77_CALL(omxunsafedgemm)(Z->majority, S->majority, &(Z->rows), &(S->cols), &(S->rows), &oned, Z->data, &(Z->leading), S->data, &(S->leading), &zerod, Ax->data, &(Ax->leading)); 	// X = ZS
	omxDGEMM(FALSE, FALSE, oned, Z, S, zerod, Ax);

	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(Ax->majority), *(Z->minority), (Ax->rows), (Z->rows), (Z->cols), oned, X->data, (X->leading), Z->data, (Z->lagging), zerod, Ax->data, (Ax->leading));}
	// F77_CALL(omxunsafedgemm)(Ax->majority, Z->minority, &(Ax->rows), &(Z->rows), &(Z->cols), &oned, Ax->data, &(Ax->leading), Z->data, &(Z->leading), &zerod, Ax->data, &(Ax->leading)); 
	omxDGEMM(FALSE, TRUE, oned, Ax, Z, zerod, Ax);
	// Ax = ZSZ' = Covariance matrix including latent variables
	
	SEXP expCovExt;
	PROTECT(expCovExt = allocMatrix(REALSXP, Ax->rows, Ax->cols));
	for(int row = 0; row < Ax->rows; row++)
		for(int col = 0; col < Ax->cols; col++)
			REAL(expCovExt)[col * Ax->rows + row] =
				omxMatrixElement(Ax, row, col);
	setAttrib(algebra, install("UnfilteredExpCov"), expCovExt);
	UNPROTECT(1);
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

void omxFastRAMInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I ) {
	
	omxMatrix* origZ = Z;
    double oned = 1, minusOned = -1.0;

	if(numIters == NA_INTEGER) {
		int ipiv[A->rows], lwork = 4 * A->rows * A->cols, k;		// TODO: Speedups can be made by preallocating this.
		double work[lwork];
		if(OMX_DEBUG_ALGEBRA) { Rprintf("RAM Algebra (I-A) inversion using standard (general) inversion.\n"); }

		/* Z = (I-A)^-1 */
		if(I->colMajor != A->colMajor) {
			omxTransposeMatrix(I);
		}
		omxCopyMatrix(Z, A, TRUE);

		/* Z = (I-A)^-1 */
		// F77_CALL(omxunsafedgemm)(I->majority, Z->majority, &(I->cols), &(I->rows), &(Z->rows), &oned, I->data, &(I->cols), I->data, &(I->cols), &minusOned, Z->data, &(Z->cols));
		omxDGEMM(FALSE, FALSE, oned, I, Z, minusOned, Z);

		// F77_CALL(dgetrf)(&(Z->rows), &(Z->cols), Z->data, &(Z->leading), ipiv, &k);
		k = omxDGETRF(Z, ipiv);
		if(OMX_DEBUG) { Rprintf("Info on LU Decomp: %d\n", k); }
		if(k > 0) {
		        char errStr[250];
		        strncpy(errStr, "(I-A) is exactly singular.", 100);
		        omxRaiseError(A->currentState, -1, errStr);                    // Raise Error
		        return;
		}
		// F77_CALL(dgetri)(&(Z->rows), Z->data, &(Z->leading), ipiv, work, &lwork, &k);
		k = omxDGETRI(Z, ipiv, work, lwork);
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
	
		omxCopyMatrix(Z, A, TRUE);
	
		/* Optimized I-A inversion: Z = (I-A)^-1 */
		// F77_CALL(omxunsafedgemm)(I->majority, A->majority, &(I->cols), &(I->rows), &(A->rows), &oned, I->data, &(I->cols), I->data, &(I->cols), &oned, Z->data, &(Z->cols));  // Z = I + A = A^0 + A^1
		omxDGEMM(FALSE, FALSE, 1.0, I, I, 1.0, Z); // Z == A + I

		for(int i = 1; i <= numIters; i++) { // TODO: Efficiently determine how many times to do this
			// The sequence goes like this: (I + A), I + (I + A) * A, I + (I + (I + A) * A) * A, ...
			// Which means only one DGEMM per iteration.
			if(OMX_DEBUG_ALGEBRA) { Rprintf("....RAM: Iteration #%d/%d\n", i, numIters);}
			omxCopyMatrix(Ax, I, TRUE);
			// F77_CALL(omxunsafedgemm)(A->majority, A->majority, &(Z->cols), &(Z->rows), &(A->rows), &oned, Z->data, &(Z->cols), A->data, &(A->cols), &oned, Ax->data, &(Ax->cols));  // Ax = Z %*% A + I
			omxDGEMM(FALSE, FALSE, oned, A, Z, oned, Ax);
			omxMatrix* m = Z; Z = Ax; Ax = m;	// Juggle to make Z equal to Ax
		}
		if(origZ != Z) { 	// Juggling has caused Ax and Z to swap
			omxCopyMatrix(Z, Ax, TRUE);
		}
	}
}

void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax) {
	
	if(OMX_DEBUG) { Rprintf("Running RAM computation."); }
		
	double oned = 1.0, zerod=0.0;
	
	if(Ax == NULL || I == NULL || Z == NULL || Y == NULL || X == NULL) {
		error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
	}
		
	if(Cov == NULL && Means == NULL) {
		return; // We're not populating anything, so why bother running the calculation?
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
	
	omxFastRAMInverse(numIters, A, Z, Ax, I );
		
	/* Cov = FZSZ'F' */
	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(F->majority), *(Z->majority), (F->rows), (Z->cols), (Z->rows), oned, F->data, (F->leading), Z->data, (Z->leading), zerod, Y->data, (Y->leading));}
	// F77_CALL(omxunsafedgemm)(F->majority, Z->majority, &(F->rows), &(Z->cols), &(Z->rows), &oned, F->data, &(F->leading), Z->data, &(Z->leading), &zerod, Y->data, &(Y->leading)); 	// Y = FZ
	omxDGEMM(FALSE, FALSE, 1.0, F, Z, 0.0, Y);

	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(Y->majority), *(S->majority), (Y->rows), (S->cols), (S->rows), oned, Y->data, (Y->leading), S->data, (S->leading), zerod, X->data, (X->leading));}
	// F77_CALL(omxunsafedgemm)(Y->majority, S->majority, &(Y->rows), &(S->cols), &(S->rows), &oned, Y->data, &(Y->leading), S->data, &(S->leading), &zerod, X->data, &(X->leading)); 	// X = FZS
	omxDGEMM(FALSE, FALSE, 1.0, Y, S, 0.0, X);

	if(OMX_DEBUG_ALGEBRA) { Rprintf("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.\n", *(X->majority), *(Y->minority), (X->rows), (Y->rows), (Y->cols), oned, X->data, (X->leading), Y->data, (Y->lagging), zerod, Cov->data, (Cov->leading));}
	// F77_CALL(omxunsafedgemm)(X->majority, Y->minority, &(X->rows), &(Y->rows), &(Y->cols), &oned, X->data, &(X->leading), Y->data, &(Y->leading), &zerod, Cov->data, &(Cov->leading));
	omxDGEMM(FALSE, TRUE, 1.0, X, Y, 0.0, Cov);
	 // Cov = FZSZ'F' (Because (FZ)' = Z'F')
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Cov, "....RAM: Model-implied Covariance Matrix:");}
	
	if(M != NULL && Means != NULL) {
		// F77_CALL(omxunsafedgemv)(Y->majority, &(Y->rows), &(Y->cols), &oned, Y->data, &(Y->leading), M->data, &onei, &zerod, Means->data, &onei);
		omxDGEMV(FALSE, 1.0, Y, M, 0.0, Means);
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....RAM: Model-implied Means Vector:");}
	}

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

	if(OMX_DEBUG) { Rprintf("Initializing RAM objective.\n"); }
	
	int l, k;

	SEXP slotValue;
	
	/* Create and register subobjective */

    omxObjective *subObjective = omxCreateSubObjective(oo);
	
	omxRAMObjective *RAMobj = (omxRAMObjective*) R_alloc(1, sizeof(omxRAMObjective));
	
	/* Set Subobjective Calls and Structures */
	subObjective->objectiveFun = omxCallRAMObjective;
	subObjective->needsUpdateFun = omxNeedsUpdateRAMObjective;
	subObjective->destructFun = omxDestroyRAMObjective;
	subObjective->setFinalReturns = NULL;
	subObjective->populateAttrFun = omxPopulateRAMAttributes;
	subObjective->argStruct = (void*) RAMobj;
	
	/* Set up objective structures */
	if(OMX_DEBUG) { Rprintf("Initializing RAM Meta Data for objective function.\n"); }

	if(OMX_DEBUG) { Rprintf("Processing M.\n"); }
	RAMobj->M = omxNewMatrixFromIndexSlot(rObj, currentState, "M");

	if(OMX_DEBUG) { Rprintf("Processing A.\n"); }
	RAMobj->A = omxNewMatrixFromIndexSlot(rObj, currentState, "A");

	if(OMX_DEBUG) { Rprintf("Processing S.\n"); }
	RAMobj->S = omxNewMatrixFromIndexSlot(rObj, currentState, "S");

	if(OMX_DEBUG) { Rprintf("Processing F.\n"); }
	RAMobj->F = omxNewMatrixFromIndexSlot(rObj, currentState, "F");

	if(OMX_DEBUG) { Rprintf("Processing usePPML.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("usePPML")));
	RAMobj->usePPML = INTEGER(slotValue)[0]; 
	UNPROTECT(1);

	if(RAMobj->usePPML) {
		PROTECT(slotValue = GET_SLOT(rObj, install("ppmlData")));
		RAMobj->ppmlData = omxNewDataFromMxData(NULL, slotValue, currentState);
		UNPROTECT(1);

		RAMobj->cov = omxDataMatrix(RAMobj->ppmlData, NULL);

		if(OMX_DEBUG) { Rprintf("Processing PPML observed means.\n"); }
		RAMobj->ppmlMeans = omxDataMeans(RAMobj->ppmlData, 0, NULL);
		if(OMX_DEBUG && RAMobj->means == NULL) { Rprintf("RAM: No PPML Observed Means.\n"); }
	} else {
		RAMobj->ppmlData  = NULL;
		RAMobj->ppmlCov   = NULL;
		RAMobj->ppmlMeans = NULL;
	}

	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { Rprintf("Generating I.\n"); }
	RAMobj->I = omxNewIdentityMatrix(RAMobj->A->rows, currentState);
	omxRecompute(RAMobj->I);

	if(OMX_DEBUG) { Rprintf("Processing expansion iteration depth.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("depth")));
	RAMobj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { Rprintf("Using %d iterations.", RAMobj->numIters); }
	UNPROTECT(1);

	l = RAMobj->F->rows;
	k = RAMobj->A->cols;

	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }

	RAMobj->Z = 	omxInitMatrix(NULL, k, k, TRUE, currentState);
	RAMobj->Ax = 	omxInitMatrix(NULL, k, k, TRUE, currentState);
	RAMobj->Y = 	omxInitMatrix(NULL, l, k, TRUE, currentState);
	RAMobj->X = 	omxInitMatrix(NULL, l, k, TRUE, currentState);

	RAMobj->cov = 		omxInitMatrix(NULL, l, l, TRUE, currentState);

	if(RAMobj->M != NULL) {
		RAMobj->means = 	omxInitMatrix(NULL, 1, l, TRUE, currentState);
	} else RAMobj->means  = 	NULL;

	/* Create parent objective */

	omxCreateMLObjective(oo, rObj, RAMobj->cov, RAMobj->means);
	
}
