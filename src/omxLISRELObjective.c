/*
 *  Copyright 2007-2012 The OpenMx Project
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
#include "omxLISRELObjective.h"

extern void omxCreateMLObjective(omxObjective* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means);

void omxCallLISRELObjective(omxObjective* oo) {
    if(OMX_DEBUG) { Rprintf("LISREL Subobjective called.\n"); }
	omxLISRELObjective* oro = (omxLISRELObjective*)(oo->argStruct);
	
	omxRecompute(oro->LX);
	omxRecompute(oro->LY);
	omxRecompute(oro->BE);
	omxRecompute(oro->GA);
	omxRecompute(oro->PH);
	omxRecompute(oro->PS);
	omxRecompute(oro->TD);
	omxRecompute(oro->TE);
	omxRecompute(oro->TH);
	/*if(oro->TX != NULL) {     // Update means?
	    omxRecompute(oro->TX);
	    omxRecompute(oro->TY);
	    omxRecompute(oro->KA);
	    omxRecompute(oro->AL);
	}*/

	omxCalculateLISRELCovarianceAndMeans(oro->LX, oro->LY, oro->BE, oro->GA, oro->PH, oro->PS, oro->TD, oro->TE, oro->TH, oro->cov, oro->means, oro->numIters, oro->I, oro->LXPH, oro->W, oro->GAPH, oro->U, oro->TOP, oro->BOT);
}

void omxDestroyLISRELObjective(omxObjective* oo) {

	if(OMX_DEBUG) { Rprintf("Destroying LISREL Objective.\n"); }
	
	omxLISRELObjective* argStruct = (omxLISRELObjective*)(oo->argStruct);

	/* We allocated 'em, so we destroy 'em. */
	if(argStruct->cov != NULL)
		omxFreeMatrixData(argStruct->cov);

	if(argStruct->means != NULL)
		omxFreeMatrixData(argStruct->means);
	
	omxFreeMatrixData(argStruct->I);
	omxFreeMatrixData(argStruct->U);
	omxFreeMatrixData(argStruct->W);
	omxFreeMatrixData(argStruct->LXPH);
	omxFreeMatrixData(argStruct->GAPH);


	if(argStruct->ppmlData != NULL) 
		omxFreeData(argStruct->ppmlData);
	
}

void omxPopulateLISRELAttributes(omxObjective *oo, SEXP algebra) {
    if(OMX_DEBUG) { Rprintf("Populating LISREL Attributes.  Currently this does very little!\n"); }
	
	/*
	omxLISRELObjective* oro = (omxLISRELObjective*) (oo->argStruct);
	omxMatrix* LX = oro->LX;
	omxMatrix* LY = oro->LY;
	omxMatrix* BE = oro->BE;
	omxMatrix* GA = oro->GA;
	omxMatrix* PH = oro->PH;
	omxMatrix* PS = oro->PS;
	omxMatrix* TD = oro->TD;
	omxMatrix* TE = oro->TE;
	omxMatrix* TH = oro->TH;
	omxMatrix* LXPH = oro->LXPH;
	omxMatrix* GAPH = oro->GAPH;
	omxMatrix* W = oro->W;
	omxMatrix* U = oro->U;
	omxMatrix* I = oro->I;
	int numIters = oro->numIters;
	double oned = 1.0, zerod = 0.0;
	
	omxRecompute(LX);
	omxRecompute(LY);
	*/ //This block of code works fine but because I do not use any of it later, it throws a huge block of warnings about unused variables.
	// In general, I do not yet understand what this function needs to do.
	
	/*
	omxFastRAMInverse(numIters, BE, Z, Ax, I ); // Z = (I-A)^-1
	
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
	*/
}

/* omxFastLISRELInverse would go here */


void omxCalculateLISRELCovarianceAndMeans(omxMatrix* LX, omxMatrix* LY, omxMatrix* BE, omxMatrix* GA, omxMatrix* PH, omxMatrix* PS,  omxMatrix* TD, omxMatrix* TE, omxMatrix* TH, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* LXPH, omxMatrix* W, omxMatrix* GAPH, omxMatrix* U, omxMatrix* TOP, omxMatrix* BOT) {
	if(OMX_DEBUG) { Rprintf("Running LISREL computation."); }
	double oned = 1.0, zerod=0.0, minusOned = -1.0;
	int ipiv[BE->rows], lwork = 4 * BE->rows * BE->cols; //This is copied from omxFastRAMInverse()
	double work[lwork];									// It lets you get the inverse of a matrix via omxDGETRI()
	omxMatrix** args = R_alloc(2, sizeof(omxMatrix*)); // Used to block construct covariance matrix
	// Note: the above give the warning message: initialization from incompatible pointer type
	
	/* Calculate the lower right quadrant: the covariance of the Xs */
	omxDGEMM(FALSE, FALSE, oned, LX, PH, zerod, LXPH);
	omxDGEMM(FALSE, TRUE, oned, LXPH, LX, oned, TD);
	
	/* Calculate (I-BE)^(-1) and LY*(I-BE)^(-1) */
	omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, BE);
	omxDGETRI(BE, ipiv, work, lwork);
	omxDGEMM(FALSE, FALSE, oned, LY, BE, zerod, LY);
	
	/* Calculate the lower left quadrant: the covariance of Xs and Ys, nX by nY */
	omxDGEMM(FALSE, TRUE, oned, LXPH, GA, zerod, W);
	omxDGEMM(FALSE, TRUE, oned, W, LY, oned, TH);
	
	/* Calculate the upper right quadrant: NOTE THIS IS MERELY THE LOWER LEFT QUADRANT TRANSPOSED. */
	//DONE as omxTranspose(TH)
	
	/* Calculate the upper left quadrant: the covariance of the Ys */
	omxDGEMM(FALSE, FALSE, oned, GA, PH, zerod, GAPH);
	omxDGEMM(FALSE, TRUE, oned, GAPH, GA, oned, PS);
	omxDGEMM(FALSE, FALSE, oned, LY, PS, zerod, U);
	omxDGEMM(FALSE, TRUE, oned, U, LY, oned, TE);
	
	/* Construct the full model-implied covariance matrix from the blocks previously calculated */
	// SigmaHat = ( TE  t(TH) )
	//            ( TH    TD  )
	args[0] = TH;
	args[1] = TD;
	omxMatrixHorizCat(args, 2, BOT);
	args[0] = TD;
	omxTransposeMatrix(TH);
	args[1] = TH;
	omxMatrixHorizCat(args, 2, TOP);
	omxTransposeMatrix(TH);
	// So that it's back where it was.
	args[0] = TOP;
	args[1] = BOT;
	omxMatrixVertCat(args, 2, Cov);
	
/*	
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
		
	// IMPORTANT: Cov = FZSZ'F'
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
*/
}

/*
void omxUpdateChildRAMObjective(omxObjective* tgt, omxObjective* src) {

	omxRAMObjective* tgtRAM = (omxRAMObjective*)(tgt->argStruct);
	omxRAMObjective* srcRAM = (omxRAMObjective*)(src->argStruct);

	omxUpdateMatrix(tgtRAM->cov, srcRAM->cov);
	if (tgtRAM->means && srcRAM->means) {
		omxUpdateMatrix(tgtRAM->means, srcRAM->means);	
	}

	if (tgt->subObjective != NULL) {
		tgt->subObjective->updateChildObjectiveFun(tgt->subObjective, src->subObjective);
	}

}
*/

/*
unsigned short int omxNeedsUpdateRAMObjective(omxObjective* oo) {
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Checking if RAM needs update.  RAM uses A:0x%x, S:0x%x, F:0x%x, M:0x%x.\n",((omxRAMObjective*)oo->argStruct)->A, ((omxRAMObjective*)oo->argStruct)->S, ((omxRAMObjective*)oo->argStruct)->F, ((omxRAMObjective*)oo->argStruct)->M); }
	return(omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->A)
	 	|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->S)
	 	|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->F)
		|| omxNeedsUpdate(((omxRAMObjective*)oo->argStruct)->M));

	// Note: cov is data, and should never need updating.
}
*/


void omxInitLISRELObjective(omxObjective* oo, SEXP rObj) {
	
	if(OMX_DEBUG) { Rprintf("Initializing LISREL objective.\n"); }

	omxState* currentState = oo->matrix->currentState;	
	
	int nx, nxi, ny, neta, ntotal;

	SEXP slotValue;
	
	/* Create and register subobjective */

    omxObjective *subObjective = omxCreateSubObjective(oo);
	
	omxLISRELObjective *LISobj = (omxLISRELObjective*) R_alloc(1, sizeof(omxLISRELObjective));
	
	/* Set Subobjective Calls and Structures */
	subObjective->objectiveFun = omxCallLISRELObjective;
	subObjective->needsUpdateFun = NULL; //omxNeedsUpdateRAMObjective;
	subObjective->destructFun = omxDestroyLISRELObjective;
	subObjective->setFinalReturns = NULL;
	subObjective->populateAttrFun = omxPopulateLISRELAttributes;
	subObjective->updateChildObjectiveFun = NULL; //omxUpdateChildRAMObjective;
	subObjective->argStruct = (void*) LISobj;
	
	/* Set up objective structures */
	if(OMX_DEBUG) { Rprintf("Initializing LISREL Meta Data for objective function.\n"); }
	
	if(OMX_DEBUG) { Rprintf("Processing LX.\n"); }
	LISobj->LX = omxNewMatrixFromIndexSlot(rObj, currentState, "LX");
	
	if(OMX_DEBUG) { Rprintf("Processing LY.\n"); }
	LISobj->LY = omxNewMatrixFromIndexSlot(rObj, currentState, "LY");
	
	if(OMX_DEBUG) { Rprintf("Processing BE.\n"); }
	LISobj->BE = omxNewMatrixFromIndexSlot(rObj, currentState, "BE");
	
	if(OMX_DEBUG) { Rprintf("Processing GA.\n"); }
	LISobj->GA = omxNewMatrixFromIndexSlot(rObj, currentState, "GA");
	
	if(OMX_DEBUG) { Rprintf("Processing PH.\n"); }
	LISobj->PH = omxNewMatrixFromIndexSlot(rObj, currentState, "PH");
	
	if(OMX_DEBUG) { Rprintf("Processing PS.\n"); }
	LISobj->PS = omxNewMatrixFromIndexSlot(rObj, currentState, "PS");
	
	if(OMX_DEBUG) { Rprintf("Processing TD.\n"); }
	LISobj->TD = omxNewMatrixFromIndexSlot(rObj, currentState, "TD");
	
	if(OMX_DEBUG) { Rprintf("Processing TE.\n"); }
	LISobj->TE = omxNewMatrixFromIndexSlot(rObj, currentState, "TE");
	
	if(OMX_DEBUG) { Rprintf("Processing TH.\n"); }
	LISobj->TH = omxNewMatrixFromIndexSlot(rObj, currentState, "TH");
	
	// TODO: Add means specification
	
	/* PPML Code: Perhaps comment out this block */
	/*
	if(OMX_DEBUG) { Rprintf("Processing usePPML.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("usePPML")));
	LISobj->usePPML = INTEGER(slotValue)[0]; 
	UNPROTECT(1);

	if(LISobj->usePPML) {
		PROTECT(slotValue = GET_SLOT(rObj, install("ppmlData")));
		LISobj->ppmlData = omxNewDataFromMxData(NULL, slotValue, currentState);
		UNPROTECT(1);

		LISobj->cov = omxDataMatrix(LISobj->ppmlData, NULL);

		if(OMX_DEBUG) { Rprintf("Processing PPML observed means.\n"); }
		LISobj->ppmlMeans = omxDataMeans(LISobj->ppmlData, 0, NULL);
		if(OMX_DEBUG && LISobj->means == NULL) { Rprintf("LISREL: No PPML Observed Means.\n"); }
	} else {
		LISobj->ppmlData  = NULL;
		LISobj->ppmlCov   = NULL;
		LISobj->ppmlMeans = NULL;
	}
	*/
	
	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { Rprintf("Generating I.\n"); }
	LISobj->I = omxNewIdentityMatrix(LISobj->BE->rows, currentState);
	omxRecompute(LISobj->I);
	
	/*
	if(OMX_DEBUG) { Rprintf("Processing expansion iteration depth.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("depth")));
	LISobj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { Rprintf("Using %d iterations.", LISobj->numIters); }
	UNPROTECT(1);
	*/
	
	/* Initialize the place holder matrices used in calculations */
	nx = LISobj->LX->rows;
	nxi = LISobj->LX->cols;
	ny = LISobj->LY->rows;
	neta = LISobj->LY->cols;
	ntotal = nx + ny;
	
	
	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }
	
	LISobj->LXPH = 	omxInitMatrix(NULL, nx, nxi, TRUE, currentState);
	LISobj->W = 	omxInitMatrix(NULL, nx, neta, TRUE, currentState);
	LISobj->GAPH = 	omxInitMatrix(NULL, neta, nxi, TRUE, currentState);
	LISobj->U = 	omxInitMatrix(NULL, ny, neta, TRUE, currentState);
	LISobj->TOP = 	omxInitMatrix(NULL, ny, ntotal, TRUE, currentState);
	LISobj->BOT = 	omxInitMatrix(NULL, nx, ntotal, TRUE, currentState);


	LISobj->cov = 		omxInitMatrix(NULL, ntotal, ntotal, TRUE, currentState);
	
	/* Uncomment this when means are implemented
	if(RAMobj->M != NULL) {
		RAMobj->means = 	omxInitMatrix(NULL, 1, l, TRUE, currentState);
	} else RAMobj->means  = 	NULL;
	*/

	/* Create parent objective */

	omxCreateMLObjective(oo, rObj, LISobj->cov, LISobj->means);
	
}

