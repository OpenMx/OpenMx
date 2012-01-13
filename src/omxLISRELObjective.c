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
	if(oro->TX != NULL) {
	    omxRecompute(oro->TX);
	    omxRecompute(oro->TY);
	    omxRecompute(oro->KA);
	    omxRecompute(oro->AL);
	}

	omxCalculateLISRELCovarianceAndMeans(oro->LX, oro->LY, oro->BE, oro->GA, oro->cov, oro->means, oro->numIters, oro->I, oro->Z, oro->Y, oro->X, oro->Ax);
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
	omxFreeMatrixData(argStruct->X);
	omxFreeMatrixData(argStruct->Y);
	omxFreeMatrixData(argStruct->Z);
	omxFreeMatrixData(argStruct->Ax);

	if(argStruct->ppmlData != NULL) 
		omxFreeData(argStruct->ppmlData);
	
}

void omxPopulateLISRELAttributes(omxObjective *oo, SEXP algebra) {
    if(OMX_DEBUG) { Rprintf("Populating LISREL Attributes.  Currently this dones nothing!\n"); }
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
	omxMatrix* X = oro->X;
	omxMatrix* Ax= oro->Ax;
	omxMatrix* Z = oro->Z;
	omxMatrix* I = oro->I;
    int numIters = oro->numIters;
    double oned = 1.0, zerod = 0.0;
    
    omxRecompute(LX);
    omxRecompute(LY);
	
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


void omxCalculateLISRELCovarianceAndMeans(omxMatrix* LX, omxMatrix* LY, omxMatrix* BE, omxMatrix* GA, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax) {
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

/*
void omxInitRAMObjective(omxObjective* oo, SEXP rObj) {
	
	omxState* currentState = oo->matrix->currentState;	

	if(OMX_DEBUG) { Rprintf("Initializing RAM objective.\n"); }
	
	int l, k;

	SEXP slotValue;
	
	// IMPORTANT: Create and register subobjective

    omxObjective *subObjective = omxCreateSubObjective(oo);
	
	omxRAMObjective *RAMobj = (omxRAMObjective*) R_alloc(1, sizeof(omxRAMObjective));
	
	// IMPORTANT: Set Subobjective Calls and Structures
	subObjective->objectiveFun = omxCallRAMObjective;
	subObjective->needsUpdateFun = omxNeedsUpdateRAMObjective;
	subObjective->destructFun = omxDestroyRAMObjective;
	subObjective->setFinalReturns = NULL;
	subObjective->populateAttrFun = omxPopulateRAMAttributes;
	subObjective->updateChildObjectiveFun = omxUpdateChildRAMObjective;
	subObjective->argStruct = (void*) RAMobj;
	
	// IMPORTANT Set up objective structures
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

	// IMPORTANT: Identity Matrix, Size Of A
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

	//IMPORTANT: Create parent objective

	omxCreateMLObjective(oo, rObj, RAMobj->cov, RAMobj->means);
	
}
*/
