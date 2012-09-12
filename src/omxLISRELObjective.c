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
#include "omxRAMObjective.h" //included for omxFastRAMInverse

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
	if(oro->TX != NULL) {     // Update means?
		omxRecompute(oro->TX);
		omxRecompute(oro->KA);
	}
	if(oro->TY != NULL) {
		omxRecompute(oro->TY);
		omxRecompute(oro->AL);
	}
	
	omxCalculateLISRELCovarianceAndMeans(oro);
}

void omxDestroyLISRELObjective(omxObjective* oo) {

	if(OMX_DEBUG) { Rprintf("Destroying LISREL Objective.\n"); }
	
	omxLISRELObjective* argStruct = (omxLISRELObjective*)(oo->argStruct);

	/* We allocated 'em, so we destroy 'em. */
	if(argStruct->cov != NULL)
		omxFreeMatrixData(argStruct->cov);

	if(argStruct->means != NULL)
		omxFreeMatrixData(argStruct->means);
	
	omxFreeMatrixData(argStruct->A);
	omxFreeMatrixData(argStruct->B);
	omxFreeMatrixData(argStruct->C);
	omxFreeMatrixData(argStruct->D);
	omxFreeMatrixData(argStruct->E);
	omxFreeMatrixData(argStruct->F);
	omxFreeMatrixData(argStruct->G);
	omxFreeMatrixData(argStruct->H);
	omxFreeMatrixData(argStruct->I);
	omxFreeMatrixData(argStruct->J);
	omxFreeMatrixData(argStruct->K);
	omxFreeMatrixData(argStruct->L);
	omxFreeMatrixData(argStruct->TOP);
	omxFreeMatrixData(argStruct->BOT);
	omxFreeMatrixData(argStruct->MUX);
	omxFreeMatrixData(argStruct->MUY);
	
	
	/* Comment out the ppml things I do not use.
	if(argStruct->ppmlData != NULL) 
		omxFreeData(argStruct->ppmlData);
	*/
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


void omxCalculateLISRELCovarianceAndMeans(omxLISRELObjective* oro) {
	omxMatrix* LX = oro->LX;
	omxMatrix* LY = oro->LY;
	omxMatrix* BE = oro->BE;
	omxMatrix* GA = oro->GA;
	omxMatrix* PH = oro->PH;
	omxMatrix* PS = oro->PS;
	omxMatrix* TD = oro->TD;
	omxMatrix* TE = oro->TE;
	omxMatrix* TH = oro->TH;
	omxMatrix* TX = oro->TX;
	omxMatrix* TY = oro->TY;
	omxMatrix* KA = oro->KA;
	omxMatrix* AL = oro->AL;
	omxMatrix* Cov = oro->cov;
	omxMatrix* Means = oro->means;
	int numIters = oro->numIters; //Used for fast RAM/LISREL inverse
	omxMatrix* A = oro->A;
	omxMatrix* B = oro->B;
	omxMatrix* C = oro->C;
	omxMatrix* D = oro->D;
	omxMatrix* E = oro->E;
	omxMatrix* F = oro->F;
	omxMatrix* G = oro->G;
	omxMatrix* H = oro->H;
	omxMatrix* I = oro->I;
	omxMatrix* J = oro->J;
	omxMatrix* K = oro->K;
	omxMatrix* L = oro->L;
	omxMatrix* TOP = oro->TOP;
	omxMatrix* BOT = oro->BOT;
	omxMatrix* MUX = oro->MUX;
	omxMatrix* MUY = oro->MUY;
	omxMatrix** args = oro->args;
	if(OMX_DEBUG) { Rprintf("Running LISREL computation in omxCalculateLISRELCovarianceAndMeans.\n"); }
	double oned = 1.0, zerod=0.0; //, minusOned = -1.0;
	//int ipiv[BE->rows], lwork = 4 * BE->rows * BE->cols; //This is copied from omxFastRAMInverse()
	//double work[lwork];									// It lets you get the inverse of a matrix via omxDGETRI()
	
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(LX, "....LISREL: LX:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(LY, "....LISREL: LY:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(BE, "....LISREL: BE:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(GA, "....LISREL: GA:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(PH, "....LISREL: PH:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(PS, "....LISREL: PS:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(TD, "....LISREL: TD:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(TE, "....LISREL: TE:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(TH, "....LISREL: TH:");}
	
	/* Calculate the lower right quadrant: the covariance of the Xs */
	if(LX->cols != 0 && LY->cols != 0) {
	//if( (LX != NULL) && (LY != NULL) ) {
		if(OMX_DEBUG) {Rprintf("Calculating Lower Right Quadrant of Expected Covariance Matrix.\n"); }
		omxDGEMM(FALSE, FALSE, oned, LX, PH, zerod, A); // A = LX*PH
		omxCopyMatrix(B, TD); // B = TD
		omxDGEMM(FALSE, TRUE, oned, A, LX, oned, B);  // B = LX*PH*LX^T + TD
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(B, "....LISREL: Lower Right Quadrant of Model-implied Covariance Matrix:");}
		
		/* Calculate (I-BE)^(-1) and LY*(I-BE)^(-1) */
		if(OMX_DEBUG) {Rprintf("Calculating Inverse of I-BE.\n"); }
		omxFastRAMInverse(numIters, BE, C, L, I ); // C = (I-BE)^-1
		//omxCopyMatrix(C, BE); // C = BE
		//omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, C); // C = I - BE
		//omxDGETRF(C, ipiv); //LU Decomp
		//omxDGETRI(C, ipiv, work, lwork); //Inverse based on LU Decomp ... C = C^(-1) = (I - BE)^(-1)
		
		
		omxDGEMM(FALSE, FALSE, oned, LY, C, zerod, D); // D = LY*C = LY * (I - BE)^(-1)
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....LISREL: LY*(I-BE)^(-1)");}
		
		/* Calculate the lower left quadrant: the covariance of Xs and Ys, nX by nY */
		if(OMX_DEBUG) {Rprintf("Calculating Lower Left Quadrant of Expected Covariance Matrix.\n"); }
		omxDGEMM(FALSE, TRUE, oned, A, GA, zerod, E); // E = A*GA^T = LX*PH*GA^T
		omxCopyMatrix(F, TH); // F = TH
		omxDGEMM(FALSE, TRUE, oned, E, D, oned, F); // F = E*D^T + F = LX*PH*GA^T * (LY * (I - BE)^(-1))^T + TH
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(F, "....LISREL: Lower Left Quadrant of Model-implied Covariance Matrix:");}
	
		
	/* Calculate the upper right quadrant: NOTE THIS IS MERELY THE LOWER LEFT QUADRANT TRANSPOSED. */
	//DONE as omxTranspose(F)
		
		
		/* Calculate the upper left quadrant: the covariance of the Ys */
		if(OMX_DEBUG) {Rprintf("Calculating Upper Left Quadrant of Expected Covariance Matrix.\n"); }
		if(OMX_DEBUG_ALGEBRA) {Rprintf("Calculating G = GA*PH.\n");}
		omxDGEMM(FALSE, FALSE, oned, GA, PH, zerod, G); // G = GA*PH
		if(OMX_DEBUG_ALGEBRA) {Rprintf("Copying C = PS.\n");}
		omxCopyMatrix(C, PS); // C = PS
		if(OMX_DEBUG_ALGEBRA) {Rprintf("Calculating G = GA*PH.\n");}
		omxDGEMM(FALSE, TRUE, oned, G, GA, oned, C); // C = G*GA^T + C = GA*PH*GA^T + PS
		omxDGEMM(FALSE, FALSE, oned, D, C, zerod, H); // H = D*C = LY * (I - BE)^(-1) * (GA*PH*GA^T + PS)
		omxCopyMatrix(J, TE); // J = TE
		omxDGEMM(FALSE, TRUE, oned, H, D, oned, J); // J = H*D^T + J = LY * (I - BE)^(-1) * (GA*PH*GA^T + PS) * (LY * (I - BE)^(-1))^T + TE
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(J, "....LISREL: Upper Left Quadrant of Model-implied Covariance Matrix:");}
	
	
	/* Construct the full model-implied covariance matrix from the blocks previously calculated */
	// SigmaHat = ( J  t(F) )
	//            ( F    B  )
		args[0] = F;
		args[1] = B;
		omxMatrixHorizCat(args, 2, BOT);
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(BOT, "....LISREL: BOT = cbind(F, B):");}
		args[0] = J;
		omxTransposeMatrix(F);
		args[1] = F;
		omxMatrixHorizCat(args, 2, TOP);
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(args[0], "....LISREL: TOP Debugging args[0] = J:");}
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(args[1], "....LISREL: TOP Debugging args[1] = F:");}
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(F, "....LISREL: TOP Debugging F (should be 2 rows):");}
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(TOP, "....LISREL: TOP = cbind(J, t(F)):");}
		omxTransposeMatrix(F); // So that it's back where it was.
		args[0] = TOP;
		args[1] = BOT;
		omxMatrixVertCat(args, 2, Cov);
	
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Cov, "....LISREL: Model-implied Covariance Matrix:");}
	
	
		/* Now Calculate the Expected Means */
		if(Means != NULL) {
			/* Mean of the Xs */
			//if(TX != NULL) {
				omxCopyMatrix(MUX, TX);
				omxDGEMV(FALSE, oned, LX, KA, oned, MUX);
			//}
			
			/* Mean of Ys */
			//if(TY != NULL) {
				omxCopyMatrix(K, AL);
				omxDGEMV(FALSE, oned, GA, KA, oned, K);
				omxCopyMatrix(MUY, TY);
				omxDGEMV(FALSE, oned, D, K, oned, MUY);
			//}
		
			/* Build means from blocks */
			args[0] = MUY;
			args[1] = MUX;
			omxMatrixVertCat(args, 2, Means);
			
			if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....LISREL: Model-implied Means Vector:");}
		}
	}
	else if(LX->cols != 0) { /* IF THE MODEL ONLY HAS Xs */
	//else if(LX != NULL) { /* IF THE MODEL ONLY HAS Xs */
		if(OMX_DEBUG) {Rprintf("Calculating Lower Right Quadrant of Expected Covariance Matrix.\n"); }
		omxDGEMM(FALSE, FALSE, oned, LX, PH, zerod, A); // A = LX*PH
		omxCopyMatrix(Cov, TD); // Cov = TD
		omxDGEMM(FALSE, TRUE, oned, A, LX, oned, Cov);  // Cov = LX*PH*LX^T + Cov
		if(Means != NULL) {
				/* Mean of the Xs */
				omxCopyMatrix(Means, TX);
				omxDGEMV(FALSE, oned, LX, KA, oned, Means);
		}
	}
	
	/* IF THE MODEL ONLY HAS Ys */
	else if(LY->cols != 0) {
	//else if(LY != NULL) {
		/* Calculate (I-BE)^(-1) and LY*(I-BE)^(-1) */
		if(OMX_DEBUG) {Rprintf("Calculating Inverse of I-BE.\n"); }
		omxFastRAMInverse(numIters, BE, C, L, I ); // C = (I-BE)^-1
		//omxCopyMatrix(C, BE); // C = BE
		//omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, C); // C = I - BE
		//omxDGETRF(C, ipiv); //LU Decomp
		//omxDGETRI(C, ipiv, work, lwork); //Inverse based on LU Decomp ... C = C^(-1) = (I - BE)^(-1)
		omxDGEMM(FALSE, FALSE, oned, LY, C, zerod, D); // D = LY*C = LY * (I - BE)^(-1)
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....LISREL: LY*(I-BE)^(-1)");}
		/* Calculate the upper left quadrant: the covariance of the Ys */
		if(OMX_DEBUG) {Rprintf("Calculating Upper Left Quadrant of Expected Covariance Matrix.\n"); }
		if(OMX_DEBUG_ALGEBRA) {Rprintf("Copying C = PS.\n");}
		omxDGEMM(FALSE, FALSE, oned, D, PS, zerod, H); // H = D*PS = LY * (I - BE)^(-1) * PS
		omxCopyMatrix(Cov, TE); // Cov = TE
		omxDGEMM(FALSE, TRUE, oned, H, D, oned, Cov); // Cov = H*D^T + Cov = LY * (I - BE)^(-1) * PS * (LY * (I - BE)^(-1))^T + TE
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(J, "....LISREL: Upper Left Quadrant of Model-implied Covariance Matrix:");}
		if(Means != NULL) {
				omxCopyMatrix(Means, TY);
				omxDGEMV(FALSE, oned, D, AL, oned, Means);
		}
	}
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

unsigned short int omxNeedsUpdateLISRELObjective(omxObjective* oo) {
	return(omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->LX)
	 	|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->LY)
	 	|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->BE)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->GA)
	 	|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->PH)
	 	|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->PS)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->TD)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->TE)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->TH)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->TX)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->TY)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->KA)
		|| omxNeedsUpdate(((omxLISRELObjective*)oo->argStruct)->AL));
	
	// Note: cov is data, and should never need updating.
}


void omxInitLISRELObjective(omxObjective* oo, SEXP rObj) {
	
	if(OMX_DEBUG) { Rprintf("Initializing LISREL objective.\n"); }
	
	omxState* currentState = oo->matrix->currentState;	
	
	int nx, nxi, ny, neta, ntotal;
	
	SEXP slotValue;   //Used to get LISREL depth, for I-BE inverse speedup
	
	/* Create and register subobjective */
	omxObjective *subObjective = omxCreateSubObjective(oo);
	strncpy(subObjective->objType, "omxLISRELObjective", 19);
	omxLISRELObjective *LISobj = (omxLISRELObjective*) R_alloc(1, sizeof(omxLISRELObjective));
	
	/* Set Subobjective Calls and Structures */
	subObjective->objectiveFun = omxCallLISRELObjective;
	subObjective->needsUpdateFun = omxNeedsUpdateLISRELObjective;
	subObjective->destructFun = omxDestroyLISRELObjective;
	subObjective->setFinalReturns = NULL;
	subObjective->populateAttrFun = NULL; //omxPopulateLISRELAttributes;
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

	if(OMX_DEBUG) { Rprintf("Processing TX.\n"); }
	LISobj->TX = omxNewMatrixFromIndexSlot(rObj, currentState, "TX");

	if(OMX_DEBUG) { Rprintf("Processing TY.\n"); }
	LISobj->TY = omxNewMatrixFromIndexSlot(rObj, currentState, "TY");

	if(OMX_DEBUG) { Rprintf("Processing KA.\n"); }
	LISobj->KA = omxNewMatrixFromIndexSlot(rObj, currentState, "KA");

	if(OMX_DEBUG) { Rprintf("Processing AL.\n"); }
	LISobj->AL = omxNewMatrixFromIndexSlot(rObj, currentState, "AL");
	
	
	if(LISobj->LY == NULL) {
		LISobj->LY = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
		LISobj->PS = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
		LISobj->BE = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
		LISobj->TE = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
	}
	
	if(LISobj->LX == NULL) {
		LISobj->LX = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
		LISobj->PH = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
		LISobj->TD = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
	}
	
	if(LISobj->LY->cols == 0 || LISobj->LX->cols == 0) {
		LISobj->GA = omxInitMatrix(NULL, LISobj->LY->cols, LISobj->LX->cols, TRUE, currentState);
		LISobj->TH = omxInitMatrix(NULL, LISobj->LX->rows, LISobj->LY->rows, TRUE, currentState);
	}
	
	
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
	
	/* Identity Matrix, Size Of BE */
	if(OMX_DEBUG) { Rprintf("Generating I.\n"); }
	LISobj->I = omxNewIdentityMatrix(LISobj->BE->rows, currentState);
	omxRecompute(LISobj->I);
	
	
	/* Get the nilpotency index of the BE matrix for I-BE inverse speedup */
	if(OMX_DEBUG) { Rprintf("Processing expansion iteration depth.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("depth")));
	LISobj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { Rprintf("Using %d iterations.", LISobj->numIters); }
	UNPROTECT(1);
	
	
	/* Initialize the place holder matrices used in calculations */
	nx = LISobj->LX->rows;
	nxi = LISobj->LX->cols;
	ny = LISobj->LY->rows;
	neta = LISobj->LY->cols;
	ntotal = nx + ny;
	
	
	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }
	
	LISobj->A = 	omxInitMatrix(NULL, nx, nxi, TRUE, currentState);
	LISobj->B = 	omxInitMatrix(NULL, nx, nx, TRUE, currentState);
	LISobj->C = 	omxInitMatrix(NULL, neta, neta, TRUE, currentState);
	LISobj->D = 	omxInitMatrix(NULL, ny, neta, TRUE, currentState);
	LISobj->E = 	omxInitMatrix(NULL, nx, neta, TRUE, currentState);
	LISobj->F = 	omxInitMatrix(NULL, nx, ny, TRUE, currentState);
	LISobj->G = 	omxInitMatrix(NULL, neta, nxi, TRUE, currentState);
	LISobj->H = 	omxInitMatrix(NULL, ny, neta, TRUE, currentState);
	LISobj->J = 	omxInitMatrix(NULL, ny, ny, TRUE, currentState);
	LISobj->K = 	omxInitMatrix(NULL, neta, 1, TRUE, currentState);
	LISobj->L = 	omxInitMatrix(NULL, neta, neta, TRUE, currentState);
	LISobj->TOP = 	omxInitMatrix(NULL, ny, ntotal, TRUE, currentState);
	LISobj->BOT = 	omxInitMatrix(NULL, nx, ntotal, TRUE, currentState);
	LISobj->MUX = 	omxInitMatrix(NULL, nx, 1, TRUE, currentState);
	LISobj->MUY = 	omxInitMatrix(NULL, ny, 1, TRUE, currentState);
	
	
	LISobj->cov = 	omxInitMatrix(NULL, ntotal, ntotal, TRUE, currentState);

	LISobj->args = (omxMatrix**) R_alloc(2, sizeof(omxMatrix*));
	
	/* Means */
	if(LISobj->TX != NULL || LISobj->TY != NULL || LISobj->KA != NULL || LISobj->AL != NULL) {
		LISobj->means = 	omxInitMatrix(NULL, 1, ntotal, TRUE, currentState);
	} else LISobj->means  = 	NULL;
	//TODO: Adjust means processing to allow only Xs or only Ys
	

	/* Create parent objective */

	omxCreateMLObjective(oo, rObj, LISobj->cov, LISobj->means);
	
}

