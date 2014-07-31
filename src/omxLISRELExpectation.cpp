/*
 *  Copyright 2007-2014 The OpenMx Project
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

#include "omxExpectation.h"
#include "omxBLAS.h"
#include "omxFIMLFitFunction.h"
#include "omxLISRELExpectation.h"

extern void omxCreateMLFitFunction(omxFitFunction* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means);
// TODO: Merge ML and FIML Fit Functions into one unit.

void omxCallLISRELExpectation(omxExpectation* oo, const char *, const char *) {
    if(OMX_DEBUG) { mxLog("LISREL Expectation Called."); }
	omxLISRELExpectation* oro = (omxLISRELExpectation*)(oo->argStruct);
	
	omxRecompute(oro->LX, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->LY, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->BE, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->GA, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->PH, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->PS, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->TD, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->TE, FF_COMPUTE_FIT, NULL);
	omxRecompute(oro->TH, FF_COMPUTE_FIT, NULL);
	if(oro->TX != NULL) {     // Update means?
		omxRecompute(oro->TX, FF_COMPUTE_FIT, NULL);
		omxRecompute(oro->KA, FF_COMPUTE_FIT, NULL);
	}
	if(oro->TY != NULL) {
		omxRecompute(oro->TY, FF_COMPUTE_FIT, NULL);
		omxRecompute(oro->AL, FF_COMPUTE_FIT, NULL);
	}
	
	omxCalculateLISRELCovarianceAndMeans(oro);
}

void omxDestroyLISRELExpectation(omxExpectation* oo) {

	if(OMX_DEBUG) { mxLog("Destroying LISREL Expectation."); }
	
	omxLISRELExpectation* argStruct = (omxLISRELExpectation*)(oo->argStruct);

	omxFreeMatrix(argStruct->cov);
	omxFreeMatrix(argStruct->means);
	
	omxFreeMatrix(argStruct->A);
	omxFreeMatrix(argStruct->B);
	omxFreeMatrix(argStruct->C);
	omxFreeMatrix(argStruct->D);
	omxFreeMatrix(argStruct->E);
	omxFreeMatrix(argStruct->F);
	omxFreeMatrix(argStruct->G);
	omxFreeMatrix(argStruct->H);
	omxFreeMatrix(argStruct->I);
	omxFreeMatrix(argStruct->J);
	omxFreeMatrix(argStruct->K);
	omxFreeMatrix(argStruct->L);
	omxFreeMatrix(argStruct->TOP);
	omxFreeMatrix(argStruct->BOT);
	omxFreeMatrix(argStruct->MUX);
	omxFreeMatrix(argStruct->MUY);
	
	if(argStruct->Lnocol) {
		omxFreeMatrix(argStruct->GA);
		omxFreeMatrix(argStruct->TH);
	}
	
	if(argStruct->noLY) {
		omxFreeMatrix(argStruct->LY);
		omxFreeMatrix(argStruct->PS);
		omxFreeMatrix(argStruct->BE);
		omxFreeMatrix(argStruct->TE);
	}
	
	if(argStruct->noLX) {
		omxFreeMatrix(argStruct->LX);
		omxFreeMatrix(argStruct->PH);
		omxFreeMatrix(argStruct->TD);
	}
}

void omxPopulateLISRELAttributes(omxExpectation *oo, SEXP algebra)
{
    Rf_setAttrib(algebra, Rf_install("numStats"), Rf_ScalarReal(omxDataDF(oo->data)));

	/*
	omxLISRELExpectation* oro = (omxLISRELExpectation*) (oo->argStruct);
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
	*/ //This block of code works fine but because I do not use any of it later, it throws a huge block of Rf_warnings about unused variables.
	// In general, I do not yet understand what this function needs to do.
	
	/*
	omxShallowInverse(numIters, BE, Z, Ax, I ); // Z = (I-A)^-1
	
	if(OMX_DEBUG_ALGEBRA) { mxLog("....DGEMM: %c %c \n %d %d %d \n %f \n %x %d %x %d \n %f %x %d.", *(Z->majority), *(S->majority), (Z->rows), (S->cols), (S->rows), oned, Z->data, (Z->leading), S->data, (S->leading), zerod, Ax->data, (Ax->leading));}
	// F77_CALL(omxunsafedgemm)(Z->majority, S->majority, &(Z->rows), &(S->cols), &(S->rows), &oned, Z->data, &(Z->leading), S->data, &(S->leading), &zerod, Ax->data, &(Ax->leading)); 	// X = ZS
	omxDGEMM(FALSE, FALSE, oned, Z, S, zerod, Ax);

	if(OMX_DEBUG_ALGEBRA) { mxLog("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.", *(Ax->majority), *(Z->minority), (Ax->rows), (Z->rows), (Z->cols), oned, X->data, (X->leading), Z->data, (Z->lagging), zerod, Ax->data, (Ax->leading));}
	// F77_CALL(omxunsafedgemm)(Ax->majority, Z->minority, &(Ax->rows), &(Z->rows), &(Z->cols), &oned, Ax->data, &(Ax->leading), Z->data, &(Z->leading), &zerod, Ax->data, &(Ax->leading)); 
	omxDGEMM(FALSE, TRUE, oned, Ax, Z, zerod, Ax);
	// Ax = ZSZ' = Covariance matrix including latent variables
	
	SEXP expCovExt;
	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, Ax->rows, Ax->cols));
	for(int row = 0; row < Ax->rows; row++)
		for(int col = 0; col < Ax->cols; col++)
			REAL(expCovExt)[col * Ax->rows + row] =
				omxMatrixElement(Ax, row, col);
	setAttrib(algebra, Rf_install("UnfilteredExpCov"), expCovExt);
	*/
}

/* omxFastLISRELInverse would go here */


void omxCalculateLISRELCovarianceAndMeans(omxLISRELExpectation* oro) {
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
	if(OMX_DEBUG) { mxLog("Running LISREL computation in omxCalculateLISRELCovarianceAndMeans."); }
	double oned = 1.0, zerod=0.0; //, minusOned = -1.0;
	//int ipiv[BE->rows], lwork = 4 * BE->rows * BE->cols; //This is copied from omxShallowInverse()
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
		if(OMX_DEBUG) {mxLog("Calculating Lower Right Quadrant of Expected Covariance Matrix."); }
		omxDGEMM(FALSE, FALSE, oned, LX, PH, zerod, A); // A = LX*PH
		omxCopyMatrix(B, TD); // B = TD
		omxDGEMM(FALSE, TRUE, oned, A, LX, oned, B);  // B = LX*PH*LX^T + TD
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(B, "....LISREL: Lower Right Quadrant of Model-implied Covariance Matrix:");}
		
		/* Calculate (I-BE)^(-1) and LY*(I-BE)^(-1) */
		if(OMX_DEBUG) {mxLog("Calculating Inverse of I-BE."); }
		omxShallowInverse(NULL, numIters, BE, C, L, I ); // C = (I-BE)^-1
		//omxCopyMatrix(C, BE); // C = BE
		//omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, C); // C = I - BE
		//omxDGETRF(C, ipiv); //LU Decomp
		//omxDGETRI(C, ipiv, work, lwork); //Inverse based on LU Decomp ... C = C^(-1) = (I - BE)^(-1)
		
		
		omxDGEMM(FALSE, FALSE, oned, LY, C, zerod, D); // D = LY*C = LY * (I - BE)^(-1)
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....LISREL: LY*(I-BE)^(-1)");}
		
		/* Calculate the lower left quadrant: the covariance of Xs and Ys, nX by nY */
		if(OMX_DEBUG) {mxLog("Calculating Lower Left Quadrant of Expected Covariance Matrix."); }
		omxDGEMM(FALSE, TRUE, oned, A, GA, zerod, E); // E = A*GA^T = LX*PH*GA^T
		omxCopyMatrix(F, TH); // F = TH
		omxDGEMM(FALSE, TRUE, oned, E, D, oned, F); // F = E*D^T + F = LX*PH*GA^T * (LY * (I - BE)^(-1))^T + TH
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(F, "....LISREL: Lower Left Quadrant of Model-implied Covariance Matrix:");}
	
		
	/* Calculate the upper right quadrant: NOTE THIS IS MERELY THE LOWER LEFT QUADRANT TRANSPOSED. */
	//DONE as omxTranspose(F)
		
		
		/* Calculate the upper left quadrant: the covariance of the Ys */
		if(OMX_DEBUG) {mxLog("Calculating Upper Left Quadrant of Expected Covariance Matrix."); }
		if(OMX_DEBUG_ALGEBRA) {mxLog("Calculating G = GA*PH.");}
		omxDGEMM(FALSE, FALSE, oned, GA, PH, zerod, G); // G = GA*PH
		if(OMX_DEBUG_ALGEBRA) {mxLog("Copying C = PS.");}
		omxCopyMatrix(C, PS); // C = PS
		if(OMX_DEBUG_ALGEBRA) {mxLog("Calculating G = GA*PH.");}
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
		if(OMX_DEBUG) {mxLog("Calculating Lower Right Quadrant of Expected Covariance Matrix."); }
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
		if(OMX_DEBUG) {mxLog("Calculating Inverse of I-BE."); }
		omxShallowInverse(NULL, numIters, BE, C, L, I ); // C = (I-BE)^-1
		//omxCopyMatrix(C, BE); // C = BE
		//omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, C); // C = I - BE
		//omxDGETRF(C, ipiv); //LU Decomp
		//omxDGETRI(C, ipiv, work, lwork); //Inverse based on LU Decomp ... C = C^(-1) = (I - BE)^(-1)
		omxDGEMM(FALSE, FALSE, oned, LY, C, zerod, D); // D = LY*C = LY * (I - BE)^(-1)
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....LISREL: LY*(I-BE)^(-1)");}
		/* Calculate the upper left quadrant: the covariance of the Ys */
		if(OMX_DEBUG) {mxLog("Calculating Upper Left Quadrant of Expected Covariance Matrix."); }
		if(OMX_DEBUG_ALGEBRA) {mxLog("Copying C = PS.");}
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
	if(OMX_DEBUG) { mxLog("Running RAM computation."); }
		
	double oned = 1.0, zerod=0.0;
	
	if(Ax == NULL || I == NULL || Z == NULL || Y == NULL || X == NULL) {
		Rf_error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
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
	// 		Rf_error("INTERNAL ERROR: Incorrectly sized matrices being passed to omxRAMObjective Calculation.\n Please report this to the OpenMx development team.");
	// }
	
	omxShallowInverse(numIters, A, Z, Ax, I );
		
	// IMPORTANT: Cov = FZSZ'F'
	if(OMX_DEBUG_ALGEBRA) { mxLog("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.", *(F->majority), *(Z->majority), (F->rows), (Z->cols), (Z->rows), oned, F->data, (F->leading), Z->data, (Z->leading), zerod, Y->data, (Y->leading));}
	// F77_CALL(omxunsafedgemm)(F->majority, Z->majority, &(F->rows), &(Z->cols), &(Z->rows), &oned, F->data, &(F->leading), Z->data, &(Z->leading), &zerod, Y->data, &(Y->leading)); 	// Y = FZ
	omxDGEMM(FALSE, FALSE, 1.0, F, Z, 0.0, Y);

	if(OMX_DEBUG_ALGEBRA) { mxLog("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.", *(Y->majority), *(S->majority), (Y->rows), (S->cols), (S->rows), oned, Y->data, (Y->leading), S->data, (S->leading), zerod, X->data, (X->leading));}
	// F77_CALL(omxunsafedgemm)(Y->majority, S->majority, &(Y->rows), &(S->cols), &(S->rows), &oned, Y->data, &(Y->leading), S->data, &(S->leading), &zerod, X->data, &(X->leading)); 	// X = FZS
	omxDGEMM(FALSE, FALSE, 1.0, Y, S, 0.0, X);

	if(OMX_DEBUG_ALGEBRA) { mxLog("....DGEMM: %c %c %d %d %d %f %x %d %x %d %f %x %d.", *(X->majority), *(Y->minority), (X->rows), (Y->rows), (Y->cols), oned, X->data, (X->leading), Y->data, (Y->lagging), zerod, Cov->data, (Cov->leading));}
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

void omxInitLISRELExpectation(omxExpectation* oo) {
	SEXP rObj = oo->rObj;
	
	if(OMX_DEBUG) { mxLog("Initializing LISREL Expectation."); }
		
	int nx, nxi, ny, neta, ntotal;
	
	SEXP slotValue;
	
	/* Create and fill expectation */
	omxLISRELExpectation *LISobj = (omxLISRELExpectation*) R_alloc(1, sizeof(omxLISRELExpectation));
	omxState* currentState = oo->currentState;
	
	/* Set Expectation Calls and Structures */
	oo->computeFun = omxCallLISRELExpectation;
	oo->destructFun = omxDestroyLISRELExpectation;
	oo->componentFun = omxGetLISRELExpectationComponent;
	oo->populateAttrFun = omxPopulateLISRELAttributes;
	oo->argStruct = (void*) LISobj;
	
	/* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Initializing LISREL Meta Data for expectation."); }
	
	if(OMX_DEBUG) { mxLog("Processing LX."); }
	LISobj->LX = omxNewMatrixFromSlot(rObj, currentState, "LX");
	
	if(OMX_DEBUG) { mxLog("Processing LY."); }
	LISobj->LY = omxNewMatrixFromSlot(rObj, currentState, "LY");
	
	if(OMX_DEBUG) { mxLog("Processing BE."); }
	LISobj->BE = omxNewMatrixFromSlot(rObj, currentState, "BE");
	
	if(OMX_DEBUG) { mxLog("Processing GA."); }
	LISobj->GA = omxNewMatrixFromSlot(rObj, currentState, "GA");
	
	if(OMX_DEBUG) { mxLog("Processing PH."); }
	LISobj->PH = omxNewMatrixFromSlot(rObj, currentState, "PH");
	
	if(OMX_DEBUG) { mxLog("Processing PS."); }
	LISobj->PS = omxNewMatrixFromSlot(rObj, currentState, "PS");
	
	if(OMX_DEBUG) { mxLog("Processing TD."); }
	LISobj->TD = omxNewMatrixFromSlot(rObj, currentState, "TD");
	
	if(OMX_DEBUG) { mxLog("Processing TE."); }
	LISobj->TE = omxNewMatrixFromSlot(rObj, currentState, "TE");
	
	if(OMX_DEBUG) { mxLog("Processing TH."); }
	LISobj->TH = omxNewMatrixFromSlot(rObj, currentState, "TH");

	if(OMX_DEBUG) { mxLog("Processing TX."); }
	LISobj->TX = omxNewMatrixFromSlot(rObj, currentState, "TX");

	if(OMX_DEBUG) { mxLog("Processing TY."); }
	LISobj->TY = omxNewMatrixFromSlot(rObj, currentState, "TY");

	if(OMX_DEBUG) { mxLog("Processing KA."); }
	LISobj->KA = omxNewMatrixFromSlot(rObj, currentState, "KA");

	if(OMX_DEBUG) { mxLog("Processing AL."); }
	LISobj->AL = omxNewMatrixFromSlot(rObj, currentState, "AL");
	
	LISobj->noLY = LISobj->LY == NULL;
	if(LISobj->noLY) {
		LISobj->LY = omxInitMatrix(0, 0, TRUE, currentState);
		LISobj->PS = omxInitMatrix(0, 0, TRUE, currentState);
		LISobj->BE = omxInitMatrix(0, 0, TRUE, currentState);
		LISobj->TE = omxInitMatrix(0, 0, TRUE, currentState);
	}
	
	LISobj->noLX = LISobj->LX == NULL;
	if(LISobj->noLX) {
		LISobj->LX = omxInitMatrix(0, 0, TRUE, currentState);
		LISobj->PH = omxInitMatrix(0, 0, TRUE, currentState);
		LISobj->TD = omxInitMatrix(0, 0, TRUE, currentState);
	}
	
	LISobj->Lnocol = LISobj->LY->cols == 0 || LISobj->LX->cols == 0;
	if(LISobj->Lnocol) {
		LISobj->GA = omxInitMatrix(LISobj->LY->cols, LISobj->LX->cols, TRUE, currentState);
		LISobj->TH = omxInitMatrix(LISobj->LX->rows, LISobj->LY->rows, TRUE, currentState);
	}
	
	
	/* Identity Matrix, Size Of BE */
	if(OMX_DEBUG) { mxLog("Generating I."); }
	LISobj->I = omxNewIdentityMatrix(LISobj->BE->rows, currentState);
	
	
	/* Get the nilpotency index of the BE matrix for I-BE inverse speedup */
	if(OMX_DEBUG) { mxLog("Processing expansion iteration depth."); }
	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("depth")));
	LISobj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { mxLog("Using %d iterations.", LISobj->numIters); }
	}
	
	/* Initialize the place holder matrices used in calculations */
	nx = LISobj->LX->rows;
	nxi = LISobj->LX->cols;
	ny = LISobj->LY->rows;
	neta = LISobj->LY->cols;
	ntotal = nx + ny;
	
	
	if(OMX_DEBUG) { mxLog("Generating internals for computation."); }
	
	LISobj->A = 	omxInitMatrix(nx, nxi, TRUE, currentState);
	LISobj->B = 	omxInitMatrix(nx, nx, TRUE, currentState);
	LISobj->C = 	omxInitMatrix(neta, neta, TRUE, currentState);
	LISobj->D = 	omxInitMatrix(ny, neta, TRUE, currentState);
	LISobj->E = 	omxInitMatrix(nx, neta, TRUE, currentState);
	LISobj->F = 	omxInitMatrix(nx, ny, TRUE, currentState);
	LISobj->G = 	omxInitMatrix(neta, nxi, TRUE, currentState);
	LISobj->H = 	omxInitMatrix(ny, neta, TRUE, currentState);
	LISobj->J = 	omxInitMatrix(ny, ny, TRUE, currentState);
	LISobj->K = 	omxInitMatrix(neta, 1, TRUE, currentState);
	LISobj->L = 	omxInitMatrix(neta, neta, TRUE, currentState);
	LISobj->TOP = 	omxInitMatrix(ny, ntotal, TRUE, currentState);
	LISobj->BOT = 	omxInitMatrix(nx, ntotal, TRUE, currentState);
	LISobj->MUX = 	omxInitMatrix(nx, 1, TRUE, currentState);
	LISobj->MUY = 	omxInitMatrix(ny, 1, TRUE, currentState);
	
	
	LISobj->cov = 	omxInitMatrix(ntotal, ntotal, TRUE, currentState);

	LISobj->args = (omxMatrix**) R_alloc(2, sizeof(omxMatrix*));
	
	/* Means */
	if(LISobj->TX != NULL || LISobj->TY != NULL || LISobj->KA != NULL || LISobj->AL != NULL) {
		LISobj->means = 	omxInitMatrix(1, ntotal, TRUE, currentState);
	} else LISobj->means  = 	NULL;
	//TODO: Adjust means processing to allow only Xs or only Ys

}

omxMatrix* omxGetLISRELExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component) {
	omxLISRELExpectation* ore = (omxLISRELExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = ore->cov;
	} else if(strEQ("means", component)) {
		retval = ore->means;
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}
	
	return retval;
}
