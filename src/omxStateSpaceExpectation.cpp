/*
 *  Copyright 2007-2015 The OpenMx Project
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
 *
 */


/***********************************************************
*
*  omxStateSpaceExpectation.c
*
*  Created: Michael D. Hunter 	Date: 2012-10-28 20:07:36
*
*  Contains code to calculate the objective function for a
*   state space model.  Currently, this is done with a 
*   Kalman filter in separate Predict and Update steps.
*   Later this could be done with one of several Kalman 
*   filter-smoothers (a forward-backward algorithm).
*
**********************************************************/

#include "omxExpectation.h"
#include "omxBLAS.h"
#include "omxFIMLFitFunction.h"
#include "omxStateSpaceExpectation.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>


void omxCallStateSpaceExpectation(omxExpectation* ox, const char *, const char *) {
    if(OMX_DEBUG) { mxLog("State Space Expectation Called."); }
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	omxRecompute(ose->A, NULL);
	omxRecompute(ose->B, NULL);
	omxRecompute(ose->C, NULL);
	omxRecompute(ose->D, NULL);
	omxRecompute(ose->Q, NULL);
	omxRecompute(ose->R, NULL);
	
	// Probably should loop through all the data here!!!
	if(ose->t == NULL){
		omxKalmanPredict(ose);
	} else {
		omxKalmanBucyPredict(ose);
	}
	omxKalmanUpdate(ose);
}



void omxDestroyStateSpaceExpectation(omxExpectation* ox) {
	
	if(OMX_DEBUG) { mxLog("Destroying State Space Expectation."); }
	
	omxStateSpaceExpectation* argStruct = (omxStateSpaceExpectation*)(ox->argStruct);
	
	omxFreeMatrix(argStruct->r);
	omxFreeMatrix(argStruct->s);
	omxFreeMatrix(argStruct->z);
	omxFreeMatrix(argStruct->x);
	omxFreeMatrix(argStruct->y);
	omxFreeMatrix(argStruct->K);
	omxFreeMatrix(argStruct->P);
	omxFreeMatrix(argStruct->S);
	omxFreeMatrix(argStruct->Y);
	omxFreeMatrix(argStruct->Z);
	omxFreeMatrix(argStruct->det);
	omxFreeMatrix(argStruct->covInfo);
	omxFreeMatrix(argStruct->cov);
	omxFreeMatrix(argStruct->means);
	omxFreeMatrix(argStruct->smallC);
	omxFreeMatrix(argStruct->smallD);
	omxFreeMatrix(argStruct->smallR);
	omxFreeMatrix(argStruct->smallr);
	omxFreeMatrix(argStruct->smallK);
	omxFreeMatrix(argStruct->smallS);
	omxFreeMatrix(argStruct->smallY);
	
	delete argStruct;
	
}


void omxPopulateSSMAttributes(omxExpectation *ox, SEXP algebra) {
	if(OMX_DEBUG) { mxLog("Populating State Space Attributes.  Currently this does very little!"); }
	
	/* Initialize */
	omxSetExpectationComponent(ox, NULL, "Reset", NULL); //maybe shoulde be on ose?  after next line?
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	if( !(ose->returnScores) ){
		if(OMX_DEBUG) { mxLog("Not asking for attributes, this is being skipped!"); }
		return;
	}
	
	SEXP xpred, ypred, ppred, spred, xupda, pupda, xsmoo, psmoo;
	// yupda, supda
	
	omxRecompute(ose->A, NULL);
	omxRecompute(ose->B, NULL);
	omxRecompute(ose->C, NULL);
	omxRecompute(ose->D, NULL);
	omxRecompute(ose->Q, NULL);
	omxRecompute(ose->R, NULL);
	
	
	// allocate matrices to be returned
	if(OMX_DEBUG_ALGEBRA) { mxLog("Allocating initial population matrices ..."); }
	int nx = ose->C->cols;
	int ny = ose->C->rows;
	if(OMX_DEBUG_ALGEBRA) { mxLog("Find number of rows of data ..."); }
	//int nt = ox->data->dataMat->rows;
	int nt = ox->data->numObs;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... numObs:\n" << nt << std::endl; }
	if(OMX_DEBUG_ALGEBRA) { mxLog("Done Finding rows of data ..."); }
	int np = ((nx+1)*nx)/2;
	int ns = ((ny+1)*ny)/2;
	Rf_protect(xpred = Rf_allocMatrix(REALSXP, nt+1, nx));
	Rf_protect(ypred = Rf_allocMatrix(REALSXP, nt+1, ny));
	Rf_protect(ppred = Rf_allocMatrix(REALSXP, nt+1, np));
	Rf_protect(spred = Rf_allocMatrix(REALSXP, nt+1, ns));
	Rf_protect(xupda = Rf_allocMatrix(REALSXP, nt+1, nx));
	//Rf_protect(yupda = Rf_allocMatrix(REALSXP, nt+1, ny));
	Rf_protect(pupda = Rf_allocMatrix(REALSXP, nt+1, np));
	//Rf_protect(supda = Rf_allocMatrix(REALSXP, nt+1, ns));
	Rf_protect(xsmoo = Rf_allocMatrix(REALSXP, nt+1, nx));
	Rf_protect(psmoo = Rf_allocMatrix(REALSXP, nt+1, np));
	
	
	if(OMX_DEBUG_ALGEBRA) { mxLog("Setting zeroth row ..."); }
	// Set first row of xpred to x0
	int row = 0;
	for(int col = 0; col < nx; col++){
		REAL(xpred)[col * (nt+1) + row] =
			omxMatrixElement(ose->x, col, 0);
		REAL(xupda)[col * (nt+1) + row] =
			omxMatrixElement(ose->x, col, 0);
	}
	
	// Set first row of ppred to vech(P0)
	int counter = 0;
	for(int i = 0; i < ose->P->cols; i++) {
		for(int j = i; j < ose->P->rows; j++) {
			REAL(ppred)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
			REAL(pupda)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
			counter++;
		}
	}
	
	Eigen::VectorXd oldDefs;
	oldDefs.resize(ox->data->defVars.size());
	oldDefs.setConstant(NA_REAL);
	
	// Probably should loop through all the data here!!!
	if(OMX_DEBUG_ALGEBRA) { mxLog("Beginning forward loop ..."); }
	for(row=1; row < (nt+1); row++){
		if(OMX_DEBUG_ALGEBRA) { mxLog("Setting first data row ..."); }
		// Set row of data
		for(int i = 0; i < ny; i++) {
			omxSetMatrixElement(ose->y, i, 0, omxDoubleDataElement(ox->data, row-1, i));
		}
		
		// handle definition variables
		int numVarsFilled = 0;
		numVarsFilled = ox->data->handleDefinitionVarList(ox->currentState, row-1, oldDefs.data());
		
		/* Run Kalman prediction */
		if(ose->t == NULL){
			omxKalmanPredict(ose);
		} else {
			omxKalmanBucyPredict(ose);
		}
		
		// Copy latent state
		for(int col = 0; col < nx; col++)
			REAL(xpred)[col * (nt+1) + row] =
				omxMatrixElement(ose->x, col, 0);
		
		// Copy latent cov
		counter = 0;
		for(int i = 0; i < ose->P->cols; i++) {
			for(int j = i; j < ose->P->rows; j++) {
				REAL(ppred)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
				counter++;
			}
		}
		
		// Create Full observed cov prediction
		if(OMX_DEBUG_ALGEBRA) { mxLog("Hand prediction of full observed cov ..."); }
		omxDSYMM(FALSE, 1.0, ose->P, ose->C, 0.0, ose->Y); // Y = C P
		omxCopyMatrix(ose->S, ose->R); // S = R
		omxDGEMM(FALSE, TRUE, 1.0, ose->Y, ose->C, 1.0, ose->S); // S = Y C^T + S THAT IS C P C^T + R
		
		// Copy observed cov
		counter = 0;
		for(int i = 0; i < ose->S->cols; i++) {
			for(int j = i; j < ose->S->rows; j++) {
				REAL(spred)[counter * (nt+1) + row] = omxMatrixElement(ose->S, j, i);
				counter++;
			}
		}
		
		/* Run Kalman update */
		omxKalmanUpdate(ose);
		
		// Copy latent state
		for(int col = 0; col < nx; col++)
			REAL(xupda)[col * (nt+1) + row] =
				omxMatrixElement(ose->x, col, 0);
		
		// Copy latent cov
		counter = 0;
		for(int i = 0; i < ose->P->cols; i++) {
			for(int j = i; j < ose->P->rows; j++) {
				REAL(pupda)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
				counter++;
			}
		}
		
		// Copy observed means prediction
		for(int col = 0; col < ny; col++)
			REAL(ypred)[col * (nt+1) + row] =
				omxMatrixElement(ose->s, col, 0);
		
		// TODO Add m2ll calculation here.
		// Probably like this
		/*m2ll = y^T S y */ // n.b. y originally is the data row but becomes the data residual!
		//omxDSYMV(1.0, S, y, 0.0, s); // s = S y
		//m2ll = omxDDOT(y, s); // m2ll = y s THAT IS y^T S y
		//m2ll += det; // m2ll = m2ll + det THAT IS m2ll = log(det(S)) + y^T S y
		// Note: this leaves off the S->cols * log(2*pi) THAT IS k*log(2*pi)
	}
	
	/* TODO Add Backward pass through data for Kalman smoother*/
	// Initialize end of smoothed latents to last updated latents
	// P = last updated P
	// x = last updated x
	
	// Copy x and P as smoothed estimates for export
	row = nt;
	// Copy latent state
	for(int col = 0; col < nx; col++)
		REAL(xsmoo)[col * (nt+1) + row] =
			omxMatrixElement(ose->x, col, 0);
	
	// Copy latent cov
	counter = 0;
	for(int i = 0; i < ose->P->cols; i++) {
		for(int j = i; j < ose->P->rows; j++) {
			REAL(psmoo)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
			counter++;
		}
	}
	
	// loop backwars through all data
	if(OMX_DEBUG_ALGEBRA) { mxLog("Beginning backward loop ..."); }
	for(row = nt-1; row > -1; row--){
		// Copy Z = updated P from pupda
		counter = 0;
		for(int i = 0; i < nx; i++) {
			for(int j = i; j < nx; j++) {
				double next = REAL(pupda)[counter * (nt+1) + row];
				omxSetMatrixElement(ose->Z, i, j, next);
				omxSetMatrixElement(ose->Z, j, i, next);
				counter++;
			}
		}
		
		// Copy z = updated x from xupda
		for(int col = 0; col < nx; col++)
			omxSetMatrixElement(ose->z, col, 0, REAL(xupda)[col * (nt+1) + row]);
		
		// Copy eigenExpA = predicted P from later row of ppred
		counter = 0;
		for(int i = 0; i < nx; i++) {
			for(int j = i; j < nx; j++) {
				double next = REAL(ppred)[counter * (nt+1) + row + 1];
				ose->eigenExpA(i, j) = next;
				ose->eigenExpA(j, i) = next;
				counter++;
			}
		}
		
		// Copy eigenPreX = predicted x from later row of xpred
		for(int col = 0; col < nx; col++)
			ose->eigenPreX(col, 0) =  REAL(xpred)[col * (nt+1) + row + 1];
		
		// Run Smoother
		omxRauchTungStriebelSmooth(ose);
		
		// Copy latent state to xsmoo
		for(int col = 0; col < nx; col++)
			REAL(xsmoo)[col * (nt+1) + row] =
				omxMatrixElement(ose->x, col, 0);
		
		// Copy latent cov to psmoo
		counter = 0;
		for(int i = 0; i < ose->P->cols; i++) {
			for(int j = i; j < ose->P->rows; j++) {
				REAL(psmoo)[counter * (nt+1) + row] = omxMatrixElement(ose->P, j, i);
				counter++;
			}
		}
	}
	
	
	// TODO check on definition variable population
	//  I suspect this does yet work properly for def vars.
	
	Rf_setAttrib(algebra, Rf_install("xPredicted"), xpred);
	Rf_setAttrib(algebra, Rf_install("yPredicted"), ypred);
	Rf_setAttrib(algebra, Rf_install("PPredicted"), ppred);
	Rf_setAttrib(algebra, Rf_install("SPredicted"), spred);
	Rf_setAttrib(algebra, Rf_install("xUpdated"), xupda);
	Rf_setAttrib(algebra, Rf_install("PUpdated"), pupda);
	Rf_setAttrib(algebra, Rf_install("xSmoothed"), xsmoo);
	Rf_setAttrib(algebra, Rf_install("PSmoothed"), psmoo);
	
	
	/*
	omxMatrix *expCovInt, *expMeanInt;
	expCovInt = argStruct->cov;
	expMeanInt = argStruct->means;
	
	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL && expMeanInt->rows > 0  && expMeanInt->cols > 0) {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, 0, 0));		
	}

	Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
	
	if(argStruct->populateRowDiagnostics){
		omxMatrix *rowLikelihoodsInt = argStruct->rowLikelihoods;
		Rf_protect(rowLikelihoodsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		for(int row = 0; row < rowLikelihoodsInt->rows; row++)
			REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);
		Rf_setAttrib(algebra, Rf_install("likelihoods"), rowLikelihoodsExt);
	}
	*/
	
}




void omxKalmanPredict(omxStateSpaceExpectation* ose) {
	if(OMX_DEBUG) { mxLog("Kalman Predict Called."); }
	/* Creat local copies of State Space Matrices */
	omxMatrix* A = ose->A;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(A, "....State Space: A"); }
	omxMatrix* B = ose->B;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(B, "....State Space: B"); }
	omxMatrix* Q = ose->Q;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Q, "....State Space: Q"); }
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	omxMatrix* z = ose->z;
	omxMatrix* P = ose->P;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P"); }
	//omxMatrix* S = ose->S;
	//omxMatrix* Y = ose->Y;
	omxMatrix* Z = ose->Z;

	/* x = A x + B u */
	omxDGEMV(FALSE, 1.0, A, x, 0.0, z); // x = A x
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(z, "....State Space: z = A x"); }
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(A, "....State Space: A"); }
	omxDGEMV(FALSE, 1.0, B, u, 1.0, z); // x = B u + x THAT IS x = A x + B u
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(z, "....State Space: z = A x + B u"); }
	omxCopyMatrix(x, z); // x = z THAT IS x = A x + B u
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x = A x + B u"); }
	
	/* P = A P A^T + Q */
	omxDSYMM(FALSE, 1.0, P, A, 0.0, Z); // Z = A P
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Z, "....State Space: Z = A P"); }
	omxCopyMatrix(P, Q); // P = Q
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P = Q"); }
	omxDGEMM(FALSE, TRUE, 1.0, Z, A, 1.0, P); // P = Z A^T + P THAT IS P = A P A^T + Q
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P = A P A^T + Q"); }
}


void omxKalmanUpdate(omxStateSpaceExpectation* ose) {
	//TODO: Clean up this hack of function.
	if(OMX_DEBUG) { mxLog("Kalman Update Called."); }
	/* Creat local copies of State Space Matrices */
	omxMatrix* C = ose->C;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(C, "....State Space: C"); }
	omxMatrix* D = ose->D;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....State Space: D"); }
	omxMatrix* r = ose->r;
	omxMatrix* s = ose->s;
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	omxMatrix* y = ose->y;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(y, "....State Space: y"); }
	omxMatrix* P = ose->P;
	// omxMatrix* Cov = ose->cov; //unused
	omxMatrix* Means = ose->means;
	omxMatrix* Det = ose->det;
	*Det->data = 0.0; // the value pointed to by Det->data is assigned to be zero
	omxMatrix* smallC = ose->smallC;
	omxMatrix* smallD = ose->smallD;
	omxMatrix* smallr = ose->smallr;
	omxMatrix* R = ose->R;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(R, "....State Space: R on entry"); }
	omxMatrix* smallR = ose->smallR;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallR, "....State Space: small R on entry"); }
	omxMatrix* smallK = ose->smallK;
	omxMatrix* smallS = ose->smallS;
	omxMatrix* smallY = ose->smallY;
	int ny = y->rows;
	int nx = x->rows;
	Eigen::VectorXi toRemoveSS(ny);
	int numRemovesSS = 0;
	Eigen::VectorXi toRemoveNoneLat(nx);
	toRemoveNoneLat.setZero();
	Eigen::VectorXi toRemoveNoneOne(1);
	toRemoveNoneOne.setZero();
	
	omxMatrix* covInfo = ose->covInfo;
	int info = 0; // Used for computing inverse for Kalman gain
	
	omxCopyMatrix(smallS, ose->S);
	
	/* Reset/Resample aliased matrices */
	omxCopyMatrix(smallC, ose->C);
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallC, "....State Space: C (Reset)"); }
	omxCopyMatrix(smallD, ose->D);
	omxCopyMatrix(smallR, ose->R);
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallR, "....State Space: small R (Reset)"); }
	omxCopyMatrix(smallr, ose->r);
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallr, "....State Space: r (Reset)"); }
	omxCopyMatrix(smallK, ose->K);
	omxCopyMatrix(smallS, ose->S);
	omxCopyMatrix(smallY, ose->Y);
	
	/* r = r - C x - D u */
	/* Alternatively, create just the expected value for the data row, x. */
	omxDGEMV(FALSE, 1.0, smallC, x, 0.0, s); // s = C x
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(s, "....State Space: s = C x"); }
	omxDGEMV(FALSE, 1.0, D, u, 1.0, s); // s = D u + s THAT IS s = C x + D u
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(s, "....State Space: s = C x + D u"); }
	omxCopyMatrix(Means, s); // Means = s THAT IS Means = C x + D u
	//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....State Space: Means"); }
	omxTransposeMatrix(Means);
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....State Space: Means"); }
	
	//If entire data vector, y, is missing, then set residual, r, to zero.
	//otherwise, compute residual.
	toRemoveSS.setZero();
	for(int j = 0; j < y->rows; j++) {
		double dataValue = omxMatrixElement(y, j, 0);
		int dataValuefpclass = std::fpclassify(dataValue);
		if(dataValuefpclass == FP_NAN || dataValuefpclass == FP_INFINITE) {
			numRemovesSS++;
			toRemoveSS[j] = 1;
			omxSetMatrixElement(r, j, 0, 0.0);
		} else {
			omxSetMatrixElement(r, j, 0, (dataValue -  omxMatrixElement(s, j, 0)));
		}
	}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(r, "....State Space: Residual (Loop)"); }
	/* Now compute the residual */
	//omxCopyMatrix(r, y); // r = y
	//omxDAXPY(-1.0, s, r); // r = r - s THAT IS r = y - (C x + D u)
	omxCopyMatrix(smallr, ose->r);
	
	/* Filter S Here */
	// N.B. if y is completely missing or completely present, leave S alone.
	// Otherwise, filter S.
	if(numRemovesSS < ny && numRemovesSS > 0) {
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: S"); }
		if(OMX_DEBUG) { mxLog("Filtering S, R, C, r, K, and Y."); }
		omxRemoveRowsAndColumns(smallS, numRemovesSS, numRemovesSS, toRemoveSS.data(), toRemoveSS.data());
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: S (Filtered)"); }
		
		omxRemoveRowsAndColumns(smallR, numRemovesSS, numRemovesSS, toRemoveSS.data(), toRemoveSS.data());
		
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallC, "....State Space: C"); }
		omxRemoveRowsAndColumns(smallC, numRemovesSS, 0, toRemoveSS.data(), toRemoveNoneLat.data());
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallC, "....State Space: C (Filtered)"); }
		
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(r, "....State Space: r"); }
		omxRemoveRowsAndColumns(smallr, numRemovesSS, 0, toRemoveSS.data(), toRemoveNoneOne.data());
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallr, "....State Space: r (Filtered)"); }
		
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallK, "....State Space: K"); }
		omxRemoveRowsAndColumns(smallK, numRemovesSS, 0, toRemoveSS.data(), toRemoveNoneLat.data());
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallK, "....State Space: K (Filtered)"); }
		
		omxRemoveRowsAndColumns(smallY, numRemovesSS, 0, toRemoveSS.data(), toRemoveNoneLat.data());
	}
	if(numRemovesSS == ny) {
		if(OMX_DEBUG_ALGEBRA) { mxLog("Completely missing row of data found."); }
		if(OMX_DEBUG_ALGEBRA) { mxLog("Skipping much of Kalman Update."); }
		return ;
	}
	
	
	/* S = C P C^T + R */
	omxDSYMM(FALSE, 1.0, P, smallC, 0.0, smallY); // Y = C P
	//omxCopyMatrix(S, smallR); // S = R
	memcpy(smallS->data, smallR->data, smallR->rows * smallR->cols * sizeof(double)); // Less safe omxCopyMatrix that keeps smallS aliased to S.
	omxDGEMM(FALSE, TRUE, 1.0, smallY, smallC, 1.0, smallS); // S = Y C^T + S THAT IS C P C^T + R
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: S = C P C^T + R"); }
	
	
	/* Now compute the Kalman Gain and update the error covariance matrix */
	/* S = S^-1 */
	omxDPOTRF(smallS, &info); // S replaced by the lower triangular matrix of the Cholesky factorization
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: Cholesky of S"); }
	covInfo->data[0] = (double) info;
	for(int i = 0; i < smallS->cols; i++) {
		*Det->data += log(fabs(omxMatrixElement(smallS, i, i)));
	}
	//det *= 2.0; //sum( log( abs( diag( chol(S) ) ) ) )*2
	omxDPOTRI(smallS, &info); // S = S^-1 via Cholesky factorization
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: Inverse of S"); }
	// If Cholesky of exp cov failed (i.e. non-positive def), Populate 1,1 element of smallS (inverse of exp cov) with NA_REAL
	if(covInfo->data[0] > 0) {
		omxSetMatrixElement(smallS, 0, 0, NA_REAL);
	}
	
	/* K = P C^T S^-1 */
	/* Computed as K^T = S^-1 C P */
	omxDSYMM(TRUE, 1.0, smallS, smallY, 0.0, smallK); // K = Y^T S THAT IS K = P C^T S^-1
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallK, "....State Space: K^T = S^-1 C P"); }
	
	/* x = x + K r */
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallr, "....State Space Check Residual: r"); }
	omxDGEMV(TRUE, 1.0, smallK, smallr, 1.0, x); // x = K r + x
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x = K r + x"); }
	
	/* P = (I - K C) P */
	/* P = P - K C P */
	omxDGEMM(TRUE, FALSE, -1.0, smallK, smallY, 1.0, P); // P = -K Y + P THAT IS P = P - K C P
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P = P - K C P"); }
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(smallS, "....State Space: Inverse of S"); }
	
	
	/*m2ll = y^T S y */ // n.b. y originally is the data row but becomes the data residual!
	//omxDSYMV(1.0, S, y, 0.0, s); // s = S y
	//m2ll = omxDDOT(y, s); // m2ll = y s THAT IS y^T S y
	//m2ll += det; // m2ll = m2ll + det THAT IS m2ll = log(det(S)) + y^T S y
	// Note: this leaves off the S->cols * log(2*pi) THAT IS k*log(2*pi)
}


void omxKalmanBucyPredict(omxStateSpaceExpectation* ose) {
	if(OMX_DEBUG) { mxLog("Kalman Bucy Predict Called."); }
	/* Creat local copies of State Space Matrices */
	omxMatrix* A = ose->A;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(A, "....State Space: A"); }
	omxMatrix* B = ose->B;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(B, "....State Space: B"); }
	omxMatrix* Q = ose->Q;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Q, "....State Space: Q"); }
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	omxMatrix* P = ose->P;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P"); }
	omxMatrix* Z = ose->Z;
	omxMatrix* t = ose->t;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(t, "....State Space: t"); }
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space oldT on entrance:\n" << ose->oldT << std::endl; }
	ose->deltaT = omxMatrixElement(t, 0, 0) - ose->oldT;
	ose->oldT = omxMatrixElement(t, 0, 0);
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space oldT on exit:\n" << ose->oldT << std::endl; }
	double deltaT = ose->deltaT;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space deltaT:\n" << deltaT << std::endl; }
	
	//EigenMatrixAdaptor eigenA(A);
	// intializes eigenA as an instance of EigenMatrixAdaptor class, initialized to omxMatrix A
	// EigenVectorAdaptor
	// eigenB = eigenA.exp(); // matrix exponential
	// for scalar multiplication
	// A * omxElement(B,1,4); A eigen matrix, B omxMatrix
	// Subtract from diagonal 1.0
	// A.diagonal() -= 1
	
	/* Eigen Matrix reference setting */
	Eigen::MatrixXd &eigenExpA = ose->eigenExpA;
	Eigen::MatrixXd &eigenIA = ose->eigenIA;
	Eigen::MatrixXd &PSI = ose->PSI;
	Eigen::MatrixXd &IP = ose->IP;
	Eigen::MatrixXd &I = ose->I;
	
	
	/*R code for the next few lines
		ldim <- nrow(x$A)
		I <- diag(1, nrow=ldim)
		expA <- as.matrix(expm(x$A * x$deltaT))
		intA <- solve(x$A) %*% (expA - I)
		x$x <- expA %*% x$x + intA %*% x$B %*% x$u
	*/
	
	
	/*  Z = A
		eigenA = Z
		eigenExpA = eigenA*deltaT  THAT IS eigenExpA = A*deltaT
	*/
	omxCopyMatrix(Z, A);
	EigenMatrixAdaptor eigenA(Z);
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenA:\n" << eigenA << std::endl; }
	eigenExpA = eigenA * deltaT;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenA:\n" << eigenA << std::endl; }
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space deltaT:\n" << deltaT << std::endl; }
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenExpA:\n" << eigenExpA << std::endl; }
	
	/* eigenExpA = expm(eigenExpA)  THAT IS eigenExpA = expm(A*deltaT) */
	eigenExpA = eigenExpA.exp();
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenExpA:\n" << eigenExpA << std::endl; }
	
	/* eigenIA = eigenExpA - I*/
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space I:\n" << I << std::endl; }
	eigenIA = eigenExpA - I;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space expA - I:\n" << eigenIA << std::endl; }
	eigenIA = eigenA.lu().solve(eigenIA);
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space A^-1 (expA - I):\n" << eigenIA << std::endl; }
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenA:\n" << eigenA << std::endl; }
	/* eigenIA = A^-1 IA  THAT IS eigenIA = A^-1 (expm(A*deltaT) - I) */
	
	/* x = expm(A*deltaT) * x + IA * B * u */
	EigenMatrixAdaptor eigenx(x); //or vector
	EigenMatrixAdaptor eigenu(u); //or vector
	EigenMatrixAdaptor eigenB(B);
	eigenx.derived() = eigenExpA * eigenx + eigenIA * eigenB * eigenu;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space XPred:\n" << eigenx << std::endl; }
	
	
	/* SUMMARY */
	// EA = expm(A*deltaT)
	// IA = A^-1 (EA - I)
	// x = EA x + IA B u
	
	
	/*R code for the next few lines
		ldim <- nrow(x$A)
		psi <- cbind(rbind(x$A, matrix(0, nrow=ldim, ncol=ldim)), rbind(x$Q, -x$A))
		epsi <- as.matrix(expm(psi * x$deltaT)) %*% rbind(x$P, I)
		x$P <- epsi[1:ldim, ] %*% solve(epsi[(ldim+1):(2*ldim), ])
	*/
	
	
	/* PSI = block matrix
	-A^T  0
	 Q    A
	*/
	EigenMatrixAdaptor eigenQ(Q);
	PSI << -1.0*eigenA.transpose(), Eigen::MatrixXd::Zero(A->rows, A->rows), eigenQ, eigenA;
	PSI = PSI * deltaT;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space PSI:\n" << PSI << std::endl; }
	/* PSI = PSI * deltaT */
	/* PSI = expm(PSI) */
	PSI = PSI.exp();
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space expPSI*deltaT:\n" << PSI << std::endl; }
	
	/* IP = block matrix
	I
	P
	*/
	EigenMatrixAdaptor eigenP(P);
	IP << I, eigenP;
	//IP << eigenP, Eigen::MatrixXd::Identity(rows, rows)
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space IP:\n" << IP << std::endl; }
	
	/* IP = PSI IP */
	IP = PSI * IP;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space Blocks:\n" << IP << std::endl; }
	
	/* eigenIA = block 1 of IP; eigenExpA = block 2 of IP */
	eigenP.derived() = IP.block(0, 0, A->rows, A->cols);
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space Block1:\n" << eigenP << std::endl; }
	eigenExpA = IP.block(A->rows, 0, A->rows, A->cols);
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space Block2:\n" << eigenExpA << std::endl; }
	eigenP.derived() = eigenP.transpose().lu().solve(eigenExpA.transpose());
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenP:\n" << eigenP << std::endl; }
	/* P = B1 B2^-1  THAT IS  solve B2^T P = B1^T for P */
	
	
	/* SUMMARY */
	// PSI = 	-A^T  0
	//			 Q    A
	// PSI = expm(PSI * deltaT)
	// IP = 	I
	// 			P
	// IP = PSI IP
	// P = IP[1] IP[2]^-1
	
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenX" << eigenx << std::endl; }
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space eigenP:\n" << eigenP << std::endl; }
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P"); }
}


void omxRauchTungStriebelSmooth(omxStateSpaceExpectation* ose) {
	if(OMX_DEBUG) { mxLog("Rauch Tung Striebel Smooth Called."); }
	/* Creat local copies of State Space Matrices */
	omxMatrix* A = ose->A;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(A, "....State Space: A"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x smoothed"); }
	omxMatrix* z = ose->z;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(z, "....State Space: x updated"); }
	omxMatrix* P = ose->P;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P smoothed"); }
	omxMatrix* Z = ose->Z;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Z, "....State Space: P updated"); }
	
	/* Eigen Matrix reference setting */
	Eigen::MatrixXd &eigenPreX = ose->eigenPreX;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space x predicted:\n" << eigenPreX << std::endl; }
	Eigen::MatrixXd &eigenExpA = ose->eigenExpA;
	if(OMX_DEBUG_ALGEBRA) {std::cout << "... State Space P predicted:\n" << eigenExpA << std::endl; }
	Eigen::MatrixXd &eigenIA = ose->eigenIA; // Storage for Sg (RTS smoother gain matrix)
	
	/* Eigen Adaptor copies (not really copies, more like wrappers) */
	EigenMatrixAdaptor eigenA(A);
	EigenMatrixAdaptor eigenx(x);
	EigenMatrixAdaptor eigenz(z);
	EigenMatrixAdaptor eigenP(P);
	EigenMatrixAdaptor eigenZ(Z);
	
	
	/* Create the RTS Gain matrix*/
	// Sg = Pui * A * Ppi+1 ^-1
	// eigenIA = Z * A * eigenExpA^-1
	// Possible typo above, A should be A^T
	eigenIA = eigenExpA.lu().solve( eigenA * eigenZ ).transpose();
	// try also
	//eigenIA = eigenExpA.ldlt().solve( eigenA.transpose() * eigenZ ).transpose();
	// with #include <Eigen/Cholesky>
	
	/* Smooth the latent state */
	// xsi = xui + Sg * (xsi+1 - xpi+1)
	// x = z + eigenIA * (x - eigenPreX)
	eigenx.derived() = eigenz + eigenIA * (eigenx - eigenPreX);
	
	/* Smooth the latent covariance */
	// Psi = Pui + Sg * (Psi+1 - Ppi+1) * Sg^T
	// P = Z + eigenIA * (P - eigenExpA) * eigenIA^T
	eigenP.derived() = eigenZ + eigenIA * (eigenP - eigenExpA) * eigenIA.transpose();
	
	// TODO add eigenPreX to struct and initialize
}


void omxInitStateSpaceExpectation(omxExpectation* ox) {
	
	SEXP rObj = ox->rObj;
	if(OMX_DEBUG) { mxLog("Initializing State Space Expectation."); }
		
	int nx, ny, nu;
	
	//SEXP slotValue;   //Used by PPML
	
	/* Create and fill expectation */
	//omxStateSpaceExpectation *SSMexp = (omxStateSpaceExpectation*) R_alloc(1, sizeof(omxStateSpaceExpectation));
	omxStateSpaceExpectation *SSMexp = new omxStateSpaceExpectation;
	
	omxState* currentState = ox->currentState;
	
	/* Set Expectation Calls and Structures */
	ox->computeFun = omxCallStateSpaceExpectation;
	ox->destructFun = omxDestroyStateSpaceExpectation;
	ox->componentFun = omxGetStateSpaceExpectationComponent;
	ox->mutateFun = omxSetStateSpaceExpectationComponent;
	ox->populateAttrFun = omxPopulateSSMAttributes;
	ox->argStruct = (void*) SSMexp;
	ox->canDuplicate = false;
	
	/* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Initializing State Space Meta Data for expectation."); }
	
	if(OMX_DEBUG) { mxLog("Processing A."); }
	SSMexp->A = omxNewMatrixFromSlot(rObj, currentState, "A");
	
	if(OMX_DEBUG) { mxLog("Processing B."); }
	SSMexp->B = omxNewMatrixFromSlot(rObj, currentState, "B");
	
	if(OMX_DEBUG) { mxLog("Processing C."); }
	SSMexp->C = omxNewMatrixFromSlot(rObj, currentState, "C");
	
	if(OMX_DEBUG) { mxLog("Processing D."); }
	SSMexp->D = omxNewMatrixFromSlot(rObj, currentState, "D");
	
	if(OMX_DEBUG) { mxLog("Processing Q."); }
	SSMexp->Q = omxNewMatrixFromSlot(rObj, currentState, "Q");
	
	if(OMX_DEBUG) { mxLog("Processing R."); }
	SSMexp->R = omxNewMatrixFromSlot(rObj, currentState, "R");
	
	if(OMX_DEBUG) { mxLog("Processing initial x."); }
	SSMexp->x0 = omxNewMatrixFromSlot(rObj, currentState, "x0");
	
	if(OMX_DEBUG) { mxLog("Processing initial P."); }
	SSMexp->P0 = omxNewMatrixFromSlot(rObj, currentState, "P0");
	
	if(OMX_DEBUG) { mxLog("Processing u."); }
	SSMexp->u = omxNewMatrixFromSlot(rObj, currentState, "u");
	
	if(OMX_DEBUG) { mxLog("Processing t."); }
	SSMexp->t = omxNewMatrixFromSlot(rObj, currentState, "t");
	
	
	/* Initialize the place holder matrices used in calculations */
	nx = SSMexp->C->cols;
	ny = SSMexp->C->rows;
	nu = SSMexp->D->cols;
	
	if(OMX_DEBUG) { mxLog("Processing first data row for y."); }
	SSMexp->y = omxInitMatrix(ny, 1, TRUE, currentState);
	for(int i = 0; i < ny; i++) {
		omxSetMatrixElement(SSMexp->y, i, 0, omxDoubleDataElement(ox->data, 0, i));
	}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(SSMexp->y, "....State Space: y"); }
	
	// TODO Make x0 and P0 static (if possible) to save memory
	// TODO Look into omxMatrix.c/h for a possible new matrix from omxMatrix function
	if(OMX_DEBUG) { mxLog("Generating static internals for resetting initial values."); }
	SSMexp->x = 	omxInitMatrix(nx, 1, TRUE, currentState);
	SSMexp->P = 	omxInitMatrix(nx, nx, TRUE, currentState);
	omxCopyMatrix(SSMexp->x, SSMexp->x0);
	omxCopyMatrix(SSMexp->P, SSMexp->P0);
	
	if(OMX_DEBUG) { mxLog("Generating internals for computation."); }
	
	SSMexp->covInfo = 	omxInitMatrix(1, 1, TRUE, currentState);
	SSMexp->det = 	omxInitMatrix(1, 1, TRUE, currentState);
	SSMexp->r = 	omxInitMatrix(ny, 1, TRUE, currentState);
	SSMexp->s = 	omxInitMatrix(ny, 1, TRUE, currentState);
	SSMexp->z = 	omxInitMatrix(nx, 1, TRUE, currentState);
	SSMexp->K = 	omxInitMatrix(ny, nx, TRUE, currentState); // Actually the tranpose of the Kalman gain
	SSMexp->S = 	omxInitMatrix(ny, ny, TRUE, currentState);
	SSMexp->Y = 	omxInitMatrix(ny, nx, TRUE, currentState);
	SSMexp->Z = 	omxInitMatrix(nx, nx, TRUE, currentState);
	
	SSMexp->cov = 		omxInitMatrix(ny, ny, TRUE, currentState);
	SSMexp->means = 	omxInitMatrix(1, ny, TRUE, currentState);
	
	/* Create alias matrices for missing data filtering */
	SSMexp->smallC = 	omxInitMatrix(ny, nx, TRUE, currentState);
	SSMexp->smallD = 	omxInitMatrix(ny, nu, TRUE, currentState);
	SSMexp->smallR = 	omxInitMatrix(ny, ny, TRUE, currentState);
	SSMexp->smallr = 	omxInitMatrix(ny, 1, TRUE, currentState);
	SSMexp->smallK = 	omxInitMatrix(ny, nx, TRUE, currentState);
	SSMexp->smallS = 	omxInitMatrix(ny, ny, TRUE, currentState);
	SSMexp->smallY = 	omxInitMatrix(ny, nx, TRUE, currentState);
	
	//SSMexp->deltaT = 	omxInitMatrix(1, 1, TRUE, currentState);
	SSMexp->oldT = 0.0;
	SSMexp->deltaT = 0.0;
	
	/* Eigen Matrix initialization */
	SSMexp->eigenExpA.resize(nx, nx);
	SSMexp->I.resize(nx, nx);
	SSMexp->I = Eigen::MatrixXd::Identity(nx, nx);
	SSMexp->eigenIA.resize(nx, nx);
	SSMexp->PSI.resize(2*nx, 2*nx);
	SSMexp->IP.resize(2*nx, nx);
	SSMexp->eigenPreX.resize(nx, 1);
	
	/* Population of Kalman scores*/
	if(OMX_DEBUG) {
		mxLog("Accessing Kalman score population option.");
	}
	SSMexp->returnScores = Rf_asInteger(R_do_slot(rObj, Rf_install("scores")));
	
	
	omxCopyMatrix(SSMexp->smallC, SSMexp->C);
	omxCopyMatrix(SSMexp->smallD, SSMexp->D);
	omxCopyMatrix(SSMexp->smallR, SSMexp->R);
	omxCopyMatrix(SSMexp->smallr, SSMexp->r);
	omxCopyMatrix(SSMexp->smallK, SSMexp->K);
	omxCopyMatrix(SSMexp->smallS, SSMexp->S);
	omxCopyMatrix(SSMexp->smallY, SSMexp->Y);

}


omxMatrix* omxGetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component) {
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = ose->cov;
	} else if(strEQ("means", component)) {
		retval = ose->means;
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	} else if(strEQ("inverse", component)) {
		retval = ose->smallS;
	} else if(strEQ("determinant", component)) {
		retval = ose->det;
	} else if(strEQ("r", component)) {
		retval = ose->r;
	} else if(strEQ("covInfo", component)) {
		retval = ose->covInfo;
	}
	
	return retval;
}

void omxSetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component, omxMatrix* om) {
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	if(!strcmp("y", component)) {
		for(int i = 0; i < ose->y->rows; i++) {
			omxSetMatrixElement(ose->y, i, 0, omxVectorElement(om, i));
		}
		//ose->y = om;
	}
	if(!strcmp("Reset", component)) {
		omxRecompute(ose->x0, NULL);
		omxRecompute(ose->P0, NULL);
		omxCopyMatrix(ose->x, ose->x0);
		omxCopyMatrix(ose->P, ose->P0);
		if(ose->t != NULL){
			ose->oldT = 0.0;
		}
	}
}



