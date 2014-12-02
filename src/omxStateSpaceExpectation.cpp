/*
 *  Copyright 2007-2014 The OpenMx Project
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


void omxCallStateSpaceExpectation(omxExpectation* ox, const char *, const char *) {
    if(OMX_DEBUG) { mxLog("State Space Expectation Called."); }
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	omxRecompute(ose->A, NULL);
	omxRecompute(ose->B, NULL);
	omxRecompute(ose->C, NULL);
	omxRecompute(ose->D, NULL);
	omxRecompute(ose->Q, NULL);
	omxRecompute(ose->R, NULL);
	//omxRecompute(ose->x0); //Don't think I need to recompute these b.c. kalmanP/U functions do not use these.  They are recomputed on resetting these components.
	//omxRecompute(ose->P0);
	//omxRecompute(ose->u);
	
	// Probably should loop through all the data here!!!
	omxKalmanPredict(ose);
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
}


void omxPopulateSSMAttributes(omxExpectation *ox, SEXP algebra) {
    if(OMX_DEBUG) { mxLog("Populating State Space Attributes.  Currently this does very little!"); }
	
}




void omxKalmanPredict(omxStateSpaceExpectation* ose) {
	if(OMX_DEBUG) { mxLog("Kalman Predict Called."); }
	/* Creat local copies of State Space Matrices */
	omxMatrix* A = ose->A;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(A, "....State Space: A"); }
	omxMatrix* B = ose->B;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(B, "....State Space: B"); }
	//omxMatrix* C = ose->C;
	//omxMatrix* D = ose->D;
	omxMatrix* Q = ose->Q;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Q, "....State Space: Q"); }
	//omxMatrix* R = ose->R;
	//omxMatrix* r = ose->r;
	//omxMatrix* s = ose->s;
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	//omxMatrix* y = ose->y;
	omxMatrix* z = ose->z;
	//omxMatrix* K = ose->K;
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
	// omxMatrix* R = ose->R; //unused
	omxMatrix* r = ose->r;
	omxMatrix* s = ose->s;
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	omxMatrix* y = ose->y;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(y, "....State Space: y"); }
	// omxMatrix* K = ose->K; //unused
	omxMatrix* P = ose->P;
	//omxMatrix* S = ose->S;  //unused
	// omxMatrix* Y = ose->Y; //unused
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



void omxInitStateSpaceExpectation(omxExpectation* ox) {
	
	SEXP rObj = ox->rObj;
	if(OMX_DEBUG) { mxLog("Initializing State Space Expectation."); }
		
	int nx, ny, nu;
	
	//SEXP slotValue;   //Used by PPML
	
	/* Create and fill expectation */
	omxStateSpaceExpectation *SSMexp = (omxStateSpaceExpectation*) R_alloc(1, sizeof(omxStateSpaceExpectation));
	omxState* currentState = ox->currentState;
	
	/* Set Expectation Calls and Structures */
	ox->computeFun = omxCallStateSpaceExpectation;
	ox->destructFun = omxDestroyStateSpaceExpectation;
	ox->componentFun = omxGetStateSpaceExpectationComponent;
	ox->mutateFun = omxSetStateSpaceExpectationComponent;
	ox->populateAttrFun = omxPopulateSSMAttributes;
	ox->argStruct = (void*) SSMexp;
	
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
	}
}



