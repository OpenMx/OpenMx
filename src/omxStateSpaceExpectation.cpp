/*
 *  Copyright 2007-2013 The OpenMx Project
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


void omxCallStateSpaceExpectation(omxExpectation* ox, enum ComputeExpectationContext ign) {
    if(OMX_DEBUG) { mxLog("State Space Expectation Called."); }
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	omxRecompute(ose->A);
	omxRecompute(ose->B);
	omxRecompute(ose->C);
	omxRecompute(ose->D);
	omxRecompute(ose->Q);
	omxRecompute(ose->R);
	
	// Probably should loop through all the data here!!!
	omxKalmanPredict(ose);
	omxKalmanUpdate(ose);
}



void omxDestroyStateSpaceExpectation(omxExpectation* ox) {
	
	if(OMX_DEBUG) { mxLog("Destroying State Space Expectation."); }
	
	omxStateSpaceExpectation* argStruct = (omxStateSpaceExpectation*)(ox->argStruct);
	
	/* We allocated 'em, so we destroy 'em. */
	omxFreeMatrixData(argStruct->r);
	omxFreeMatrixData(argStruct->s);
	omxFreeMatrixData(argStruct->z);
	//omxFreeMatrixData(argStruct->u); // This is data, destroy it?
	//omxFreeMatrixData(argStruct->x); // This is latent data, destroy it?
	omxFreeMatrixData(argStruct->y); // This is data, destroy it?
	omxFreeMatrixData(argStruct->K); // This is the Kalman gain, destroy it?
	//omxFreeMatrixData(argStruct->P); // This is latent cov, destroy it?
	omxFreeMatrixData(argStruct->S); // This is data error cov, destroy it?
	omxFreeMatrixData(argStruct->Y);
	omxFreeMatrixData(argStruct->Z);
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
    if(OMX_DEBUG) { mxLog("Kalman Update Called."); }
	/* Creat local copies of State Space Matrices */
	//omxMatrix* A = ose->A;
	//omxMatrix* B = ose->B;
	omxMatrix* C = ose->C;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(C, "....State Space: C"); }
	omxMatrix* D = ose->D;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(D, "....State Space: D"); }
	//omxMatrix* Q = ose->Q;
	omxMatrix* R = ose->R;
	omxMatrix* r = ose->r;
	omxMatrix* s = ose->s;
	omxMatrix* u = ose->u;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(u, "....State Space: u"); }
	omxMatrix* x = ose->x;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x"); }
	omxMatrix* y = ose->y;
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(y, "....State Space: y"); }
	omxMatrix* K = ose->K;
	omxMatrix* P = ose->P;
	omxMatrix* S = ose->S;
	omxMatrix* Y = ose->Y;
	//omxMatrix* Z = ose->Z;
	omxMatrix* Cov = ose->cov;
	omxMatrix* Means = ose->means;
	omxMatrix* Det = ose->det;
	*Det->data = 0.0; // the value pointed to by Det->data is assigned to be zero
	
	int info = 0; // Used for computing inverse for Kalman gain
	
	/* r = r - C x - D u */
	/* Alternatively, create just the expected value for the data row, x. */
	omxDGEMV(FALSE, 1.0, C, x, 0.0, s); // s = C x
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(s, "....State Space: s = C x"); }
	omxDGEMV(FALSE, 1.0, D, u, 1.0, s); // s = D u + s THAT IS s = C x + D u
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(s, "....State Space: s = C x + D u"); }
	omxCopyMatrix(Means, s); // Means = s THAT IS Means = C x + D u
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....State Space: Means"); }
	omxTransposeMatrix(Means);
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....State Space: Means"); }
	/* Now compute the residual */
	omxCopyMatrix(r, y); // r = y
	omxDAXPY(-1.0, s, r); // r = r - s THAT IS r = y - (C x + D u)
	
	/* S = C P C^T + R */
	omxDSYMM(FALSE, 1.0, P, C, 0.0, Y); // Y = C P
	omxCopyMatrix(S, R); // S = R
	omxDGEMM(FALSE, TRUE, 1.0, Y, C, 1.0, S); // S = Y C^T + S THAT IS C P C^T + R
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(S, "....State Space: S = C P C^T + R"); }
	
	omxCopyMatrix(Cov, S); //Note: I know this is inefficient memory use, but for now it is more clear.-MDH
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Cov, "....State Space: Cov"); }
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....State Space: Means"); }
	
	/* Now compute the Kalman Gain and update the error covariance matrix */
	/* S = S^-1 */
	omxDPOTRF(S, &info); // S replaced by the lower triangular matrix of the Cholesky factorization
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(S, "....State Space: Cholesky of S"); }
	for(int i = 0; i < S->cols; i++) {
		*Det->data += log(fabs(omxMatrixElement(S, i, i)));
	}
	//det *= 2.0; //sum( log( abs( diag( chol(S) ) ) ) )*2
	omxDPOTRI(S, &info); // S = S^-1 via Cholesky factorization
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(S, "....State Space: Inverse of S"); }
	
	/* K = P C^T S^-1 */
	/* Computed as K^T = S^-1 C P */
	omxDSYMM(TRUE, 1.0, S, Y, 0.0, K); // K = Y^T S THAT IS K = P C^T S^-1
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(K, "....State Space: K^T = S^-1 C P"); }
	
	/* x = x + K r */
	omxDGEMV(TRUE, 1.0, K, r, 1.0, x); // x = K r + x
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(x, "....State Space: x = K r + x"); }
	
	/* P = (I - K C) P */
	/* P = P - K C P */
	omxDGEMM(TRUE, FALSE, -1.0, K, Y, 1.0, P); // P = -K Y + P THAT IS P = P - K C P
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(P, "....State Space: P = P - K C P"); }
	
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
	SSMexp->x = omxNewMatrixFromSlot(rObj, currentState, "x");
	
	if(OMX_DEBUG) { mxLog("Processing initial P."); }
	SSMexp->P = omxNewMatrixFromSlot(rObj, currentState, "P");
	
	
	/* Initialize the place holder matrices used in calculations */
	nx = SSMexp->C->cols;
	ny = SSMexp->C->rows;
	nu = SSMexp->D->cols;
	
	if(OMX_DEBUG) { mxLog("Processing first data row for y."); }
	SSMexp->y = omxInitMatrix(NULL, ny, 1, TRUE, currentState);
	for(int i = 0; i < ny; i++) {
		omxSetMatrixElement(SSMexp->y, i, 0, omxMatrixElement(ox->data->dataMat, 0, i));
	}
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(SSMexp->y, "....State Space: y"); }
	
	// TODO Make x0 and P0 static (if possible) to save memory
	// TODO Look into omxMatrix.c/h for a possible new matrix from omxMatrix function
	if(OMX_DEBUG) { mxLog("Generating static internals for resetting initial values."); }
	SSMexp->x0 = 	omxInitMatrix(NULL, nx, 1, TRUE, currentState);
	SSMexp->P0 = 	omxInitMatrix(NULL, nx, nx, TRUE, currentState);
	omxCopyMatrix(SSMexp->x0, SSMexp->x);
	omxCopyMatrix(SSMexp->P0, SSMexp->P);
	
	if(OMX_DEBUG) { mxLog("Generating internals for computation."); }
	
	SSMexp->det = 	omxInitMatrix(NULL, 1, 1, TRUE, currentState);
	SSMexp->r = 	omxInitMatrix(NULL, ny, 1, TRUE, currentState);
	SSMexp->s = 	omxInitMatrix(NULL, ny, 1, TRUE, currentState);
	SSMexp->u = 	omxInitMatrix(NULL, nu, 1, TRUE, currentState);
	SSMexp->z = 	omxInitMatrix(NULL, nx, 1, TRUE, currentState);
	SSMexp->K = 	omxInitMatrix(NULL, ny, nx, TRUE, currentState); // Actually the tranpose of the Kalman gain
	SSMexp->S = 	omxInitMatrix(NULL, ny, ny, TRUE, currentState);
	SSMexp->Y = 	omxInitMatrix(NULL, ny, nx, TRUE, currentState);
	SSMexp->Z = 	omxInitMatrix(NULL, nx, nx, TRUE, currentState);
	
	SSMexp->cov = 		omxInitMatrix(NULL, ny, ny, TRUE, currentState);
	SSMexp->means = 	omxInitMatrix(NULL, 1, nx, TRUE, currentState);
}


omxMatrix* omxGetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component) {
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(!strncmp("cov", component, 3)) {
		retval = ose->cov;
	} else if(!strncmp("mean", component, 4)) {
		retval = ose->means;
	} else if(!strncmp("pvec", component, 4)) {
		// Once implemented, change compute function and return pvec
	} else if(!strncmp("inverse", component, 7)) {
		retval = ose->S;
	} else if(!strncmp("determinant", component, 11)) {
		retval = ose->det;
	}
	
	if(OMX_DEBUG) { mxLog("Returning 0x%x.", retval); }

	return retval;
}

void omxSetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component, omxMatrix* om) {
	omxStateSpaceExpectation* ose = (omxStateSpaceExpectation*)(ox->argStruct);
	
	if(!strcmp("y", component)) {
		ose->y = om;
	}
	if(!strcmp("Reset", component)) {
		omxCopyMatrix(ose->x, ose->x0);
		omxCopyMatrix(ose->P, ose->P0);
	}
}



