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
 
#ifndef _OMXSTATESPACEEXPECTATION_H_
#define _OMXSTATESPACEEXPECTATION_H_

typedef struct {
	omxMatrix *cov, *means;
	omxMatrix *A, *B, *C, *D, *Q, *R, *t; // State Space model Matrices
	omxMatrix *r, *s, *u, *x, *y, *z; // Data and place holder vectors
	omxMatrix *K, *P, *S, *Y, *Z; // Behind the scenes state space matrices (P, S, and K) and place holder matrices
	omxMatrix *x0, *P0; // Placeholders for initial state and initial Rf_error cov
	omxMatrix *det; // Determinant of expected covariance matrix S
	omxMatrix *smallC, *smallD, *smallr, *smallR, *smallK, *smallS, *smallY; //aliases of C, D, r, R, K, and S for missing data handling
	omxMatrix *covInfo; //info from Cholesky decomp of small expected cov (smallS) to be passed to FIML single iteration
	double oldT;
	double deltaT;
	int returnScores; // Whether or not to populated the row-by-row expected means, covariances, Kalman scores, and likelihoods back to R
	
	/* Eigen Matrix initialization */
	Eigen::MatrixXd eigenExpA;
	Eigen::MatrixXd I;
	Eigen::MatrixXd eigenIA;
	Eigen::MatrixXd PSI;
	Eigen::MatrixXd IP;
	
} omxStateSpaceExpectation;


void omxKalmanPredict(omxStateSpaceExpectation* ose);
void omxKalmanUpdate(omxStateSpaceExpectation* ose);
void omxKalmanBucyPredict(omxStateSpaceExpectation* ose);

void omxInitStateSpaceExpectation(omxExpectation* ox);

omxMatrix* omxGetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component);

void omxSetStateSpaceExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component, omxMatrix* om);


#endif /* _OMXSTATESPAACEEXPECTATION_H_ */
