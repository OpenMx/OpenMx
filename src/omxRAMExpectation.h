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
 
#ifndef _OMXRAMEXPECTATION_H_
#define _OMXRAMEXPECTATION_H_

typedef struct {

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *C, *EF, *ZM, *U, *V, *W, *X, *Y, *Z, *Ax, *P, *bCB, *ZSBC, *Mns, *lilI, *eCov, *beCov;
    omxMatrix *dA, *dS, *dM, *b, *D, *paramVec;
    // For gradients
    omxMatrix ***eqnOuts;
    omxMatrix **dAdts, **dSdts, **dMdts;

    // Each vector keeps track of the number
    // of non-zero elements in the 
    // corresponding dXdts matrix.
    int *dAdtsCount, *dSdtsCount, *dMdtsCount;

    int *dAdtsRowCache, *dAdtsColCache;
    int *dSdtsRowCache, *dSdtsColCache;
    int *dMdtsRowCache, *dMdtsColCache;

    omxMatrix *tempVec, *bigSum, *lilSum;

    int *pNums, nParam; // For Fast Gradient/Hessian Computation
	int numIters;
	double logDetObserved;
	double n;
	double *work;
	int lwork;

} omxRAMExpectation;

void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, 
    omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, 
    omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);

void omxInitRAMExpectation(omxExpectation* oo, SEXP rObj);
omxMatrix* omxGetRAMExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component);
void fastRAMGradientML(omxExpectation* oo, omxFitFunction* off, double* result);

#endif /* _OMXRAMEXPECTATION_H_ */
