/*
 *  Copyright 2007-2012 The OpenMx Project
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
 
#ifndef _OMXRAMOBJECTIVE_H_
#define _OMXRAMOBJECTIVE_H_

typedef struct {

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *C, *EF, *ZM, *U, *V, *W, *X, *Y, *Z, *Ax, *P, *Mns;
    omxMatrix *dA, *dS, *dM;

	int numIters;
	double logDetObserved;
	double n;
	double* work;
	int lwork;

	int usePPML;
	omxData *ppmlData;
	omxMatrix *ppmlCov, *ppmlMeans;

} omxRAMObjective;

void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, 
    omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, 
    omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);

void omxInitRAMObjective(omxObjective* oo, SEXP rObj);

void omxFastRAMInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I );

#endif /* _OMXRAMOBJECTIVE_H_ */
