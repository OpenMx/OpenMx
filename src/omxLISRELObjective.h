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
 
#ifndef _OMXLISRELOBJECTIVE_H_
#define _OMXLISRELOBJECTIVE_H_

typedef struct {

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *LX, *LY, *BE, *GA, *PH, *PS, *TD, *TE, *TH, *TX, *TY, *KA, *AL, *I; // LISREL model Matrices
	omxMatrix *C, *X, *Y, *Z, *Ax, *P, *V, *Mns; // Other Matrices used in computations.  NOTE: THESE HAVE NOT BEEN UPDATED YET

	int numIters;
	double logDetObserved;
	double n;
	double* work;
	int lwork;

	int usePPML;
	omxData *ppmlData;
	omxMatrix *ppmlCov, *ppmlMeans;

} omxLISRELObjective;

void omxCalculateLISRELCovarianceAndMeans(omxMatrix* LX, omxMatrix* LY, omxMatrix* BE, 
    omxMatrix* GA, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, 
    omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);

void omxInitLISRELObjective(omxObjective* oo, SEXP rObj);

/*
void omxFastLISRELInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I ); // same as RAM inverse
*/

#endif /* _OMXRAMOBJECTIVE_H_ */
