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

	omxMatrix *cov, *means; // expected covariance and means
	omxMatrix *LX, *LY, *BE, *GA, *PH, *PS, *TD, *TE, *TH; // LISREL model Matrices
	//omxMatrix *TX, *TY, *KA, *AL; //LISREL Means Matrices
	omxMatrix *A, *B, *C, *D, *E, *F, *G, *H, *I, *J; // Place holder matrices used in computations
	omxMatrix *TOP, *BOT; // Place holder matrices for building covariance matrix from blocks
	//omxMatrix *C, *P, *V, *Mns; // Other Matrices, not sure what these are for.

	int numIters; // used by omxFastRAM/LISRELInverse
	double logDetObserved;
	double n;
	double* work; // used by omxFastRAM/LISRELInverse
	int lwork; // used by omxFastRAM/LISRELInverse

	int usePPML;
	omxData *ppmlData;
	omxMatrix *ppmlCov, *ppmlMeans;

	omxMatrix **args;

} omxLISRELObjective;

void omxCalculateLISRELCovarianceAndMeans(omxLISRELObjective* oro);

void omxInitLISRELObjective(omxObjective* oo, SEXP rObj);

/*
void omxFastLISRELInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I ); // same as RAM inverse
*/

#endif /* _OMXRAMOBJECTIVE_H_ */
