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
 
#ifndef _OMXLISRELEXPECTATION_H_
#define _OMXLISRELEXPECTATION_H_

typedef struct {

	omxMatrix *cov, *means; // expected covariance and means
	omxMatrix *LX, *LY, *BE, *GA, *PH, *PS, *TD, *TE, *TH; // LISREL model Matrices
	omxMatrix *TX, *TY, *KA, *AL; //LISREL Means Matrices
	omxMatrix *A, *B, *C, *D, *E, *F, *G, *H, *I, *J, *K, *L; // Place holder matrices used in computations.  Note: L is analogous to Ax in RAM and is used in I-BE inverse
	omxMatrix *TOP, *BOT; // Place holder matrices for building covariance matrix from blocks
	omxMatrix *MUX, *MUY; //Place holder matrices for building means from blocks
	//omxMatrix *C, *P, *V, *Mns; // Other Matrices, not sure what these are for.

	int numIters; // used by omxFastRAM/LISRELInverse
	double logDetObserved;
	double n;
	//double* work; // used by omxFastRAM/LISRELInverse
	//int lwork; // used by omxFastRAM/LISRELInverse

	omxMatrix **args;

	bool noLX;
	bool noLY;
	bool Lnocol;

} omxLISRELExpectation;

void omxCalculateLISRELCovarianceAndMeans(omxLISRELExpectation* oro);

void omxInitLISRELExpectation(omxExpectation* oo);

omxMatrix* omxGetLISRELExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component);

/*
void omxFastLISRELInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I ); // same as RAM inverse
*/

#endif /* _OMXLISRELEXPECTATION_H_ */
