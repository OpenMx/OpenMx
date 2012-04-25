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
 */


#ifndef _OMX_ML_OBJECTIVE_
#define _OMX_ML_OBJECTIVE_ TRUE

typedef struct omxMLObjective {

	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* localCov;
	omxMatrix* localProd;
	omxMatrix* P;
	omxMatrix* C;
	omxMatrix* I;
    
	double n;
	double logDetObserved;

	double* work;
	int lwork;
	
    // Subobjective Storage;
    omxMatrix** dSigma;         // dSigma/dTheta
    omxMatrix** dMu;            // dMu/dTheta
    omxMatrix* Mu;
    omxMatrix* Ms;
    omxMatrix* X;
    omxMatrix* Y;
    
    // TODO: Add space for second derivatives
    
    // Allows the subobjective function to calculate its own derivates in an optimal way.
    void (*derivativeFun)(omxObjective*, omxMatrix**, omxMatrix**, int*);

} omxMLObjective;

void omxCreateMLObjective(omxObjective* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means);

void omxSetMLObjectiveGradient(omxObjective* oo, void (*)(omxObjective*, double*));

void omxSetMLObjectiveGradientComponents(omxObjective* oo, void (*)(omxObjective*, omxMatrix**, omxMatrix**, int*));

void omxCalculateMLGradient(omxObjective* oo, double* gradient);

#endif /* _OMX_ML_OBJECTIVE_ */
