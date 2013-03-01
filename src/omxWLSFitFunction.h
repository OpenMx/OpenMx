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
 */


#ifndef _OMX_WLS_FITFUNCTION_
#define _OMX_WLS_FITFUNCTION_ TRUE

typedef struct omxWLSFitFunction {

	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* observedFlattened;
	omxMatrix* expectedFlattened;
	omxMatrix* weights;
	omxMatrix* P;
	omxMatrix* B;

} omxWLSFitFunction;

void omxCreateWLSFitFunction(omxFitFunction* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means, omxMatrix* weights);

#endif /* _OMX_WLS_FITFUNCTION_ */
