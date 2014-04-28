/*
 *  Copyright 2007-2014 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxBLAS.h"
#include "omxDefines.h"
#include "omxNormalExpectation.h"

void omxComputeNormalExpectation(omxExpectation* ox, const char *, const char *) {
    if(OMX_DEBUG) { mxLog("Normal Expectation calculating."); }

	omxNormalExpectation* one = (omxNormalExpectation*) (ox->argStruct);

	omxRecompute(one->cov);
	if(one->means != NULL)
	    omxRecompute(one->means);
}

void omxDestroyNormalExpectation(omxExpectation* ox) {

	if(OMX_DEBUG) { mxLog("Destroying Normal Expectation."); }

}

void omxPopulateNormalAttributes(omxExpectation *ox, SEXP algebra) {
    if(OMX_DEBUG) { mxLog("Populating Normal Attributes."); }

	omxNormalExpectation* one = (omxNormalExpectation*) (ox->argStruct);
    
	omxMatrix *cov = one->cov;
	omxMatrix *means = one->means;

	SEXP expMeanExt, expCovExt;

    omxRecompute(cov);
	if(means != NULL)
    	omxRecompute(means);

	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, cov->rows, cov->cols));
	for(int row = 0; row < cov->rows; row++)
		for(int col = 0; col < cov->cols; col++)
			REAL(expCovExt)[col * cov->rows + row] =
				omxMatrixElement(cov, row, col);

	
	if (means != NULL) {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, means->rows, means->cols));
		for(int row = 0; row < means->rows; row++)
			for(int col = 0; col < means->cols; col++)
				REAL(expMeanExt)[col * means->rows + row] =
					omxMatrixElement(means, row, col);
	} else {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, 0, 0));
	}

	Rf_setAttrib(algebra, Rf_install("ExpCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("ExpMean"), expMeanExt);
	Rf_setAttrib(algebra, Rf_install("numStats"), Rf_ScalarReal(omxDataDF(ox->data)));
	Rf_unprotect(2);
}

void omxInitNormalExpectation(omxExpectation* ox) {
	
	SEXP rObj = ox->rObj;
	omxState* currentState = ox->currentState;

    if(OMX_DEBUG) { mxLog("Initializing Normal expectation."); }

	omxNormalExpectation *one = (omxNormalExpectation*) R_alloc(1, sizeof(omxNormalExpectation));
	
	/* Set Expectation Calls and Structures */
	ox->computeFun = omxComputeNormalExpectation;
	ox->destructFun = omxDestroyNormalExpectation;
	ox->componentFun = omxGetNormalExpectationComponent;
	ox->populateAttrFun = omxPopulateNormalAttributes;
	ox->argStruct = (void*) one;
	
	/* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Processing cov."); }
	one->cov = omxNewMatrixFromSlot(rObj, currentState, "covariance");

	if(OMX_DEBUG) { mxLog("Processing Means."); }
	one->means = omxNewMatrixFromSlot(rObj, currentState, "means");
}

omxMatrix* omxGetNormalExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component){
/* Return appropriate parts of Expectation to the Fit Function */
	if(OMX_DEBUG) { mxLog("Normal expectation: %s requested--", component); }

	omxNormalExpectation* one = (omxNormalExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(!strncmp("cov", component, 3)) {
		retval = one->cov;
	} else if(!strncmp("mean", component, 4)) {
		retval = one->means;
	} else if(!strncmp("pvec", component, 4)) {
		// Once implemented, change compute function and return pvec
	}
	if (retval) omxRecompute(retval);
	
	return retval;
}
