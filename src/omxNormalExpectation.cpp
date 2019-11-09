/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
#include "omxDefines.h"
#include "EnableWarnings.h"

struct omxNormalExpectation : public omxExpectation {
	typedef omxExpectation super;

	omxMatrix *cov, *means; // observed covariance and means

	double logDetObserved;
	double n;

	omxNormalExpectation(omxState *st, int num) : super(st, num) {}
	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual void populateAttr(SEXP expectation);
	virtual omxMatrix *getComponent(const char*);
};

void omxNormalExpectation::compute(FitContext *fc, const char *, const char *) {
	omxNormalExpectation* one = this;

	omxRecompute(one->cov, fc);
	if(one->means != NULL)
	    omxRecompute(one->means, fc);
	if (one->thresholdsMat) omxRecompute(one->thresholdsMat, fc);
}

void omxNormalExpectation::populateAttr(SEXP algebra) {
    if(OMX_DEBUG) { mxLog("Populating Normal Attributes."); }

	omxRecompute(cov, NULL);
	if(means != NULL) omxRecompute(means, NULL);

	{
		SEXP expCovExt;
	ScopedProtect p1(expCovExt, Rf_allocMatrix(REALSXP, cov->rows, cov->cols));
	for(int row = 0; row < cov->rows; row++)
		for(int col = 0; col < cov->cols; col++)
			REAL(expCovExt)[col * cov->rows + row] =
				omxMatrixElement(cov, row, col);
	Rf_setAttrib(algebra, Rf_install("ExpCov"), expCovExt);
	}

	
	if (means != NULL) {
		SEXP expMeanExt;
		ScopedProtect p1(expMeanExt, Rf_allocMatrix(REALSXP, means->rows, means->cols));
		for(int row = 0; row < means->rows; row++)
			for(int col = 0; col < means->cols; col++)
				REAL(expMeanExt)[col * means->rows + row] =
					omxMatrixElement(means, row, col);
		Rf_setAttrib(algebra, Rf_install("ExpMean"), expMeanExt);
	} else {
		SEXP expMeanExt;
		ScopedProtect p1(expMeanExt, Rf_allocMatrix(REALSXP, 0, 0));
		Rf_setAttrib(algebra, Rf_install("ExpMean"), expMeanExt);
	}

	ProtectedSEXP RnumStats(Rf_ScalarReal(omxDataDF(data)));
	Rf_setAttrib(algebra, Rf_install("numStats"), RnumStats);
}

omxExpectation *omxInitNormalExpectation(omxState *st, int num)
{ return new omxNormalExpectation(st, num); }

void omxNormalExpectation::init()
{
    if(OMX_DEBUG) { mxLog("Initializing Normal expectation."); }

    omxNormalExpectation *one = this;
	
	/* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Processing cov."); }
	one->cov = omxNewMatrixFromSlot(rObj, currentState, "covariance");

	if(OMX_DEBUG) { mxLog("Processing Means."); }
	one->means = omxNewMatrixFromSlot(rObj, currentState, "means");
}

omxMatrix* omxNormalExpectation::getComponent(const char* component){
/* Return appropriate parts of Expectation to the Fit Function */
	if(OMX_DEBUG) { mxLog("Normal expectation: %s requested--", component); }

	omxNormalExpectation* one = this;
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = one->cov;
	} else if(strEQ("means", component)) {
		retval = one->means;
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}
	if (retval) omxRecompute(retval, NULL);
	
	return retval;
}
