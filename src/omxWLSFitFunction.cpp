 /*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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

#include "omxData.h"
#include <Eigen/Core>
// #include <Eigen/Dense>
#include "EnableWarnings.h"

struct omxWLSFitFunction : omxFitFunction {

	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* expectedSlope;
	omxMatrix* observedFlattened;
	omxMatrix* expectedFlattened;
	omxMatrix* P;
	omxMatrix* B;

	const char *type;
	const char *continuousType;
	bool fullWeight;

	omxWLSFitFunction() :
		type("WLS"), continuousType("cumulants"), fullWeight(true) {};
	virtual ~omxWLSFitFunction();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
	virtual void populateAttr(SEXP algebra);

	void prepData();
	virtual void invalidateCache()
	{
		omxFreeMatrix(observedFlattened);
		observedFlattened = 0;
	}
};

omxWLSFitFunction::~omxWLSFitFunction()
{
	if(OMX_DEBUG) {mxLog("Freeing WLS FitFunction.");}

	invalidateCache();
	omxWLSFitFunction* owo = this;
	omxFreeMatrix(owo->expectedFlattened);
	omxFreeMatrix(owo->B);
	omxFreeMatrix(owo->P);
}

void omxWLSFitFunction::compute(int want, FitContext *fc)
{
	auto *oo = this;
	auto *owo = this;
	if (want & FF_COMPUTE_INITIAL_FIT) return;

	/* Recompute and recopy */
	if(OMX_DEBUG) { mxLog("WLSFitFunction Computing expectation"); }
	omxExpectationCompute(fc, expectation, NULL);

	if ((want & FF_COMPUTE_PREOPTIMIZE) && !observedFlattened) {
		prepData();
		return;
	}

	if(OMX_DEBUG) { mxLog("Beginning WLS Evaluation.");}
	// Requires: Data, means, covariances.

	double sum = 0.0;

	omxMatrix *eCov, *eMeans, *oFlat, *eFlat;

	/* Locals for readability.  Compiler should cut through this. */
	eCov		= owo->expectedCov;
	eMeans 		= owo->expectedMeans;
	auto &eThresh   = oo->expectation->getThresholdInfo();
	if (!observedFlattened) return;
	oFlat		= observedFlattened;
	eFlat		= owo->expectedFlattened;

	EigenVectorAdaptor EeFlat(eFlat);
	omxExpectation *ex = expectation;
	normalToStdVector(eCov, eMeans, expectedSlope,
										[ex](int r, int c)->double{ return ex->getThreshold(r,c); },
										eThresh, EeFlat);

	omxCopyMatrix(B, oFlat);

	if(OMX_DEBUG) {omxPrintMatrix(eFlat, "....WLS Expected Vector: "); }
	omxDAXPY(-1.0, eFlat, B);
	//if(OMX_DEBUG) {omxPrintMatrix(B, "....WLS Observed - Expected Vector: "); }

	omxMatrix *weights = expectation->data->getSingleObsSummaryStats().acovMat;
	if(weights != NULL) {
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(weights, "....WLS Weight Matrix: "); }
		omxDGEMV(TRUE, 1.0, weights, B, 0.0, P);
	} else {
		// ULS Case: Memcpy faster than dgemv.
		omxCopyMatrix(P, B);
		omxTransposeMatrix(P);
	}

	sum = omxDDOT(P, B);

	oo->matrix->data[0] = sum;

	if(OMX_DEBUG) { mxLog("WLSFitFunction value comes to: %f.", oo->matrix->data[0]); }

}

void omxWLSFitFunction::populateAttr(SEXP algebra)
{
	if(OMX_DEBUG) { mxLog("Populating WLS Attributes."); }

	omxWLSFitFunction *argStruct = this;
	omxMatrix *expCovInt = argStruct->expectedCov;	    		// Expected covariance
	omxMatrix *expMeanInt = argStruct->expectedMeans;			// Expected means
	omxMatrix *weightInt = expectation->data->getSingleObsSummaryStats().acovMat;

	SEXP expCovExt, expMeanExt, gradients;
	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);

	if (expMeanInt != NULL) {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, 0, 0));
	}

	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(weightInt, "...WLS Weight Matrix: W"); }
	SEXP weightExt = NULL;
	if (weightInt) {
		Rf_protect(weightExt = Rf_allocMatrix(REALSXP, weightInt->rows, weightInt->cols));
		for(int row = 0; row < weightInt->rows; row++)
			for(int col = 0; col < weightInt->cols; col++)
				REAL(weightExt)[col * weightInt->rows + row] = weightInt->data[col * weightInt->rows + row];
	}


	if(0) {  /* TODO fix for new internal API
		int nLocs = Global->numFreeParams;
		double gradient[Global->numFreeParams];
		for(int loc = 0; loc < nLocs; loc++) {
			gradient[loc] = NA_REAL;
		}
		//oo->gradientFun(oo, gradient);
		Rf_protect(gradients = Rf_allocMatrix(REALSXP, 1, nLocs));

		for(int loc = 0; loc < nLocs; loc++)
			REAL(gradients)[loc] = gradient[loc];
		 */
	} else {
		Rf_protect(gradients = Rf_allocMatrix(REALSXP, 0, 0));
	}

	if(OMX_DEBUG) { mxLog("Installing populated WLS Attributes."); }
	Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
	if (weightExt) Rf_setAttrib(algebra, Rf_install("weights"), weightExt);
	Rf_setAttrib(algebra, Rf_install("gradients"), gradients);

	ProtectedSEXP Rsat(Rf_ScalarReal(0));
	Rf_setAttrib(algebra, Rf_install("SaturatedLikelihood"), Rsat);
	//Rf_setAttrib(algebra, Rf_install("IndependenceLikelihood"), Rf_ScalarReal(0));
	ProtectedSEXP Rmisfit(Rf_ScalarReal(omxMatrixElement(matrix, 0, 0)));
	Rf_setAttrib(algebra, Rf_install("ADFMisfit"), Rmisfit);
}

omxFitFunction *omxInitWLSFitFunction()
{ return new omxWLSFitFunction; }

void omxWLSFitFunction::prepData()
{
	auto *oo = this;
	auto *newObj = this;

	omxData* dataMat = oo->expectation->data;

	if (!matrix->currentState->isClone()) {
		std::vector<int> exoPred;
		expectation->getExogenousPredictors(exoPred);
		// how to prohibit def vars? TODO

		dataMat->prepObsStats(matrix->currentState, expectation->getDataColumnNames(),
				      exoPred, type, continuousType, fullWeight);
		if (isErrorRaised()) return;
	}

	auto &obsStat = dataMat->getSingleObsSummaryStats();
	//obsStat.log();
	omxMatrix *cov = obsStat.covMat;
	if (!cov) {
		omxRaiseErrorf("%s: an observed covariance matrix is required", name());
		return;
	}

	omxMatrix *means = obsStat.meansMat;
	omxMatrix *obsThresholdsMat = obsStat.thresholdMat;
	omxMatrix *weights = obsStat.acovMat;
	std::vector< omxThresholdColumn > &oThresh = obsStat.thresholdCols;

	auto &eThresh = oo->expectation->getThresholdInfo();

	// Error Checking: Observed/Expected means must agree.
	// ^ is XOR: true when one is false and the other is not.
	if((newObj->expectedMeans == NULL) ^ (means == NULL)) {
		if(newObj->expectedMeans != NULL) {
			if (eThresh.size() == 0) {
				omxRaiseError("Observed means not detected, but expected means specified.\n"
				   "To model means with all continuous data, you need to set allContinuousMethod='marginals'");
			} else {
				omxRaiseError("Means are required when the data include ordinal measurements");
			}
			return;
		} else {
			omxRaiseError("Observed means were provided, but an expected means matrix was not specified.\n  If you provide observed means, you must specify a model for the means.\n");
			return;
		}
	}

	if((eThresh.size()==0) ^ (obsThresholdsMat==0)) {
		if (eThresh.size()) {
			omxRaiseError("Observed thresholds not detected, but an expected thresholds matrix was specified.\n   If you wish to model the thresholds, you must provide observed thresholds.\n ");
			return;
		} else {
			omxRaiseError("Observed thresholds were provided, but an expected thresholds matrix was not specified.\nIf you provide observed thresholds, you must specify a model for the thresholds.\n");
			return;
		}
	}

	if (obsThresholdsMat) {
		for(int i = 0; i < int(oThresh.size()); i++) {
			eThresh[i].numThresholds = oThresh[i].numThresholds;
		}
	}

	if (OMX_DEBUG) {
		mxLog("expected thresholds:");
		for (auto &th : eThresh) { th.log(); }
		mxLog("observed thresholds:");
		for (auto &th : oThresh) { th.log(); }
	}

	int vectorSize = expectation->numSummaryStats();

	if(weights) {
		if (weights->rows != weights->cols || weights->cols != vectorSize) {
			mxThrow("Developer Error in WLS-based FitFunction object: WLS-based expectation specified an incorrectly-sized weight matrix (%d != %d).\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.", vectorSize, weights->cols);
		}

		EigenMatrixAdaptor Eweight(weights);
		Eigen::MatrixXd offDiagW = Eweight.triangularView<Eigen::StrictlyUpper>();
		double offDiag = offDiagW.array().abs().sum();
		oo->units = offDiag > 0.0? FIT_UNITS_SQUARED_RESIDUAL_CHISQ : FIT_UNITS_SQUARED_RESIDUAL;
	} else {
		oo->units = FIT_UNITS_SQUARED_RESIDUAL;
	}

	observedFlattened = omxInitMatrix(vectorSize, 1, TRUE, matrix->currentState);
	EigenVectorAdaptor oFlat(observedFlattened);
	if (obsThresholdsMat) {
		EigenMatrixAdaptor oTh(obsThresholdsMat);
		normalToStdVector(cov, means, obsStat.slopeMat,
                      [&oThresh, &oTh](int r,int c)->double{ return oTh(r, oThresh[c].column); },
                      oThresh, oFlat);
	} else {
		normalToStdVector(cov, means, obsStat.slopeMat, [](int r,int c)->double{ return 0; },
											oThresh, oFlat);
	}
	if(OMX_DEBUG) {omxPrintMatrix(observedFlattened, "....WLS Observed Vector: "); }

	/* Temporary storage for calculation */
	expectedFlattened = omxInitMatrix(vectorSize, 1, TRUE, matrix->currentState);
	P = omxInitMatrix(1, vectorSize, TRUE, matrix->currentState);
	B = omxInitMatrix(vectorSize, 1, TRUE, matrix->currentState);
}

void omxWLSFitFunction::init()
{
	auto *oo = this;

	if (!oo->expectation) { mxThrow("%s requires an expectation", name()); }

	if (R_has_slot(rObj, Rf_install("type"))) {
		ProtectedSEXP RwlsType(R_do_slot(rObj, Rf_install("type")));
		type = CHAR(STRING_ELT(RwlsType,0));
	}
	if (R_has_slot(rObj, Rf_install("continuousType"))) {
		ProtectedSEXP RwlsContType(R_do_slot(rObj, Rf_install("continuousType")));
		continuousType = CHAR(STRING_ELT(RwlsContType,0));
	}
	if (R_has_slot(rObj, Rf_install("fullWeight"))) {
		ProtectedSEXP RwlsFullWeight(R_do_slot(rObj, Rf_install("fullWeight")));
		fullWeight = Rf_asLogical(RwlsFullWeight);
	}

	if (!fullWeight && !strEQ(type, "ULS")) {
		mxThrow("%s: !fullWeight && !strEQ(type, ULS)", name());
	}

	expectedCov = omxGetExpectationComponent(oo->expectation, "cov");
	expectedMeans = omxGetExpectationComponent(oo->expectation, "means");
	expectedSlope = omxGetExpectationComponent(expectation, "slope");

  if (expectedSlope) {
    expectation->invalidateCache();
    expectation->connectToData();
  }

	observedFlattened = 0;
	expectedFlattened = 0;
	B = 0;
	P = 0;
  canDuplicate = true;
}
