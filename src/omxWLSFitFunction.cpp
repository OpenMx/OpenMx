 /*
 *  Copyright 2007-2018 by the individuals mentioned in the source code history
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
	omxMatrix* weights;
	omxMatrix* P;
	omxMatrix* B;
	int numOrdinal;
	int vectorSize;
	
	omxWLSFitFunction() {};
	virtual ~omxWLSFitFunction();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
	
	void flattenDataToVector(omxMatrix* cov, omxMatrix* means, omxMatrix *slope, omxMatrix *thresholdMat,
				 std::vector< omxThresholdColumn > &thresholds, omxMatrix* vector);
	void prepData();
	virtual void invalidateCache() { prepData(); }
};

void omxWLSFitFunction::flattenDataToVector(omxMatrix* cov, omxMatrix* means, omxMatrix *slope,
					    omxMatrix *thresholdMat,
			 std::vector< omxThresholdColumn > &thresholds, omxMatrix* vector)
{
	EigenVectorAdaptor vec1(vector);
	normalToStdVector(cov, means, slope, thresholdMat, numOrdinal, thresholds, vec1);
}

omxWLSFitFunction::~omxWLSFitFunction()
{
	if(OMX_DEBUG) {mxLog("Freeing WLS FitFunction.");}
	
	omxWLSFitFunction* owo = this;
	omxFreeMatrix(owo->observedFlattened);
	omxFreeMatrix(owo->expectedFlattened);
	omxFreeMatrix(owo->B);
	omxFreeMatrix(owo->P);
}


void omxWLSFitFunction::compute(int want, FitContext *fc)
{
	auto *oo = this;
	auto *owo = this;
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;
	
	if(OMX_DEBUG) { mxLog("Beginning WLS Evaluation.");}
	// Requires: Data, means, covariances.
	
	double sum = 0.0;
	
	omxMatrix *eCov, *eMeans, *oFlat, *eFlat;
	
	/* Locals for readability.  Compiler should cut through this. */
	eCov		= owo->expectedCov;
	eMeans 		= owo->expectedMeans;
	auto &eThresh   = oo->expectation->getThresholdInfo();
	oFlat		= owo->observedFlattened;
	eFlat		= owo->expectedFlattened;
	int onei	= 1;
	
	/* Recompute and recopy */
	if(OMX_DEBUG) { mxLog("WLSFitFunction Computing expectation"); }
	omxExpectationCompute(fc, expectation, NULL);
	
	omxMatrix *expThresholdsMat = expectation->thresholdsMat;
	
	flattenDataToVector(eCov, eMeans, expectedSlope, expThresholdsMat, eThresh, eFlat);
	
	omxCopyMatrix(B, oFlat);
	
	if(OMX_DEBUG) {omxPrintMatrix(eFlat, "....WLS Expected Vector: "); }
	omxDAXPY(-1.0, eFlat, B);
	//if(OMX_DEBUG) {omxPrintMatrix(B, "....WLS Observed - Expected Vector: "); }
	
	if(weights != NULL) {
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(weights, "....WLS Weight Matrix: "); }
		omxDGEMV(TRUE, 1.0, weights, B, 0.0, P);
	} else {
		// ULS Case: Memcpy faster than dgemv.
		omxCopyMatrix(P, B);
		omxTransposeMatrix(P);
	}
	
	sum = F77_CALL(ddot)(&(P->cols), P->data, &onei, B->data, &onei);
	
	oo->matrix->data[0] = sum;
	
	if(OMX_DEBUG) { mxLog("WLSFitFunction value comes to: %f.", oo->matrix->data[0]); }
	
}

void omxWLSFitFunction::populateAttr(SEXP algebra)
{
	if(OMX_DEBUG) { mxLog("Populating WLS Attributes."); }
	
	omxWLSFitFunction *argStruct = this;
	omxMatrix *expCovInt = argStruct->expectedCov;	    		// Expected covariance
	omxMatrix *expMeanInt = argStruct->expectedMeans;			// Expected means
	omxMatrix *weightInt = argStruct->weights;			// Expected means
	
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

	if (vectorSize != expectation->numSummaryStats())
		Rf_error("%s: vectorSize changed from %d -> %d",
			 vectorSize, expectation->numSummaryStats());

	std::vector<int> exoPred;
	expectation->getExogenousPredictors(exoPred);

	omxData* dataMat = oo->expectation->data;

	if (dataMat->defVars.size() == exoPred.size()) {
		// OK
	} else if (dataMat->hasDefinitionVariables()) Rf_error("%s: def vars not implemented", oo->name());
	
	// For multiple threads, need to grab parent's info TODO
	dataMat->prepObsStats(matrix->currentState, expectation->getDataColumnNames(), exoPred);

	auto &obsStat = dataMat->getSingleObsSummaryStats();
	//obsStat.log();
	omxMatrix *cov = obsStat.covMat;
	if (!cov) {
		omxRaiseErrorf("%s: an observed covariance matrix is required", name());
		return;
	}

	omxMatrix *means = obsStat.meansMat;
	omxMatrix *obsThresholdsMat = obsStat.thresholdMat;
	weights = obsStat.acovMat;
	std::vector< omxThresholdColumn > &oThresh = obsStat.thresholdCols;

	numOrdinal = oo->expectation->numOrdinal;
	auto &eThresh = oo->expectation->getThresholdInfo();

	if (eThresh.size() && !means) {
		omxRaiseError("Means are required when the data include ordinal measurements");
		return;
	}

	// Error Checking: Observed/Expected means must agree.  
	// ^ is XOR: true when one is false and the other is not.
	if((newObj->expectedMeans == NULL) ^ (means == NULL)) {
		if(newObj->expectedMeans != NULL) {
			omxRaiseError("Observed means not detected, but an expected means matrix was specified.\n  If you  wish to model the means, you must provide observed means.\n");
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
	
	if(weights) {
		if (weights->rows != weights->cols || weights->cols != vectorSize) {
			omxRaiseErrorf("Developer Error in WLS-based FitFunction object: WLS-based expectation specified an incorrectly-sized weight matrix (%d != %d).\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.", vectorSize, weights->cols);
			return;
		}
	
		EigenMatrixAdaptor Eweight(weights);
		Eigen::MatrixXd offDiagW = Eweight.triangularView<Eigen::StrictlyUpper>();
		double offDiag = offDiagW.array().abs().sum();
		oo->units = offDiag > 0.0? FIT_UNITS_SQUARED_RESIDUAL_CHISQ : FIT_UNITS_SQUARED_RESIDUAL;
	} else {
		oo->units = FIT_UNITS_SQUARED_RESIDUAL;
	}
	
	if (obsThresholdsMat && oo->expectation->thresholdsMat) {
		if (obsThresholdsMat->rows != oo->expectation->thresholdsMat->rows ||
		    obsThresholdsMat->cols != oo->expectation->thresholdsMat->cols) {
			omxRaiseError("Observed and expected threshold matrices must have the same number of rows and columns");
		}
	}
	
	flattenDataToVector(cov, means, obsStat.slopeMat, obsThresholdsMat,
			    oThresh, newObj->observedFlattened);
	if(OMX_DEBUG) {omxPrintMatrix(newObj->observedFlattened, "....WLS Observed Vector: "); }
	flattenDataToVector(newObj->expectedCov, newObj->expectedMeans,
			    newObj->expectedSlope, oo->expectation->thresholdsMat,
				eThresh, newObj->expectedFlattened);
}

void omxWLSFitFunction::init()
{
	auto *oo = this;
	
	omxState *currentState = oo->matrix->currentState;
	
	if (!oo->expectation) { Rf_error("%s requires an expectation", name()); }
	
	expectedCov = omxGetExpectationComponent(oo->expectation, "cov");
	expectedMeans = omxGetExpectationComponent(oo->expectation, "means");
	expectedSlope = omxGetExpectationComponent(expectation, "slope");

	/* Error check weight matrix size */
	vectorSize = expectation->numSummaryStats();
	if(OMX_DEBUG) { mxLog("Intial WLSFitFunction vectorSize comes to: %d.", vectorSize); }
	
	/* Temporary storage for calculation */
	observedFlattened = omxInitMatrix(vectorSize, 1, TRUE, currentState);
	expectedFlattened = omxInitMatrix(vectorSize, 1, TRUE, currentState);
	P = omxInitMatrix(1, vectorSize, TRUE, currentState);
	B = omxInitMatrix(vectorSize, 1, TRUE, currentState);

	prepData();
}
