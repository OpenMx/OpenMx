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
	omxMatrix* observedFlattened;
	omxMatrix* expectedFlattened;
	omxMatrix* weights;
	omxMatrix* P;
	omxMatrix* B;
	int n;
	int fullWls;
	int numOrdinal;
	
	omxWLSFitFunction() :standardMeans(0), standardThresholds(0) {};
	virtual ~omxWLSFitFunction();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
	
	// 'standard' prefix variables are temp space used by flattenDataToVector
	omxMatrix* standardCov;
	omxMatrix* standardMeans;
	omxMatrix* standardThresholds;
	void flattenDataToVector(omxMatrix* cov, omxMatrix* means, omxMatrix *thresholdMat,
				 std::vector< omxThresholdColumn > &thresholds, omxMatrix* vector);
};

void omxWLSFitFunction::flattenDataToVector(omxMatrix* cov, omxMatrix* means, omxMatrix *thresholdMat,
			 std::vector< omxThresholdColumn > &thresholds, omxMatrix* vector)
{
	EigenVectorAdaptor vec1(vector);
	normalToStdVector(cov, means, thresholdMat, numOrdinal, thresholds, vec1);
}

omxWLSFitFunction::~omxWLSFitFunction()
{
	if(OMX_DEBUG) {mxLog("Freeing WLS FitFunction.");}
	
	omxWLSFitFunction* owo = this;
	omxFreeMatrix(owo->observedFlattened);
	omxFreeMatrix(owo->expectedFlattened);
	omxFreeMatrix(owo->B);
	omxFreeMatrix(owo->P);
	omxFreeMatrix(owo->standardCov);
	omxFreeMatrix(owo->standardMeans);
	omxFreeMatrix(owo->standardThresholds);
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
	
	flattenDataToVector(eCov, eMeans, expThresholdsMat, eThresh, eFlat);
	
	omxCopyMatrix(B, oFlat);
	
	//if(OMX_DEBUG) {omxPrintMatrix(B, "....WLS Observed Vector: "); }
	if(OMX_DEBUG) {omxPrintMatrix(eFlat, "....WLS Expected Vector: "); }
	omxDAXPY(-1.0, eFlat, B);
	//if(OMX_DEBUG) {omxPrintMatrix(B, "....WLS Observed - Expected Vector: "); }
	
	if(weights != NULL) {
		//if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(weights, "....WLS Weight Matrix: "); }
		omxDGEMV(TRUE, 1.0, weights, B, 0.0, P);
	} else {
		// ULS Case: Memcpy faster than dgemv.
		omxCopyMatrix(P, B);
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

void omxWLSFitFunction::init()
{
	auto *oo = this;
	auto *newObj = this;
	
	omxState *currentState = oo->matrix->currentState;
	omxMatrix *cov, *means;
	
	if(OMX_DEBUG) { mxLog("Initializing WLS FitFunction function."); }
	
	if(OMX_DEBUG) { mxLog("Retrieving expectation.\n"); }
	if (!oo->expectation) { Rf_error("%s requires an expectation", oo->fitType); }
	
	if(OMX_DEBUG) { mxLog("Retrieving data.\n"); }
	omxData* dataMat = oo->expectation->data;
	if (dataMat->hasDefinitionVariables()) Rf_error("%s: def vars not implemented", oo->name());
	
	if(!strEQ(omxDataType(dataMat), "acov") && !strEQ(omxDataType(dataMat), "cov")) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "WLS FitFunction unable to handle data type %s.  Data must be of type 'acov'.\n", omxDataType(dataMat));
		omxRaiseError(errstr);
		free(errstr);
		if(OMX_DEBUG) { mxLog("WLS FitFunction unable to handle data type %s.  Aborting.", omxDataType(dataMat)); }
		return;
	}
	
	fullWls = strEQ("WLS", CHAR(Rf_asChar(R_do_slot(rObj, Rf_install("weights")))));

	oo->units = fullWls? FIT_UNITS_SQUARED_RESIDUAL_CHISQ : FIT_UNITS_SQUARED_RESIDUAL;
	
	/* Get Expectation Elements */
	newObj->expectedCov = omxGetExpectationComponent(oo->expectation, "cov");
	newObj->expectedMeans = omxGetExpectationComponent(oo->expectation, "means");
	
	/* Read and set expected means, variances, and weights */
	dataMat->permute(oo->expectation->getDataColumns());

	cov = omxDataCovariance(dataMat);
	means = omxDataMeans(dataMat);
	weights = omxDataAcov(dataMat);

	std::vector< omxThresholdColumn > &oThresh = omxDataThresholds(oo->expectation->data);

	newObj->n = omxDataNumObs(dataMat);
	
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
	
	if((eThresh.size()==0) ^ (oThresh.size()==0)) {
		if (eThresh.size()) {
			omxRaiseError("Observed thresholds not detected, but an expected thresholds matrix was specified.\n   If you wish to model the thresholds, you must provide observed thresholds.\n ");
			return;
		} else {
			omxRaiseError("Observed thresholds were provided, but an expected thresholds matrix was not specified.\nIf you provide observed thresholds, you must specify a model for the thresholds.\n");
			return;
		}
	}

	for(int i = 0, ei=0; i < int(oThresh.size()); i++) {
		while (ei < int(eThresh.size()) && eThresh[ei].dColumn != oThresh[i].dColumn) ++ei;
		eThresh[ei].numThresholds = oThresh[i].numThresholds;
	}

	if (OMX_DEBUG) {
		mxLog("expected thresholds:");
		for (auto &th : eThresh) { th.log(); }
		mxLog("observed thresholds:");
		for (auto &th : oThresh) { th.log(); }
	}
	
	/* Error check weight matrix size */
	int ncol = cov->cols;
	int vectorSize = expectation->numSummaryStats();
	if(OMX_DEBUG) { mxLog("Intial WLSFitFunction vectorSize comes to: %d.", vectorSize); }
	
	if(weights != NULL && (weights->rows != weights->cols || weights->cols != vectorSize)) {
		omxRaiseError("Developer Error in WLS-based FitFunction object: WLS-based expectation specified an incorrectly-sized weight matrix.\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.");
		return;
	}
	
	
	// FIXME: More Rf_error checking for incoming Fit Functions
	
	/* Temporary storage for calculation */
	newObj->observedFlattened = omxInitMatrix(vectorSize, 1, TRUE, currentState);
	newObj->expectedFlattened = omxInitMatrix(vectorSize, 1, TRUE, currentState);
	newObj->P = omxInitMatrix(1, vectorSize, TRUE, currentState);
	newObj->B = omxInitMatrix(vectorSize, 1, TRUE, currentState);
	newObj->standardCov = omxInitMatrix(ncol, ncol, TRUE, currentState);
	if (oo->expectation->thresholdsMat) {
		newObj->standardThresholds = omxInitMatrix(oo->expectation->thresholdsMat->rows, oo->expectation->thresholdsMat->cols, TRUE, currentState);
	}
	if(means){
		newObj->standardMeans = omxInitMatrix(1, ncol, TRUE, currentState);
	}
	omxMatrix *obsThresholdsMat = oo->expectation->data->obsThresholdsMat;
	if (obsThresholdsMat && oo->expectation->thresholdsMat) {
		if (obsThresholdsMat->rows != oo->expectation->thresholdsMat->rows ||
		    obsThresholdsMat->cols != oo->expectation->thresholdsMat->cols) {
			omxRaiseError("Observed and expected threshold matrices must have the same number of rows and columns");
		}
	}
	
	flattenDataToVector(cov, means, obsThresholdsMat, oThresh, newObj->observedFlattened);
	flattenDataToVector(newObj->expectedCov, newObj->expectedMeans, oo->expectation->thresholdsMat,
				eThresh, newObj->expectedFlattened);
}
