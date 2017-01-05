 /*
 *  Copyright 2007-2017 The OpenMx Project
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

struct omxWLSFitFunction {

	std::vector< omxThresholdColumn > eThresh; // with column mapped by expectation's dataColumns
	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCovPermuted;
	omxMatrix* expectedMeansPermuted;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* observedFlattened;
	omxMatrix* expectedFlattened;
	omxMatrix* weights;
	omxMatrix* P;
	omxMatrix* B;
	omxMatrix* standardExpectedCov;
	omxMatrix* standardExpectedMeans;
	omxMatrix* standardExpectedThresholds;
	int n;

	omxWLSFitFunction() :standardExpectedMeans(0), standardExpectedThresholds(0) {};

	void permuteExpectation(omxExpectation *expectation) {
		auto dc = expectation->getDataColumns();
		if (!dc.size()) {
			omxCopyMatrix(expectedCov, expectedCovPermuted);
			if (expectedMeans) omxCopyMatrix(expectedMeans, expectedMeansPermuted);
			return;
		}
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm(dc);
		EigenMatrixAdaptor Ecovp(expectedCovPermuted);
		EigenMatrixAdaptor Ecov(expectedCov);
		Ecov.derived() = perm * Ecovp * perm.transpose();
		if (expectedMeans) {
			EigenVectorAdaptor Emeanp(expectedMeansPermuted);
			EigenVectorAdaptor Emean(expectedMeans);
			Emean.derived() = perm * Emeanp;
		}
	}
};

#ifdef SHADOW_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#endif

static void flattenDataToVector(omxMatrix* cov, omxMatrix* means, omxMatrix *thresholdMat,
			 std::vector< omxThresholdColumn > &thresholds, omxMatrix* vector) {
	// TODO: vectorize data flattening
	// if(OMX_DEBUG) { mxLog("Flattening out data vectors: cov 0x%x, mean 0x%x, thresh 0x%x[n=%d] ==> 0x%x", 
	//         cov, means, thresholds, nThresholds, vector); }
	
	int nextLoc = 0;
	for(int j = 0; j < cov->rows; j++) {
		for(int k = j; k < cov->rows; k++) {
			omxSetVectorElement(vector, nextLoc, omxMatrixElement(cov, j, k)); // Use upper triangle in case of SYMM-style mat.
			nextLoc++;
		}
	}
	if (means != NULL) {
		for(int j = 0; j < cov->rows; j++) {
			omxSetVectorElement(vector, nextLoc, omxVectorElement(means, j));
			nextLoc++;
		}
	}
	for(int j = 0; j < int(thresholds.size()); j++) {
		omxThresholdColumn* thresh = &thresholds[j];
		for(int k = 0; k < thresh->numThresholds; k++) {
			omxSetVectorElement(vector, nextLoc, omxMatrixElement(thresholdMat, k, thresh->column));
			nextLoc++;
		}
	}
}

static void omxDestroyWLSFitFunction(omxFitFunction *oo) {
	
	if(OMX_DEBUG) {mxLog("Freeing WLS FitFunction.");}
	if(oo->argStruct == NULL) return;
	
	omxWLSFitFunction* owo = ((omxWLSFitFunction*)oo->argStruct);
	omxFreeMatrix(owo->expectedCov);
	omxFreeMatrix(owo->expectedMeans);
	omxFreeMatrix(owo->observedFlattened);
	omxFreeMatrix(owo->expectedFlattened);
	omxFreeMatrix(owo->B);
	omxFreeMatrix(owo->P);
	omxFreeMatrix(owo->standardExpectedCov);
	omxFreeMatrix(owo->standardExpectedMeans);
	omxFreeMatrix(owo->standardExpectedThresholds);
	delete owo;
}


static void standardizeCovMeansThresholds(omxMatrix* inCov, omxMatrix* inMeans,
			omxMatrix* inThresholdsMat, std::vector< omxThresholdColumn > &thresholds,
			omxMatrix* outCov, omxMatrix* outMeans, omxMatrix* outThresholdsMat) {
	//omxMatrix* pass[1];
	//pass[0] = inCov;
	//omxCovToCor(fc, pass, 1, outCov);
	
	Eigen::ArrayXd stddev;
	EigenMatrixAdaptor egInCov(inCov);
	EigenMatrixAdaptor egOutCov(outCov);
	
	stddev = egInCov.diagonal().array().sqrt();
	Eigen::ArrayXd stddevUse = Eigen::ArrayXd::Ones(egInCov.rows());
	
	// standardize mean and thresholds
	if(inMeans != NULL) {
		EigenMatrixAdaptor egInM(inMeans);
		EigenMatrixAdaptor egOutM(outMeans);
		// means
		egOutM = egInM;
		
		// thresholds
		if(inThresholdsMat != NULL){
			EigenMatrixAdaptor egInThr(inThresholdsMat);
			EigenMatrixAdaptor egOutThr(outThresholdsMat);
			if(OMX_DEBUG) {
				mxLog("Expected thresholds have size %d.", int(thresholds.size()));
			}
			for(int j = 0; j < int(thresholds.size()); j++) {
				omxThresholdColumn* thresh = &thresholds[j];
				if(OMX_DEBUG) {
					mxLog("Data column %d has %d thresholds stored in threshold column %d",
					      thresh->dColumn, thresh->numThresholds, thresh->column);
				}
				int dcol = thresh->dColumn;
				for(int k = 0; k < thresh->numThresholds; k++) {
					egOutThr(k, thresh->column) =
						( egInThr(k, thresh->column) - egInM(0, dcol) ) / stddev[dcol];
				}
				if(thresh->numThresholds > 0){
					egOutM(0, dcol) = 0.0;
					stddevUse[dcol] = stddev[dcol];
				}
			}
		}
	} else {
		// means
		//omxCopyMatrix(outMeans, inMeans); //egOutM = egInM;
		
		// thresholds
		if(inThresholdsMat != NULL){
			EigenMatrixAdaptor egInThr(inThresholdsMat);
			EigenMatrixAdaptor egOutThr(outThresholdsMat);
			for(int j = 0; j < int(thresholds.size()); j++) {
				omxThresholdColumn* thresh = &thresholds[j];
				// thresh->column should probably be thresh->dColumn as above
				// This code path is untested and maybe should be deleted. TODO
				for(int k = 0; k < thresh->numThresholds; k++) {
					egOutThr(k, thresh->column) = egInThr(k, thresh->column) / stddev[thresh->column];
				}
				stddevUse[thresh->column] = stddev[thresh->column];
			}
		}
	}
	
	// standardize covariance
	for(int i = 0; i < egInCov.rows(); i++) {
		for(int j = 0; j <= i; j++) {
			egOutCov(i,j) = egInCov(i, j) / (stddevUse[i] * stddevUse[j]);
			egOutCov(j,i) = egOutCov(i,j);
		}
	}
	
}


static void omxCallWLSFitFunction(omxFitFunction *oo, int want, FitContext *fc) {
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;
	
	if(OMX_DEBUG) { mxLog("Beginning WLS Evaluation.");}
	// Requires: Data, means, covariances.
	
	double sum = 0.0;
	
	omxMatrix *eCov, *eMeans, *P, *B, *weights, *oFlat, *eFlat;
	omxMatrix *seCov, *seMeans, *seThresholdsMat;
	
	omxWLSFitFunction *owo = ((omxWLSFitFunction*)oo->argStruct);
	
	/* Locals for readability.  Compiler should cut through this. */
	eCov		= owo->expectedCov;
	eMeans 		= owo->expectedMeans;
	auto &eThresh   = owo->eThresh;
	oFlat		= owo->observedFlattened;
	eFlat		= owo->expectedFlattened;
	weights		= owo->weights;
	B			= owo->B;
	P			= owo->P;
	seCov		= owo->standardExpectedCov;
	seMeans		= owo->standardExpectedMeans;
	seThresholdsMat = owo->standardExpectedThresholds;
	int onei	= 1;
	
	omxExpectation* expectation = oo->expectation;
	
	/* Recompute and recopy */
	if(OMX_DEBUG) { mxLog("WLSFitFunction Computing expectation"); }
	omxExpectationCompute(fc, expectation, NULL);
	owo->permuteExpectation(expectation);
	
	omxMatrix *expThresholdsMat = expectation->thresholdsMat;
	
	standardizeCovMeansThresholds(eCov, eMeans, expThresholdsMat, eThresh,
			seCov, seMeans, seThresholdsMat);
	if(expThresholdsMat != NULL){
		flattenDataToVector(seCov, seMeans, seThresholdsMat, eThresh, eFlat);
	} else {
		flattenDataToVector(eCov, eMeans, expThresholdsMat, eThresh, eFlat);
	}
	
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

static void omxPopulateWLSAttributes(omxFitFunction *oo, SEXP algebra) {
	if(OMX_DEBUG) { mxLog("Populating WLS Attributes."); }
	
	omxWLSFitFunction *argStruct = ((omxWLSFitFunction*)oo->argStruct);
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
	
	Rf_setAttrib(algebra, Rf_install("SaturatedLikelihood"), Rf_ScalarReal(0));
	//Rf_setAttrib(algebra, Rf_install("IndependenceLikelihood"), Rf_ScalarReal(0));
	Rf_setAttrib(algebra, Rf_install("ADFMisfit"), Rf_ScalarReal(omxMatrixElement(oo->matrix, 0, 0)));
}

static void omxSetWLSFitFunctionCalls(omxFitFunction* oo) {
	
	/* Set FitFunction Calls to WLS FitFunction Calls */
	oo->computeFun = omxCallWLSFitFunction;
	oo->destructFun = omxDestroyWLSFitFunction;
	oo->populateAttrFun = omxPopulateWLSAttributes;
}

void omxInitWLSFitFunction(omxFitFunction* oo) {
	
	omxMatrix *cov, *means, *weights;
	
	if(OMX_DEBUG) { mxLog("Initializing WLS FitFunction function."); }
	
	int vectorSize = 0;
	
	omxSetWLSFitFunctionCalls(oo);
	
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
	
	omxWLSFitFunction *newObj = new omxWLSFitFunction;
	oo->argStruct = (void*)newObj;
	oo->units = FIT_UNITS_SQUARED_RESIDUAL;
	
	/* Get Expectation Elements */
	newObj->expectedCovPermuted = omxGetExpectationComponent(oo->expectation, "cov");
	int esize = newObj->expectedCovPermuted->rows;
	newObj->expectedCov = omxInitMatrix(esize, esize, TRUE, oo->matrix->currentState);
	newObj->expectedMeansPermuted = omxGetExpectationComponent(oo->expectation, "means");
	newObj->expectedMeans = 0;
	if (newObj->expectedMeansPermuted) {
		newObj->expectedMeans = omxInitMatrix(1, esize, TRUE, oo->matrix->currentState);
	}
	newObj->permuteExpectation(oo->expectation);
	
	// FIXME: threshold structure should be asked for by omxGetExpectationComponent
	
	/* Read and set expected means, variances, and weights */
	cov = omxDataCovariance(dataMat);
	means = omxDataMeans(dataMat);
	weights = omxDataAcov(dataMat);
	
	newObj->observedCov = cov;
	newObj->observedMeans = means;
	newObj->weights = weights;
	newObj->n = omxDataNumObs(dataMat);
	
	auto dc = oo->expectation->getDataColumns();
	auto &origThresh = oo->expectation->thresholds;
	auto &eThresh = newObj->eThresh;
	for (auto &th : origThresh) {
		if (th.numThresholds == 0) continue;
		omxThresholdColumn adj = th;
		adj.dColumn = dc[th.dColumn];
		eThresh.push_back(adj);
	}
	std::sort(eThresh.begin(), eThresh.end(),
		  [](omxThresholdColumn &a, omxThresholdColumn &b) -> bool
		  { return a.dColumn < b.dColumn; });

	std::vector< omxThresholdColumn > &oThresh = omxDataThresholds(oo->expectation->data);
	
	// Error Checking: Observed/Expected means must agree.  
	// ^ is XOR: true when one is false and the other is not.
	if((newObj->expectedMeans == NULL) ^ (newObj->observedMeans == NULL)) {
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
	
	/* Error check weight matrix size */
	int ncol = newObj->observedCov->cols;
	vectorSize = (ncol * (ncol + 1) ) / 2;
	if(newObj->expectedMeans != NULL) {
		vectorSize = vectorSize + ncol;
	}
	for(int i = 0; i < int(oThresh.size()); i++) {
		vectorSize = vectorSize + oThresh[i].numThresholds;
	}
	if(OMX_DEBUG) { mxLog("Intial WLSFitFunction vectorSize comes to: %d.", vectorSize); }
	
	if(weights != NULL && (weights->rows != weights->cols || weights->cols != vectorSize)) {
		omxRaiseError("Developer Error in WLS-based FitFunction object: WLS-based expectation specified an incorrectly-sized weight matrix.\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.");
		return;
	}
	
	
	// FIXME: More Rf_error checking for incoming Fit Functions
	
	/* Temporary storage for calculation */
	newObj->observedFlattened = omxInitMatrix(vectorSize, 1, TRUE, oo->matrix->currentState);
	newObj->expectedFlattened = omxInitMatrix(vectorSize, 1, TRUE, oo->matrix->currentState);
	newObj->P = omxInitMatrix(1, vectorSize, TRUE, oo->matrix->currentState);
	newObj->B = omxInitMatrix(vectorSize, 1, TRUE, oo->matrix->currentState);
	newObj->standardExpectedCov = omxInitMatrix(ncol, ncol, TRUE, oo->matrix->currentState);
	if (oo->expectation->thresholdsMat) {
		newObj->standardExpectedThresholds = omxInitMatrix(oo->expectation->thresholdsMat->rows, oo->expectation->thresholdsMat->cols, TRUE, oo->matrix->currentState);
	}
	if(means){
		newObj->standardExpectedMeans = omxInitMatrix(1, ncol, TRUE, oo->matrix->currentState);
	}
	omxMatrix *obsThresholdsMat = oo->expectation->data->obsThresholdsMat;
	
	flattenDataToVector(newObj->observedCov, newObj->observedMeans, obsThresholdsMat, oThresh, newObj->observedFlattened);
	flattenDataToVector(newObj->expectedCov, newObj->expectedMeans, oo->expectation->thresholdsMat,
				eThresh, newObj->expectedFlattened);

}
