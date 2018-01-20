 /*
 *  Copyright 2007-2018 The OpenMx Project
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

	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;
	omxMatrix* observedFlattened;
	omxMatrix* expectedFlattened;
	omxMatrix* weights;
	omxMatrix* P;
	omxMatrix* B;
	int n;
	int fullWls;
	
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
	// TODO: vectorize data flattening

	Eigen::ArrayXd stddev;
	EigenMatrixAdaptor egInCov(cov);
	EigenMatrixAdaptor egOutCov(standardCov);

	stddev = egInCov.diagonal().array().sqrt();
	Eigen::ArrayXd stddevUse = Eigen::ArrayXd::Ones(egInCov.rows());

	// standardize thresholds with corresponding transform to mean and covarinace
	if (means) {
		EigenMatrixAdaptor egInM(means);
		EigenMatrixAdaptor egOutM(standardMeans);
		egOutM = egInM;

		if (thresholdMat) {
			EigenMatrixAdaptor egInThr(thresholdMat);
			EigenMatrixAdaptor egOutThr(standardThresholds);
			for(int j = 0; j < int(thresholds.size()); j++) {
				omxThresholdColumn* thresh = &thresholds[j];
				if (thresh->numThresholds == 0) continue;
				int dcol = thresh->dColumn;
				for(int k = 0; k < thresh->numThresholds; k++) {
					egOutThr(k, thresh->column) =
						( egInThr(k, thresh->column) - egInM(0, dcol) ) / stddev[dcol];
				}
				egOutM(0, dcol) = 0.0;
				stddevUse[dcol] = stddev[dcol];
			}
		}
	}

	// standardize covariance of ordinal indicators
	for(int i = 0; i < egInCov.rows(); i++) {
		for(int j = 0; j <= i; j++) {
			egOutCov(i,j) = egInCov(i, j) / (stddevUse[i] * stddevUse[j]);
			egOutCov(j,i) = egOutCov(i,j);
		}
	}
	
	int nextLoc = 0;
	for(int j = 0; j < cov->rows; j++) {
		for(int k = j; k < cov->rows; k++) {
			omxSetVectorElement(vector, nextLoc, egOutCov(j, k));
			nextLoc++;
		}
	}
	if (means) {
		EigenMatrixAdaptor egOutM(standardMeans);
		for(int j = 0; j < cov->rows; j++) {
			omxSetVectorElement(vector, nextLoc, egOutM(j));
			nextLoc++;
		}
	}
	for(int j = 0; j < int(thresholds.size()); j++) {
		EigenMatrixAdaptor egOutThr(standardThresholds);
		omxThresholdColumn* thresh = &thresholds[j];
		for(int k = 0; k < thresh->numThresholds; k++) {
			omxSetVectorElement(vector, nextLoc, egOutThr(k, thresh->column));
			nextLoc++;
		}
	}
}

omxWLSFitFunction::~omxWLSFitFunction()
{
	if(OMX_DEBUG) {mxLog("Freeing WLS FitFunction.");}
	
	omxWLSFitFunction* owo = this;
	omxFreeMatrix(owo->observedCov);
	omxFreeMatrix(owo->observedMeans);
	omxFreeMatrix(owo->weights);
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
	
	int vectorSize = 0;
	
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
	cov = omxCreateCopyOfMatrix(omxDataCovariance(dataMat), currentState);
	means = omxCreateCopyOfMatrix(omxDataMeans(dataMat), currentState);
	weights = omxCreateCopyOfMatrix(omxDataAcov(dataMat), currentState);

	std::vector< omxThresholdColumn > &origThresh = omxDataThresholds(oo->expectation->data);
	std::vector< omxThresholdColumn > oThresh = origThresh;
	
	auto dc = oo->expectation->getDataColumns();
	if (dc.size()) {
		Eigen::VectorXi invDataColumns(dc.size()); // data -> expectation order
		for (int cx=0; cx < int(dc.size()); ++cx) {
		 	invDataColumns[dc[cx]] = cx;
		}
		//mxPrintMat("invDataColumns", invDataColumns);
		//Eigen::VectorXi invDataColumns = dc;
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pm(invDataColumns);
		EigenMatrixAdaptor Ecov(cov);
		Ecov.derived() = (pm * Ecov * pm.transpose()).eval();
		if (means) {
			EigenVectorAdaptor Emean(means);
			Emean.derived() = (pm * Emean).eval();
		}

		Eigen::MatrixXi mm(dc.size(), dc.size());
		for (int cx=0, en=0; cx < dc.size(); ++cx) {
			for (int rx=cx; rx < dc.size(); ++rx) {
				mm(rx,cx) = en;
				en += 1;
			}
		}
		mm = mm.selfadjointView<Eigen::Lower>();
		mm = (pm * mm * pm.transpose()).eval();
		//mxPrintMat("mm", mm);

		Eigen::VectorXi tstart(origThresh.size() + 1);
		tstart[0] = 0;
		int totalThresholds = 0;
		for (int tx=0; tx < int(origThresh.size()); ++tx) {
			totalThresholds += origThresh[tx].numThresholds;
			tstart[tx+1] = totalThresholds;
		}

		int wpermSize = triangleLoc1(dc.size()) + totalThresholds;
		if (means) wpermSize += dc.size();
		Eigen::VectorXi wperm(wpermSize);

		for (int cx=0, en=0; cx < dc.size(); ++cx) {
			for (int rx=cx; rx < dc.size(); ++rx) {
				wperm[en] = mm(rx,cx);
				en += 1;
			}
		}

		if (means) {
			wperm.segment(triangleLoc1(dc.size()), dc.size()) = dc.array() + triangleLoc1(dc.size());
		}

		std::vector<int> newOrder;
		newOrder.reserve(origThresh.size());
		for (int xx=0; xx < int(origThresh.size()); ++xx) newOrder.push_back(xx);

		std::sort(newOrder.begin(), newOrder.end(),
			  [&](const int &a, const int &b) -> bool
			  { return invDataColumns[origThresh[a].dColumn] < invDataColumns[origThresh[b].dColumn]; });

		//for (auto &order : newOrder) mxLog("new order %d lev %d", order, origThresh[order].numThresholds);

		int thStart = triangleLoc1(dc.size());
		if (means) thStart += dc.size();
		for (int t1=0, dest=0; t1 < int(newOrder.size()); ++t1) {
			int oldIndex = newOrder[t1];
			auto &th = oThresh[oldIndex];
			for (int t2=0; t2 < th.numThresholds; ++t2) {
				wperm[thStart + dest] = thStart + tstart[oldIndex] + t2;
				dest += 1;
			}
		}

		for (auto &th : oThresh) th.dColumn = invDataColumns[th.dColumn];
		std::sort(oThresh.begin(), oThresh.end(),
			  [](const omxThresholdColumn &a, const omxThresholdColumn &b) -> bool
			  { return a.dColumn < b.dColumn; });

		//mxPrintMat("wperm", wperm);
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> wpm(wperm);
		EigenMatrixAdaptor Eweights(weights);
		Eweights.derived() = (wpm.transpose() * Eweights * wpm).eval();
		//mxPrintMat("ew", Eweights);
	}

	newObj->observedCov = cov;
	newObj->observedMeans = means;
	newObj->n = omxDataNumObs(dataMat);
	
	auto &eThresh = oo->expectation->getThresholdInfo();

	if (eThresh.size() && !means) {
		omxRaiseError("Means are required when the data include ordinal measurements");
		return;
	}

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
	if (OMX_DEBUG) {
		mxLog("expected thresholds:");
		for (auto &th : eThresh) { th.log(); }
		mxLog("observed thresholds:");
		for (auto &th : oThresh) { th.log(); }
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
	
	flattenDataToVector(newObj->observedCov, newObj->observedMeans, obsThresholdsMat, oThresh, newObj->observedFlattened);
	flattenDataToVector(newObj->expectedCov, newObj->expectedMeans, oo->expectation->thresholdsMat,
				eThresh, newObj->expectedFlattened);
}
