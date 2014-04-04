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
 */

#include "omxAlgebraFunctions.h"
#include "omxExpectation.h"
#include "omxFIMLFitFunction.h"
#include "omxRAMExpectation.h"
#include "omxBuffer.h"
#include "matrix.h"

static const double MIN_VARIANCE = 1e-6;

struct MLFitState {

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

    // Expectation Storage;
    omxMatrix** dSigma;         // dSigma/dTheta
    omxMatrix** dMu;            // dMu/dTheta
    omxMatrix* Mu;
    omxMatrix* Ms;
    omxMatrix* X;
    omxMatrix* Y;

	// fisher information
	int haveLatentMap;
	std::vector<int> latentMap;
	std::vector<HessianBlock> lhBlocks;
};

static void omxDestroyMLFitFunction(omxFitFunction *oo) {

	if(OMX_DEBUG) {mxLog("Freeing ML Fit Function.");}
	MLFitState* omlo = ((MLFitState*)oo->argStruct);

	omxFreeMatrix(omlo->localCov);
	omxFreeMatrix(omlo->localProd);
	omxFreeMatrix(omlo->P);
	omxFreeMatrix(omlo->C);
	omxFreeMatrix(omlo->I);

	delete omlo;
}

static void calcExtraLikelihoods(omxFitFunction *oo, double *saturated_out, double *independence_out)
{
	MLFitState *state = (MLFitState*) oo->argStruct;
	double det = 0.0;
	omxMatrix* cov = state->observedCov;
	int ncols = state->observedCov->cols;
    
	*saturated_out = (state->logDetObserved + ncols) * (state->n - 1);

	// Independence model assumes all-zero manifest covariances.
	// (det(expected) + tr(observed * expected^-1)) * (n - 1);
	// expected is the diagonal of the observed.  Inverse expected is 1/each diagonal value.
	// Therefore the diagonal elements of observed * expected^-1 are each 1.
	// So the trace of the matrix is the same as the number of columns.
	// The determinant of a diagonal matrix is the product of the diagonal elements.
	// Since these are the same in the expected as in the observed, we can get 'em direct.

	for(int i = 0; i < ncols; i++) {
		// We sum logs instead of logging the product.
		det += log(omxMatrixElement(cov, i, i));
	}
	if(OMX_DEBUG) { omxPrint(cov, "Observed:"); }

	*independence_out = (ncols + det) * (state->n - 1);
}

static void addOutput(omxFitFunction *oo, MxRList *out)
{
	// DEPRECATED, use omxPopulateMLAttributes
	double saturated_out;
	double independence_out;
	calcExtraLikelihoods(oo, &saturated_out, &independence_out);
	out->add("SaturatedLikelihood", Rf_ScalarReal(saturated_out));
	out->add("IndependenceLikelihood", Rf_ScalarReal(independence_out));
}

static void mvnFit(omxFitFunction *oo, FitContext *fc)
{
	double sum = 0.0, det = 0.0;
	char u = 'U';
	char r = 'R';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	double minusoned = -1.0;
	int onei = 1;
	double fmean = 0.0;

	MLFitState *omo = ((MLFitState*)oo->argStruct);
	omxMatrix *scov, *smeans, *cov, *means, *localCov, *localProd, *P, *C;
	scov 		= omo->observedCov;
	smeans		= omo->observedMeans;
	cov			= omo->expectedCov;
	means 		= omo->expectedMeans;
	localCov 	= omo->localCov;
	localProd 	= omo->localProd;
	P		 	= omo->P;
	C		 	= omo->C;
	double n 	= omo->n;
	double Q	= omo->logDetObserved;

	omxCopyMatrix(localCov, cov);				// But expected cov is destroyed in inversion

	if(OMX_DEBUG_ALGEBRA) {
		omxPrint(scov, "Observed Covariance is");
		omxPrint(localCov, "Implied Covariance Is");
		omxPrint(cov, "Original Covariance Is");
	}

	/* Calculate |expected| */

//	F77_CALL(dgetrf)(&(localCov->cols), &(localCov->rows), localCov->data, &(localCov->cols), ipiv, &info);
	F77_CALL(dpotrf)(&u, &(localCov->cols), localCov->data, &(localCov->cols), &info);

	if(OMX_DEBUG_ALGEBRA) { mxLog("Info on LU Decomp: %d", info);}
	if(info > 0) {
		oo->matrix->data[0] = NA_REAL;
		if (fc) fc->recordIterationError("Expected covariance matrix is non-positive-definite");
		return;
	}

	//det = log(det)	// TVO: changed multiplication of det to first log and the summing up; this line should just set det to zero.
	for(info = 0; info < localCov->cols; info++) { 	    	// |cov| is the square of the product of the diagonal elements of U from the LU factorization.
		det += log(fabs(localCov->data[info+localCov->rows*info])); // TVO: changed * to + and added fabs command
	}
	det *= 2.0;		// TVO: instead of det *= det;

	if(OMX_DEBUG_ALGEBRA) { mxLog("Determinant of Expected Cov: %f", exp(det)); }
	// TVO: removed det = log(fabs(det))
	if(OMX_DEBUG_ALGEBRA) { mxLog("Log of Determinant of Expected Cov: %f", det); }

	/* Calculate Expected^(-1) */
//	F77_CALL(dgetri)(&(localCov->rows), localCov->data, &(localCov->cols), ipiv, work, lwork, &info);
	F77_CALL(dpotri)(&u, &(localCov->rows), localCov->data, &(localCov->cols), &info);
	if(OMX_DEBUG_ALGEBRA) { mxLog("Info on Invert: %d", info); }

	if(OMX_DEBUG_ALGEBRA) {omxPrint(cov, "Expected Covariance Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localCov, "Inverted Matrix:");}

	/* Calculate C = Observed * expected^(-1) */

	// Stop gcc from issuing a Rf_warning
	int majority = *(scov->majority) == 'n' ? scov->rows : scov->cols;

	/*  TODO:  Make sure leading edges are being appropriately calculated, and sub them back into this */
	F77_CALL(dsymm)(&r, &u, &(localCov->rows), &(scov->cols),
					&oned, localCov->data, &(majority),
 					scov->data, &(majority),
					&zerod, localProd->data, &(localProd->leading));

    /* And get the trace of the result */

	for(info = 0; info < localCov->cols; info++) {
		sum += localProd->data[info*localCov->cols + info];
	}

//	for(info = 0; info < (localCov->cols * localCov->rows); info++) {
//		sum += localCov->data[info] * scov->data[info];
//	}

	if(OMX_DEBUG_ALGEBRA) {omxPrint(scov, "Observed Covariance Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localCov, "Inverse Matrix:");}
	if(OMX_DEBUG_ALGEBRA) {omxPrint(localProd, "Product Matrix:");}

	if(means != NULL) {
		if(OMX_DEBUG_ALGEBRA) { mxLog("Means Likelihood Calculation"); }
		omxRecompute(means);
		omxCopyMatrix(P, means);
		// P = means - smeans
		if(OMX_DEBUG_ALGEBRA) {omxPrint(means, "means");}
		if(OMX_DEBUG_ALGEBRA) {omxPrint(smeans, "smeans");}
		F77_CALL(daxpy)(&(smeans->cols), &minusoned, smeans->data, &onei, P->data, &onei);
		if(OMX_DEBUG_ALGEBRA) {omxPrint(P, "means - smeans");}
		// C = P * Cov
		F77_CALL(dsymv)(&u, &(localCov->rows), &oned, localCov->data, &(localCov->leading), P->data, &onei, &zerod, C->data, &onei);
		// P = C * P'
		fmean = F77_CALL(ddot)(&(C->cols), P->data, &onei, C->data, &onei);

		if(OMX_DEBUG_ALGEBRA) { mxLog("Mean contribution to likelihood is %f per row.", fmean); }
		if(fmean < 0.0) fmean = 0.0;
	}

	oo->matrix->data[0] = (sum + det) * (n - 1) + fmean * (n);

	if(OMX_DEBUG) { mxLog("MLFitFunction value comes to: %f (Chisq: %f).", oo->matrix->data[0], (sum + det) - Q - cov->cols); }
}

static void buildLatentParamMap(omxFitFunction* oo, FitContext *fc)
{
	FreeVarGroup *fvg = fc->varGroup;
	MLFitState *state = (MLFitState*) oo->argStruct;
	std::vector<int> &latentMap = state->latentMap;
	int meanNum = ~state->expectedMeans->matrixNumber;
	int covNum = ~state->expectedCov->matrixNumber;
	int maxAbilities = state->expectedCov->rows;
	int numLatents = maxAbilities + triangleLoc1(maxAbilities);

	if (state->haveLatentMap == fvg->id[0]) return;
	if (0) mxLog("%s: rebuild latent parameter map for var group %d",
		     oo->matrix->name, fvg->id[0]); // TODO add runtime verbose setting

	latentMap.assign(numLatents + triangleLoc1(numLatents), -1);

	state->lhBlocks.clear();
	state->lhBlocks.resize(2);  // 0=mean, 1=cov

	int numParam = int(fvg->vars.size());
	for (int px=0; px < numParam; px++) {
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = loc->matrix;
			if (matNum == meanNum) {
				latentMap[loc->row + loc->col] = px;
				state->lhBlocks[0].vars.push_back(px);
			} else if (matNum == covNum) {
				int a1 = loc->row;
				int a2 = loc->col;
				if (a1 < a2) std::swap(a1, a2);
				int cell = maxAbilities + triangleLoc1(a1) + a2;
				if (latentMap[cell] == -1) {
					latentMap[cell] = px;
					state->lhBlocks[1].vars.push_back(px);

					if (a1 == a2 && fv->lbound == NEG_INF) {
						fv->lbound = MIN_VARIANCE;  // variance must be positive
						if (fc->est[px] < fv->lbound) {
							Rf_error("Starting value for variance %s is not positive", fv->name);
						}
					}
				} else if (latentMap[cell] != px) {
					// doesn't detect similar problems in multigroup constraints TODO
					Rf_error("Covariance matrix must be constrained to preserve symmetry");
				}
			}
		}
	}
	state->haveLatentMap = fvg->id[0];

	for (int p1=0; p1 < maxAbilities; p1++) {
		HessianBlock &hb = state->lhBlocks[0];
		int at1 = latentMap[p1];
		if (at1 < 0) continue;
		int hb1 = std::lower_bound(hb.vars.begin(), hb.vars.end(), at1) - hb.vars.begin();

		for (int p2=0; p2 <= p1; p2++) {
			int at2 = latentMap[p2];
			if (at2 < 0) continue;
			int hb2 = std::lower_bound(hb.vars.begin(), hb.vars.end(), at2) - hb.vars.begin();

			if (hb1 < hb2) std::swap(hb1, hb2);
			int at = numLatents + triangleLoc1(p1) + p2;
			latentMap[at] = hb1 * hb.vars.size() + hb2;
		}
	}

	for (int p1=maxAbilities; p1 < numLatents; p1++) {
		HessianBlock &hb = state->lhBlocks[1];
		int at1 = latentMap[p1];
		if (at1 < 0) continue;
		int hb1 = std::lower_bound(hb.vars.begin(), hb.vars.end(), at1) - hb.vars.begin();

		for (int p2=maxAbilities; p2 <= p1; p2++) {
			int at2 = latentMap[p2];
			if (at2 < 0) continue;
			int hb2 = std::lower_bound(hb.vars.begin(), hb.vars.end(), at2) - hb.vars.begin();

			if (hb1 < hb2) std::swap(hb1, hb2);
			int at = numLatents + triangleLoc1(p1) + p2;
			latentMap[at] = hb1 * hb.vars.size() + hb2;
		}
	}
}

static void mvnInfo(omxFitFunction *oo, FitContext *fc)
{
	buildLatentParamMap(oo, fc);

	const double Scale = Global->llScale;
	MLFitState *state = (MLFitState*) oo->argStruct;
	const int maxAbilities = state->expectedCov->rows;
	omxMatrix *cov = state->expectedCov;
	const double numObs = state->n;
	int numLatents = maxAbilities + triangleLoc1(maxAbilities);
	std::vector<int> &latentMap = state->latentMap;

	if (0) mxLog("%s: latentHessian", oo->matrix->name);

	omxBuffer<double> icovBuffer(maxAbilities * maxAbilities);
	memcpy(icovBuffer.data(), cov->data, sizeof(double) * maxAbilities * maxAbilities);
	Matrix icovMat(icovBuffer.data(), maxAbilities, maxAbilities);
	int info = InvertSymmetricPosDef(icovMat, 'U');
	if (info) return;

	for (int m1=0; m1 < maxAbilities; ++m1) {
		for (int m2=0; m2 < m1; ++m2) {
			icovBuffer[m2 * maxAbilities + m1] = icovBuffer[m1 * maxAbilities + m2];
		}
	}

	{
		HessianBlock *hb = state->lhBlocks[0].clone();

		int px=numLatents;
		for (int m1=0; m1 < maxAbilities; ++m1) {
			for (int m2=0; m2 <= m1; ++m2) {
				int to = latentMap[px];
				++px;
				if (to < 0) continue;
				hb->mat.data()[to] = -Scale * numObs * icovBuffer[m1 * maxAbilities + m2];
			}
		}
		fc->queue(hb);
	}

	HessianBlock *hb = state->lhBlocks[1].clone();

	std::vector<double> term1(maxAbilities * maxAbilities);
	std::vector<double> term2(maxAbilities * maxAbilities);

	int f1=0;
	for (int r1=0; r1 < maxAbilities; ++r1) {
		for (int c1=0; c1 <= r1; ++c1) {
			memcpy(term1.data()      + c1 * maxAbilities,
			       icovBuffer.data() + r1 * maxAbilities, maxAbilities * sizeof(double));
			if (r1 != c1) {
				memcpy(term1.data()      + r1 * maxAbilities,
				       icovBuffer.data() + c1 * maxAbilities, maxAbilities * sizeof(double));
			}
			int f2 = f1;
			for (int r2=r1; r2 < maxAbilities; ++r2) {
				for (int c2 = (r1==r2? c1 : 0); c2 <= r2; ++c2) {
					int to = latentMap[numLatents + triangleLoc1(f2 + maxAbilities) + f1 + maxAbilities];
					++f2;
					if (to < 0) continue;

					memcpy(term2.data()      + c2 * maxAbilities,
					       icovBuffer.data() + r2 * maxAbilities, maxAbilities * sizeof(double));
					if (r2 != c2) {
						memcpy(term2.data()      + r2 * maxAbilities,
						       icovBuffer.data() + c2 * maxAbilities, maxAbilities * sizeof(double));
					}

					double tr = 0;
					for (int d1=0; d1 < maxAbilities; ++d1) {
						for (int d2=0; d2 < maxAbilities; ++d2) {
							tr += term1[d2 * maxAbilities + d1] * term2[d1 * maxAbilities + d2];
						}
					}

					// Simulation suggests the sample size should be
					// numObs-2 but this is tedious to accomodate
					// when there are parameter equality constraints. Whether
					// the sample size is adjusted or not seems to make
					// no detectable difference in tests.
					hb->mat.data()[to] = Scale * numObs * -.5 * tr;

					OMXZERO(term2.data() + c2 * maxAbilities, maxAbilities);
					if (c2 != r2) OMXZERO(term2.data() + r2 * maxAbilities, maxAbilities);
				}
			}
			OMXZERO(term1.data() + c1 * maxAbilities, maxAbilities);
			if (c1 != r1) OMXZERO(term1.data() + r1 * maxAbilities, maxAbilities);
			++f1;
		}
	}
	fc->queue(hb);
}

static void omxCallMLFitFunction(omxFitFunction *oo, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_PREOPTIMIZE)) return;

	omxExpectation* expectation = oo->expectation;
	omxExpectationCompute(expectation, NULL);

	if ((want & FF_COMPUTE_INFO) && strcmp(expectation->expType, "MxExpectationNormal")==0) {
		if (fc->infoMethod != INFO_METHOD_HESSIAN) {
			omxRaiseErrorf(globalState, "Information matrix approximation method %d is not available",
				       fc->infoMethod);
			return;
		}
		mvnInfo(oo, fc);
	}

	if (want & FF_COMPUTE_FIT) {
		mvnFit(oo, fc);
	}
}

static void omxPopulateMLAttributes(omxFitFunction *oo, SEXP algebra) {
    if(OMX_DEBUG) { mxLog("Populating ML Attributes."); }

	MLFitState *argStruct = ((MLFitState*)oo->argStruct);
	omxMatrix *expCovInt = argStruct->expectedCov;	    		// Expected covariance
	omxMatrix *expMeanInt = argStruct->expectedMeans;			// Expected means

	SEXP expCovExt, expMeanExt;
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

	Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
	
	double saturated_out;
	double independence_out;
	calcExtraLikelihoods(oo, &saturated_out, &independence_out);
	Rf_setAttrib(algebra, Rf_install("SaturatedLikelihood"), Rf_ScalarReal(saturated_out));
	Rf_setAttrib(algebra, Rf_install("IndependenceLikelihood"), Rf_ScalarReal(independence_out));
}

void omxInitMLFitFunction(omxFitFunction* oo)
{
	if (!oo->expectation) { Rf_error("%s requires an expectation", oo->fitType); }

	omxExpectation *expectation = oo->expectation;
	if (strcmp(expectation->expType, "MxExpectationBA81")==0) {
		omxInitFitFunctionBA81(oo);
		return;
	}

	if(OMX_DEBUG) { mxLog("Initializing ML fit function."); }

	int info = 0;
	double det = 1.0;
	char u = 'U';
	
	oo->computeFun = omxCallMLFitFunction;
	oo->destructFun = omxDestroyMLFitFunction;
	oo->addOutput = addOutput;
	oo->populateAttrFun = omxPopulateMLAttributes;

	omxData* dataMat = oo->expectation->data;

	if(!(dataMat == NULL) && strncmp(omxDataType(dataMat), "cov", 3) != 0 && strncmp(omxDataType(dataMat), "cor", 3) != 0) {
		if(strncmp(omxDataType(dataMat), "raw", 3) == 0) {
			if(OMX_DEBUG) { mxLog("Raw Data: Converting to FIML."); }
			omxInitFIMLFitFunction(oo);
			return;
		}
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "ML FitFunction unable to handle data type %s.\n", omxDataType(dataMat));
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
		if(OMX_DEBUG) { mxLog("ML FitFunction unable to handle data type %s.  Aborting.", omxDataType(dataMat)); }
		return;
	}

	MLFitState *newObj = new MLFitState;
	newObj->haveLatentMap = FREEVARGROUP_INVALID;
	oo->argStruct = (void*)newObj;

	if(OMX_DEBUG) { mxLog("Processing Observed Covariance."); }
	newObj->observedCov = omxDataMatrix(dataMat, NULL);
	if(OMX_DEBUG) { mxLog("Processing Observed Means."); }
	newObj->observedMeans = omxDataMeans(dataMat, NULL, NULL);
	if(OMX_DEBUG && newObj->observedMeans == NULL) { mxLog("ML: No Observed Means."); }
	if(OMX_DEBUG) { mxLog("Processing n."); }
	newObj->n = omxDataNumObs(dataMat);

	newObj->expectedCov = omxGetExpectationComponent(oo->expectation, oo, "cov");
	newObj->expectedMeans = omxGetExpectationComponent(oo->expectation, oo, "means");

	if(newObj->expectedCov == NULL) {
		omxRaiseError(0, 0,
			"Developer Error in ML-based fit function object: ML's expectation must specify a model-implied covariance matrix.\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.");
		return;
	}

	// Error Checking: Observed/Expected means must agree.  
	// ^ is XOR: true when one is false and the other is not.
	if((newObj->expectedMeans == NULL) ^ (newObj->observedMeans == NULL)) {
		if(newObj->expectedMeans != NULL) {
			omxRaiseError(0,0,
				"Observed means not detected, but an expected means matrix was specified.\n  If you provide observed means, you must specify a model for the means.\n");
			return;
		} else {
			omxRaiseError(0,0,
				"Observed means were provided, but an expected means matrix was not specified.\n  If you  wish to model the means, you must provide observed means.\n");
			return;	        
		}
	}

	/* Temporary storage for calculation */
	int rows = newObj->observedCov->rows;
	int cols = newObj->observedCov->cols;
	newObj->localCov = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->localProd = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->P = omxInitMatrix(NULL, 1, cols, TRUE, oo->matrix->currentState);
	newObj->C = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
	newObj->I = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);

	for(int i = 0; i < rows; i++) omxSetMatrixElement(newObj->I, i, i, 1.0);

	omxCopyMatrix(newObj->localCov, newObj->observedCov);

	newObj->lwork = newObj->expectedCov->rows;
	newObj->work = (double*)R_alloc(newObj->lwork, sizeof(double));

	// TODO: Determine where the saturated model computation should go.

	F77_CALL(dpotrf)(&u, &(newObj->localCov->cols), newObj->localCov->data, &(newObj->localCov->cols), &info);

	if(OMX_DEBUG) { mxLog("Info on LU Decomp: %d", info); }
	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Observed Covariance Matrix is non-positive-definite.\n");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
		return;
	}
	for(info = 0; info < newObj->localCov->cols; info++) {
		det *= omxMatrixElement(newObj->localCov, info, info);
	}
	det *= det;					// Product of squares.

	if(OMX_DEBUG) { mxLog("Determinant of Observed Cov: %f", det); }
	newObj->logDetObserved = log(det);
	if(OMX_DEBUG) { mxLog("Log Determinant of Observed Cov: %f", newObj->logDetObserved); }

	omxCopyMatrix(newObj->localCov, newObj->expectedCov);
}

static void omxSetMLFitFunctionGradientComponents(omxFitFunction* oo, void (*derivativeFun)(omxFitFunction*, omxMatrix**, omxMatrix**, int*)) {
    if(OMX_DEBUG) { mxLog("Setting up gradient component function for ML FitFunction."); }
    if(!strncmp("omxFIMLFitFunction", oo->fitType, 16)) {
        if(OMX_DEBUG) { mxLog("FIML FitFunction gradients not yet implemented. Skipping."); }
        return; // ERROR:NYI.
    }
    
    if(derivativeFun == NULL) {
        char Rf_errorstr[250];
        sprintf(Rf_errorstr, "Programmer Rf_error: ML gradient components given NULL gradient function.");
        omxRaiseError(oo->matrix->currentState, -2, Rf_errorstr);
        return;
    }
    
    MLFitState *omo = ((MLFitState*) oo->argStruct);
    int rows = omo->observedCov->rows;
    int cols = omo->observedCov->cols;
    size_t nFreeVars = oo->freeVarGroup->vars.size();
            
    omo->X  = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
    omo->Y  = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
    omo->Ms = omxInitMatrix(NULL, 1, cols, TRUE, oo->matrix->currentState);
    omo->Mu = omxInitMatrix(NULL, 1, cols, TRUE, oo->matrix->currentState);
    omo->dSigma = (omxMatrix**) R_alloc(nFreeVars, sizeof(omxMatrix*));
    omo->dMu = (omxMatrix**) R_alloc(nFreeVars, sizeof(omxMatrix*));
    for(size_t i = 0; i < nFreeVars; i++) {
        omo->dSigma[i] = omxInitMatrix(NULL, rows, cols, TRUE, oo->matrix->currentState);
        omo->dMu[i] = omxInitMatrix(NULL, rows, 1, TRUE, oo->matrix->currentState);
    }
    //oo->gradientFun = omxCalculateMLGradient; TODO
}

static void omxCalculateMLGradient(omxFitFunction* oo, double* gradient) {

    if(OMX_DEBUG) { mxLog("Beginning ML Gradient Calculation."); }
    // mxLog("Beginning ML Gradient Calculation, Iteration %d.%d (%d)\n", 
        // oo->matrix->currentState->majorIteration, oo->matrix->currentState->minorIteration,
        // oo->matrix->currentState->computeCount); //:::DEBUG:::
    // 1) Calculate current Expected Covariance
    // 2) Calculate eCov, the Inverse Expected Covariance matrix 
    // 3) Calculate C = I - eCov D, where D is the observed covariance matrix
    // 4) Calculate b = M - [observed M]
    // 5) For each location in locs:
    //   gradient[loc] = tr(eCov^-1 %*% dEdt %*% C) - (b^T %*% eCov^-1 %*% dEdt + 2 dMdt^T))eCov^-1 b)

    MLFitState *omo = ((MLFitState*)oo->argStruct);
    
    /* Locals for readability.  Compiler should cut through this. */
    omxMatrix *scov         = omo->observedCov;
    omxMatrix *smeans       = omo->observedMeans;
    omxMatrix *cov          = omo->expectedCov;
    omxMatrix *M            = omo->expectedMeans;
    omxMatrix *eCov         = omo->localCov;        // TODO: Maybe need to avoid reusing these
    omxMatrix *I            = omo->I;
    omxMatrix *C            = omo->C;
    omxMatrix *X            = omo->X;
    omxMatrix *Y            = omo->Y;
    omxMatrix *Mu           = omo->Mu;
    omxMatrix *Ms           = omo->Ms;
    omxMatrix *P            = omo->P;
    double n                = omo->n;
    omxMatrix** dSigmas     = omo->dSigma;
    omxMatrix** dMus        = omo->dMu;
    
    size_t gradientSize = oo->freeVarGroup->vars.size();
    
    char u = 'U';
    int info;
    double minusoned = -1.0;
    int onei = 1;
    int status[gradientSize];
    int nLocs = gradientSize;
    
    // Calculate current FitFunction values
    // We can safely assume this has been done
    // omxFitFunctionCompute(oo);
    
    // Calculate current eCov
    
    omxCopyMatrix(eCov, cov);				// But expected cov is destroyed in inversion
    
    F77_CALL(dpotrf)(&u, &(eCov->cols), eCov->data, &(eCov->cols), &info);

    if(OMX_DEBUG_ALGEBRA) { mxLog("Info on LU Decomp: %d", info);}
    if(info > 0) {
	    char *errstr = (char*) calloc(250, sizeof(char));
        sprintf(errstr, "Expected covariance matrix is non-positive-definite");
        if(oo->matrix->currentState->computeCount <= 0) {
            strncat(errstr, " at starting values", 20);
        }
        strncat(errstr, ".\n", 3);
        omxRaiseError(oo->matrix->currentState, -1, errstr);                        // Raise Rf_error
        free(errstr);
        return;                                                                     // Leave output untouched
    }
    
    F77_CALL(dpotri)(&u, &(eCov->rows), eCov->data, &(eCov->cols), &info);
    if(info > 0) {
	    char *errstr = (char*) calloc(250, sizeof(char));
        sprintf(errstr, "Expected covariance matrix is not invertible");
        if(oo->matrix->currentState->computeCount <= 0) {
            strncat(errstr, " at starting values", 20);
        }
        strncat(errstr, ".\n", 3);
        omxRaiseError(oo->matrix->currentState, -1, errstr);                        // Raise Rf_error
        free(errstr);
        return;
    }
    // Calculate P = expected means - observed means
    if(M != NULL) {
        omxCopyMatrix(P, smeans);
    	F77_CALL(daxpy)(&(smeans->cols), &minusoned, M->data, &onei, P->data, &onei);
    }
	
	// Reset C and Calculate C = I - eCov * oCov
    omxCopyMatrix(C, I);
    omxDSYMM(TRUE, -1.0, eCov, scov, 1.0, C);
    
    // For means, calculate Ms = eCov-1 %*% P
    if(M != NULL)
        omxDSYMM(FALSE, 1.0, eCov, P, 0.0, Ms);
    
    // Calculate parameter-level derivatives
    // TODO: Parallelize Here.

    if(OMX_DEBUG)  { mxLog("Calling component function."); }
    // omo->derivativeFun(oo, dSigmas, dMus, status);
    
    for(int currentLoc = 0; currentLoc < nLocs; currentLoc++) {
        double meanInfluence, covInfluence;
        if(status[currentLoc] < 0) continue;  // Failure in computation--skip.
        //   gradient[loc] = tr(eCov^-1 %*% dEdt %*% C) - 
        //    (b^T %*% eCov^-1 %*% dEdt + 2 dMdt^T))eCov^-1 b)
        // omxDGEMM(FALSE, FALSE, 1.0, dSigmas[currentLoc], C, 0.0, Y);
        omxDSYMM(TRUE, 1.0, dSigmas[currentLoc], C, 0.0, Y);
        omxDSYMM(TRUE, 1.0, eCov, Y, 0.0, X);
        gradient[currentLoc] = 0;
        covInfluence = 0.0;
        for(int i = 0; i < eCov->cols; i++) 
            covInfluence += omxMatrixElement(X, i, i);
        if(M != NULL) {
            omxCopyMatrix(Mu, dMus[currentLoc]);
            omxDSYMV(1.0, dSigmas[currentLoc], Ms, 2.0, Mu);
            meanInfluence = F77_CALL(ddot)(&(eCov->cols), Mu->data, &onei, Ms->data, &onei);
        } else {
            meanInfluence = 0;
        }
        gradient[currentLoc] = (covInfluence * (n-1)) - (meanInfluence * n);
        if(OMX_DEBUG) { 
            mxLog("Calculation for Gradient value %d: Cov: %3.9f, Mean: %3.9f, total: %3.9f",
            currentLoc, covInfluence, meanInfluence, gradient[currentLoc]); 
        }
    }
}
