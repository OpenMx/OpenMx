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

#include "omxDefines.h"

#include <utility>
#include <stan/math/mix/mat.hpp>
#include "multi_normal_sufficient.hpp"

#include "omxExpectation.h"
#include "RAMInternal.h"
#include "matrix.h"
#include "Compute.h"
#include "EnableWarnings.h"

static const double MIN_VARIANCE = 1e-6;

struct MLFitState : omxFitFunction {
	bool copiedData;
	omxMatrix* observedCov;
	omxMatrix* observedMeans;
	omxMatrix* expectedCov;
	omxMatrix* expectedMeans;

	double n;
	double logDetObserved;

	MLFitState() : copiedData(false) {};
	virtual ~MLFitState();
	virtual omxFitFunction *initMorph();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
	virtual void addOutput(MxRList *out);
};

MLFitState::~MLFitState()
{
	if(OMX_DEBUG) {mxLog("Freeing ML Fit Function.");}
	MLFitState* omlo = this;
	if (omlo->copiedData) {
		omxFreeMatrix(omlo->observedCov);
		omxFreeMatrix(omlo->observedMeans);
	}
}

static void calcExtraLikelihoods(omxFitFunction *oo, double *saturated_out, double *independence_out)
{
	MLFitState *state = (MLFitState*) oo;
	double det = 0.0;
	omxMatrix* cov = state->observedCov;
	int ncols = state->observedCov->cols;

	*saturated_out = (state->logDetObserved * state->n) + ncols * (state->n - 1);

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

	*independence_out = ncols * (state->n - 1) + det * state->n;
}

void MLFitState::addOutput(MxRList *out)
{
	auto *oo = this;
	// DEPRECATED, use omxPopulateMLAttributes
	if(OMX_DEBUG) { mxLog("Deprecated ML Attribute Population Code Running."); }
	double saturated_out;
	double independence_out;
	omxData* dataMat = oo->expectation->data;
	if (!strEQ(omxDataType(dataMat), "raw")) {
		calcExtraLikelihoods(oo, &saturated_out, &independence_out);
		out->add("SaturatedLikelihood", Rf_ScalarReal(saturated_out));
		out->add("IndependenceLikelihood", Rf_ScalarReal(independence_out));
	}
}

struct multi_normal_deriv {
	FitContext *fc;
	std::vector<bool> &fvMask;
	MLFitState *omo;

	multi_normal_deriv(FitContext *_fc, std::vector<bool> &_fvMask, MLFitState *_omo) :
		fc(_fc), fvMask(_fvMask), omo(_omo) {};

	template <typename T>
	T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
		EigenMatrixAdaptor obCovAdapter(omo->observedCov);
		Eigen::MatrixXd obCov = obCovAdapter;
		EigenMatrixAdaptor exCovAdapter(omo->expectedCov);
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> exCov = exCovAdapter.cast<T>();
		
		Eigen::VectorXd obMeans(obCov.rows());
		Eigen::Matrix<T, Eigen::Dynamic, 1> exMeans;
		if (omo->observedMeans) {
			EigenVectorAdaptor obMeansAdapter(omo->observedMeans);
			obMeans = obMeansAdapter;
			EigenVectorAdaptor exMeansAdapter(omo->expectedMeans);
			exMeans = exMeansAdapter.cast<T>();
		} else {
			obMeans = Eigen::VectorXd::Zero(obCov.rows());
			exMeans.resize(obMeans.size());
			exMeans.setZero();
		}

		int xx=0;
		for (size_t fx=0; fx < fvMask.size(); ++fx) {
			if (!fvMask[fx]) continue;
			omxFreeVar *fv = fc->varGroup->vars[fx];
			for (size_t lx=0; lx < fv->locations.size(); ++lx) {
				omxFreeVarLocation &loc = fv->locations[lx];
				int matNum = ~loc.matrix;
				if (matNum == omo->expectedCov->matrixNumber) {
					exCov(loc.row, loc.col) = x[xx];
				} else if (omo->expectedMeans && matNum == omo->expectedMeans->matrixNumber) {
					exMeans(loc.row + loc.col) = x[xx];
				}
			}
			++xx;
		}

		return stan::prob::multi_normal_sufficient_log<true>(omo->n, obMeans, obCov, exMeans, exCov);
	}
};

void MLFitState::compute(int want, FitContext *fc)
{
	auto *oo = this;
	const double Scale = Global->llScale;
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;

	omxExpectationCompute(fc, expectation, NULL);

	if ((want & FF_COMPUTE_FIT) &&
	    !(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN | FF_COMPUTE_INFO))) {
		// works for any multivariate normal expectation (e.g. vanilla, RAM, LISREL, etc)

		MLFitState *omo = (MLFitState*) oo;
		EigenMatrixAdaptor obCovAdapter(omo->observedCov);
		Eigen::MatrixXd obCov = obCovAdapter;
		EigenMatrixAdaptor exCovAdapter(omo->expectedCov);
		Eigen::MatrixXd exCov = exCovAdapter;
		double fit;
		try {
			if (omo->observedMeans) {
				EigenVectorAdaptor obMeansAdapter(omo->observedMeans);
				Eigen::VectorXd obMeans = obMeansAdapter;
				EigenVectorAdaptor exMeansAdapter(omo->expectedMeans);
				Eigen::VectorXd exMeans = exMeansAdapter;
				fit = stan::prob::multi_normal_sufficient_log<false>(omo->n, obMeans, obCov, exMeans, exCov);
			} else {
				Eigen::VectorXd means(obCov.rows());
				means.setZero();
				fit = stan::prob::multi_normal_sufficient_log<false>(omo->n, means, obCov, means, exCov);
			}
			using stan::math::LOG_TWO_PI;
			fit += .5 * omo->n * obCov.rows() * LOG_TWO_PI;
		} catch (const std::exception& e) {
			fit = NA_REAL;
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		} catch (...) {
			fit = NA_REAL;
			if (fc) fc->recordIterationError("%s: unknown error", oo->name());
		}
		oo->matrix->data[0] = Scale * fit;
	} else if (strEQ(expectation->name, "MxExpectationNormal") &&
		   (want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN | FF_COMPUTE_INFO))) {
		if ((want & FF_COMPUTE_INFO) && fc->infoMethod != INFO_METHOD_HESSIAN) {
			omxRaiseErrorf("Information matrix approximation method %d is not available",
				       fc->infoMethod);
			return;
		}
		MLFitState *omo = (MLFitState*) oo;
		// should forward computation to the expectation TODO

		int num_param = 0;
		std::vector<bool> fvMask(fc->numParam);
		for (int fx=0; fx < int(fc->numParam); ++fx) {
			omxFreeVar *fv = fc->varGroup->vars[fx];
			bool relevant = fv->getLocation(omo->expectedCov) != NULL;
			if (omo->expectedMeans) {
				relevant |= fv->getLocation(omo->expectedMeans) != NULL;
			}
			fvMask[fx] = relevant;
			if (relevant) ++num_param;
		}

		if (num_param == 0) {
			// if num_param == 0 then hessian() doesn't compute the fit
			want &= ~(FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN | FF_COMPUTE_INFO);
			if (want) compute(want, fc);
			return;
		}

		double init_log_prob = 0.0;
		Eigen::VectorXd init_grad = Eigen::VectorXd::Zero(num_param);
		HessianBlock *hb = new HessianBlock;

		Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(num_param);
		int cpx=0;
		for (int fx=0; fx < int(fc->numParam); ++fx) {
			if (!fvMask[fx]) continue;
			hb->vars.push_back(fx);
			cont_params[cpx++] = fc->est[fx];
		}
		hb->mat.resize(num_param, num_param);
		hb->mat.setZero();
		
		try {
			multi_normal_deriv model(fc, fvMask, omo);
			stan::math::hessian(model, cont_params, init_log_prob, init_grad, hb->mat);
		} catch (const std::exception& e) {
			init_log_prob = NA_REAL;
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		} catch (...) {
			init_log_prob = NA_REAL;
			if (fc) fc->recordIterationError("%s: unknown error", oo->name());
		}

		if (want & FF_COMPUTE_FIT) {
			oo->matrix->data[0] = Scale * init_log_prob;
		}
		if (want & FF_COMPUTE_GRADIENT) {
			int px=0;
			for (int fx=0; fx < int(fc->numParam); ++fx) {
				if (!fvMask[fx]) continue;
				fc->haveGrad[fx] = true;
				fc->gradZ[fx] += Scale * init_grad[px];
				++px;
			}
		}
		if (want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN | FF_COMPUTE_INFO)) {
			const double HScale = (want & FF_COMPUTE_INFO)? -fabs(Scale) : Scale;
			hb->mat *= HScale;
			fc->queue(hb);
		} else {
			delete hb;
		}
	} else {
		mxThrow("Not implemented");
	}
}

void MLFitState::populateAttr(SEXP algebra) {
    if(OMX_DEBUG) { mxLog("Populating ML Attributes."); }

    auto *oo = this;
    MLFitState *argStruct = this;
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
	ProtectedSEXP Rsat(Rf_ScalarReal(saturated_out));
	Rf_setAttrib(algebra, Rf_install("SaturatedLikelihood"), Rsat);
	ProtectedSEXP Rind(Rf_ScalarReal(independence_out));
	Rf_setAttrib(algebra, Rf_install("IndependenceLikelihood"), Rind);
}

omxFitFunction *omxInitMLFitFunction()
{ return new MLFitState; }

omxFitFunction *MLFitState::initMorph()
{
	auto *oo = this;

	if (!oo->expectation) { mxThrow("%s requires an expectation", oo->fitType); }
	oo->units = FIT_UNITS_MINUS2LL;

	if (strcmp(expectation->name, "MxExpectationBA81")==0) {
		return omxChangeFitType(oo, "imxFitFunctionBA81");
	}
	
	if (strEQ(expectation->name, "MxExpectationGREML")) {
		return omxChangeFitType(oo, "imxFitFunciontGRMFIML");
	}
	
	if (strEQ(expectation->name, "MxExpectationStateSpace")) {
		return omxChangeFitType(oo, "imxFitFunciontStateSpace");
	}

	if (strEQ(expectation->name, "MxExpectationHiddenMarkov") ||
	    strEQ(expectation->name, "MxExpectationMixture")) {
		return omxChangeFitType(oo, "imxFitFunciontHiddenMarkov");
	}

	if(OMX_DEBUG) { mxLog("Initializing ML fit function."); }

	omxData* dataMat = oo->expectation->data;

	ProtectedSEXP Rfellner(R_do_slot(oo->rObj, Rf_install("fellner")));
	int wantRowwiseLikelihood = Rf_asInteger(R_do_slot(oo->rObj, Rf_install("vector")));

	bool fellnerPossible = (strEQ(omxDataType(dataMat), "raw") && expectation->numOrdinal == 0 &&
				strEQ(oo->expectation->name, "MxExpectationRAM") && !wantRowwiseLikelihood);

	if (Rf_asLogical(Rfellner) == 1 && !fellnerPossible) {
		mxThrow("%s: fellner requires raw data (have %s), "
			 "all continuous indicators (%d are ordinal), "
			 "MxExpectationRAM (have %s), and no row-wise likelihoods (want %d)",
			 oo->name(), dataMat->getType(), expectation->numOrdinal,
			expectation->name, wantRowwiseLikelihood);
	}

	if (strEQ(omxDataType(dataMat), "raw")) {
		int useFellner = Rf_asLogical(Rfellner);
		if (strEQ(oo->expectation->name, "MxExpectationRAM")) {
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			if (ram->between.size()) {
				if (useFellner == 0) {
					mxThrow("%s: fellner=TRUE is required for %s",
						 oo->name(), expectation->name);
				}
				useFellner = 1;
			}
		}

		// Cannot enable unconditionally because performance
		// suffers with models that make heavy use of defvars
		// such as continuous time models.
		//
		//if (useFellner == NA_LOGICAL && fellnerPossible) useFellner = 1;

		const char *to;
		if (useFellner == 1) {
			to = "imxFitFunctionFellner";
		} else {
			to = "imxFitFunctionFIML";
		}
		if(OMX_DEBUG) { mxLog("Raw Data: Converting from %s to %s", oo->fitType, to); }
		return omxChangeFitType(oo, to);
	}

	init();
	return this;
}

void MLFitState::init()
{
	auto *oo = this;
	auto *newObj = this;
	omxData* dataMat = expectation->data;

	if(!strEQ(omxDataType(dataMat), "cov") && !strEQ(omxDataType(dataMat), "cor")) {
		omxRaiseErrorf("ML FitFunction unable to handle data type %s", omxDataType(dataMat));
		return;
	}

	oo->canDuplicate = true;

	newObj->observedCov = omxDataCovariance(dataMat);
	newObj->observedMeans = omxDataMeans(dataMat);
	newObj->copiedData = false;

	auto dc = oo->expectation->getDataColumns();
	if (dc.size()) {
		if (dataMat->isDynamic()) mxThrow("%s: dynamic data & column reordering"
						   " is not implemented yet", name());
		newObj->copiedData = true;
		newObj->observedCov = omxCreateCopyOfMatrix(newObj->observedCov, oo->matrix->currentState);
		newObj->observedMeans = omxCreateCopyOfMatrix(newObj->observedMeans, oo->matrix->currentState);
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pm(dc);
		EigenMatrixAdaptor Ecov(newObj->observedCov);
		Ecov.derived() = (pm.transpose() * Ecov * pm).eval();
		if (newObj->observedMeans) {
			EigenVectorAdaptor Emean(newObj->observedMeans);
			Emean.derived() = (pm.transpose() * Emean).eval();
		}
	}

	if(OMX_DEBUG && newObj->observedMeans == NULL) { mxLog("ML: No Observed Means."); }

	newObj->n = omxDataNumObs(dataMat);

	newObj->expectedCov = omxGetExpectationComponent(oo->expectation, "cov");
	newObj->expectedMeans = omxGetExpectationComponent(oo->expectation, "means");

	if(newObj->expectedCov == NULL) {
		omxRaiseError("Developer Error in ML-based fit function object: ML's expectation must specify a model-implied covariance matrix.\nIf you are not developing a new expectation type, you should probably post this to the OpenMx forums.");
		return;
	}

	// Error Checking: Observed/Expected means must agree.  
	// ^ is XOR: true when one is false and the other is not.
	if((newObj->expectedMeans == NULL) ^ (newObj->observedMeans == NULL)) {
		if(newObj->expectedMeans != NULL) {
			omxRaiseError("Observed means not detected, but an expected means matrix was specified.\n  If you provide observed means, you must specify a model for the means.\n");
			return;
		} else {
			omxRaiseErrorf("%s: Observed means were provided, but an expected means matrix was not specified.\n  If you  wish to model the means, you must provide observed means.\n", oo->name());
			return;
		}
	}

	// add expectation API for derivs TODO
	if (strEQ(expectation->name, "MxExpectationNormal") &&
	    newObj->expectedCov->isSimple() &&
	    (!newObj->expectedMeans || newObj->expectedMeans->isSimple())) {
		oo->gradientAvailable = true;
		oo->hessianAvailable = true;
	}

	EigenMatrixAdaptor obCov(newObj->observedCov);
	stan::math::LDLT_factor<double,Eigen::Dynamic,Eigen::Dynamic> ldlt_obCov(obCov);
	if (!ldlt_obCov.success()) {
		omxRaiseErrorf("Observed Covariance Matrix is non-positive-definite.");
		return;
	}
	newObj->logDetObserved = log_determinant_ldlt(ldlt_obCov);
	if(OMX_DEBUG) { mxLog("Log Determinant of Observed Cov: %f", newObj->logDetObserved); }
}
