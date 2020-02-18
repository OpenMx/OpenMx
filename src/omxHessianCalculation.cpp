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

/**
 * Based on:
 *
 * Paul Gilbert and Ravi Varadhan (2012). numDeriv: Accurate Numerical Derivatives. R package
 * version 2012.9-1. http://CRAN.R-project.org/package=numDeriv
 *
 **/

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>

#include "omxDefines.h"
#include "glue.h"
#include "omxState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxNPSOLSpecific.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "CovEntrywisePar.h"
#include "EnableWarnings.h"
#include "finiteDifferences.h"

struct hess_struct {
	int probeCount;
	double* Haprox;
	double *Gcentral;
	double *Gforward;
	double *Gbackward;
	FitContext *fc;
	omxMatrix* fitMatrix;

	~hess_struct() {
		delete [] Haprox;
		delete [] Gcentral;
		delete [] Gforward;
		delete [] Gbackward;
	}
};

class omxComputeNumericDeriv : public omxCompute {
	typedef omxCompute super;
	double stepSize;
	int numIter;
	bool parallel;
	int totalProbeCount;
	int verbose;
	bool wantHessian;
	bool checkGradient;
	double *knownHessian;
	std::vector<int> khMap;

	omxMatrix *fitMat;
	double minimum;
	Eigen::ArrayXd optima;
	int numParams;
	std::vector<int> paramMap;
	double *gcentral, *gforward, *gbackward;
	double *hessian;
	struct hess_struct *hessWorkVector;
	bool recordDetail;
	SEXP detail;

	void omxPopulateHessianWork(struct hess_struct *hess_work, FitContext* fc);
	void omxEstimateHessianOnDiagonal(int i, struct hess_struct* hess_work);
	void omxEstimateHessianOffDiagonal(int i, int l, struct hess_struct* hess_work);
	void omxCalcFinalConstraintJacobian(FitContext* fc, int npar);

	struct calcHessianEntry {
		omxComputeNumericDeriv &cnd;
		const int numParams;
		double *hessian;
		struct hess_struct *hessWorkVector;
		calcHessianEntry(omxComputeNumericDeriv *_cnd) :
			cnd(*_cnd), numParams(cnd.numParams), hessian(cnd.hessian),
			hessWorkVector(cnd.hessWorkVector)
		{};
		int getNumCols() { return numParams; };
		bool isDone(int rx, int cx) {
			if (hessian) {
				return std::isfinite(hessian[cx*numParams + rx]);
			} else {
				return rx != cx;  // Only interested in gradient
			}
		}
		void onDiag(int ii) {
			int threadId = omx_absolute_thread_num();
			cnd.omxEstimateHessianOnDiagonal(ii, hessWorkVector + threadId);
		}
		void offDiag(int rx, int cx) {
			int threadId = omx_absolute_thread_num();
			cnd.omxEstimateHessianOffDiagonal(rx, cx, hessWorkVector + threadId);
		}
		void reportProgress(int numDone) {
			std::string detail = std::to_string(numDone) + "/" +
				std::to_string(triangleLoc1(numParams));
			Global->reportProgress1(cnd.name, detail);
		}
	};
	
 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

void omxComputeNumericDeriv::omxPopulateHessianWork(struct hess_struct *hess_work, FitContext* fc)
{
	hess_work->probeCount = 0;
	hess_work->Haprox = new double[numIter];
	hess_work->Gcentral = new double[numIter];
	hess_work->Gforward = new double[numIter];
	hess_work->Gbackward = new double[numIter];
	hess_work->fitMatrix = fc->lookupDuplicate(fitMat);
	hess_work->fc = fc;
}

/**
  @params i              parameter number
  @params hess_work      local copy
  @params optima         shared read-only variable
  @params gradient       shared write-only variable
  @params hessian        shared write-only variable
 */
void omxComputeNumericDeriv::omxEstimateHessianOnDiagonal(int i, struct hess_struct* hess_work)
{
	int ix = paramMap[i];
	static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm

	double *Haprox             = hess_work->Haprox;
	double *Gcentral             = hess_work->Gcentral;
	double *Gforward             = hess_work->Gforward;
	double *Gbackward            = hess_work->Gbackward;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 
	FitContext* fc = hess_work->fc; 
	double *freeParams         = fc->est;

	/* Part the first: Gradient and diagonal */
	double iOffset = std::max(fabs(stepSize * optima[i]), stepSize);
	for(int k = 0; k < numIter; k++) {			// Decreasing step size, starting at k == 0
		freeParams[ix] = optima[i] + iOffset;
		
		fc->copyParamToModel();

		++hess_work->probeCount;
		omxRecompute(fitMatrix, fc);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		freeParams[ix] = optima[i] - iOffset;

		fc->copyParamToModel();

		++hess_work->probeCount;
		omxRecompute(fitMatrix, fc);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Gcentral[k] = (f1 - f2) / (2.0*iOffset); 						// This is for the gradient
		Gforward[k] = (minimum - f2) / iOffset;
		Gbackward[k] = (f1 - minimum) / iOffset;
		Haprox[k] = (f1 - 2.0 * minimum + f2) / (iOffset * iOffset);		// This is second derivative
		freeParams[ix] = optima[i];									// Reset parameter value
		iOffset /= v;
		if(verbose >= 2) {
			mxLog("Hessian: diag[%s] Δ%.16g (#%d) F1 %.15f F2 %.15f grad %f hess %f",
			      fc->varGroup->vars[i]->name, iOffset, k, f1, f2, Gcentral[k], Haprox[k]);
		}
	}

	for(int m = 1; m < numIter; m++) {						// Richardson Step
		for(int k = 0; k < (numIter - m); k++) {
			// NumDeriv Hard-wires 4s for r here. Why?
			Gcentral[k] = (Gcentral[k+1] * pow(4.0, m) - Gcentral[k])/(pow(4.0, m)-1);
			Gforward[k] = (Gforward[k+1] * pow(4.0, m) - Gforward[k])/(pow(4.0, m)-1);
			Gbackward[k] = (Gbackward[k+1] * pow(4.0, m) - Gbackward[k])/(pow(4.0, m)-1);
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k])/(pow(4.0, m)-1);
		}
	}

	if(verbose >= 2) {
		mxLog("Hessian: diag[%s] final grad %f hess %f", fc->varGroup->vars[i]->name, Gcentral[0], Haprox[0]);
	}
	gcentral[i]  = Gcentral[0];
	gforward[i]  = Gforward[0];
	gbackward[i] = Gbackward[0];
	if (hessian) hessian[i*numParams + i] = Haprox[0];
}

void omxComputeNumericDeriv::omxEstimateHessianOffDiagonal(int i, int l, struct hess_struct* hess_work)
{
	int ix = paramMap[i];
	int lx = paramMap[l];
    static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm

	double *Haprox             = hess_work->Haprox;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 
	FitContext* fc = hess_work->fc; 
	double *freeParams         = fc->est;

	double iOffset = std::max(fabs(stepSize*optima[i]), stepSize);
	double lOffset = std::max(fabs(stepSize*optima[l]), stepSize);

	for(int k = 0; k < numIter; k++) {
		freeParams[ix] = optima[i] + iOffset;
		freeParams[lx] = optima[l] + lOffset;

		fc->copyParamToModel();

		++hess_work->probeCount;
		omxRecompute(fitMatrix, fc);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		freeParams[ix] = optima[i] - iOffset;
		freeParams[lx] = optima[l] - lOffset;

		fc->copyParamToModel();

		++hess_work->probeCount;
		omxRecompute(fitMatrix, fc);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Haprox[k] = (f1 - 2.0 * minimum + f2 - hessian[i*numParams+i]*iOffset*iOffset -
						hessian[l*numParams+l]*lOffset*lOffset)/(2.0*iOffset*lOffset);
		if(verbose >= 2) {
			mxLog("Hessian first off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f)",
			      Haprox[k], iOffset, lOffset, f1, hessian[i*numParams+i], hessian[l*numParams+l],
			      v, k, pow(v, k), stepSize*optima[i], stepSize*optima[l]);
		}

		freeParams[ix] = optima[i];				// Reset parameter values
		freeParams[lx] = optima[l];

		iOffset = iOffset / v;					//  And shrink step
		lOffset = lOffset / v;
	}

	for(int m = 1; m < numIter; m++) {						// Richardson Step
		for(int k = 0; k < (numIter - m); k++) {
			//if(OMX_DEBUG) {mxLog("Hessian off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f, and %f).", Haprox[k], iOffset, lOffset, stepSize, optima[i], optima[l], v, m, pow(4.0, m), stepSize*optima[i], stepSize*optima[l], k);}
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k]) / (pow(4.0, m)-1);
		}
	}

	if(verbose >= 2) {
		mxLog("Hessian estimation: Populating Hessian"
		      " ([%d, %d] = %d and %d) with value %f...",
		      i, l, i*numParams+l, l*numParams+i, Haprox[0]);
	}
	hessian[i*numParams+l] = Haprox[0];
	hessian[l*numParams+i] = Haprox[0];
}

void omxComputeNumericDeriv::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	/*if (state->conListX.size()) {
		mxThrow("%s: cannot proceed with constraints (%d constraints found)",
			name, int(state->conListX.size()));
	}*/

	fitMat = omxNewMatrixFromSlot(rObj, state, "fitfunction");

	SEXP slotValue;

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("iterations")));
	numIter = INTEGER(slotValue)[0];
	if (numIter < 2) mxThrow("stepSize must be 2 or greater");

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("parallel")));
	parallel = Rf_asLogical(slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("checkGradient")));
	checkGradient = Rf_asLogical(slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);

	{
		ProtectedSEXP Rhessian(R_do_slot(rObj, Rf_install("hessian")));
		wantHessian = Rf_asLogical(Rhessian);
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("stepSize")));
	stepSize = GRADIENT_FUDGE_FACTOR(3.0) * REAL(slotValue)[0];
	if (stepSize <= 0) mxThrow("stepSize must be positive");

	knownHessian = NULL;
	{
		ScopedProtect(slotValue, R_do_slot(rObj, Rf_install("knownHessian")));
		if (!Rf_isNull(slotValue)) {
			knownHessian = REAL(slotValue);
			SEXP dimnames;
			ScopedProtect pdn(dimnames, Rf_getAttrib(slotValue, R_DimNamesSymbol));
			{
				SEXP names;
				ScopedProtect p1(names, VECTOR_ELT(dimnames, 0));
				{
					int nlen = Rf_length(names);
					khMap.assign(nlen, -1);
					for (int nx=0; nx < nlen; ++nx) {
						const char *vname = CHAR(STRING_ELT(names, nx));
						for (int vx=0; vx < int(varGroup->vars.size()); ++vx) {
							if (strEQ(vname, varGroup->vars[vx]->name)) {
								khMap[nx] = vx;
								if (verbose >= 1) mxLog("%s: knownHessian[%d] '%s' mapped to %d",
											name, nx, vname, vx);
								break;
							}
						}
					}
				}
			}
		}
	}

	numParams = 0;
	totalProbeCount = 0;
	numParams = 0;
	recordDetail = true;
	detail = 0;
}

void omxComputeNumericDeriv::omxCalcFinalConstraintJacobian(FitContext* fc, int npar){
	allconstraints_functional acf(*fc, verbose);
	Eigen::MatrixWrapper< Eigen::ArrayXd > optimaM(optima);
	Eigen::VectorXd resulttmp(fc->state->numEqC + fc->state->numIneqC);
	Eigen::MatrixXd jactmp(fc->state->numEqC + fc->state->numIneqC, npar);
	acf(optimaM, resulttmp, jactmp);
	/*Gradient algorithm, iterations, and stepsize are hardcoded as they are for two reasons.
	 * 1.  Differentiating the constraint functions should not take long, expecially compared to 
	 * twice-differentiating the fitfunction, so it might as well be done carefully.
	 * 2.  The default behavior during the ComputeNumericDeriv step uses different values of 
	 * gradient stepsize and iterations depending on whether or not the MxModel contains thresholds,
	 * since the numerical accuracy of the -2logL is worse when multivariate-normal integration is involved.
	 * But, that has no bearing on the constraint functions, so it doesn't really make sense to use
	 * the stepsize and iterations stored in the omxComputeNumericDeriv object.
	 */
	fd_jacobian<true>(
		GradientAlgorithm_Central, 4, 1.0e-7,
    acf, resulttmp, optimaM, jactmp);
	
	fc->constraintFunVals = resulttmp;
	fc->constraintJacobian = jactmp;
	/* Subsequent code assumes that fc->est is set to the optimium
	  point. If there is more than 1 thread then the top-level fc->est
	  is unchanged, but if there is only 1 thread then the last
	  parameter in fc->est is offset due to numeric differentiation.
	  Hence, we set it back to its optimum value: */
	fc->est[npar-1L] = optima.coeff(npar-1L);
	return;
}

void omxComputeNumericDeriv::computeImpl(FitContext *fc)
{
	if (fc->fitUnits == FIT_UNITS_SQUARED_RESIDUAL ||
	    fc->fitUnits == FIT_UNITS_SQUARED_RESIDUAL_CHISQ) {  // refactor TODO
		numParams = 0;
		if (verbose >= 1) mxLog("%s: derivatives %s units are meaningless",
					name, fitUnitsToName(fc->fitUnits));
		return; //Possible TODO: calculate Hessian anyway?
	}

	int newWanted = fc->wanted | FF_COMPUTE_GRADIENT;
	if (wantHessian) newWanted |= FF_COMPUTE_HESSIAN;

	int nf = fc->calcNumFree();
	if (numParams != 0 && numParams != nf) {
		mxThrow("%s: number of parameters changed from %d to %d",
			 name, numParams, nf);
	}

	numParams = nf;
	if (numParams <= 0) { complainNoFreeParam(); return; }

	optima.resize(numParams);
	fc->copyEstToOptimizer(optima);
	paramMap.resize(numParams);
	for (int px=0,ex=0; px < numParams; ++ex) {
		if (fc->profiledOut[ex]) continue;
		paramMap[px++] = ex;
	}

	omxAlgebraPreeval(fitMat, fc);
	fc->createChildren(fitMat); // allow FIML rowwiseParallel even when parallel=false

	fc->state->countNonlinearConstraints(fc->state->numEqC, fc->state->numIneqC, false);
	int c_n = fc->state->numEqC + fc->state->numIneqC;
	fc->constraintFunVals.resize(c_n);
	fc->constraintJacobian.resize(c_n, numParams);
	if(c_n){
		omxCalcFinalConstraintJacobian(fc, numParams);
	}
	// TODO: Allow more than one hessian value for calculation

	int numChildren = 1;
	if (parallel && !fc->openmpUser && fc->childList.size()) numChildren = fc->childList.size();

	if (!fc->haveReferenceFit(fitMat)) return;

	minimum = fc->fit;

	hessWorkVector = new hess_struct[numChildren];
	if (numChildren == 1) {
		omxPopulateHessianWork(hessWorkVector, fc);
	} else {
		for(int i = 0; i < numChildren; i++) {
			omxPopulateHessianWork(hessWorkVector + i, fc->childList[i]);
		}
	}
	if(verbose >= 1) mxLog("Numerical Hessian approximation (%d children, ref fit %.16g)",
			       numChildren, minimum);

	hessian = NULL;
	if (wantHessian) {
		hessian = fc->getDenseHessUninitialized();
		Eigen::Map< Eigen::MatrixXd > eH(hessian, numParams, numParams);
		eH.setConstant(NA_REAL);

		if (knownHessian) {
			int khSize = int(khMap.size());
			Eigen::Map< Eigen::MatrixXd > kh(knownHessian, khSize, khMap.size());
			int warns=0;
			for (int rx=0; rx < khSize; ++rx) {
				for (int cx=0; cx < khSize; ++cx) {
					if (khMap[rx] < 0 || khMap[cx] < 0) continue;
					if (!std::isfinite(kh(rx, cx))) continue;
					if (rx < cx && fabs(kh(rx, cx) - kh(cx, rx)) > 1e-9) {
						if (++warns < 3) {
							Rf_warning("%s: knownHessian[%d,%d] is not symmetric",
								   name, 1+rx, 1+cx);
						}
					}
					eH(khMap[rx], khMap[cx]) = (kh(rx, cx) + kh(cx, rx))/2.0;
				}
			}
		}
	}

	if (detail) {
		recordDetail = false; // already done it once
	} else {
		Rf_protect(detail = Rf_allocVector(VECSXP, 4));
		SET_VECTOR_ELT(detail, 0, Rf_allocVector(LGLSXP, numParams));
		for (int gx=0; gx < 3; ++gx) {
			SET_VECTOR_ELT(detail, 1+gx, Rf_allocVector(REALSXP, numParams));
		}
		SEXP detailCols;
		Rf_protect(detailCols = Rf_allocVector(STRSXP, 4));
		Rf_setAttrib(detail, R_NamesSymbol, detailCols);
		SET_STRING_ELT(detailCols, 0, Rf_mkChar("symmetric"));
		SET_STRING_ELT(detailCols, 1, Rf_mkChar("forward"));
		SET_STRING_ELT(detailCols, 2, Rf_mkChar("central"));
		SET_STRING_ELT(detailCols, 3, Rf_mkChar("backward"));

		SEXP detailRowNames;
		Rf_protect(detailRowNames = Rf_allocVector(STRSXP, numParams));
		Rf_setAttrib(detail, R_RowNamesSymbol, detailRowNames);
		for (int nx=0; nx < int(numParams); ++nx) {
			SET_STRING_ELT(detailRowNames, nx, Rf_mkChar(fc->varGroup->vars[nx]->name));
		}
		markAsDataFrame(detail);
	}

	gforward = REAL(VECTOR_ELT(detail, 1));
	gcentral = REAL(VECTOR_ELT(detail, 2));
	gbackward = REAL(VECTOR_ELT(detail, 3));
	Eigen::Map< Eigen::ArrayXd > Gf(gforward, numParams);
	Eigen::Map< Eigen::ArrayXd > Gc(gcentral, numParams);
	Eigen::Map< Eigen::ArrayXd > Gb(gbackward, numParams);
	Gf.setConstant(NA_REAL);
	Gc.setConstant(NA_REAL);
	Gb.setConstant(NA_REAL);

	calcHessianEntry che(this);
	CovEntrywiseParallel(numChildren, che);

	for(int i = 0; i < numChildren; i++) {
		struct hess_struct *hw = hessWorkVector + i;
		totalProbeCount += hw->probeCount;
	}
	delete [] hessWorkVector;
	if (isErrorRaised()) return;

	Eigen::Map< Eigen::ArrayXi > Gsymmetric(LOGICAL(VECTOR_ELT(detail, 0)), numParams);
	double gradNorm = 0.0;
	
	double feasibilityTolerance = Global->feasibilityTolerance;
	for (int px=0; px < numParams; ++px) {
		// factor out simliar code in ComputeNR
		Gsymmetric[px] = true;
		omxFreeVar &fv = *fc->varGroup->vars[ paramMap[px] ];
		if ((fabs(optima[px] - fv.lbound) < feasibilityTolerance && Gc[px] > 0) ||
		    (fabs(optima[px] - fv.ubound) < feasibilityTolerance && Gc[px] < 0)) {
			Gsymmetric[px] = false;
			continue;
		}
		gradNorm += Gc[px] * Gc[px];
		double relsym = 2 * fabs(Gf[px] + Gb[px]) / (Gb[px] - Gf[px]);
		Gsymmetric[px] &= (Gf[px] < 0 && 0 < Gb[px] && relsym < 1.5);
		if (checkGradient && verbose >= 2 && !Gsymmetric[px]) {
			mxLog("%s: param[%d] %d %f", name, px, Gsymmetric[px], relsym);
		}
	}
	
	fc->haveGrad.assign(fc->numParam, true);
	fc->gradZ.resize(fc->numParam);
	fc->gradZ.setZero();
	fc->copyGradFromOptimizer(Gc);
	
	if(c_n){
		fc->inequality.resize(fc->state->numIneqC);
		fc->analyticIneqJacTmp.resize(fc->state->numIneqC, numParams);
		fc->myineqFun(true, verbose, omxConstraint::LESS_THAN, false);
	}

	gradNorm = sqrt(gradNorm);
	double gradThresh = Global->getGradientThreshold(minimum);
	//The gradient will generally not be near zero at a local minimum if there are equality constraints 
	//or active inequality constraints:
	if ( checkGradient && gradNorm > gradThresh && !(fc->state->numEqC || fc->inequality.array().sum()) ) {
		if (verbose >= 1) {
			mxLog("Some gradient entries are too large, norm %f", gradNorm);
		}
		if (fc->getInform() < INFORM_NOT_AT_OPTIMUM) fc->setInform(INFORM_NOT_AT_OPTIMUM);
	}

	fc->setEstFromOptimizer(optima);
	// auxillary information like per-row likelihoods need a refresh
	ComputeFit(name, fitMat, FF_COMPUTE_FIT, fc);
	fc->wanted = newWanted;
}

void omxComputeNumericDeriv::reportResults(FitContext *fc, MxRList *slots, MxRList *result)
{
	if (!numParams || !(fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN | FF_COMPUTE_GRADIENT))) return;

	if (wantHessian) {
		SEXP calculatedHessian;
		Rf_protect(calculatedHessian = Rf_allocMatrix(REALSXP, numParams, numParams));
		fc->copyDenseHess(REAL(calculatedHessian));
		result->add("calculatedHessian", calculatedHessian);
	}

	MxRList out;
	out.add("probeCount", Rf_ScalarInteger(totalProbeCount));
	if (detail && recordDetail) {
		Eigen::Map< Eigen::ArrayXi > Gsymmetric(LOGICAL(VECTOR_ELT(detail, 0)), numParams);
		out.add("gradient", detail);
	}
	slots->add("output", out.asR());
}

omxCompute *newComputeNumericDeriv()
{
	return new omxComputeNumericDeriv;
}

