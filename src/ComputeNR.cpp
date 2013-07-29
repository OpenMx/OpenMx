/*
 *  Copyright 2013 The OpenMx Project
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

#include <valarray>

#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"

class ComputeNR : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMatrix;
	bool adjustStart;

	int maxIter;
	double tolerance;
	int inform, iter;
	int verbose;

public:
	ComputeNR();
	virtual void initFromFrontend(SEXP rObj);
	virtual void compute(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus() { return inform; }  // backward compatibility
};

class omxCompute *newComputeNewtonRaphson()
{
	return new ComputeNR();
}

ComputeNR::ComputeNR()
{
	inform = 0;
	iter = 0;
}

void ComputeNR::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	if (!fitMatrix->fitFunction->hessianAvailable ||
	    !fitMatrix->fitFunction->gradientAvailable) {
		error("Newton-Raphson requires derivatives");
	}

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("adjustStart")));
	adjustStart = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	PROTECT(slotValue = GET_SLOT(rObj, install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) error("tolerance must be positive");

	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);
}

void ComputeNR::compute(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		error("Model has no free parameters");
		return;
	}

	if (adjustStart) {
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
		fc->copyParamToModel(globalState);
	}

	iter = 0;
	//double prevLL = nan("unset");
	//bool decreasing = TRUE;

	std::valarray<double> prevAdj1(numParam);
	std::valarray<double> prevAdj2(numParam);
	double caution = 0.5;

	while (1) {
		if (verbose >= 2) {
			mxLog("Begin %d/%d iter of Newton-Raphson with tolerance %.02g",
			      iter+1, maxIter, tolerance);
		}

		const int want = FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN;

		OMXZERO(fc->grad, numParam);
		OMXZERO(fc->hess, numParam * numParam);

		omxFitFunctionCompute(fitMatrix->fitFunction, want, fc);

		if (verbose >= 4) {
			fc->log("Newton-Raphson", FF_COMPUTE_ESTIMATE);
		}

		// Only need LL for diagnostics; Can avoid computing it? TODO
		//double LL = fitMatrix->data[0];
		//if (isfinite(prevLL) && prevLL < LL - tolerance) decreasing = FALSE;
		//prevLL = LL;

		//		fc->log(FF_COMPUTE_ESTIMATE|FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN);

		if (verbose >= 3) {
			for (size_t h1=1; h1 < numParam; h1++) {
				for (size_t h2=0; h2 < h1; h2++) {
					if (fc->hess[h2 * numParam + h1] == 0) continue;
					error("Hessian is not lower triangular");
				}
			}
		}

		std::vector<double> ihess(numParam * numParam);
		int dim = int(numParam);
		const char uplo = 'L';
		int info;
		int retries = 0;
		const double diag_adj = -38;
		do {
			memcpy(ihess.data(), fc->hess, sizeof(double) * numParam * numParam);

			if (retries >= 1) {
				double adj = exp(diag_adj + retries);
				for (size_t px=0; px < numParam; ++px) {
					ihess[px * numParam + px] += adj;
				}
			}

			F77_CALL(dpotrf)(&uplo, &dim, ihess.data(), &dim, &info);
			if (info < 0) error("Arg %d is invalid", -info);
			if (info == 0) break;
		} while (++retries < -diag_adj/2);

		if (info > 0) {
			omxRaiseErrorf(globalState, "Hessian is not even close to positive definite (order %d)", info);
			//fc->log("Hessian", FF_COMPUTE_HESSIAN);
			break;
		}
		if (verbose >= 2 && retries >= 1) {
			mxLog("Forced to add diag(%.3g) to Hessian (this is bad)", exp(diag_adj + retries));
		}

		F77_CALL(dpotri)(&uplo, &dim, ihess.data(), &dim, &info);
		if (info < 0) error("Arg %d is invalid", -info);
		if (info > 0) {
			// Impossible to fail if dpotrf worked?
			omxRaiseErrorf(globalState, "Hessian is not of full rank");
			break;
		}

		if ((iter+1) % 3 == 0) {
			// Ramsay, J. O. (1975). Solving Implicit
			// Equations in Psychometric Data Analysis.
			// Psychometrika, 40(3), 337-360.

			std::valarray<double> adjDiff(numParam);
			adjDiff = prevAdj1 - prevAdj2;
			//double normPrevAdj1 = sqrt((prevAdj1 * prevAdj1).sum());
			double normPrevAdj2 = (prevAdj2 * prevAdj2).sum();
			double normAdjDiff = (adjDiff * adjDiff).sum();
			double ratio = sqrt(normPrevAdj2 / normAdjDiff);
			//mxLog("normPrevAdj %.2f %.2f normAdjDiff %.2f ratio %.2f",
			//normPrevAdj1, normPrevAdj2, normAdjDiff, ratio);
			caution = 1 - (1-caution) * ratio;
			if (caution < 0) caution = 0;  // doesn't make sense to go faster than full speed
			//mxLog("ratio %.2f so %.2f", ratio, caution);
		}

		std::vector<double> newEst(numParam);
		memcpy(newEst.data(), fc->est, sizeof(double) * numParam);
		int incx = 1;
		double alpha = -(1 - caution);
		double beta = 1;
		F77_CALL(dsymv)(&uplo, &dim, &alpha, ihess.data(), &dim,
				fc->grad, &incx, &beta, newEst.data(), &incx);

		double maxAdj = 0;
		for (size_t px=0; px < numParam; ++px) {
			double param = newEst[px];
			omxFreeVar *fv = fc->varGroup->vars[px];

			bool hitBound=false;
			if (param < fv->lbound) {
				hitBound=true;
				param = fc->est[px] - (fc->est[px] - fv->lbound) / 2;
			}
			if (param > fv->ubound) {
				hitBound=true;
				param = fc->est[px] + (fv->ubound - fc->est[px]) / 2;
			}

			double badj = fabs(param - fc->est[px]);
			if (maxAdj < badj)
				maxAdj = badj;

			if (verbose >= 3) {
				std::string buf(fv->name);
				buf += string_snprintf(": %.4f -> %.4f", fc->est[px], param);
				if (hitBound) {
					buf += string_snprintf(" wanted %.4f but hit bound", newEst[px]);
				}
				bool osc = false;
				if (iter) osc = (param - fc->est[px]) * prevAdj1[px] < 0;
				if (osc) {
					buf += " *OSC*";
				}
				buf += "\n";
				mxLogBig(buf);
			}

			prevAdj2[px] = prevAdj1[px];
			prevAdj1[px] = param - fc->est[px];
			fc->est[px] = param;
		}

		fc->copyParamToModel(globalState);

		R_CheckUserInterrupt();

		if (maxAdj < tolerance) {
			if (verbose >= 1) {
				mxLog("Newton-Raphson converged in %d cycles", iter);
			}
			break;
		} else if (++iter >= maxIter) {
			if (verbose >= 1) {
				mxLog("Newton-Raphson failed to converge after %d cycles", iter);
			}
			break;
		}
	}


	// The check is too dependent on numerical precision to enable by default.
	// Anyway, it's just a tool for developers.
	//if (0 && !decreasing) warning("Newton-Raphson iterations did not converge");

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_POSTOPTIMIZE, fc);
}

void ComputeNR::reportResults(FitContext *fc, MxRList *out)
{
	if (Global->numIntervals) {
		warning("Confidence intervals are not implemented for Newton-Raphson");
	}  

	omxPopulateFitFunction(fitMatrix, out);

	size_t numFree = varGroup->vars.size();

	SEXP estimate;
	PROTECT(estimate = allocVector(REALSXP, numFree));
	memcpy(REAL(estimate), fc->est, sizeof(double) * numFree);

	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fc->fit)));
	out->push_back(std::make_pair(mkChar("Minus2LogLikelihood"), ScalarReal(fc->fit)));
	out->push_back(std::make_pair(mkChar("estimate"), estimate));

	SEXP iterations;

	// SEXP code;
	// PROTECT(code = NEW_NUMERIC(1));
	// REAL(code)[0] = inform;
	// out->push_back(std::make_pair(mkChar("nr.code"), code));

	PROTECT(iterations = NEW_NUMERIC(1));
	REAL(iterations)[0] = iter;
	out->push_back(std::make_pair(mkChar("nr.iterations"), iterations));
}
