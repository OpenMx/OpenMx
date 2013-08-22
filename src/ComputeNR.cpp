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

	int maxIter;
	double tolerance;
	int inform, iter;
	int verbose;
	bool carefully;

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
	PROTECT(slotValue = GET_SLOT(rObj, install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	PROTECT(slotValue = GET_SLOT(rObj, install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) error("tolerance must be positive");

	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("carefully")));
	carefully = asLogical(slotValue);
}

void omxApproxInvertPosDefTriangular(int dim, double *hess, double *ihess, double *stress)
{
	const char uplo = 'L';
	int info;
	int retries = 0;
	const int maxRetries = 31; // assume >=32 bit integers
	double adj = 0;
	do {
		memcpy(ihess, hess, sizeof(double) * dim * dim);

		if (retries >= 1) {
			int th = maxRetries - retries;
			if (th > 0) {
				adj = 1.0/(1 << th);
			} else {
				adj = (1 << -th);
			}
			for (int px=0; px < dim; ++px) {
				ihess[px * dim + px] += adj;
			}
		}

		F77_CALL(dpotrf)(&uplo, &dim, ihess, &dim, &info);
		if (info < 0) error("Arg %d is invalid", -info);
		if (info == 0) break;
	} while (++retries < maxRetries * 1.5);

	if (info > 0) {
		omxRaiseErrorf(globalState, "Hessian is not even close to positive definite (order %d)", info);
		return;
	}
	F77_CALL(dpotri)(&uplo, &dim, ihess, &dim, &info);
	if (info < 0) error("Arg %d is invalid", -info);
	if (info > 0) {
		// Impossible to fail if dpotrf worked?
		omxRaiseErrorf(globalState, "Hessian is not of full rank");
	}

	if (stress) *stress = adj;
}

void omxApproxInvertPackedPosDefTriangular(int dim, int *mask, double *packedHess, double *stress)
{
	int mdim = 0;
	for (int dx=0; dx < dim; ++dx) if (mask[dx]) mdim += 1;

	std::vector<double> hess(mdim * mdim, 0.0);
	for (int d1=0, px=0, m1=-1; d1 < dim; ++d1) {
		if (mask[d1]) ++m1;
		for (int d2=0, m2=-1; d2 <= d1; ++d2) {
			if (mask[d2]) ++m2;
			if (mask[d1] && mask[d2]) {
				hess[m2 * mdim + m1] = packedHess[px];
			}
			++px;
		}
	}

	std::vector<double> ihess(mdim * mdim);
	omxApproxInvertPosDefTriangular(mdim, hess.data(), ihess.data(), stress);

	for (int d1=0, px=0, m1=-1; d1 < dim; ++d1) {
		if (mask[d1]) ++m1;
		for (int d2=0, m2=-1; d2 <= d1; ++d2) {
			if (mask[d2]) ++m2;
			if (mask[d1] && mask[d2]) {
				packedHess[px] = *stress? 0 : ihess[m2 * mdim + m1];
			}
			++px;
		}
	}
}

class Ramsay1975 {
	// Ramsay, J. O. (1975). Solving Implicit Equations in
	// Psychometric Data Analysis.  Psychometrika, 40(3), 337-360.

	FitContext *fc;
	size_t numParam;
	int flavor;
	int verbose;
	double highWatermark;
	std::vector<double> prevAdj1;
	std::vector<double> prevAdj2;

public:
	double maxCaution;
	double caution;

	Ramsay1975(FitContext *fc, int flavor, double caution, int verbose);
	void recordEstimate(int px, double newEst);
	void recalibrate(bool *restart);
	void restart();
};

Ramsay1975::Ramsay1975(FitContext *fc, int flavor, double caution, int verbose)
{
	this->fc = fc;
	this->flavor = flavor;
	this->verbose = verbose;
	this->caution = caution;
	maxCaution = 0.0;
	highWatermark = std::max(0.5, caution);  // arbitrary guess

	numParam = fc->varGroup->vars.size();
	prevAdj1.assign(numParam, 0);
	prevAdj2.resize(numParam);
}

void Ramsay1975::recordEstimate(int px, double newEst)
{
	omxFreeVar *fv = fc->varGroup->vars[px];
	bool hitBound=false;
	double param = newEst;
	if (param < fv->lbound) {
		hitBound=true;
		param = fc->est[px] - (fc->est[px] - fv->lbound) / 2;
	}
	if (param > fv->ubound) {
		hitBound=true;
		param = fc->est[px] + (fv->ubound - fc->est[px]) / 2;
	}
	
	prevAdj2[px] = prevAdj1[px];
	prevAdj1[px] = param - fc->est[px];
	
	if (verbose >= 4) {
		std::string buf;
		buf += string_snprintf("~%d~%s: %.4f -> %.4f", px, fv->name, fc->est[px], param);
		if (hitBound) {
			buf += string_snprintf(" wanted %.4f but hit bound", newEst);
		}
		if (prevAdj1[px] * prevAdj2[px] < 0) {
			buf += " *OSC*";
		}
		buf += "\n";
		mxLogBig(buf);
	}

	fc->est[px] = param;
}

void Ramsay1975::recalibrate(bool *restart)
{
	double normPrevAdj2 = 0;
	double normAdjDiff = 0;
	std::vector<double> adjDiff(numParam);

	// The choice of norm is also arbitrary. Other norms might work better.
	for (size_t px=0; px < numParam; ++px) {
		if (fc->flavor[px] != flavor) continue;
		adjDiff[px] = prevAdj1[px] - prevAdj2[px];
		normPrevAdj2 += prevAdj2[px] * prevAdj2[px];
	}

	for (size_t px=0; px < numParam; ++px) {
		if (fc->flavor[px] != flavor) continue;
		normAdjDiff += adjDiff[px] * adjDiff[px];
	}
	double ratio = sqrt(normPrevAdj2 / normAdjDiff);
	//if (verbose >= 3) mxLog("Ramsay[%d]: sqrt(%.5f/%.5f) = %.5f",
	// flavor, normPrevAdj2, normAdjDiff, ratio);

	double newCaution = 1 - (1-caution) * ratio;
	if (newCaution < 0) newCaution = 0;  // doesn't make sense to go faster than full speed
	if (newCaution > .95) newCaution = .95;  // arbitrary guess
	if (newCaution < caution) {
		caution = newCaution/3 + 2*caution/3;  // don't speed up too fast, arbitrary ratio
	} else {
		caution = newCaution;
	}
	maxCaution = std::max(maxCaution, caution);
	if (caution < highWatermark || (normPrevAdj2 < 1e-3 && normAdjDiff < 1e-3)) {
		if (verbose >= 3) mxLog("Ramsay[%d]: %.2f caution", flavor, caution);
	} else {
		if (verbose >= 3) mxLog("Ramsay[%d]: caution %.2f > %.2f, extreme oscillation, restart recommended",
					flavor, caution, highWatermark);
		*restart = TRUE;
	}
	highWatermark += .02; // arbitrary guess
}

void Ramsay1975::restart()
{
	prevAdj1.assign(numParam, 0);
	prevAdj2.assign(numParam, 0);
	highWatermark = 1 - (1 - highWatermark) * .5; // arbitrary guess
	caution = std::max(caution, highWatermark);   // arbitrary guess
	maxCaution = std::max(maxCaution, caution);
	highWatermark = caution;
	if (verbose >= 3) {
		mxLog("Ramsay[%d]: restart with %.2f caution %.2f highWatermark",
		      flavor, caution, highWatermark);
	}
}

void pda(const double *ar, int rows, int cols);

void ComputeNR::compute(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		error("Model has no free parameters");
		return;
	}

	if (fitMatrix->fitFunction->parametersHaveFlavor) {
		for (size_t px=0; px < numParam; ++px) {
			fc->flavor[px] = -1;
		}
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PARAMFLAVOR, fc);
	} else {
		OMXZERO(fc->flavor, numParam);
	}

	if (1) { // add conditions to disable TODO
		fc->hgProd.resize(0);
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_HGPROD, fc);
		std::sort(fc->hgProd.begin(), fc->hgProd.end());
		std::vector<matrixVectorProdTerm>::iterator iter = std::unique(fc->hgProd.begin(), fc->hgProd.end());
		fc->hgProd.resize( std::distance(fc->hgProd.begin(), iter) );
		//fc->log("NR", FF_COMPUTE_HGPROD);
	}

	std::vector<Ramsay1975*> ramsay;
	for (size_t px=0; px < numParam; ++px) {
		// global namespace for flavors? TODO
		if (fc->flavor[px] < 0 || fc->flavor[px] > 100) {  // max flavor? TODO
			error("Invalid parameter flavor %d", fc->flavor[px]);
		}
		while (int(ramsay.size()) < fc->flavor[px]+1) {
			Ramsay1975 *ram = new Ramsay1975(fc, fc->flavor[px], fc->caution, verbose);
			ramsay.push_back(ram);
		}
	}

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
	fc->maybeCopyParamToModel(globalState);

	iter = 0;
	int sinceRestart = 0;
	bool converged = false;
	bool approaching = false;
	bool restarted = carefully || fc->caution >= .5;
	double maxAdj = 0;
	double maxAdjSigned = 0;
	int maxAdjFlavor = 0;
	int maxAdjParam = -1;
	double bestLL = 0;

	std::valarray<double> startBackup(fc->est, numParam);
	std::vector<double> bestBackup(numParam);

	fitMatrix->data[0] = 0;  // may not recompute it, don't leave stale data

	if (verbose >= 2) {
		mxLog("Welcome to Newton-Raphson (tolerance %.3g, max iter %d, %ld flavors)",
		      tolerance, maxIter, ramsay.size());
	}
	while (1) {
		if (verbose >= 2) {
			const char *pname = "none";
			if (maxAdjParam >= 0) pname = fc->varGroup->vars[maxAdjParam]->name;
			mxLog("Begin %d/%d iter of Newton-Raphson (prev maxAdj %.4f for %s, flavor %d)",
			      iter+1, maxIter, maxAdjSigned, pname, maxAdjFlavor);
		}

		int want = FF_COMPUTE_GRADIENT|FF_COMPUTE_IHESSIAN;
		if (restarted) want |= FF_COMPUTE_FIT;

		OMXZERO(fc->grad, numParam);
		OMXZERO(fc->ihess, numParam * numParam);

		omxFitFunctionCompute(fitMatrix->fitFunction, want, fc);
		if (isErrorRaised(globalState)) break;

		if (verbose >= 5) {
			int show = FF_COMPUTE_ESTIMATE;
			if (verbose >= 6) show |= FF_COMPUTE_GRADIENT|FF_COMPUTE_IHESSIAN;
			fc->log("Newton-Raphson", show);
		}

		double LL = fitMatrix->data[0];
		if (bestLL == 0 || bestLL > LL) {
			bestLL = LL;
			memcpy(bestBackup.data(), fc->est, sizeof(double) * numParam);
		}
		if (bestLL > 0) {
			if (verbose >= 3) mxLog("Newton-Raphson cantankerous fit %.5f (best %.5f)", LL, bestLL);
			if (approaching && (LL > bestLL + 0.5 || !std::isfinite(LL))) {  // arbitrary threshold
				memcpy(fc->est, bestBackup.data(), sizeof(double) * numParam);
				fc->copyParamToModel(globalState);
				break;
			}
		}

		if (verbose >= 1 || OMX_DEBUG) {
			for (size_t h1=1; h1 < numParam; h1++) {
				for (size_t h2=0; h2 < h1; h2++) {
					if (fc->ihess[h2 * numParam + h1] == 0) continue;
					omxRaiseErrorf(globalState, "Inverse Hessian is not upper triangular");
				}
			}
		}

		if (0) {
			const double maxDiag = 4;
			for (size_t dx=0; dx < numParam; dx++) {
				double mag = fabs(fc->ihess[dx * numParam + dx]);
				if (maxDiag < mag) {
					double old = fc->grad[dx];
					double logBad = log(1 + mag - maxDiag);
					fc->grad[dx] /= (1 + logBad);  // arbitrary guess
					mxLog("ihess bad at diag %lu grad %.8g -> %.8g", dx, old, fc->grad[dx]);
				}
			}
			//fc->log("bad", FF_COMPUTE_IHESSIAN);
		}

		bool restart = false;
		if ((++sinceRestart) % 3 == 0) {
			for (size_t rx=0; rx < ramsay.size(); ++rx) {
				ramsay[rx]->recalibrate(&restart);
			}
		}
		for (size_t px=0; px < numParam; ++px) {
			if (!std::isfinite(fc->grad[px])) {
				if (!restart) {
					if (verbose >= 3) {
						mxLog("Newton-Raphson: grad[%lu] not finite, restart recommended", px);
					}
					restart = true;
					break;
				}
			}
		}
		if ((sinceRestart % 3) == 0 && !restart) {
			// 3 iterations without extreme oscillation or NaN gradients
			if (!approaching && verbose >= 2) {
				mxLog("Newton-Raphson: Probably approaching solution");
			}
			approaching = true;
		}

		maxAdjParam = -1;
		maxAdj = 0;
		if (restart && (!approaching || !restarted)) {
			approaching = false;
			restarted = true;
			bestLL = 0;
			sinceRestart = 0;

			if (verbose >= 2) {
				mxLog("Newton-Raphson: Increasing damping ...");
			}
			for (size_t px=0; px < numParam; ++px) {
				fc->est[px] = startBackup[px];
			}
			for (size_t rx=0; rx < ramsay.size(); ++rx) {
				ramsay[rx]->restart();
			}
		} else {
			double *grad = fc->grad;
			double *ihess = fc->ihess;

			std::vector<double> move(numParam, 0.0);
			for (size_t px=0; px < fc->hgProd.size(); ++px) {
				matrixVectorProdTerm &mvp = fc->hgProd[px];
				move[mvp.dest] += ihess[mvp.hentry] * grad[mvp.gentry];
			}

			if (0) {
				std::vector<double> blasMove(numParam); // uninit
				const char uplo = 'U';
				int dim = int(numParam);
				int incx = 1;
				double alpha = 1;
				double beta = 0;
				F77_CALL(dsymv)(&uplo, &dim, &alpha, fc->ihess, &dim,
						fc->grad, &incx, &beta, blasMove.data(), &incx);

				double blasDiff = 0;
				for (size_t px=0; px < numParam; ++px) {
					blasDiff += fabs(move[px] - blasMove[px]);
				}
				mxLog("blasdiff %.10g", blasDiff);
				if (blasDiff > 1) {
					fc->log("N-R", FF_COMPUTE_GRADIENT|FF_COMPUTE_IHESSIAN|FF_COMPUTE_HGPROD);
					pda(move.data(), 1, numParam);
					pda(blasMove.data(), 1, numParam);
					abort();
				}
			}

			for (size_t px=0; px < numParam; ++px) {
				Ramsay1975 *ramsay1 = ramsay[ fc->flavor[px] ];
				double oldEst = fc->est[px];
				double speed = 1 - ramsay1->caution;

				ramsay1->recordEstimate(px, oldEst - speed * move[px]);

				double badj = fabs(oldEst - fc->est[px]);
				if (maxAdj < badj) {
					maxAdj = badj;
					maxAdjSigned = oldEst - fc->est[px];
					maxAdjParam = px;
					maxAdjFlavor = fc->flavor[px];
				}
			}
			converged = maxAdj < tolerance;
		}

		fc->copyParamToModel(globalState);

		R_CheckUserInterrupt();

		if (converged || ++iter >= maxIter) break;
	}

	fc->caution = 0;
	for (size_t rx=0; rx < ramsay.size(); ++rx) {
		fc->caution = std::max(fc->caution, ramsay[rx]->maxCaution);
		delete ramsay[rx];
	}

	if (verbose >= 1) {
		if (converged) {
			mxLog("Newton-Raphson converged in %d cycles (max change %.12f, max caution %.4f)",
			      iter, maxAdj, fc->caution);
		} else if (iter < maxIter) {
			mxLog("Newton-Raphson not improving on %.6f after %d cycles (max caution %.4f)",
			      bestLL, iter, fc->caution);
		} else if (iter == maxIter) {
			mxLog("Newton-Raphson failed to converge after %d cycles (max caution %.4f)",
			      iter, fc->caution);
		}
	}

	// better status reporting TODO
	if (!converged && iter == maxIter) {
		if (bestLL > 0) {
			memcpy(fc->est, bestBackup.data(), sizeof(double) * numParam);
			fc->copyParamToModel(globalState);
		} else {
			omxRaiseErrorf(globalState, "Newton-Raphson failed to converge after %d cycles", iter);
		}
	}
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
