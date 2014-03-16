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
#include "matrix.h"

class ComputeNR : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMatrix;

	int maxIter;
	double tolerance;
	int inform, iter;
	int verbose;
	bool carefully;
	std::vector<Ramsay1975*> ramsay;
	double getMaxCaution() const;
	void clearRamsay();

public:
	ComputeNR();
	virtual ~ComputeNR();
	virtual void initFromFrontend(SEXP rObj);
	virtual omxFitFunction *getFitFunction();
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
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

void ComputeNR::clearRamsay()
{
	for (size_t rx=0; rx < ramsay.size(); ++rx) {
		delete ramsay[rx];
	}
	ramsay.clear();
}

ComputeNR::~ComputeNR()
{
	clearRamsay();
}

void ComputeNR::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	if (!fitMatrix->fitFunction->hessianAvailable ||
	    !fitMatrix->fitFunction->gradientAvailable) {
		Rf_error("Newton-Raphson requires derivatives");
	}

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) Rf_error("tolerance must be positive");

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("carefully")));
	carefully = Rf_asLogical(slotValue);
}

omxFitFunction *ComputeNR::getFitFunction()
{ return fitMatrix->fitFunction; }

void omxApproxInvertPosDefTriangular(int dim, double *hess, double *ihess, double *stress)
{
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

		Matrix ihessMat(ihess, dim, dim);
		info = InvertSymmetricPosDef(ihessMat, 'L');
		if (info == 0) break;
	} while (++retries < maxRetries * 1.5);

	if (info > 0) {
		// or just set stress to something high and return? TODO
		omxRaiseErrorf(globalState, "Hessian is not even close to positive definite (order %d)", info);
		return;
	}

	if (stress) *stress = adj;
}

void omxApproxInvertPackedPosDefTriangular(int dim, int *mask, double *packedHess, double *stress)
{
	int mdim = 0;
	for (int dx=0; dx < dim; ++dx) if (mask[dx]) mdim += 1;
	if (mdim == 0) {
		*stress = 0;
		return;
	}

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

void pda(const double *ar, int rows, int cols);

double ComputeNR::getMaxCaution() const
{
	double caution = 0;
	for (size_t rx=0; rx < ramsay.size(); ++rx) {
		caution = std::max(caution, ramsay[rx]->maxCaution);
	}
	return caution;
}

void ComputeNR::computeImpl(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		Rf_error("Model has no free parameters");
		return;
	}

	OMXZERO(fc->flavor, numParam);
	if (fitMatrix->fitFunction->parametersHaveFlavor) {
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PARAMFLAVOR, fc);
	}

	for (size_t px=0; px < numParam; ++px) {
		// global namespace for flavors? TODO
		if (fc->flavor[px] < 0 || fc->flavor[px] > 100) {  // max flavor? TODO
			Rf_error("Invalid parameter flavor %d", fc->flavor[px]);
		}
		while (int(ramsay.size()) < fc->flavor[px]+1) {
			const double minCaution = 0; // doesn't make sense to go faster
			Ramsay1975 *ram = new Ramsay1975(fc, int(ramsay.size()), fc->caution,
							 verbose, minCaution);
			ramsay.push_back(ram);
		}
	}

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

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
		++iter;
		if (verbose >= 2) {
			const char *pname = "none";
			if (maxAdjParam >= 0) pname = fc->varGroup->vars[maxAdjParam]->name;
			mxLog("Begin %d/%d iter of Newton-Raphson (prev maxAdj %.4f for %s, flavor %d)",
			      iter, maxIter, maxAdjSigned, pname, maxAdjFlavor);
		}

		int want = FF_COMPUTE_GRADIENT|FF_COMPUTE_IHESSIAN;
		if (restarted) want |= FF_COMPUTE_FIT;

		fc->grad = Eigen::VectorXd::Zero(fc->numParam);
		fc->clearHessian();

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

		bool restart = false;
		if ((++sinceRestart) % 3 == 0) {
			for (size_t rx=0; rx < ramsay.size(); ++rx) {
				ramsay[rx]->recalibrate(&restart);
			}
		}
		for (size_t px=0; px < numParam; ++px) {
			if (!std::isfinite(fc->grad(px))) {
				if (!restart) {
					if (verbose >= 3) {
						mxLog("Newton-Raphson: grad[%d] not finite, restart recommended", int(px));
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
			Eigen::VectorXd move(fc->ihessGradProd());

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
			converged = maxAdj < tolerance * (1-getMaxCaution());
			//converged = maxAdj < tolerance;  // makes practically no difference
		}

		fc->copyParamToModel(globalState);

		R_CheckUserInterrupt();

		if (converged || iter >= maxIter) break;
	}

	fc->caution = getMaxCaution();
	clearRamsay();

	if (converged) {
		fc->inform = INFORM_CONVERGED_OPTIMUM;
		fc->wanted |= FF_COMPUTE_BESTFIT;
		if (verbose >= 1) {
			mxLog("Newton-Raphson converged in %d cycles (max change %.12f, max caution %.4f)",
			      iter, maxAdj, fc->caution);
		}
	} else if (iter < maxIter) {
		fc->inform = INFORM_UNCONVERGED_OPTIMUM;
		if (verbose >= 1) {
			mxLog("Newton-Raphson not improving on %.6f after %d cycles (max caution %.4f)",
			      bestLL, iter, fc->caution);
		}
	} else if (iter == maxIter) {
		if (bestLL > 0) {
			fc->inform = INFORM_NOT_AT_OPTIMUM;
			memcpy(fc->est, bestBackup.data(), sizeof(double) * numParam);
			fc->copyParamToModel(globalState);
			if (verbose >= 1) {
				mxLog("Newton-Raphson improved but fail to converge after %d cycles (max caution %.4f)",
				      iter, fc->caution);
			}
		} else {
			fc->inform = INFORM_ITERATION_LIMIT;
			if (verbose >= 1) {
				mxLog("Newton-Raphson failed to converge after %d cycles (max caution %.4f)",
				      iter, fc->caution);
			}
		}
	}

	fc->iterations = iter;
}

void ComputeNR::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	if (Global->numIntervals) {
		Rf_warning("Confidence intervals are not implemented for Newton-Raphson");
	}  

	omxPopulateFitFunction(fitMatrix, out);

	slots->add("inform", Rf_ScalarInteger(inform));
	slots->add("iterations", Rf_ScalarInteger(iter));
}
