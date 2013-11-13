/*
  Copyright 2012-2013 Joshua Nathaniel Pritikin and contributors

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>

#include "omxFitFunction.h"
#include "omxExpectationBA81.h"
#include "omxOpenmpWrap.h"
#include "libifa-rpf.h"

struct BA81FitState {

	bool haveLatentMap;
	std::vector<int> latentMap;
	bool freeLatents;

	bool haveItemMap;
	int itemDerivPadSize;     // maxParam + maxParam*(1+maxParam)/2
	std::vector<int> paramFlavor;        // freeParam
	std::vector<int> paramMap;           // itemParam->cols * itemDerivPadSize -> index of free parameter
	std::vector<int> paramLocations;     // param# -> count of appearances in ItemParam
	std::vector<int> itemParamFree;      // itemParam->cols * itemParam->rows
	std::vector<int> ihessDivisor;       // freeParam * freeParam
	std::vector< matrixVectorProdTerm > hgProd;

	std::vector< FreeVarGroup* > varGroups;
	size_t numItemParam;

	omxMatrix *itemParam;
	omxMatrix *latentMean;
	omxMatrix *latentCov;

	BA81FitState();
	~BA81FitState();
	void copyEstimates(BA81Expect *estate);
};

BA81FitState::BA81FitState()
{
	haveItemMap = false;
	haveLatentMap = false;
	freeLatents = false;
}

void BA81FitState::copyEstimates(BA81Expect *estate)
{
	omxCopyMatrix(itemParam, estate->itemParam);
	omxCopyMatrix(latentMean, estate->latentMeanOut);
	omxCopyMatrix(latentCov, estate->latentCovOut);
}

static void buildLatentParamMap(omxFitFunction* oo, FitContext *fc)
{
	FreeVarGroup *fvg = fc->varGroup;
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	std::vector<int> &latentMap = state->latentMap;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	int meanNum = estate->latentMeanOut->matrixNumber;
	int covNum = estate->latentCovOut->matrixNumber;
	int itemNum = estate->itemParam->matrixNumber;
	int maxAbilities = estate->maxAbilities;
	int numLatents = maxAbilities + triangleLoc1(maxAbilities);

	latentMap.assign(numLatents, -1);

	int numParam = int(fvg->vars.size());
	for (int px=0; px < numParam; px++) {
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			if (matNum == meanNum) {
				latentMap[loc->row + loc->col] = px;
				state->freeLatents = true;
			} else if (matNum == covNum) {
				int a1 = loc->row;
				int a2 = loc->col;
				if (a1 < a2) std::swap(a1, a2);
				int cell = maxAbilities + triangleLoc1(a1) + a2;
				if (latentMap[cell] == -1) {
					latentMap[cell] = px;

					if (a1 == a2 && fv->lbound == NEG_INF) {
						fv->lbound = 1e-6;  // variance must be positive
						if (fc->est[px] < fv->lbound) {
							error("Starting value for variance %s is negative", fv->name);
						}
					}
				} else if (latentMap[cell] != px) {
					// doesn't work for multigroup constraints TODO
					error("In covariance matrix, %s and %s must be constrained equal to preserve symmetry",
					      fvg->vars[latentMap[cell]]->name, fv->name);
				}
				state->freeLatents = true;
			} else if (matNum == itemNum) {
				omxRaiseErrorf(globalState, "The fitfunction free.set should consist of "
					       "latent distribution parameters, excluding item parameters");
			}
		}
	}
	state->haveLatentMap = TRUE;
}

static void buildItemParamMap(omxFitFunction* oo, FitContext *fc)
{
	FreeVarGroup *fvg = fc->varGroup;
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemParam = estate->itemParam;
	int size = itemParam->cols * state->itemDerivPadSize;
	state->paramMap.assign(size, -1);  // matrix location to free param index
	state->itemParamFree.assign(itemParam->rows * itemParam->cols, FALSE);

	size_t numFreeParams = state->numItemParam = fvg->vars.size();
	state->paramLocations.assign(numFreeParams, 0);
	state->paramFlavor.assign(numFreeParams, -1);

	for (size_t px=0; px < numFreeParams; px++) {
		omxFreeVar *fv = fvg->vars[px];
		state->paramLocations[px] = int(fv->locations.size());
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			// prohibit mean & cov TODO
			if (matNum == itemParam->matrixNumber) {
				int at = loc->col * state->itemDerivPadSize + loc->row;
				state->paramMap[at] = px;
				state->itemParamFree[loc->col * itemParam->rows + loc->row] = TRUE;

				const double *spec = estate->itemSpec[loc->col];
				int id = spec[RPF_ISpecID];
				int flavor;
				double upper, lower;
				(*rpf_model[id].paramInfo)(spec, loc->row, &flavor, &upper, &lower);
				if (state->paramFlavor[px] < 0) {
					state->paramFlavor[px] = flavor;
				} else if (state->paramFlavor[px] != flavor) {
					error("Cannot equate %s with %s[%d,%d]", fv->name,
					      itemParam->name, loc->row, loc->col);
				}
				if (fv->lbound == NEG_INF && isfinite(lower)) {
					fv->lbound = lower;
					if (fc->est[px] < fv->lbound) {
						error("Starting value %s %f less than lower bound %f",
						      fv->name, fc->est[px], lower);
					}
				}
				if (fv->ubound == INF && isfinite(upper)) {
					fv->ubound = upper;
					if (fc->est[px] > fv->ubound) {
						error("Starting value %s %f greater than upper bound %f",
						      fv->name, fc->est[px], upper);
					}
				}
			}
		}
	}

	state->ihessDivisor.resize(size);

	for (int cx=0; cx < itemParam->cols; ++cx) {
		const double *spec = estate->itemSpec[cx];
		int id = spec[RPF_ISpecID];
		int numParam = (*rpf_model[id].numParam)(spec);

		for (int p1=0; p1 < numParam; p1++) {
			int at1 = state->paramMap[cx * state->itemDerivPadSize + p1];
			if (at1 < 0) continue;

			for (int p2=0; p2 <= p1; p2++) {
				int at2 = state->paramMap[cx * state->itemDerivPadSize + p2];
				if (at2 < 0) continue;

				if (at1 < at2) std::swap(at1, at2);  // lower triangle

				//mxLog("Item %d param(%d,%d) -> H[%d,%d]", cx, p1, p2, at1, at2);
				int at = cx * state->itemDerivPadSize + numParam + triangleLoc1(p1) + p2;
				int hoffset = at1 * numFreeParams + at2;

				//mxLog("? H %d * g %d = p %d", hoffset, at2, at1);
				matrixVectorProdTerm mvpt(hoffset, at2, at1);
				state->hgProd.push_back(mvpt);

				if (at1 != at2) {
					matrixVectorProdTerm mvpt(hoffset, at1, at2);
					state->hgProd.push_back(mvpt);
				}

				state->paramMap[at] = numFreeParams + hoffset;

				state->ihessDivisor[at] =
					state->paramLocations[at1] * state->paramLocations[at2];
			}
		}
	}

	state->haveItemMap = TRUE;
	//pia(state->paramMap.data(), state->itemDerivPadSize, itemParam->cols);
}

static double
ba81ComputeEMFit(omxFitFunction* oo, int want, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *customPrior = estate->customPrior;
	omxMatrix *itemParam = estate->itemParam;
	std::vector<const double*> &itemSpec = estate->itemSpec;
        std::vector<int> &cumItemOutcomes = estate->cumItemOutcomes;
	const int maxDims = estate->maxDims;
	const size_t numItems = estate->itemSpec.size();
	const int do_fit = want & FF_COMPUTE_FIT;
	const int do_deriv = want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN);

	if (estate->verbose) mxLog("%s: em.fit(want fit=%d deriv=%d)", oo->matrix->name, do_fit, do_deriv);

	double ll = 0;
	if (customPrior) {
		omxRecompute(customPrior);
		ll = customPrior->data[0];
		// need deriv adjustment TODO
	}

	if (!isfinite(ll)) {
		omxPrint(itemParam, "item param");
		error("Bayesian prior returned %g; do you need to add a lbound/ubound?", ll);
	}

	if (do_fit) ba81OutcomeProb(estate, FALSE, TRUE);

	const int thrDerivSize = itemParam->cols * state->itemDerivPadSize;
	std::vector<double> thrDeriv(thrDerivSize * Global->numThreads);
	double *wherePrep = estate->wherePrep.data();

#pragma omp parallel for num_threads(Global->numThreads) reduction(+:ll)
	for (size_t ix=0; ix < numItems; ix++) {
		const int thrId = omx_absolute_thread_num();
		const double *spec = estate->itemSpec[ix];
		const int id = spec[RPF_ISpecID];
		const rpf_dLL1_t dLL1 = rpf_model[id].dLL1;
		const int iOutcomes = estate->itemOutcomes[ix];
		const int outcomeBase = cumItemOutcomes[ix] * estate->totalQuadPoints;
		const double *weight = estate->expected + outcomeBase;
                const double *oProb = estate->outcomeProb + outcomeBase;
		const double *iparam = omxMatrixColumn(itemParam, ix);
		double *myDeriv = thrDeriv.data() + thrDerivSize * thrId + ix * state->itemDerivPadSize;

		for (long qx=0; qx < estate->totalQuadPoints; qx++) {
			if (do_fit) {
				for (int ox=0; ox < iOutcomes; ox++) {
					ll += weight[ox] * oProb[ox];
				}
			}
			if (do_deriv) {
				(*dLL1)(spec, iparam, wherePrep + qx * maxDims, weight, myDeriv);
			}
			weight += iOutcomes;
			oProb += iOutcomes;
		}
	}

	size_t excluded = 0;

	if (do_deriv) {
		double *deriv0 = thrDeriv.data();

		int perThread = itemParam->cols * state->itemDerivPadSize;
		for (int th=1; th < Global->numThreads; th++) {
			double *thrD = thrDeriv.data() + th * perThread;
			for (int ox=0; ox < perThread; ox++) deriv0[ox] += thrD[ox];
		}

		for (size_t ix=0; ix < numItems; ix++) {
			const double *spec = itemSpec[ix];
			int id = spec[RPF_ISpecID];
			double *iparam = omxMatrixColumn(itemParam, ix);
			double *pad = deriv0 + ix * state->itemDerivPadSize;
			(*rpf_model[id].dLL2)(spec, iparam, pad);
		}

		int numFreeParams = int(state->numItemParam);
		int numParams = itemParam->cols * state->itemDerivPadSize;
		for (int ox=0; ox < numParams; ox++) {
			int to = state->paramMap[ox];
			if (to == -1) continue;

			// Need to check because this can happen if
			// lbounds/ubounds are not set appropriately.
			if (0 && !isfinite(deriv0[ox])) {
				int item = ox / itemParam->rows;
				mxLog("item parameters:\n");
				const double *spec = itemSpec[item];
				int id = spec[RPF_ISpecID];
				int numParam = (*rpf_model[id].numParam)(spec);
				double *iparam = omxMatrixColumn(itemParam, item);
				pda(iparam, numParam, 1);
				// Perhaps bounds can be pulled in from librpf? TODO
				error("Deriv %d for item %d is %f; are you missing a lbound/ubound?",
				      ox, item, deriv0[ox]);
			}

			if (to < numFreeParams) {
				if (want & FF_COMPUTE_GRADIENT) {
					fc->grad[to] += deriv0[ox];
				}
			} else {
				if (want & FF_COMPUTE_HESSIAN) {
					int Hto = to - numFreeParams;
					fc->hess[Hto] += deriv0[ox];
				}
			}
		}

		if (want & FF_COMPUTE_IHESSIAN) {
			for (size_t ix=0; ix < numItems; ix++) {
				const double *spec = itemSpec[ix];
				int id = spec[RPF_ISpecID];
				int iParams = (*rpf_model[id].numParam)(spec);
				double *pad = deriv0 + ix * state->itemDerivPadSize + iParams;
				int *mask = state->itemParamFree.data() + ix * itemParam->rows;
				double stress;
				omxApproxInvertPackedPosDefTriangular(iParams, mask, pad, &stress);
				// If items excluded then ihessDivisor is wrong TODO
				if (stress) ++excluded;
			}
			for (int ox=0; ox < numParams; ox++) {
				int to = state->paramMap[ox];
				if (to == -1) continue;
				if (to >= numFreeParams) {
					int Hto = to - numFreeParams;
					fc->ihess[Hto] += deriv0[ox] / state->ihessDivisor[ox];
				}
			}
		}
	}

	if (excluded && estate->verbose >= 1) {
		mxLog("%s: Hessian not positive definite for %lu/%lu items",
		      oo->matrix->name, excluded, numItems);
	}
	if (excluded == numItems) {
		omxRaiseErrorf(globalState, "Hessian not positive definite for %lu/%lu items",
			       excluded, numItems);
	}

	return -ll;
}

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg)
{}

static void setLatentStartingValues(omxFitFunction *oo, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	std::vector<int> &latentMap = state->latentMap;
	std::vector<double> &ElatentMean = estate->ElatentMean;
	std::vector<double> &ElatentCov = estate->ElatentCov;
	int maxAbilities = estate->maxAbilities;

	if (!estate->Qpoint.size()) return; // if evaluating fit without estimating model

	// can use Ramsay for these parameters without messing up SEM SEs? TODO
	fc->changedEstimates = true;

	for (int a1 = 0; a1 < maxAbilities; ++a1) {
		if (latentMap[a1] >= 0) {
			int to = latentMap[a1];
			fc->est[to] = ElatentMean[a1];
		}

		for (int a2 = 0; a2 <= a1; ++a2) {
			int to = latentMap[maxAbilities + triangleLoc1(a1) + a2];
			if (to < 0) continue;
			fc->est[to] = ElatentCov[a1 * maxAbilities + a2];
		}
	}

	//fc->log("setLatentStartingValues", FF_COMPUTE_ESTIMATE);
}

static double
ba81ComputeFit(omxFitFunction* oo, int want, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;

	if (estate->type == EXPECTATION_AUGMENTED) {
		if (!state->haveItemMap) buildItemParamMap(oo, fc);

		if (state->numItemParam != fc->varGroup->vars.size()) error("mismatch"); // remove TODO

		if (want & FF_COMPUTE_PARAMFLAVOR) {
			for (size_t px=0; px < state->numItemParam; ++px) {
				if (state->paramFlavor[px] < 0) continue;
				fc->flavor[px] = state->paramFlavor[px];
			}
			return 0;
		}

		if (want & FF_COMPUTE_HGPROD) {
			for (size_t px=0; px < state->hgProd.size(); ++px) {
				fc->hgProd.push_back(state->hgProd[px]);
			}
			return 0;
		}

		if (want & FF_COMPUTE_PREOPTIMIZE) {
			omxExpectationCompute(oo->expectation, NULL);
			// schilling_bock_2005_rescale(oo, fc); seems counterproductive
			return 0;
		}

		double got = ba81ComputeEMFit(oo, want, fc);
		return got;
	} else if (estate->type == EXPECTATION_OBSERVED) {
		if (!state->haveLatentMap) buildLatentParamMap(oo, fc);

		if (want & FF_COMPUTE_PREOPTIMIZE) {
			if (state->freeLatents) {
				setLatentStartingValues(oo, fc);
			}
			return 0;
		}

		if (want & (FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN)) {
			warning("%s: Derivs are not available for latent distribution parameters", oo->matrix->name);
		}

		if (want & FF_COMPUTE_MAXABSCHANGE) {
			double mac = std::max(omxMaxAbsDiff(state->itemParam, estate->itemParam),
					      omxMaxAbsDiff(state->latentMean, estate->latentMeanOut));
			fc->mac = std::max(mac, omxMaxAbsDiff(state->latentCov, estate->latentCovOut));
			state->copyEstimates(estate);
		}

		if (want & FF_COMPUTE_FIT) {
			omxExpectationCompute(oo->expectation, NULL);

			double *patternLik = estate->patternLik;
			int *numIdentical = estate->numIdentical;
			int numUnique = estate->numUnique;
			estate->excludedPatterns = 0;
			const double LogLargest = estate->LogLargestDouble;
			double got = 0;
#pragma omp parallel for num_threads(Global->numThreads) reduction(+:got)
			for (int ux=0; ux < numUnique; ux++) {
				if (!validPatternLik(estate, patternLik[ux])) {
#pragma omp atomic
					++estate->excludedPatterns;
					// somehow indicate that this -2LL is provisional TODO
					continue;
				}
				got += numIdentical[ux] * (log(patternLik[ux]) - LogLargest);
			}
			if (estate->verbose) mxLog("%s: fit (%d/%d excluded)",
						   oo->matrix->name, estate->excludedPatterns, numUnique);
			//mxLog("fit %.4f", -2 * got);
			return -2 * got;
		}

		return 0;
	} else {
		error("Confused");
	}
}

static void ba81Compute(omxFitFunction *oo, int want, FitContext *fc)
{
	if (!want) return;
	double got = ba81ComputeFit(oo, want, fc);
	if (got) oo->matrix->data[0] = got;
}

BA81FitState::~BA81FitState()
{
	omxFreeAllMatrixData(itemParam);
	omxFreeAllMatrixData(latentMean);
	omxFreeAllMatrixData(latentCov);
}

static void ba81Destroy(omxFitFunction *oo) {
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	delete state;
}

void omxInitFitFunctionBA81(omxFitFunction* oo)
{
	if (!oo->argStruct) { // ugh!
		BA81FitState *state = new BA81FitState;
		oo->argStruct = state;
	}

	BA81FitState *state = (BA81FitState*) oo->argStruct;

	omxExpectation *expectation = oo->expectation;
	BA81Expect *estate = (BA81Expect*) expectation->argStruct;

	//newObj->data = oo->expectation->data;

	oo->computeFun = ba81Compute;
	oo->setVarGroup = ba81SetFreeVarGroup;
	oo->destructFun = ba81Destroy;
	oo->gradientAvailable = TRUE;
	oo->hessianAvailable = TRUE;
	oo->parametersHaveFlavor = TRUE;

	int maxParam = estate->itemParam->rows;
	state->itemDerivPadSize = maxParam + triangleLoc1(maxParam);

	int numItems = estate->itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = estate->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= rpf_numModels) {
			error("ItemSpec %d has unknown item model %d", ix, id);
		}
	}

	

	state->itemParam = omxInitMatrix(NULL, 0, 0, TRUE, globalState);
	state->latentMean = omxInitMatrix(NULL, 0, 0, TRUE, globalState);
	state->latentCov = omxInitMatrix(NULL, 0, 0, TRUE, globalState);
	state->copyEstimates(estate);
}
