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

#include "omxFitFunction.h"
#include "omxExpectationBA81.h"
#include "omxOpenmpWrap.h"
#include "libifa-rpf.h"

static const char *NAME = "FitFunctionBA81";

struct BA81FitState {

	omxMatrix *itemParam;     // M step version
	int derivPadSize;         // maxParam + maxParam*(1+maxParam)/2
	double *thrDeriv;         // itemParam->cols * derivPadSize * thread
	int *paramMap;            // itemParam->cols * derivPadSize -> index of free parameter
	std::vector<int> latentMeanMap;
	std::vector<int> latentCovMap;
	std::vector<int> NAtriangle;
	bool rescale;
	omxMatrix *customPrior;
	int choleskyError;
	double *tmpLatentMean;    // maxDims
	double *tmpLatentCov;     // maxDims * maxDims ; only lower triangle is used
	int fitCount;
	int gradientCount;

	std::vector< FreeVarGroup* > varGroups;
	FreeVarGroup *latentFVG;

	BA81FitState();
};

BA81FitState::BA81FitState()
{
	itemParam = NULL;
	thrDeriv = NULL;
	paramMap = NULL;
	latentFVG = NULL;
	customPrior = NULL;
	fitCount = 0;
	gradientCount = 0;
	tmpLatentMean = NULL;
	tmpLatentCov = NULL;
}

static void buildParamMap(omxFitFunction* oo)
{
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemParam = state->itemParam;
	int size = itemParam->cols * state->derivPadSize;
	int maxAbilities = estate->maxAbilities;
	int meanNum = estate->latentMeanOut->matrixNumber;
	int covNum = estate->latentCovOut->matrixNumber;
	FreeVarGroup *fvg = state->latentFVG;

	state->latentMeanMap.assign(maxAbilities, -1);
	state->latentCovMap.assign(maxAbilities * maxAbilities, -1);

	for (size_t px=0; px < fvg->vars.size(); px++) {
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			if (matNum == meanNum) {
				state->latentMeanMap[loc->row + loc->col] = px;
			} else if (matNum == covNum) {
				state->latentCovMap[loc->col * maxAbilities + loc->row] = px;
			}
		}
	}

	state->paramMap = Realloc(NULL, size, int);  // matrix location to free param index
	for (int px=0; px < size; px++) {
		state->paramMap[px] = -1;
	}

	size_t numFreeParams = oo->freeVarGroup->vars.size();
	int *pRow = Realloc(NULL, numFreeParams, int);
	int *pCol = Realloc(NULL, numFreeParams, int);

	for (size_t px=0; px < numFreeParams; px++) {
		pRow[px] = -1;
		pCol[px] = -1;
		omxFreeVar *fv = oo->freeVarGroup->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			if (matNum == itemParam->matrixNumber) {
				pRow[px] = loc->row;
				pCol[px] = loc->col;
				int at = pCol[px] * state->derivPadSize + pRow[px];
				state->paramMap[at] = px;
			}
		}
	}

	for (size_t p1=0; p1 < numFreeParams; p1++) {
		for (size_t p2=p1; p2 < numFreeParams; p2++) {
			if (pCol[p1] == -1 || pCol[p1] != pCol[p2]) continue;
			const double *spec = omxMatrixColumn(estate->itemSpec, pCol[p1]);
			int id = spec[RPF_ISpecID];
			int numParam = (*rpf_model[id].numParam)(spec);
			int r1 = pRow[p1];
			int r2 = pRow[p2];
			if (r1 > r2) { int tmp=r1; r1=r2; r2=tmp; }
			int rowOffset = 0;
			for (int rx=1; rx <= r2; rx++) rowOffset += rx;
			int at = pCol[p1] * state->derivPadSize + numParam + rowOffset + r1;
			state->paramMap[at] = numFreeParams + p1 * numFreeParams + p2;
			if (p2 != p1) state->NAtriangle.push_back(p2 * numFreeParams + p1);
		}
	}

	Free(pRow);
	Free(pCol);

	state->thrDeriv = Realloc(NULL, itemParam->cols * state->derivPadSize * Global->numThreads, double);
}

OMXINLINE static double
ba81Fit1Ordinate(omxFitFunction* oo, const int *quad, const double *weight, int want)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemSpec = estate->itemSpec;
	omxMatrix *itemParam = state->itemParam;
	int numItems = itemParam->cols;
	int maxOutcomes = estate->maxOutcomes;
	int maxDims = estate->maxDims;
	double *myDeriv = state->thrDeriv + itemParam->cols * state->derivPadSize * omx_absolute_thread_num();
	int do_deriv = want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN);

	double where[maxDims];
	pointToWhere(estate->Qpoint, quad, where, maxDims);

	double *outcomeProb = computeRPF(estate->itemSpec, estate->design, itemParam, estate->maxDims,
					 estate->maxOutcomes, quad, estate->Qpoint); // avoid malloc/free? TODO
	if (!outcomeProb) return 0;

	double thr_ll = 0;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		int id = spec[RPF_ISpecID];
		int iOutcomes = spec[RPF_ISpecOutcomes];

		double area = exp(logAreaProduct(estate, quad, estate->Sgroup[ix]));   // avoid exp() here? TODO
		for (int ox=0; ox < iOutcomes; ox++) {
#if 0
#pragma omp critical(ba81Fit1OrdinateDebug1)
			if (!isfinite(outcomeProb[ix * maxOutcomes + ox])) {
				pda(itemParam->data, itemParam->rows, itemParam->cols);
				pda(outcomeProb, outcomes, numItems);
				error("RPF produced NAs");
			}
#endif
			double got = weight[ox] * outcomeProb[ix * maxOutcomes + ox];
			thr_ll += got * area;
		}

		if (do_deriv) {
			double *iparam = omxMatrixColumn(itemParam, ix);
			double *pad = myDeriv + ix * state->derivPadSize;
			(*rpf_model[id].dLL1)(spec, iparam, where, area, weight, pad);
		}
		weight += iOutcomes;
	}

	Free(outcomeProb);

	return thr_ll;
}

static double
ba81ComputeMFit1(omxFitFunction* oo, int want, double *gradient, double *hessian)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *customPrior = state->customPrior;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *itemSpec = estate->itemSpec;
	int maxDims = estate->maxDims;
	const int totalOutcomes = estate->totalOutcomes;

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

#pragma omp parallel for num_threads(Global->numThreads)
	for (long qx=0; qx < estate->totalQuadPoints; qx++) {
		//double area = exp(state->priLogQarea[qx]);  // avoid exp() here? TODO
		int quad[maxDims];
		decodeLocation(qx, maxDims, estate->quadGridSize, quad);
		double *weight = estate->expected + qx * totalOutcomes;
		double thr_ll = ba81Fit1Ordinate(oo, quad, weight, want);
		
#pragma omp atomic
		ll += thr_ll;
	}

	if (gradient) {
		double *deriv0 = state->thrDeriv;

		int perThread = itemParam->cols * state->derivPadSize;
		for (int th=1; th < Global->numThreads; th++) {
			double *thrD = state->thrDeriv + th * perThread;
			for (int ox=0; ox < perThread; ox++) deriv0[ox] += thrD[ox];
		}

		int numItems = itemParam->cols;
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = omxMatrixColumn(itemSpec, ix);
			int id = spec[RPF_ISpecID];
			double *iparam = omxMatrixColumn(itemParam, ix);
			double *pad = deriv0 + ix * state->derivPadSize;
			(*rpf_model[id].dLL2)(spec, iparam, pad);
		}

		int numFreeParams = int(oo->freeVarGroup->vars.size());
		int numParams = itemParam->cols * state->derivPadSize;
		for (int ox=0; ox < numParams; ox++) {
			int to = state->paramMap[ox];
			if (to == -1) continue;

			// Need to check because this can happen if
			// lbounds/ubounds are not set appropriately.
			if (0 && !isfinite(deriv0[ox])) {
				int item = ox / itemParam->rows;
				mxLog("item parameters:\n");
				const double *spec = omxMatrixColumn(itemSpec, item);
				int id = spec[RPF_ISpecID];
				int numParam = (*rpf_model[id].numParam)(spec);
				double *iparam = omxMatrixColumn(itemParam, item);
				pda(iparam, numParam, 1);
				// Perhaps bounds can be pulled in from librpf? TODO
				error("Deriv %d for item %d is %f; are you missing a lbound/ubound?",
				      ox, item, deriv0[ox]);
			}

			if (to < numFreeParams) {
				gradient[to] -= deriv0[ox];
			} else {
				hessian[to - numFreeParams] -= deriv0[ox];
			}
		}
	}

	return -ll;
}

static void
moveLatentDistribution(omxFitFunction *oo, FitContext *fc,
		       double *ElatentMean, double *ElatentCov)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemSpec = estate->itemSpec;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *design = estate->design;
	double *tmpLatentMean = state->tmpLatentMean;
	double *tmpLatentCov = state->tmpLatentCov;
	int maxDims = estate->maxDims;
	int maxAbilities = estate->maxAbilities;

	int numItems = itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		int id = spec[RPF_ISpecID];
		const double *rawDesign = omxMatrixColumn(design, ix);
		int idesign[design->rows];
		int idx = 0;
		for (int dx=0; dx < design->rows; dx++) {
			if (isfinite(rawDesign[dx])) {
				idesign[idx++] = rawDesign[dx]-1;
			} else {
				idesign[idx++] = -1;
			}
		}
		for (int d1=0; d1 < idx; d1++) {
			if (idesign[d1] == -1) {
				tmpLatentMean[d1] = 0;
			} else {
				tmpLatentMean[d1] = ElatentMean[idesign[d1]];
			}
			for (int d2=0; d2 <= d1; d2++) {
				int cell = idesign[d2] * maxAbilities + idesign[d1];
				if (idesign[d1] == -1 || idesign[d2] == -1) {
					tmpLatentCov[d2 * maxDims + d1] = d1==d2? 1 : 0;
				} else {
					tmpLatentCov[d2 * maxDims + d1] = ElatentCov[cell];
				}
			}
		}
		if (1) {  // ease debugging, make optional TODO
			for (int d1=idx; d1 < maxDims; d1++) tmpLatentMean[d1] = nan("");
			for (int d1=0; d1 < maxDims; d1++) {
				for (int d2=0; d2 < maxDims; d2++) {
					if (d1 < idx && d2 < idx) continue;
					tmpLatentCov[d2 * maxDims + d1] = nan("");
				}
			}
		}
		double *iparam = omxMatrixColumn(itemParam, ix);
		int *mask = state->paramMap + state->derivPadSize * ix;
		rpf_model[id].rescale(spec, iparam, mask, tmpLatentMean, tmpLatentCov);
	}

	int numFreeParams = int(oo->freeVarGroup->vars.size());
	for (int rx=0; rx < itemParam->rows; rx++) {
		for (int cx=0; cx < itemParam->cols; cx++) {
			int vx = state->paramMap[cx * state->derivPadSize + rx];
			if (vx >= 0 && vx < numFreeParams) {
				fc->est[vx] = omxMatrixElement(itemParam, rx, cx);
			}
		}
	}
}

static void
schilling_bock_2005_rescale(omxFitFunction *oo, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	double *ElatentMean = estate->ElatentMean;
	double *ElatentCov = estate->ElatentCov;
	int maxAbilities = estate->maxAbilities;

	//mxLog("schilling bock\n");
	//pda(ElatentMean, maxAbilities, 1);
	//pda(ElatentCov, maxAbilities, maxAbilities);
	//omxPrint(design, "design");

	// use omxDPOTRF instead? TODO
	const char triangle = 'L';
	F77_CALL(dpotrf)(&triangle, &maxAbilities, ElatentCov, &maxAbilities, &state->choleskyError);
	if (state->choleskyError != 0) {
		warning("Cholesky failed with %d; rescaling disabled", state->choleskyError); // make error TODO?
		return;
	}

	moveLatentDistribution(oo, fc, ElatentMean, ElatentCov);
	fc->copyParamToModel(globalState);
}

OMXINLINE static void
updateLatentParam(omxFitFunction* oo, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	int maxAbilities = estate->maxAbilities;

	// TODO need denom for multigroup
	for (int a1=0; a1 < maxAbilities; ++a1) {
		if (state->latentMeanMap[a1] >= 0) {
			fc->est[ state->latentMeanMap[a1] ] = estate->ElatentMean[a1];
		}
		for (int a2=0; a2 < maxAbilities; ++a2) {
			int cell = a2 * maxAbilities + a1;
			if (state->latentCovMap[cell] < 0) continue;
			fc->est[ state->latentCovMap[cell] ] = estate->ElatentCov[cell];
		}
	}

	fc->copyParamToModel(globalState);
}

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg) // too ad hoc? TODO
{
	if (!oo->argStruct) { // ugh!
		BA81FitState *state = new BA81FitState;
		oo->argStruct = state;
	}

	BA81FitState *state = (BA81FitState*) oo->argStruct;

	state->varGroups.push_back(fvg);
	if (state->varGroups.size() == 2) {
		int small = 0;
		if (state->varGroups[0] == state->varGroups[1])
			warning("Cannot recognize correct free parameter groups");
		if (state->varGroups[0]->vars.size() > state->varGroups[1]->vars.size())
			small = 1;
		oo->freeVarGroup = state->varGroups[small];
		state->latentFVG = state->varGroups[!small];
	} else if (state->varGroups.size() > 2) {
		// ignore
	}
}

static double
ba81ComputeFit(omxFitFunction* oo, int want, FitContext *fc)
{
	if (!want) return 0;

	BA81FitState *state = (BA81FitState*) oo->argStruct;

	if (!state->paramMap) buildParamMap(oo);

	if (want & FF_COMPUTE_PREOPTIMIZE) {
		if (state->rescale) schilling_bock_2005_rescale(oo, fc); // how does this work in multigroup? TODO
		return 0;
	}

	if (want & FF_COMPUTE_FIT) {
		++state->fitCount;
	}

	if (want & (FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN)) {
		// M-step

		if (fc->varGroup != oo->freeVarGroup) error("FreeVarGroup mismatch");

		++state->gradientCount;

		omxMatrix *itemParam = state->itemParam;
		OMXZERO(state->thrDeriv, state->derivPadSize * itemParam->cols * Global->numThreads);

		for (size_t nx=0; nx < state->NAtriangle.size(); ++nx) {
			fc->hess[ state->NAtriangle[nx] ] = nan("symmetric");
		}

		double got = ba81ComputeMFit1(oo, want, fc->grad, fc->hess);
		return got;
	} else {
		// Major EM iteration, note completely different LL calculation

		if (fc->varGroup != state->latentFVG) error("FreeVarGroup mismatch");

		updateLatentParam(oo, fc);

		BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
		double *patternLik = estate->patternLik;
		int *numIdentical = estate->numIdentical;
		int numUnique = estate->numUnique;
		double got = 0;
		for (int ux=0; ux < numUnique; ux++) {
			got += numIdentical[ux] * patternLik[ux];
		}
		return -2 * got;
	}
}

static void ba81Compute(omxFitFunction *oo, int want, FitContext *fc)
{
	oo->matrix->data[0] = ba81ComputeFit(oo, want, fc);
}

static void ba81Destroy(omxFitFunction *oo) {
	BA81FitState *state = (BA81FitState *) oo->argStruct;

	omxFreeAllMatrixData(state->customPrior);
	Free(state->paramMap);
	Free(state->thrDeriv);
	Free(state->tmpLatentMean);
	Free(state->tmpLatentCov);
	omxFreeAllMatrixData(state->itemParam);
	delete state;
}

void omxInitFitFunctionBA81(omxFitFunction* oo)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	SEXP rObj = oo->rObj;

	omxExpectation *expectation = oo->expectation;
	BA81Expect *estate = (BA81Expect*) expectation->argStruct;

	//newObj->data = oo->expectation->data;

	oo->computeFun = ba81Compute;
	oo->setVarGroup = ba81SetFreeVarGroup;
	oo->destructFun = ba81Destroy;
	oo->gradientAvailable = TRUE;
	oo->hessianAvailable = TRUE;

	SEXP tmp;
	PROTECT(tmp = GET_SLOT(rObj, install("rescale")));
	state->rescale = asLogical(tmp);

	state->itemParam =
		omxNewMatrixFromSlot(rObj, globalState, "ItemParam");

	if (estate->EitemParam->rows != state->itemParam->rows ||
	    estate->EitemParam->cols != state->itemParam->cols) {
		error("ItemParam and EItemParam matrices must be the same dimension");
	}

	state->customPrior =
		omxNewMatrixFromSlot(rObj, globalState, "CustomPrior");
	
	int maxParam = state->itemParam->rows;
	state->derivPadSize = maxParam + maxParam*(1+maxParam)/2;

	state->tmpLatentMean = Realloc(NULL, estate->maxDims, double);
	state->tmpLatentCov = Realloc(NULL, estate->maxDims * estate->maxDims, double);

	int numItems = state->itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		double *spec = omxMatrixColumn(estate->itemSpec, ix);
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= rpf_numModels) {
			error("ItemSpec column %d has unknown item model %d", ix, id);
		}
	}
}
