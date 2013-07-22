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

struct BA81FitState {

	bool haveLatentMap;
	std::vector<int> latentMap;

	bool haveItemMap;
	int itemDerivPadSize;     // maxParam + maxParam*(1+maxParam)/2
	std::vector<int> paramMap;            // itemParam->cols * itemDerivPadSize -> index of free parameter
	std::vector<size_t> paramLocations;   // itemParam->cols * itemDerivPadSize -> # of locations

	std::vector<int> NAtriangle; // TODO remove
	omxMatrix *cholCov;
	int choleskyError;
	double *tmpLatentMean;    // maxDims
	double *tmpLatentCov;     // maxDims * maxDims ; only lower triangle is used
	omxMatrix *icov;          // inverse covariance matrix
	int fitCount;             // dubious, remove? TODO
	int gradientCount;        // dubious, remove? TODO

	std::vector< FreeVarGroup* > varGroups;
	size_t numItemParam;

	BA81FitState();
	~BA81FitState();
};

BA81FitState::BA81FitState()
{
	fitCount = 0;
	gradientCount = 0;
	tmpLatentMean = NULL;
	tmpLatentCov = NULL;
	haveItemMap = false;
	haveLatentMap = false;
}

static void buildLatentParamMap(omxFitFunction* oo, FreeVarGroup *fvg)
{
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	std::vector<int> &latentMap = state->latentMap;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	int meanNum = estate->latentMeanOut->matrixNumber;
	int covNum = estate->latentCovOut->matrixNumber;
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
			} else if (matNum == covNum) {
				int a1 = loc->row;
				int a2 = loc->col;
				if (a1 < a2) std::swap(a1, a2);
				int cell = maxAbilities + triangleLoc1(a1) + a2;
				if (latentMap[cell] == -1)
					latentMap[cell] = px;
				else if (latentMap[cell] != px) {
					// doesn't work for multigroup constraints TODO
					error("In covariance matrix, %s and %s must be constrained equal to preserve symmetry",
					      fvg->vars[latentMap[cell]]->name, fv->name);
				}
				if (a1 == a2 && fv->lbound == NEG_INF) {
					fv->lbound = 1e-6;  // variance must be positive
				}
			}
		}
	}
	state->haveLatentMap = TRUE;
}

static void buildItemParamMap(omxFitFunction* oo, FreeVarGroup *fvg)
{
	BA81FitState *state = (BA81FitState *) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemParam = estate->itemParam;
	int size = itemParam->cols * state->itemDerivPadSize;
	state->paramMap.assign(size, -1);  // matrix location to free param index
	state->paramLocations.assign(size, 0);

	size_t numFreeParams = state->numItemParam = fvg->vars.size();
	int *pRow = Realloc(NULL, numFreeParams, int);
	int *pCol = Realloc(NULL, numFreeParams, int);

	for (size_t px=0; px < numFreeParams; px++) {
		pRow[px] = -1;
		pCol[px] = -1;
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			if (matNum == itemParam->matrixNumber) {
				pRow[px] = loc->row;
				pCol[px] = loc->col;
				int at = pCol[px] * state->itemDerivPadSize + pRow[px];
				state->paramMap[at] = px;
				state->paramLocations[at] = fv->locations.size();

				const double *spec = estate->itemSpec[loc->col];
				int id = spec[RPF_ISpecID];
				double upper, lower;
				(*rpf_model[id].paramBound)(spec, loc->row, &upper, &lower);
				if (fv->lbound == NEG_INF && isfinite(lower)) fv->lbound = lower;
				if (fv->ubound == INF && isfinite(upper)) fv->ubound = upper;
			}
		}
	}

	for (size_t p1=0; p1 < numFreeParams; p1++) {
		for (size_t p2=p1; p2 < numFreeParams; p2++) {
			if (pCol[p1] == -1 || pCol[p1] != pCol[p2]) continue;
			const double *spec = estate->itemSpec[pCol[p1]];
			int id = spec[RPF_ISpecID];
			int numParam = (*rpf_model[id].numParam)(spec);
			int r1 = pRow[p1];
			int r2 = pRow[p2];
			if (r1 > r2) { int tmp=r1; r1=r2; r2=tmp; }
			int rowOffset = 0;
			for (int rx=1; rx <= r2; rx++) rowOffset += rx;
			int at = pCol[p1] * state->itemDerivPadSize + numParam + rowOffset + r1;
			state->paramMap[at] = numFreeParams + p1 * numFreeParams + p2;
			if (p2 != p1) state->NAtriangle.push_back(p2 * numFreeParams + p1);
		}
	}

	Free(pRow);
	Free(pCol);

	state->haveItemMap = TRUE;
}

OMXINLINE static double
ba81Fit1Ordinate(omxFitFunction* oo, const int *quad, const double *weight, int want, double *myDeriv)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *itemParam = estate->itemParam;
	int numItems = itemParam->cols;
	int maxOutcomes = estate->maxOutcomes;
	int maxDims = estate->maxDims;
	int do_fit = want & FF_COMPUTE_FIT;
	int do_deriv = want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN);

	double where[maxDims];
	pointToWhere(estate, quad, where, maxDims);

	double *outcomeProb = NULL;
	if (do_fit) {
		outcomeProb = computeRPF(estate, itemParam, quad, TRUE); // avoid malloc/free? TODO
		if (!outcomeProb) return 0;
	}

	double thr_ll = 0;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = estate->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		int iOutcomes = spec[RPF_ISpecOutcomes];

		double area = areaProduct(estate, quad, estate->Sgroup[ix]);
		if (do_fit) {
			for (int ox=0; ox < iOutcomes; ox++) {
#if 0
#pragma omp critical(ba81Fit1OrdinateDebug1)
				if (!std::isfinite(outcomeProb[ix * maxOutcomes + ox])) {
					pda(itemParam->data, itemParam->rows, itemParam->cols);
					pda(outcomeProb, outcomes, numItems);
					error("RPF produced NAs");
				}
#endif
				double got = weight[ox] * outcomeProb[ix * maxOutcomes + ox];
				thr_ll += got * area;
			}
		}

		if (do_deriv) {
			double *iparam = omxMatrixColumn(itemParam, ix);
			double *pad = myDeriv + ix * state->itemDerivPadSize;
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
	if (estate->verbose) mxLog("%s: fit-%d", oo->matrix->name, want);

	omxMatrix *customPrior = estate->customPrior;
	omxMatrix *itemParam = estate->itemParam;
	std::vector<const double*> &itemSpec = estate->itemSpec;
	int maxDims = estate->maxDims;
	const int totalOutcomes = estate->totalOutcomes;

	double *thrDeriv = Calloc(itemParam->cols * state->itemDerivPadSize * Global->numThreads, double);

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
		int quad[maxDims];
		decodeLocation(qx, maxDims, estate->quadGridSize, quad);
		double *weight = estate->expected + qx * totalOutcomes;
		double *myDeriv = thrDeriv + itemParam->cols * state->itemDerivPadSize * omx_absolute_thread_num();
		double thr_ll = ba81Fit1Ordinate(oo, quad, weight, want, myDeriv);
		
#pragma omp atomic
		ll += thr_ll;
	}

	if (gradient) {
		double *deriv0 = thrDeriv;

		int perThread = itemParam->cols * state->itemDerivPadSize;
		for (int th=1; th < Global->numThreads; th++) {
			double *thrD = thrDeriv + th * perThread;
			for (int ox=0; ox < perThread; ox++) deriv0[ox] += thrD[ox];
		}

		int numItems = itemParam->cols;
		for (int ix=0; ix < numItems; ix++) {
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
				gradient[to] -= deriv0[ox];
			} else {
				hessian[to - numFreeParams] -= deriv0[ox];
			}
		}
	}

	Free(thrDeriv);

	return -ll;
}

static int
moveLatentDistribution(omxFitFunction *oo, FitContext *fc,
		       double *ElatentMean, double *ElatentCov)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	std::vector<const double*> &itemSpec = estate->itemSpec;
	omxMatrix *itemParam = estate->itemParam;
	omxMatrix *design = estate->design;
	double *tmpLatentMean = state->tmpLatentMean;
	double *tmpLatentCov = state->tmpLatentCov;
	int maxDims = estate->maxDims;
	int maxAbilities = estate->maxAbilities;
	int moveCount = 0;

	int numItems = itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = itemSpec[ix];
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
		int *mask = state->paramMap.data() + state->itemDerivPadSize * ix;
		rpf_model[id].rescale(spec, iparam, mask, tmpLatentMean, tmpLatentCov);
	}

	int numFreeParams = int(fc->varGroup->vars.size());
	for (int rx=0; rx < itemParam->rows; rx++) {
		for (int cx=0; cx < itemParam->cols; cx++) {
			int px = cx * state->itemDerivPadSize + rx;
			int vx = state->paramMap[px];
			if (vx >= 0 && vx < numFreeParams && state->paramLocations[px] == 1) {
				fc->est[vx] = omxMatrixElement(itemParam, rx, cx);
				++moveCount;
			}
		}
	}
	return moveCount;
}

static void
schilling_bock_2005_rescale(omxFitFunction *oo, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	omxMatrix *cholCov = state->cholCov;
	int maxAbilities = estate->maxAbilities;

	//pda(ElatentMean, maxAbilities, 1);
	//pda(estate->ElatentCov.data(), maxAbilities, maxAbilities);
	//omxPrint(design, "design");

	memcpy(cholCov->data, estate->ElatentCov.data(), sizeof(double) * maxAbilities * maxAbilities);

	const char triangle = 'L';
	F77_CALL(dpotrf)(&triangle, &maxAbilities, cholCov->data, &maxAbilities, &state->choleskyError);
	if (state->choleskyError != 0) {
		warning("Cholesky failed with %d; rescaling disabled", state->choleskyError); // make error TODO?
		return;
	}

	//pda(cholCov->data, maxAbilities, maxAbilities);

	int moveCount = moveLatentDistribution(oo, fc, estate->ElatentMean.data(), cholCov->data);
	if (estate->verbose) mxLog("%s: schilling-bock (%d param)", oo->matrix->name, moveCount);
}

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg)
{}

// can use same reorganization to avoid the gram product of where for every pattern TODO
static void mapLatentDeriv(BA81FitState *state, BA81Expect *estate, int sgroup, double piece,
			   const std::vector<double> &derivCoef,
			   double *derivOut)
{
	int maxAbilities = estate->maxAbilities;
	int maxDims = estate->maxDims;
	int pmax = maxDims;
	if (estate->numSpecific) pmax -= 1;

	if (sgroup == 0) {
		int cx = 0;
		for (int d1=0; d1 < pmax; ++d1) {
			double amt1 = piece * derivCoef[d1];
#pragma omp atomic
			derivOut[d1] += amt1;
			for (int d2=0; d2 <= d1; ++d2) {
				int to = maxAbilities + cx;
				double amt2 = piece * derivCoef[maxDims + cx];
#pragma omp atomic
				derivOut[to] += amt2;
				++cx;
			}
		}
	}

	if (estate->numSpecific) {
		int sdim = pmax + sgroup;
		double amt3 = piece * derivCoef[pmax];
#pragma omp atomic
		derivOut[sdim] += amt3;

		double amt4 = piece * derivCoef[maxDims + triangleLoc0(pmax)];
		int to = maxAbilities + triangleLoc0(sdim);
#pragma omp atomic
		derivOut[to] += amt4;
	}
}

static double *reduceForSpecific(omxMatrix *mat, int maxDims, int sgroup)
{
	double *out = Calloc(maxDims * maxDims, double);
	for (int d1=0; d1 < maxDims-1; ++d1) {
		int cell = d1 * maxDims;
		for (int d2=0; d2 < maxDims-1; ++d2) {
			out[cell + d2] = omxMatrixElement(mat, d1, d2);
		}
	}
	int sloc = maxDims-1 + sgroup;
	out[maxDims*maxDims - 1] = omxMatrixElement(mat, sloc, sloc);
	return out;
}

static void calcDerivCoef(BA81FitState *state, BA81Expect *estate,
			  double *where, int sgroup, std::vector<double> *derivCoef)
{
	omxMatrix *mean = estate->latentMeanOut;
	omxMatrix *cov = estate->latentCovOut;
	omxMatrix *icov = state->icov;
	double *covData = cov->data;
	double *icovData = icov->data;
	int maxDims = estate->maxDims;
	const char R='R';
	const char L='L';
	const char U='U';
	const double alpha = 1;
	const double beta = 0;
	const int one = 1;

	double *scov = NULL;
	double *sicov = NULL;
	if (estate->numSpecific) {
		scov = reduceForSpecific(cov, maxDims, sgroup);
		covData = scov;

		sicov = reduceForSpecific(icov, maxDims, sgroup);
		icovData = sicov;
	}

	std::vector<double> whereDiff(maxDims);
	std::vector<double> whereGram(triangleLoc1(maxDims));
	for (int d1=0; d1 < maxDims; ++d1) {
		whereDiff[d1] = where[d1] - omxVectorElement(mean, d1);
	}
	gramProduct(whereDiff.data(), whereDiff.size(), whereGram.data());

	F77_CALL(dsymv)(&U, &maxDims, &alpha, icovData, &maxDims, whereDiff.data(), &one,
			&beta, derivCoef->data(), &one);

	std::vector<double> covGrad1(maxDims * maxDims);
	std::vector<double> covGrad2(maxDims * maxDims);

	int cx=0;
	for (int d1=0; d1 < maxDims; ++d1) {
		for (int d2=0; d2 <= d1; ++d2) {
			covGrad1[d2 * maxDims + d1] = covData[d2 * maxDims + d1] - whereGram[cx];
			++cx;
		}
	}

	F77_CALL(dsymm)(&R, &L, &maxDims, &maxDims, &alpha, covGrad1.data(), &maxDims, icovData,
			&maxDims, &beta, covGrad2.data(), &maxDims);
	F77_CALL(dsymm)(&R, &L, &maxDims, &maxDims, &alpha, icovData, &maxDims, covGrad2.data(),
			&maxDims, &beta, covGrad1.data(), &maxDims);

	for (int d1=0; d1 < maxDims; ++d1) {
		covGrad1[d1 * maxDims + d1] /= 2.0;
	}

	cx = maxDims;
	for (int d1=0; d1 < maxDims; ++d1) {
		int cell = d1 * maxDims;
		for (int d2=0; d2 <= d1; ++d2) {
			(*derivCoef)[cx] = -covGrad1[cell + d2];
			++cx;
		}
	}

	Free(scov);
	Free(sicov);
}

static bool latentDeriv(omxFitFunction *oo, double *gradient)
{
	omxExpectation *expectation = oo->expectation;
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) expectation->argStruct;
	if (estate->verbose) mxLog("%s: latentDeriv", oo->matrix->name);

	int numUnique = estate->numUnique;
	int numSpecific = estate->numSpecific;
	int maxDims = estate->maxDims;
	int maxAbilities = estate->maxAbilities;
	int primaryDims = maxDims;
	omxMatrix *cov = estate->latentCovOut;
	int *numIdentical = estate->numIdentical;
	double *patternLik = estate->patternLik;

	OMXZERO(patternLik, numUnique);

	omxCopyMatrix(state->icov, cov);

	int info;
	omxDPOTRF(state->icov, &info);
	if (info != 0) {
		if (info < 0) error("dpotrf invalid argument %d", -info);
		return FALSE;
	}
	omxDPOTRI(state->icov, &info);
	if (info != 0) {
		if (info < 0) error("dpotri invalid argument %d", -info);
		return FALSE;
	}
	// fill in rest from upper triangle
	for (int rx=1; rx < maxAbilities; ++rx) {
		for (int cx=0; cx < rx; ++cx) {
			omxSetMatrixElement(state->icov, rx, cx, omxMatrixElement(state->icov, cx, rx));
		}
	}

	int maxDerivCoef = maxDims + triangleLoc1(maxDims);
	int numLatents = maxAbilities + triangleLoc1(maxAbilities);
	double *uniqueDeriv = Calloc(numUnique * numLatents, double);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < estate->totalQuadPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, maxDims, estate->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(estate, quad, where, maxDims);
			std::vector<double> derivCoef(maxDerivCoef);
			calcDerivCoef(state, estate, where, 0, &derivCoef);
			double area = estate->priQarea[qx];
			double *lxk = ba81LikelihoodFast(expectation, thrId, 0, quad);

			for (int px=0; px < numUnique; px++) {
				double tmp = (lxk[px] * area);
#pragma omp atomic
				patternLik[px] += tmp;
				mapLatentDeriv(state, estate, 0, tmp, derivCoef,
					       uniqueDeriv + px * numLatents);
			}
		}
	} else {
		primaryDims -= 1;
		int sDim = primaryDims;
		long specificPoints = estate->quadGridSize;

#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < estate->totalPrimaryPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, primaryDims, estate->quadGridSize, quad);

			cai2010(expectation, thrId, FALSE, quad);
			double *allElxk = eBase(estate, thrId);
			double *Eslxk = esBase(estate, thrId);

			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				double where[maxDims];
				pointToWhere(estate, quad, where, maxDims);
				for (int sgroup=0; sgroup < numSpecific; sgroup++) {
					std::vector<double> derivCoef(maxDerivCoef);
					calcDerivCoef(state, estate, where, sgroup, &derivCoef);
					double area = areaProduct(estate, quad, sgroup);
					double *lxk = ba81LikelihoodFast(expectation, thrId, sgroup, quad);
					for (int px=0; px < numUnique; px++) {
						double Ei = allElxk[px];
						double Eis = Eslxk[sgroup * numUnique + px];
						double tmp = ((Ei / Eis) * lxk[px] * area);
						mapLatentDeriv(state, estate, sgroup, tmp, derivCoef,
							       uniqueDeriv + px * numLatents);
					}
				}
			}

			double priArea = estate->priQarea[qx];
			for (int px=0; px < numUnique; px++) {
				double Ei = allElxk[px];
				double tmp = (Ei * priArea);
#pragma omp atomic
				patternLik[px] += tmp;
			}
		}
	}

	/*
	std::vector<double> hess1(triangleLoc1(numLatents));
	std::vector<double> hessSum(triangleLoc1(numLatents));

	// could run nicely in parallel with numUnique * triangleLoc(numLatents) buffer
	for (int px=0; px < numUnique; ++px) {
		gramProduct(uniqueDeriv + px * numLatents, numLatents, hess1.data());
		double dups = numIdentical[px];
		for (int rx=0; rx < triangleLoc1(numLatents); ++rx) {
			hessSum[rx] += hess1[rx] * dups;
		}
	}
	*/

#pragma omp parallel for num_threads(Global->numThreads)
	for (int px=0; px < numUnique; ++px) {
		double weight = numIdentical[px] / patternLik[px];
		for (int rx=0; rx < numLatents; ++rx) {
			uniqueDeriv[px * numLatents + rx] *= weight;
		}
	}

#pragma omp parallel for num_threads(Global->numThreads)
	for (int rx=0; rx < numLatents; ++rx) {
		for (int px=1; px < numUnique; ++px) {
			uniqueDeriv[rx] += uniqueDeriv[px * numLatents + rx];
                }
        }

	for (int l1=0; l1 < numLatents; ++l1) {
		int t1 = state->latentMap[l1];
		if (t1 < 0) continue;
		gradient[t1] -= 2 * uniqueDeriv[l1];

		/*
		for (int l2=0; l2 <= l1; ++l2) {
			int t2 = state->latentMap[l2];
			if (t2 < 0) continue;
			hessian[numLatents * t1 + t2] -= 2 * hessSum[triangleLoc1(l1) + l2];
		}
		*/
	}

	Free(uniqueDeriv);

	return TRUE;
}

static void setLatentStartingValues(omxFitFunction *oo, FitContext *fc)
{
	BA81FitState *state = (BA81FitState*) oo->argStruct;
	BA81Expect *estate = (BA81Expect*) oo->expectation->argStruct;
	std::vector<int> &latentMap = state->latentMap;
	std::vector<double> &ElatentMean = estate->ElatentMean;
	std::vector<double> &ElatentCov = estate->ElatentCov;
	int maxAbilities = estate->maxAbilities;

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

	++state->fitCount;

	if (want & FF_COMPUTE_POSTOPTIMIZE) return 0;

	if (estate->type == EXPECTATION_AUGMENTED) {
		if (!state->haveItemMap) buildItemParamMap(oo, fc->varGroup);

		if (want & FF_COMPUTE_PREOPTIMIZE) {
			schilling_bock_2005_rescale(oo, fc); // how does this work in multigroup? TODO
			return 0;
		}

		if (want & FF_COMPUTE_GRADIENT) ++state->gradientCount;

		for (size_t nx=0; nx < state->NAtriangle.size(); ++nx) {
			fc->hess[ state->NAtriangle[nx] ] = nan("symmetric");
		}

		if (state->numItemParam != fc->varGroup->vars.size()) error("mismatch"); // remove TODO
		double got = ba81ComputeMFit1(oo, want, fc->grad, fc->hess);
		return got;
	} else if (estate->type == EXPECTATION_OBSERVED) {
		if (!state->haveLatentMap) buildLatentParamMap(oo, fc->varGroup);

		if (want & FF_COMPUTE_PREOPTIMIZE) {
			setLatentStartingValues(oo, fc);
			return 0;
		}

		omxExpectationCompute(oo->expectation, NULL);

		if (want & FF_COMPUTE_GRADIENT) {
			if (!latentDeriv(oo, fc->grad)) {
				return INFINITY;
			}
		}

		if (want & FF_COMPUTE_HESSIAN) {
			warning("%s: Hessian is not available for latent distribution parameters", oo->matrix->name);
		}

		if (want & FF_COMPUTE_FIT) {
			if (estate->verbose) mxLog("%s: fit", oo->matrix->name);
			double *patternLik = estate->patternLik;
			int *numIdentical = estate->numIdentical;
			int numUnique = estate->numUnique;
			double got = 0;
#pragma omp parallel for num_threads(Global->numThreads) schedule(static,64) reduction(+:got)
			for (int ux=0; ux < numUnique; ux++) {
				got += numIdentical[ux] * log(patternLik[ux]);
			}
			//mxLog("fit %.2f", -2 * got);
			return -2 * got;
		}

		// if (want & FF_COMPUTE_POSTOPTIMIZE)  discard lxk cache? TODO

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
	Free(tmpLatentMean);
	Free(tmpLatentCov);
	omxFreeAllMatrixData(icov);
	omxFreeAllMatrixData(cholCov);
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

	int maxParam = estate->itemParam->rows;
	state->itemDerivPadSize = maxParam + triangleLoc1(maxParam);

	int maxAbilities = estate->maxAbilities;

	state->tmpLatentMean = Realloc(NULL, estate->maxDims, double);
	state->tmpLatentCov = Realloc(NULL, estate->maxDims * estate->maxDims, double);

	int numItems = estate->itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = estate->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= rpf_numModels) {
			error("ItemSpec %d has unknown item model %d", ix, id);
		}
	}

	state->icov = omxInitMatrix(NULL, maxAbilities, maxAbilities, TRUE, globalState);
	state->cholCov = omxInitMatrix(NULL, maxAbilities, maxAbilities, TRUE, globalState);
}
