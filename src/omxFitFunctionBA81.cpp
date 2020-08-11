/*
 * Copyright 2012-2017 Joshua Nathaniel Pritikin and contributors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>

#include "omxFitFunction.h"
#include "omxExpectationBA81.h"
#include "libifa-rpf.h"
#include "matrix.h"
#include "Compute.h"
#include "EnableWarnings.h"

struct BA81FitState : omxFitFunction {
	//private:
	void copyEstimates(BA81Expect *estate);
	//public:

	// numFreeParam is used by both the item and latent parameter
	// maps (for Hessian offsets) and is initialized first to
	// allow either the item or latent parameter map to be built
	// first. We cache it here to avoid passing the FitContext
	// around everywhere.
	size_t numFreeParam;                 // the current var group's

	int haveLatentMap;
	std::vector<int> latentMap;
	bool freeLatents;                    // only to support old style direct latents, remove TODO
	int ElatentVersion;

	int haveItemMap;
	int itemDerivPadSize;                // maxParam + maxParam*(1+maxParam)/2
	bool freeItemParams;
	std::vector<HessianBlock> hBlocks;
	std::vector<int> paramPerItem;       // itemParam->cols
	std::vector<const char *> paramFlavor;        // numFreeParam
	// gradient: itemParam->cols * itemDerivPadSize -> index of free parameter
	// Hessian:  itemParam->cols * itemDerivPadSize -> full Hessian offset in current varGroup
	std::vector<int> paramMap;
	std::vector<int> hbMap;              // itemParam->cols * itemDerivPadSize -> per-item HessianBlock offset
	std::vector<int> itemGradMap;        // index of gradient -> index of free parameter
	std::vector<int> itemParamFree;      // itemParam->cols * itemParam->rows : remove TODO

	// The following are only used to compute FF_COMPUTE_MAXABSCHANGE
	omxMatrix *itemParam;
	omxMatrix *latentMean;
	omxMatrix *latentCov;

	bool returnRowLikelihoods;

	BA81FitState();
	virtual ~BA81FitState();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
};

// writes to upper triangle of full matrix
static void addSymOuterProd(const double weight, const double *vec, const int len, double *out)
{
	for (int d1=0; d1 < len; ++d1) {
		for (int d2=0; d2 <= d1; ++d2) {
			out[d1 * len + d2] += weight * vec[d1] * vec[d2];
		}
	}
}

BA81FitState::BA81FitState()
{
	haveItemMap = FREEVARGROUP_INVALID;
	haveLatentMap = FREEVARGROUP_INVALID;
}

void BA81FitState::copyEstimates(BA81Expect *estate)
{
	omxCopyMatrix(itemParam, estate->itemParam);
	if (estate->_latentMeanOut) omxCopyMatrix(latentMean, estate->_latentMeanOut);
	if (estate->_latentCovOut)  omxCopyMatrix(latentCov, estate->_latentCovOut);
}

static void buildLatentParamMap(omxFitFunction* oo, FitContext *fc)
{
	FreeVarGroup *fvg = fc->varGroup;
	BA81FitState *state = (BA81FitState *) oo;
	std::vector<int> &latentMap = state->latentMap;
	BA81Expect *estate = (BA81Expect*) oo->expectation;
	int maxAbilities = estate->grp.itemDims;

	if (state->haveLatentMap == fc->varGroup->id[0]) return;
	if (estate->verbose >= 1) mxLog("%s: rebuild latent parameter map for var group %d",
					oo->name(), fc->varGroup->id[0]);

	state->freeLatents = false;

	int numLatents = maxAbilities + triangleLoc1(maxAbilities);
	latentMap.assign(numLatents, -1);

	int meanNum = 0;
	if (estate->_latentMeanOut) meanNum = ~estate->_latentMeanOut->matrixNumber;
	int covNum = 0;
	if (estate->_latentCovOut) covNum = ~estate->_latentCovOut->matrixNumber;

	int numParam = int(fvg->vars.size());
	for (int px=0; px < numParam; px++) {
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = loc->matrix;
			if (matNum == meanNum && estate->_latentMeanOut) {
				latentMap[loc->row + loc->col] = px;
				state->freeLatents = true;
			} else if (matNum == covNum && estate->_latentCovOut) {
				int a1 = loc->row;
				int a2 = loc->col;
				if (a1 < a2) std::swap(a1, a2);
				int cell = maxAbilities + triangleLoc1(a1) + a2;
				if (latentMap[cell] == -1) {
					latentMap[cell] = px;

					if (a1 == a2 && fv->lbound == NEG_INF) {
						fv->lbound = estate->grp.quad.MIN_VARIANCE;  // variance must be positive
						Global->boundsUpdated = true;
						if (fc->est[px] < fv->lbound) {
							mxThrow("Starting value for variance %s is not positive", fv->name);
						}
					}
				} else if (latentMap[cell] != px) {
					// doesn't detect similar problems in multigroup constraints TODO
					mxThrow("Covariance matrix must be constrained to preserve symmetry");
				}
				state->freeLatents = true;
			}
		}
	}
	state->haveLatentMap = fc->varGroup->id[0];
}

static void buildItemParamMap(omxFitFunction* oo, FitContext *fc)
{
	FreeVarGroup *fvg = fc->varGroup;
	BA81FitState *state = (BA81FitState *) oo;
	BA81Expect *estate = (BA81Expect*) oo->expectation;
	std::vector<const double*> &itemSpec = estate->grp.spec;

	if (state->haveItemMap == fc->varGroup->id[0]) return;
	if (estate->verbose >= 1) mxLog("%s: rebuild item parameter map for var group %d",
					oo->name(), fc->varGroup->id[0]);

	omxMatrix *itemParam = estate->itemParam;
	int size = itemParam->cols * state->itemDerivPadSize;
	state->freeItemParams = false;
	state->hBlocks.clear();
	state->hBlocks.resize(itemParam->cols);
	state->hbMap.assign(size, -1);     // matrix location to HessianBlock offset
	state->paramMap.assign(size, -1);  // matrix location to free param index
	state->itemParamFree.assign(itemParam->rows * itemParam->cols, FALSE);

	const size_t numFreeParam = state->numFreeParam;
	state->paramFlavor.assign(numFreeParam, NULL);

	int totalParam = 0;
	state->paramPerItem.resize(itemParam->cols);
	for (int cx=0; cx < itemParam->cols; ++cx) {
		const double *spec = itemSpec[cx];
		const int id = spec[RPF_ISpecID];
		const int numParam = (*Glibrpf_model[id].numParam)(spec);
		state->paramPerItem[cx] = numParam;
		totalParam += numParam;
	}
	state->itemGradMap.assign(totalParam, -1);

	for (size_t px=0; px < numFreeParam; px++) {
		omxFreeVar *fv = fvg->vars[px];
		for (size_t lx=0; lx < fv->locations.size(); lx++) {
			omxFreeVarLocation *loc = &fv->locations[lx];
			int matNum = ~loc->matrix;
			if (matNum != itemParam->matrixNumber) continue;

			int at = loc->col * state->itemDerivPadSize + loc->row;
			state->paramMap[at] = px;
			std::vector<int> &varMap = state->hBlocks[loc->col].vars;
			if (std::find(varMap.begin(), varMap.end(), px) == varMap.end()) {
				varMap.push_back(px);
			}
			state->itemParamFree[loc->col * itemParam->rows + loc->row] = TRUE;
			state->freeItemParams = true;

			const double *spec = itemSpec[loc->col];
			int id = spec[RPF_ISpecID];
			const char *flavor;
			double upper, lower;
			(*Glibrpf_model[id].paramInfo)(spec, loc->row, &flavor, &upper, &lower);
			if (state->paramFlavor[px] == 0) {
				state->paramFlavor[px] = flavor;
			} else if (strcmp(state->paramFlavor[px], flavor) != 0) {
				mxThrow("Cannot equate %s with %s[%d,%d]", fv->name,
					 itemParam->name(), loc->row, loc->col);
			}
			if (fv->lbound == NEG_INF && std::isfinite(lower)) {
				fv->lbound = lower;
				Global->boundsUpdated = true;
				if (fc->est[px] < fv->lbound) {
					mxThrow("Starting value %s %f less than lower bound %f",
					      fv->name, fc->est[px], lower);
				}
			}
			if (fv->ubound == INF && std::isfinite(upper)) {
				fv->ubound = upper;
				Global->boundsUpdated = true;
				if (fc->est[px] > fv->ubound) {
					mxThrow("Starting value %s %f greater than upper bound %f",
					      fv->name, fc->est[px], upper);
				}
			}
		}
	}

	int gradOffset = 0;
	for (int cx=0; cx < itemParam->cols; ++cx) {
		for (int rx=0; rx < state->paramPerItem[cx]; ++rx) {
			int at = cx * state->itemDerivPadSize + rx;
			int px = state->paramMap[at];
			if (px >= 0) state->itemGradMap[gradOffset] = px;
			++gradOffset;
		}
	}

	for (int cx=0; cx < itemParam->cols; ++cx) {
		HessianBlock &hb = state->hBlocks[cx];
		int numParam = state->paramPerItem[cx];
		for (int p1=0; p1 < numParam; p1++) {
			const int outer_at1 = state->paramMap[cx * state->itemDerivPadSize + p1];
			if (outer_at1 < 0) continue;
			const int outer_hb1 = std::lower_bound(hb.vars.begin(), hb.vars.end(), outer_at1) - hb.vars.begin();
			if (hb.vars[outer_hb1] != outer_at1) mxThrow("oops");

			for (int p2=0; p2 <= p1; p2++) {
				int at1 = outer_at1;
				int hb1 = outer_hb1;
				int at2 = state->paramMap[cx * state->itemDerivPadSize + p2];
				if (at2 < 0) continue;
				if (p1 == p2 && at1 != at2) mxThrow("oops");
				int hb2 = std::lower_bound(hb.vars.begin(), hb.vars.end(), at2) - hb.vars.begin();
				if (hb.vars[hb2] != at2) mxThrow("oops");

				if (at1 < at2) std::swap(at1, at2); // outer_at1 unaffected
				if (hb1 < hb2) std::swap(hb1, hb2); // outer_hb1 unaffected

				//				mxLog("Item %d param(%d,%d) -> H[%d,%d] B[%d,%d]",
				//				      cx, p1, p2, at1, at2, hb1, hb2);
				int at = cx * state->itemDerivPadSize + numParam + triangleLoc1(p1) + p2;
				int hoffset = at1 * numFreeParam + at2;

				state->paramMap[at] = numFreeParam + hoffset;
				state->hbMap[at] = hb1 * hb.vars.size() + hb2;
			}
		}
	}

	state->haveItemMap = fc->varGroup->id[0];
	//pia(state->paramMap.data(), state->itemDerivPadSize, itemParam->cols);
}

struct ba81mstepEval {
	const int ix;
	const double *spec;
	const int id;
	const rpf_dLL1_t dLL1;
	const double *iparam;
	double *myDeriv;
	ba81mstepEval(int _ix, const double *_spec, BA81Expect *_estate,
		      double *_myDeriv) :
		ix(_ix), spec(_spec),
		id(_spec[RPF_ISpecID]), dLL1(Glibrpf_model[id].dLL1),
		iparam(omxMatrixColumn(_estate->itemParam, ix)),
		myDeriv(_myDeriv)
	{};
	void operator()(double *abscissa, double *outcomeCol, double *iexp)
	{
		(*dLL1)(spec, iparam, abscissa, iexp, myDeriv);
	};
};

static double
ba81ComputeEMFit(omxFitFunction* oo, int want, FitContext *fc)
{
	const double Scale = Global->llScale;
	BA81FitState *state = (BA81FitState*) oo;
	BA81Expect *estate = (BA81Expect*) oo->expectation;
	omxMatrix *itemParam = estate->itemParam;
	std::vector<const double*> &itemSpec = estate->grp.spec;
	ba81NormalQuad &quad = estate->getQuad();
	const int numItems = (int) itemSpec.size();
	const int do_fit = want & FF_COMPUTE_FIT;
	const int do_deriv = want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN);

	if (do_deriv && !state->freeItemParams) {
		omxRaiseErrorf("%s: no free parameters", oo->name());
		return NA_REAL;
	}

	if (state->returnRowLikelihoods) {
		omxRaiseErrorf("%s: vector=TRUE not implemented", oo->name());
		return NA_REAL;
	}

	if (estate->verbose >= 3) mxLog("%s: complete data fit(want fit=%d deriv=%d)", oo->name(), do_fit, do_deriv);

	if (do_fit) quad.cacheOutcomeProb(itemParam->data, TRUE);

	const int thrDerivSize = itemParam->cols * state->itemDerivPadSize;
	std::vector<double> thrDeriv;

	double ll = 0;
	if (do_fit) ll = quad.mstepFit();

	if (do_deriv) {
		thrDeriv.resize(thrDerivSize * Global->numThreads);

#pragma omp parallel for num_threads(Global->numThreads)
		for (int ix=0; ix < numItems; ix++) {
			int thrId = omx_absolute_thread_num();
			double *myDeriv = thrDeriv.data() + thrDerivSize * thrId + ix * state->itemDerivPadSize;
			ba81mstepEval op(ix, itemSpec[ix], estate, myDeriv);
			quad.mstepIter(ix, op);
		}
	}

	int excluded = 0;

	if (do_deriv) {
		double *deriv0 = thrDeriv.data();

		int perThread = itemParam->cols * state->itemDerivPadSize;
		for (int th=1; th < Global->numThreads; th++) {
			double *thrD = thrDeriv.data() + th * perThread;
			for (int ox=0; ox < perThread; ox++) deriv0[ox] += thrD[ox];
		}

		int numFreeParams = int(state->numFreeParam);
		int ox=-1;
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = itemSpec[ix];
			int id = spec[RPF_ISpecID];
			double *iparam = omxMatrixColumn(itemParam, ix);
			double *pad = deriv0 + ix * state->itemDerivPadSize;
			(*Glibrpf_model[id].dLL2)(spec, iparam, pad);

			HessianBlock *hb = state->hBlocks[ix].clone();
			hb->mat.triangularView<Eigen::Upper>().setZero();

			for (int dx=0; dx < state->itemDerivPadSize; ++dx) {
				int to = state->paramMap[++ox];
				if (to == -1) continue;

				// Need to check because this can happen if
				// lbounds/ubounds are not set appropriately.
				if (0 && !std::isfinite(deriv0[ox])) {
					int item = ox / itemParam->rows;
					mxLog("item parameters:\n");
					const double *spec2 = itemSpec[item];
					int id2 = spec2[RPF_ISpecID];
					int numParam = (*Glibrpf_model[id2].numParam)(spec2);
					double *iparam2 = omxMatrixColumn(itemParam, item);
					pda(iparam2, numParam, 1);
					// Perhaps bounds can be pulled in from librpf? TODO
					mxThrow("Deriv %d for item %d is %f; are you missing a lbound/ubound?",
						 ox, item, deriv0[ox]);
				}

				if (to < numFreeParams) {
					if (want & FF_COMPUTE_GRADIENT) {
						fc->haveGrad[to] = true;
						fc->gradZ(to) -= Scale * deriv0[ox];
					}
				} else {
					if (want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
						int Hto = state->hbMap[ox];
						if (Hto >= 0) hb->mat.data()[Hto] -= Scale * deriv0[ox];
					}
				}
			}
			fc->queue(hb);
		}
	}

	fc->skippedRows += excluded;
	if (excluded && estate->verbose >= 1) {
		mxLog("%s: Hessian not positive definite for %d/%d items",
		      oo->name(), (int) excluded, (int) numItems);
	}
	if (excluded == numItems) {
		omxRaiseErrorf("Hessian not positive definite for %d/%d items",
			       (int) excluded, (int) numItems);
	}

	return Scale * ll;
}

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg)
{}

struct ba81sandwichOp {
	const int numItems;
	const int numParam;
	BA81FitState *state;
	std::vector<const int*> &dataColumns;
	std::vector<int> &itemOutcomes;
	std::vector<int> &rowMap;
	std::vector<const double*> &spec;
	omxMatrix *itemParam;
	const int itemDerivPadSize;
	const double abScale;
	Eigen::ArrayXd &rowWeight;

	Eigen::ArrayXXd gradBuf;
	Eigen::ArrayXd patternLik1;
	Eigen::ArrayXi px;
	Eigen::ArrayXi gradOffset;
	Eigen::ArrayXXd expected;
	Eigen::ArrayXXd itemDeriv;
	Eigen::ArrayXXd patGrad;
	Eigen::ArrayXXd breadG;
	Eigen::ArrayXXd breadH;

	ba81sandwichOp(int numThreads, BA81Expect *estate, int _numParam, BA81FitState *_state,
		       omxMatrix *_itemParam, double _abScale) :
		numItems(estate->grp.numItems()), numParam(_numParam), state(_state),
		dataColumns(estate->grp.dataColumns), itemOutcomes(estate->grp.itemOutcomes),
		rowMap(estate->grp.rowMap), spec(estate->grp.spec), itemParam(_itemParam),
		itemDerivPadSize(_state->itemDerivPadSize), abScale(_abScale),
		rowWeight(estate->grp.rowMult)
	{
		gradBuf.resize(numParam, numThreads);
		patternLik1.resize(numThreads);
		px.resize(numThreads);
		gradOffset.resize(numThreads);
		expected.resize(estate->grp.maxOutcomes, numThreads);
		itemDeriv.resize(itemDerivPadSize, numThreads);
		patGrad.resize(numParam, numThreads);
		breadG.resize(numParam * numParam, numThreads);
		breadH.resize(numParam * numParam, numThreads);
		breadG.setZero();
		breadH.setZero();
	};

	int getNumItems() const { return numItems; };

	void beginQuadPoint(int thrId)
	{
		gradOffset[thrId] = 0;
		gradBuf.col(thrId).setZero();
	}

	template <typename T1>
	void operator()(int thrId, Eigen::MatrixBase<T1> &abscissa, double weight, int ix)
	{
		double tmp = weight * patternLik1[thrId];
		double sqrtTmp = sqrt(tmp);

		if (ix) gradOffset(thrId) += state->paramPerItem[ix-1];
		int pick = dataColumns[ix][ rowMap[px(thrId)] ];
		if (pick == NA_INTEGER) return;
		expected.col(thrId).setZero();
		expected(pick, thrId) = 1.0;
		const double *spec1 = spec[ix];
		double *iparam = omxMatrixColumn(itemParam, ix);
		const int id = spec1[RPF_ISpecID];
		itemDeriv.col(thrId).setZero();
		(*Glibrpf_model[id].dLL1)(spec1, iparam, abscissa.derived().data(),
					  &expected.coeffRef(0, thrId),
					  &itemDeriv.coeffRef(0,thrId));
		(*Glibrpf_model[id].dLL2)(spec1, iparam, &itemDeriv.coeffRef(0,thrId));

		for (int par = 0; par < state->paramPerItem[ix]; ++par) {
			int to = state->itemGradMap[gradOffset(thrId) + par];
			if (to >= 0) {
				gradBuf(to,thrId) -= itemDeriv(par,thrId) * sqrtTmp;
				patGrad(to,thrId) -= itemDeriv(par,thrId) * tmp;
			}
		}

		int derivBase = ix * itemDerivPadSize;
		for (int ox=0; ox < itemDerivPadSize; ox++) {
			int to = state->paramMap[derivBase + ox];
			if (to >= int(numParam)) {
				int Hto = to - numParam;
				breadH(Hto, thrId) += abScale * itemDeriv(ox,thrId) * tmp * rowWeight[px(thrId)];
			}
		}
	}
	void endQuadPoint(int thrId)
	{
		addSymOuterProd(abScale * rowWeight[px(thrId)], &gradBuf.coeffRef(0,thrId),
				numParam, &breadG.coeffRef(0,thrId));
	}

};

static void sandwich(omxFitFunction *oo, FitContext *fc)
{
	const double abScale = fabs(Global->llScale);
	omxExpectation *expectation = oo->expectation;
	BA81FitState *state = (BA81FitState*) oo;
	BA81Expect *estate = (BA81Expect*) expectation;
	if (estate->verbose >= 1) mxLog("%s: sandwich", oo->name());

	omxMatrix *itemParam = estate->itemParam;
	estate->grp.quad.cacheOutcomeProb(itemParam->data, FALSE);

	const int numThreads = Global->numThreads;
	const int numUnique = estate->getNumUnique();
	ba81NormalQuad &quad = estate->getQuad();
	std::vector<bool> &rowSkip = estate->grp.rowSkip;
	const int numParam = (int) fc->varGroup->vars.size();
	std::vector<double> thrMeat(numThreads * numParam * numParam);
	Eigen::ArrayXd &rowWeight = estate->grp.rowMult;

	if (quad.hasBifactorStructure) {
		mxThrow("Sandwich information matrix method is not implemented for bifactor-optimized models");
	}

	ba81sandwichOp op(numThreads, estate, numParam, state, itemParam, abScale);

	quad.allocBuffers();

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) continue;
		int thrId = omx_absolute_thread_num();
		double *meat = thrMeat.data() + thrId * numParam * numParam;   //b

		double patternLik1 =
			quad.computePatternLik(thrId, estate->grp.rowMap[px]);

		if (!ifaGroup::validPatternLik(patternLik1)) {
			omxRaiseErrorf("%s: pattern %d has an invalid probability %g",
				       oo->name(), estate->grp.rowMap[px], patternLik1);
			continue;
		}

		op.px[thrId] = px;
		op.patternLik1[thrId] = 1 / patternLik1;
		op.patGrad.col(thrId).setZero();

		quad.computeRowDeriv(thrId, op);

		addSymOuterProd(abScale * rowWeight[px], &op.patGrad.coeffRef(0,thrId), numParam, meat);
	}

	quad.releaseBuffers();

	// only need upper triangle TODO
	for (int tx=1; tx < numThreads; ++tx) {
		op.breadG.col(0) += op.breadG.col(tx);
	}
	for (int tx=1; tx < numThreads; ++tx) {
		op.breadH.col(0) += op.breadH.col(tx);
	}
	for (int tx=1; tx < numThreads; ++tx) {
		double *th = thrMeat.data() + tx * numParam * numParam;
		for (int en=0; en < numParam * numParam; ++en) {
			thrMeat[en] += th[en];
		}
	}
	//pda(op.breadG.data(), numParam, numParam);
	//pda(op.breadH.data(), numParam, numParam);
	//pda(thrMeat.data(), numParam, numParam);
	if (fc->infoA) {
		for (int d1=0; d1 < numParam; ++d1) {
			for (int d2=0; d2 < numParam; ++d2) {
				int cell = d1 * numParam + d2;
				fc->infoA[cell] += op.breadH(cell,0) - op.breadG(cell,0) + thrMeat[cell];
			}
		}
	}
	if (fc->infoB) {
		for (int d1=0; d1 < numParam; ++d1) {
			for (int d2=0; d2 < numParam; ++d2) {
				int cell = d1 * numParam + d2;
				fc->infoB[cell] += thrMeat[cell];
			}
		}
	}
}

static void setLatentStartingValues(omxFitFunction *oo, FitContext *fc) //remove? TODO
{
	BA81FitState *state = (BA81FitState*) oo;
	BA81Expect *estate = (BA81Expect*) oo->expectation;
	std::vector<int> &latentMap = state->latentMap;
	ba81NormalQuad &quad = estate->getQuad();
	int maxAbilities = quad.abilities();
	omxMatrix *estMean = estate->estLatentMean;
	omxMatrix *estCov = estate->estLatentCov;

	for (int a1 = 0; a1 < maxAbilities; ++a1) {
		if (latentMap[a1] >= 0) {
			int to = latentMap[a1];
			fc->est[to] = omxVectorElement(estMean, a1);
		}

		for (int a2 = 0; a2 <= a1; ++a2) {
			int to = latentMap[maxAbilities + triangleLoc1(a1) + a2];
			if (to < 0) continue;
			fc->est[to] = omxMatrixElement(estCov, a1, a2);
		}
	}

	if (estate->verbose >= 1) {
		mxLog("%s: set latent parameters for version %d",
		      oo->name(), estate->ElatentVersion);
	}
}

struct ba81gradCovOp {
	const int numItems;
	std::vector<const int*> &dataColumns;
	std::vector<int> &rowMap;
	std::vector<const double*> &spec;
	omxMatrix *itemParam;
	Eigen::ArrayXXd expected;
	Eigen::ArrayXXd ideriv;
	const int itemDerivPadSize;
	Eigen::ArrayXi px;

	ba81gradCovOp(int _numItems, BA81Expect *estate,
		      int _itemDerivPadSize, omxMatrix *_itemParam,
		      int numThreads, int itemDerivSize) :
		numItems(_numItems), dataColumns(estate->grp.dataColumns),
		rowMap(estate->grp.rowMap),
		spec(estate->grp.spec), itemParam(_itemParam),
		itemDerivPadSize(_itemDerivPadSize)
	{
		px.resize(numThreads);
		expected.resize(estate->grp.maxOutcomes, numThreads);
		ideriv.resize(itemDerivSize, numThreads);
	};

	int getNumItems() const { return numItems; };

	void beginQuadPoint(int thrId) {};

	template <typename T1>
	void operator()(int thrId, Eigen::MatrixBase<T1> &abscissa, double weight, int ix)
	{
		int pick = dataColumns[ix][rowMap[px(thrId)]];
		if (pick == NA_INTEGER) return;
		expected.col(thrId).setZero();
		expected(pick, thrId) = weight;
		const double *spec1 = spec[ix];
		double *iparam = omxMatrixColumn(itemParam, ix);
		const int id = spec1[RPF_ISpecID];
		double *myDeriv = &ideriv.coeffRef(ix * itemDerivPadSize, thrId);
		(*Glibrpf_model[id].dLL1)(spec1, iparam, abscissa.derived().data(),
					  &expected.coeffRef(0, thrId), myDeriv);
	}

	void endQuadPoint(int thrId) {};
};

static void gradCov(omxFitFunction *oo, FitContext *fc)
{
	const double Scale = Global->llScale;
	omxExpectation *expectation = oo->expectation;
	BA81FitState *state = (BA81FitState*) oo;
	BA81Expect *estate = (BA81Expect*) expectation;
	if (estate->verbose >= 1) mxLog("%s: cross product approximation", oo->name());

	omxMatrix *itemParam = estate->itemParam;
	estate->grp.quad.cacheOutcomeProb(itemParam->data, FALSE);

	const int numThreads = Global->numThreads;
	const int numUnique = estate->getNumUnique();
	ba81NormalQuad &quad = estate->getQuad();
	Eigen::ArrayXd &rowWeight = estate->grp.rowMult;
	std::vector<bool> &rowSkip = estate->grp.rowSkip;
	const int itemDerivSize = itemParam->cols * state->itemDerivPadSize;
	const size_t numParam = fc->varGroup->vars.size();
	std::vector<double> thrGrad(numThreads * numParam);
	std::vector<double> thrMeat(numThreads * numParam * numParam);

	int numLatents = 0;
	bool freeLatents = state->freeLatents;
	if (freeLatents) {
		Eigen::VectorXd meanVec;
		Eigen::MatrixXd srcMat;
		estate->getLatentDistribution(fc, meanVec, srcMat);
		int info = quad.cacheDerivCoef(meanVec, srcMat);
		if (info) {
			omxRaiseErrorf("%s: latent covariance matrix is not positive definite", oo->name());
			return;
		}
		numLatents = quad.abilities() + triangleLoc1(quad.abilities());
	}

	ba81gradCovOp op(state->freeItemParams? estate->numItems() : 0,
			 estate, state->itemDerivPadSize, itemParam,
			 numThreads, itemDerivSize);

	quad.allocBuffers();

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) continue;
		int thrId = omx_absolute_thread_num();
		Eigen::ArrayXd latentGrad(numLatents);
		latentGrad.setZero();
		op.ideriv.col(thrId).setZero();

		double patternLik1 =
			quad.computePatternLik(thrId, estate->grp.rowMap[px]);

		if (!ifaGroup::validPatternLik(patternLik1)) {
			omxRaiseErrorf("%s: pattern %d has an invalid probability %g",
			 	       oo->name(), estate->grp.rowMap[px], patternLik1);
			continue;
		}

		quad.prepLatentDist(thrId);

		op.px[thrId] = px;
		quad.computeRowDeriv(thrId, op, freeLatents, latentGrad);

		Eigen::VectorXd patGrad(numParam);
		patGrad.setZero();
		double *grad = thrGrad.data() + thrId * numParam;
		double *meat = thrMeat.data() + thrId * numParam * numParam;
		double weight = 1 / patternLik1;
		int gradOffset = 0;

		for (int ix=0; ix < op.getNumItems(); ++ix) {
			const double *spec = estate->itemSpec(ix);
			double *iparam = omxMatrixColumn(itemParam, ix);
			const int id = spec[RPF_ISpecID];
			double *myDeriv = &op.ideriv.coeffRef(ix * state->itemDerivPadSize, thrId);
			(*Glibrpf_model[id].dLL2)(spec, iparam, myDeriv);

			for (int par = 0; par < state->paramPerItem[ix]; ++par) {
				int to = state->itemGradMap[gradOffset];
				if (to >= 0) patGrad[to] -= weight * myDeriv[par];
				++gradOffset;
			}
		}

		for (int lx=0; lx < numLatents; ++lx) {
			int to = state->latentMap[lx];
			if (to >= 0) patGrad[to] += weight * latentGrad[lx];
		}
		for (size_t par=0; par < numParam; ++par) {
			grad[par] += patGrad[par] * Scale * rowWeight[px];
		}
		addSymOuterProd(fabs(Scale) * rowWeight[px], patGrad.data(), numParam, meat);
	}

	quad.releaseBuffers();
	quad.releaseDerivCoefCache();

	for (int tx=1; tx < numThreads; ++tx) {
		double *th = thrGrad.data() + tx * numParam;
		for (size_t en=0; en < numParam; ++en) {
			thrGrad[en] += th[en];
		}
	}
	for (int tx=1; tx < numThreads; ++tx) {
		double *th = thrMeat.data() + tx * numParam * numParam;
		for (size_t en=0; en < numParam * numParam; ++en) {
			thrMeat[en] += th[en];
		}
	}
	for (size_t d1=0; d1 < numParam; ++d1) {
		fc->haveGrad[d1] = true;
		fc->gradZ(d1) += thrGrad[d1];
	}
	if (fc->infoB) {
		for (size_t d1=0; d1 < numParam; ++d1) {
			for (size_t d2=0; d2 < numParam; ++d2) {
				int cell = d1 * numParam + d2;
				fc->infoB[cell] += thrMeat[cell];
			}
		}
	}
}

void BA81FitState::compute(int want, FitContext *fc)
{
	auto *oo = this;
	BA81FitState *state = (BA81FitState*) this;
	BA81Expect *estate = (BA81Expect*) oo->expectation;

	if (fc) state->numFreeParam = fc->varGroup->vars.size();

	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_FINAL_FIT)) return;

	if (estate->type == EXPECTATION_AUGMENTED) {
		buildItemParamMap(oo, fc);

		if (want & FF_COMPUTE_PREOPTIMIZE) {
			omxExpectationCompute(fc, oo->expectation, NULL);
			return;
		}

		if (want & FF_COMPUTE_INFO) {
			buildLatentParamMap(oo, fc);
			if (!state->freeItemParams) {
				omxRaiseErrorf("%s: no free parameters", oo->name());
				return;
			}
			ba81RefreshQuadrature(oo->expectation);

			if (fc->infoMethod == INFO_METHOD_HESSIAN) {
				ba81ComputeEMFit(oo, FF_COMPUTE_HESSIAN, fc);
			} else {
				omxRaiseErrorf("Information matrix approximation method %d is not available",
					       fc->infoMethod);
				return;
			}
			return;
		}

		double got = ba81ComputeEMFit(oo, want, fc);
		oo->matrix->data[0] = got;
		return;
	} else if (estate->type == EXPECTATION_OBSERVED) {

		if (want == FF_COMPUTE_STARTING) {
			buildLatentParamMap(oo, fc);
			if (state->freeLatents) setLatentStartingValues(oo, fc);
			return;
		}

		if (want & (FF_COMPUTE_INFO | FF_COMPUTE_GRADIENT)) {
			buildLatentParamMap(oo, fc); // only to check state->freeLatents
			buildItemParamMap(oo, fc);
			if (!state->freeItemParams && !state->freeLatents) {
				omxRaiseErrorf("%s: no free parameters", oo->name());
				return;
			}
			ba81RefreshQuadrature(oo->expectation);

			if (want & FF_COMPUTE_GRADIENT ||
			    (want & FF_COMPUTE_INFO && fc->infoMethod == INFO_METHOD_MEAT)) {
				gradCov(oo, fc);
			} else {
				if (state->freeLatents) {
					omxRaiseErrorf("Information matrix approximation method %d is not available",
						       fc->infoMethod);
					return;
				}
				if (!state->freeItemParams) {
					omxRaiseErrorf("%s: no free parameters", oo->name());
					return;
				}
				sandwich(oo, fc);
			}
		}
		if (want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
			omxRaiseErrorf("%s: Hessian is not available for observed data", oo->name());
		}

		if (want & FF_COMPUTE_MAXABSCHANGE) {
			double mac = std::max(omxMaxAbsDiff(state->itemParam, estate->itemParam),
					      omxMaxAbsDiff(state->latentMean, estate->_latentMeanOut));
			fc->mac = std::max(mac, omxMaxAbsDiff(state->latentCov, estate->_latentCovOut));
			state->copyEstimates(estate);
		}

		if (want & FF_COMPUTE_FIT) {
			omxExpectationCompute(fc, oo->expectation, NULL);

			Eigen::ArrayXd &patternLik = estate->grp.patternLik;
			const int numUnique = estate->getNumUnique();
			if (state->returnRowLikelihoods) {
				const double OneOverLargest = estate->grp.quad.getReciprocalOfOne();
				for (int rx=0; rx < numUnique; rx++) {
					int dest = estate->grp.rowMap[rx];
					oo->matrix->data[dest] = patternLik[rx] * OneOverLargest;
				}
			} else {
				Eigen::ArrayXd &rowWeight = estate->grp.rowMult;
				const double LogLargest = estate->LogLargestDouble;
				double got = 0;
#pragma omp parallel for num_threads(Global->numThreads) reduction(+:got)
				for (int ux=0; ux < numUnique; ux++) {
					if (patternLik[ux] == 0) continue;
					got += rowWeight[ux] * (log(patternLik[ux]) - LogLargest);
				}
				double fit = nan("infeasible");
				if (estate->grp.excludedPatterns < numUnique) {
					fit = addSkippedRowPenalty(got, estate->grp.excludedPatterns);
					fit *= Global->llScale;
				}
				fc->skippedRows += estate->grp.excludedPatterns;
				if (estate->verbose >= 1) mxLog("%s: observed fit %.4f (%d/%d excluded)",
								oo->name(), fit, estate->grp.excludedPatterns, numUnique);
				oo->matrix->data[0] = fit;
			}
		}
	} else {
		mxThrow("%s: Predict nothing or scores before computing %d", oo->name(), want);
	}
}

BA81FitState::~BA81FitState()
{
	omxFreeMatrix(itemParam);
	omxFreeMatrix(latentMean);
	omxFreeMatrix(latentCov);
}

omxFitFunction *omxInitFitFunctionBA81()
{ return new BA81FitState; }

void BA81FitState::init()
{
	auto *oo = this;
	auto *state = this;
	omxState *currentState = oo->matrix->currentState;

	BA81Expect *estate = (BA81Expect*) expectation;
	estate->fit = oo;

	if (!estate->itemParam->isSimple()) {
		omxRaiseErrorf("%s: non-simple item param matrices are not supported yet",
			       name());
	}

	oo->gradientAvailable = TRUE;
	oo->hessianAvailable = TRUE;

	int maxParam = estate->itemParam->rows;
	state->itemDerivPadSize = maxParam + triangleLoc1(maxParam);

	int numItems = estate->itemParam->cols;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = estate->itemSpec(ix);
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= Glibrpf_numModels) {
			mxThrow("ItemSpec %d has unknown item model %d", ix, id);
		}
	}

	state->itemParam = omxInitMatrix(0, 0, TRUE, currentState);
	state->latentMean = omxInitMatrix(0, 0, TRUE, currentState);
	state->latentCov = omxInitMatrix(0, 0, TRUE, currentState);
	state->copyEstimates(estate);

	state->returnRowLikelihoods = Rf_asInteger(R_do_slot(oo->rObj, Rf_install("vector")));

	units = returnRowLikelihoods? FIT_UNITS_PROBABILITY : FIT_UNITS_MINUS2LL;
}
