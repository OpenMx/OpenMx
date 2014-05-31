/*
  Copyright 2012-2014 Joshua Nathaniel Pritikin and contributors

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

#include <limits>
#include <typeinfo>
#include <Rmath.h>

#include "omxExpectationBA81.h"
#include "glue.h"
#include "libifa-rpf.h"
#include "dmvnorm.h"
#include "omxBuffer.h"
#include "matrix.h"

const struct rpf *rpf_model = NULL;
int rpf_numModels;

void pda(const double *ar, int rows, int cols)
{
	if (rows == 0 || cols == 0) return;
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%.6g, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

void pia(const int *ar, int rows, int cols)
{
	if (rows == 0 || cols == 0) return;
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%d, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

void ba81LikelihoodSlow2(BA81Expect *state, const int px, double *out)
{
	const long totalQuadPoints = state->totalQuadPoints;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	const size_t numItems = state->itemSpec.size();
	omxData *data = state->data;
	const int *colMap = state->colMap;
	std::vector<int> &rowMap = state->rowMap;
	double *oProb = state->outcomeProb;
	std::vector<double> &priQarea = state->priQarea;

	for (long qx=0; qx < totalQuadPoints; ++qx) {
		out[qx] = priQarea[qx];
	}

	const int row = rowMap[px];
	for (size_t ix=0; ix < numItems; ix++) {
		int pick = omxIntDataElementUnsafe(data, row, colMap[ix]);
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		for (long qx=0; qx < totalQuadPoints; ++qx) {
			out[qx] *= oProb[pick];
			oProb += itemOutcomes[ix];
		}
	}
}

void cai2010EiEis(BA81Expect *state, const int px, double *lxk, double *Eis, double *Ei)
{
	const int numSpecific = state->numSpecific;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	double *oProb = state->outcomeProb;
	const long totalQuadPoints = state->totalQuadPoints;
	const long totalPrimaryPoints = state->totalPrimaryPoints;
	const long specificPoints = state->quadGridSize;
	const size_t numItems = state->itemSpec.size();
	const double OneOverLargest = state->OneOverLargestDouble;
	omxData *data = state->data;
	const int *colMap = state->colMap;
	std::vector<int> &rowMap = state->rowMap;
	std::vector<double> &speQarea = state->speQarea;
	std::vector<double> &priQarea = state->priQarea;

	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (long sx=0; sx < specificPoints * numSpecific; sx++) {
			lxk[qloc] = speQarea[sx];
			++qloc;
		}
	}

	const int row = rowMap[px];
	for (size_t ix=0; ix < numItems; ix++) {
		int pick = omxIntDataElementUnsafe(data, row, colMap[ix]);
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;
		int Sgroup = state->Sgroup[ix];
		double *out1 = lxk;
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			out1[Sgroup] *= oProb[pick];
			oProb += itemOutcomes[ix];
			out1 += numSpecific;
		}
	}

	for (long qx=0; qx < totalPrimaryPoints * numSpecific; ++qx) Eis[qx] = 0;
	for (long qx=0; qx < totalPrimaryPoints; ++qx) Ei[qx] = priQarea[qx];

	long eisloc = 0;
	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (long sx=0; sx < specificPoints; sx++) {
			for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
				double piece = lxk[qloc];
				Eis[eisloc + sgroup] += piece;
				++qloc;
			}
		}
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Ei[qx] *= Eis[eisloc + sgroup] * OneOverLargest;
		}
		eisloc += numSpecific;
	}

	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Eis[qloc] = Ei[qx] / Eis[qloc];
			++qloc;
		}
	}
}

void BA81LatentEstimate::mapDenseSpace(struct BA81Expect *state, double piece, const double *where,
				      const double *whereGram, double *latentDist)
{
	const int pmax = state->primaryDims;
	int gx = 0;
	int cx = state->maxAbilities;
	for (int d1=0; d1 < pmax; d1++) {
		double piece_w1 = piece * where[d1];
		latentDist[d1] += piece_w1;
		for (int d2=0; d2 <= d1; d2++) {
			double piece_cov = piece * whereGram[gx];
			latentDist[cx] += piece_cov;
			++cx; ++gx;
		}
	}
}

void BA81LatentEstimate::mapSpecificSpace(struct BA81Expect *state, int sgroup, double piece, const double *where,
					  const double *whereGram, double *latentDist)
{
	int pmax = state->primaryDims;

	int sdim = pmax + sgroup;
	double piece_w1 = piece * where[pmax];
	latentDist[sdim] += piece_w1;

	double piece_var = piece * whereGram[triangleLoc0(pmax)];
	int to = state->maxAbilities + triangleLoc0(sdim);
	latentDist[to] += piece_var;
}

// Depends on item parameters, but not latent distribution
void ba81OutcomeProb(BA81Expect *state, bool estep, bool wantLog)
{
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	std::vector<int> &cumItemOutcomes = state->cumItemOutcomes;
	omxMatrix *itemParam = state->itemParam;
	Eigen::MatrixXi &design = state->design;
	const int maxDims = state->maxDims;
	const size_t numItems = state->itemSpec.size();
	state->outcomeProb = Realloc(state->outcomeProb, state->totalOutcomes * state->totalQuadPoints, double);
	double *param = (estep && state->EitemParam)? state->EitemParam : itemParam->data;

#pragma omp parallel for num_threads(Global->numThreads)
	for (size_t ix=0; ix < numItems; ix++) {
		double *qProb = state->outcomeProb + cumItemOutcomes[ix] * state->totalQuadPoints;
		const double *spec = state->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double *iparam = param + ix * itemParam->rows;
		rpf_prob_t prob_fn = wantLog? rpf_model[id].logprob : rpf_model[id].prob;

		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			double *where = state->wherePrep.data() + qx * maxDims;

			double ptheta[dims];
			for (int dx=0; dx < dims; dx++) {
				int ability = design(dx, ix) - 1; // remove -1 here TODO
				if (ability >= maxDims) ability = maxDims-1;
				ptheta[dx] = where[ability];
			}

			(*prob_fn)(spec, iparam, ptheta, qProb);

			qProb += itemOutcomes[ix];
		}
	}
}

void BA81LatentFixed::normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double patternLik1, int thrId)
{
	double weight = state->rowWeight[px] / patternLik1;
	for (long qx=0; qx < state->ptsPerThread; ++qx) {
		Qweight[qx] *= weight;
	}
}

void BA81LatentScores::begin(struct BA81Expect *state)
{
	numLatents = state->maxAbilities + triangleLoc1(state->maxAbilities);
	thrScore.resize(numLatents * Global->numThreads);
}

void BA81LatentScores::normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double patternLik1, int thrId)
{
	const int maxAbilities = state->maxAbilities;
	const int primaryDims = state->primaryDims;
	omxData *data = state->data;

	// NOTE: Qweight remains unnormalized

	double *scorePad = thrScore.data() + numLatents * thrId;
	OMXZERO(scorePad, numLatents);
	mapSpace(state, Qweight, scorePad);

	double OneOverPatternLik = 1/patternLik1;
	for (int lx=0; lx < numLatents; ++lx) {
		scorePad[lx] *= OneOverPatternLik;
	}

	int cx = maxAbilities;
	for (int a1=0; a1 < primaryDims; ++a1) {
		for (int a2=0; a2 <= a1; ++a2) {
			double ma1 = scorePad[a1];
			double ma2 = scorePad[a2];
			scorePad[cx] -= ma1 * ma2;
			++cx;
		}
	}
	for (int sx=0; sx < state->numSpecific; sx++) {
		int sdim = primaryDims + sx;
		double ma1 = scorePad[sdim];
		scorePad[maxAbilities + triangleLoc0(sdim)] -= ma1 * ma1;
	}

	std::vector<double*> &out = state->scoresOut;
	int dups = omxDataNumIdenticalRows(data, state->rowMap[px]); // should == rowWeight[px]
	for (int dup=0; dup < dups; dup++) {
		int dest = omxDataIndex(data, state->rowMap[px]+dup);

		for (int ax=0; ax < maxAbilities; ++ax) {
			out[ax][dest] = scorePad[ax];
		}
		for (int ax=0; ax < maxAbilities; ++ax) {
			out[maxAbilities + ax][dest] = sqrt(scorePad[maxAbilities + triangleLoc0(ax)]);
		}
		for (int ax=0; ax < triangleLoc1(maxAbilities); ++ax) {
			out[2*maxAbilities + ax][dest] = scorePad[maxAbilities + ax];
		}
	}
}

void BA81LatentScores::end(struct BA81Expect *state)
{
	const int numUnique = (int) state->rowMap.size();
	std::vector<double*> &out = state->scoresOut;
	omxData *data = state->data;
	double *patternLik = state->patternLik;

	for (int px=0; px < numUnique; px++) {
		if (patternLik[px]) continue;
		int dups = omxDataNumIdenticalRows(data, state->rowMap[px]);
		for (int dup=0; dup < dups; dup++) {
			int dest = omxDataIndex(data, state->rowMap[px]+dup);
			for (int ax=0; ax < int(out.size()); ++ax) {
				out[ax][dest] = NA_REAL;
			}
		}
	}
}

void BA81LatentSummary::begin(struct BA81Expect *state)
{
	thrDweight.assign(state->ptsPerThread * Global->numThreads, 0.0);
	numLatents = state->maxAbilities + triangleLoc1(state->maxAbilities);
	latentDist.assign(numLatents, 0.0);
}

void BA81LatentSummary::normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double patternLik1, int thrId)
{
	double weight = state->rowWeight[px] / patternLik1;
	double *Dweight = thrDweight.data() + state->ptsPerThread * thrId;
	for (long qx=0; qx < state->ptsPerThread; ++qx) {
		double tmp = Qweight[qx] * weight;
		Dweight[qx] += tmp;
		Qweight[qx] = tmp;
	}
}

void BA81LatentEstimate::mapSpace(struct BA81Expect *state, double *thrDweight, double *latentDist)
{
	const double *wherePrep = state->wherePrep.data();
	const double *whereGram = state->whereGram.data();
	const int maxDims = state->maxDims;
	const int whereGramSize = triangleLoc1(maxDims);
	const int numSpecific = state->numSpecific;

	if (numSpecific == 0) { // use template to handle this branch at compile time TODO
		for (long qx=0; qx < state->totalQuadPoints; ++qx) {
			mapDenseSpace(state, thrDweight[qx], wherePrep + qx * maxDims,
				      whereGram + qx * whereGramSize, latentDist);
		}
	} else {
		long qloc=0;
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			const double *whPrep = wherePrep + qx * maxDims;
			const double *whGram = whereGram + qx * whereGramSize;
			mapDenseSpace(state, thrDweight[qloc], whPrep, whGram, latentDist);
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapSpecificSpace(state, Sgroup, thrDweight[qloc], whPrep, whGram, latentDist);
				++qloc;
			}
		}
	}
}

void BA81LatentSummary::end(struct BA81Expect *state)
{
	for (int tx=1; tx < Global->numThreads; ++tx) {
		double *Dweight = thrDweight.data() + state->ptsPerThread * tx;
		double *dest = thrDweight.data();
		for (long qx=0; qx < state->ptsPerThread; ++qx) {
			dest[qx] += Dweight[qx];
		}
	}

	mapSpace(state, thrDweight.data(), latentDist.data());

	omxMatrix *meanOut = state->estLatentMean;
	omxMatrix *covOut = state->estLatentCov;
	const double weightSum = state->weightSum;
	const int maxAbilities = state->maxAbilities;
	const int primaryDims = state->primaryDims;

	double *latentDist1 = latentDist.data();
	for (int d1=0; d1 < maxAbilities; d1++) {
		omxSetVectorElement(meanOut, d1, latentDist1[d1] / weightSum);
	}

	for (int d1=0; d1 < primaryDims; d1++) {
		int cx = maxAbilities + triangleLoc1(d1);
		for (int d2=0; d2 <= d1; d2++) {
			double cov = (latentDist1[cx] / weightSum -
				      omxVectorElement(meanOut, d1) * omxVectorElement(meanOut, d2));
			omxSetMatrixElement(covOut, d1, d2, cov);
			if (d1 != d2) omxSetMatrixElement(covOut, d2, d1, cov);
			++cx;
		}
	}
	for (int d1=primaryDims; d1 < maxAbilities; d1++) {
		int loc = maxAbilities + triangleLoc0(d1);
		double cov = (latentDist1[loc] / weightSum -
			      omxVectorElement(meanOut, d1) * omxVectorElement(meanOut, d1));
		omxSetMatrixElement(covOut, d1, d1, cov);
	}

	++state->ElatentVersion;
}

template <typename CovType>
void BA81RefreshPatLik<CovType>::begin(struct BA81Expect *state)
{
}

template <>
int BA81Config<BA81Dense>::getPrimaryPoints(struct BA81Expect *state)
{
	return state->totalQuadPoints;
}

template <>
int BA81Config<BA81TwoTier>::getPrimaryPoints(struct BA81Expect *state)
{
	return state->totalPrimaryPoints;
}

template <typename CovType>
void BA81Estep<CovType>::begin(BA81Expect *state)
{
	totalQuadPoints = state->totalQuadPoints;
	numItems = int(state->itemSpec.size());
	colMap = state->colMap;
	data = state->data;
	thrExpected.assign(state->totalOutcomes * state->totalQuadPoints * Global->numThreads, 0.0);
}

template<>
void BA81Estep<BA81Dense>::addRow(struct BA81Expect *state, int px, double *Qweight, int thrId)
{
	double *out = thrExpected.data() + thrId * state->totalOutcomes * totalQuadPoints;
	std::vector<int> &rowMap = state->rowMap;
	std::vector<int> &itemOutcomes = state->itemOutcomes;

	for (int ix=0; ix < numItems; ++ix) {
		int pick = omxIntDataElementUnsafe(data, rowMap[px], colMap[ix]);
		if (pick == NA_INTEGER) {
			out += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		for (long qx=0; qx < totalQuadPoints; ++qx) {
			out[pick] += Qweight[qx];
			out += itemOutcomes[ix];
		}
	}
}

template<>
void BA81Estep<BA81TwoTier>::addRow(struct BA81Expect *state, int px, double *Qweight, int thrId)
{
	double *out = thrExpected.data() + thrId * state->totalOutcomes * totalQuadPoints;
	std::vector<int> &rowMap = state->rowMap;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	const int numSpecific = state->numSpecific;

	for (int ix=0; ix < numItems; ++ix) {
		int pick = omxIntDataElementUnsafe(data, rowMap[px], colMap[ix]);
		if (pick == NA_INTEGER) {
			out += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		int Sgroup = state->Sgroup[ix];
		double *Qw = Qweight;
		for (long qx=0; qx < totalQuadPoints; ++qx) {
			out[pick] += Qw[Sgroup];
			out += itemOutcomes[ix];
			Qw += numSpecific;
		}
	}
}

template <typename CovType>
void BA81Estep<CovType>::recordTable(struct BA81Expect *state)
{
	const int numThreads = Global->numThreads;
	const long expectedSize = state->totalQuadPoints * state->totalOutcomes;
	double *e1 = thrExpected.data();
	memcpy(state->expected, e1, sizeof(double) * expectedSize);
	e1 += expectedSize;

	for (int tx=1; tx < numThreads; ++tx) {
		for (long ex=0; ex < expectedSize; ++ex) {
			state->expected[ex] += *e1;
			++e1;
		}
	}
}

template <
  typename CovTypePar,
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
void BA81EngineBase<CovTypePar, LatentPolicy, EstepPolicy>::verboseLog(struct BA81Expect *state)
{
	if (state->verbose >= 1) {
		const int numUnique = (int) state->rowMap.size();
		mxLog("%s: estep(item version %d)<%s, %s, %s> %d/%d rows excluded",
		      state->name, omxGetMatrixVersion(state->itemParam),
		      typeid(CovTypePar).name(),
		      typeid(LatentPolicy).name(), typeid(EstepPolicy<CovTypePar>).name(),
		      state->excludedPatterns, numUnique);
		if (state->verbose >= 2) {
			omxPrint(state->estLatentMean, "mean");
			omxPrint(state->estLatentCov, "cov");
		}
	}
}

template <
  typename CovType,
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
double BA81EngineBase<CovType, LatentPolicy, EstepPolicy>::getPatLik(struct BA81Expect *state, int px, double *lxk)
{
	const int pts = BA81Config<CovType>::getPrimaryPoints(state);
	double *patternLik = state->patternLik;
	double patternLik1 = 0;

	for (int qx=0; qx < pts; qx++) {
		patternLik1 += lxk[qx];
	}

	// This uses the previous iteration's latent distribution.
	// If we recompute patternLikelihood to get the current
	// iteration's expected scores then it speeds up convergence.
	// However, recomputing patternLikelihood and dependent
	// math takes much longer than simply using the data
	// we have available here. This is even more true for the
	// two-tier model.
	if (!validPatternLik(state, patternLik1)) {
#pragma omp atomic
		state->excludedPatterns += 1;
		patternLik[px] = 0;
		return 0;
	}

	patternLik[px] = patternLik1;
	return patternLik1;
}

template <
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
struct BA81Engine<BA81Dense, LatentPolicy, EstepPolicy> :
	LatentPolicy, EstepPolicy<BA81Dense>, BA81EngineBase<BA81Dense, LatentPolicy, EstepPolicy> {
	typedef BA81Dense CovType;
	void ba81Estep1(struct BA81Expect *state);
};

template <
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
void BA81Engine<BA81Dense, LatentPolicy, EstepPolicy>::ba81Estep1(struct BA81Expect *state)
{
	const int numThreads = Global->numThreads;
	state->ptsPerThread = state->totalQuadPoints;
	state->primaryDims  = state->maxDims;
	const int numUnique = (int) state->rowMap.size();
	Eigen::VectorXd thrQweight;
	thrQweight.resize(state->ptsPerThread * numThreads);
	state->excludedPatterns = 0;
	state->patternLik = Realloc(state->patternLik, numUnique, double);
	double *patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<CovType>::begin(state);
	LatentPolicy::begin(state);

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omx_absolute_thread_num();
		double *Qweight = thrQweight.data() + state->ptsPerThread * thrId;
		ba81LikelihoodSlow2(state, px, Qweight);

		double patternLik1 = BA81Engine<BA81Dense, LatentPolicy, EstepPolicy>::getPatLik(state, px, Qweight);
		if (patternLik1 == 0) continue;

		LatentPolicy::normalizeWeights(state, px, Qweight, patternLik1, thrId);
		EstepPolicy<CovType>::addRow(state, px, Qweight, thrId);
	}

	// Can do these last steps in parallel, but it only saves 2%
	// in one test. Plus, this optimization is counterproductive
	// when EstepPolicy does nothing.

	EstepPolicy<CovType>::recordTable(state);
	LatentPolicy::end(state);
	BA81Engine<CovType, LatentPolicy, EstepPolicy>::verboseLog(state);
}

template <
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
struct BA81Engine<BA81TwoTier, LatentPolicy, EstepPolicy> :
	LatentPolicy, EstepPolicy<BA81TwoTier>, BA81EngineBase<BA81TwoTier, LatentPolicy, EstepPolicy> {
	typedef BA81TwoTier CovType;
	void ba81Estep1(struct BA81Expect *state);
};

template <
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
void BA81Engine<BA81TwoTier, LatentPolicy, EstepPolicy>::ba81Estep1(struct BA81Expect *state)
{
	const int numSpecific = state->numSpecific;
	const int numThreads = Global->numThreads;
	state->ptsPerThread = state->totalQuadPoints * numSpecific;
	state->primaryDims  = state->maxDims - 1;
	const int numUnique = (int) state->rowMap.size();
	Eigen::VectorXd thrQweight;
	thrQweight.resize(state->ptsPerThread * numThreads);
	state->excludedPatterns = 0;
	state->patternLik = Realloc(state->patternLik, numUnique, double);
	double *patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<CovType>::begin(state);
	LatentPolicy::begin(state);

	const long totalPrimaryPoints = state->totalPrimaryPoints;
	const long specificPoints = state->quadGridSize;
	omxBuffer<double> thrEi(totalPrimaryPoints * numThreads);
	omxBuffer<double> thrEis(totalPrimaryPoints * numSpecific * numThreads);

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omx_absolute_thread_num();
		double *Qweight = thrQweight.data() + state->ptsPerThread * thrId;
		double *Ei = thrEi.data() + totalPrimaryPoints * thrId;
		double *Eis = thrEis.data() + totalPrimaryPoints * numSpecific * thrId;
		cai2010EiEis(state, px, Qweight, Eis, Ei);

		double patternLik1 = BA81Engine<BA81TwoTier, LatentPolicy, EstepPolicy>::getPatLik(state, px, Ei);
		if (patternLik1 == 0) continue;

		// Can omit rest if we only want BA81RefreshPatLik TODO
		// Move Eis normalization loop here TODO

		for (long qloc=0, eisloc=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
			for (long sx=0; sx < specificPoints; sx++) {
				for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
					Qweight[qloc] *= Eis[eisloc + Sgroup];
					++qloc;
				}
			}
		}

		LatentPolicy::normalizeWeights(state, px, Qweight, patternLik1, thrId);
		EstepPolicy<CovType>::addRow(state, px, Qweight, thrId);
	}

	EstepPolicy<CovType>::recordTable(state);
	LatentPolicy::end(state);
	BA81Engine<CovType, LatentPolicy, EstepPolicy>::verboseLog(state);
}

static int getLatentVersion(BA81Expect *state)
{
	return omxGetMatrixVersion(state->latentMeanOut) + omxGetMatrixVersion(state->latentCovOut);
}

OMXINLINE static void
pointToWhere(BA81Expect *state, const int *quad, double *where, int upto)
{
	for (int dx=0; dx < upto; dx++) {
		where[dx] = state->Qpoint[quad[dx]];
	}
}

OMXINLINE static void
decodeLocation(long qx, const int dims, const long grid, int *quad)
{
	for (int dx=dims-1; dx >= 0; --dx) {
		quad[dx] = qx % grid;
		qx = qx / grid;
	}
}

struct sortAreaHelper {  // could be generalized with a template
	std::vector<double> &target;
	bool operator() (int i,int j) { return target[i] > target[j]; }
	sortAreaHelper(std::vector<double> &tgt) : target(tgt) {}
};

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
void ba81SetupQuadrature(omxExpectation* oo)
{
	// This is required because the EM acceleration can push the
	// covariance matrix to be slightly non-pd when predictors
	// are highly correlated.
	const bool forcePositiveSemiDefinite = true;

	BA81Expect *state = (BA81Expect *) oo->argStruct;
	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	if (state->Qpoint.size() == 0 && latentClean) return;

	if (state->verbose >= 1) {
		mxLog("%s: quadrature(%d)", oo->name, getLatentVersion(state));
		if (state->verbose >= 2) {
			pda(state->latentMeanOut->data, 1, state->maxAbilities);
			pda(state->latentCovOut->data, state->maxAbilities, state->maxAbilities);
		}
	}

	int gridsize = state->targetQpoints;
	const int maxDims = state->maxDims;
	double Qwidth = state->Qwidth;
	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);

	// try starting small and increasing to the cap? TODO
	state->quadGridSize = gridsize;

	state->totalQuadPoints = 1;
	for (int dx=0; dx < maxDims; dx++) {
		state->totalQuadPoints *= gridsize;
	}
	const long totalQuadPoints = state->totalQuadPoints;

	if (int(state->Qpoint.size()) != gridsize) {
		state->Qpoint.clear();
		state->Qpoint.reserve(gridsize);
		double qgs = gridsize-1;
		for (int px=0; px < gridsize; ++px) {
			state->Qpoint.push_back(Qwidth - px * 2 * Qwidth / qgs);
		}
	}

	std::vector<double> wherePrep(totalQuadPoints * maxDims);

	for (long qx=0; qx < totalQuadPoints; qx++) {
		double *wh = wherePrep.data() + qx * maxDims;
		int quad[maxDims];
		decodeLocation(qx, maxDims, gridsize, quad);
		pointToWhere(state, quad, wh, maxDims);
	}

	state->totalPrimaryPoints = totalQuadPoints;
	if (numSpecific) {
		state->totalPrimaryPoints /= gridsize;
		state->speQarea.resize(gridsize * numSpecific);
	}
	const int totalPrimaryPoints = state->totalPrimaryPoints;

	omxBuffer<double> priCovData(priDims * priDims);
	for (int d1=0; d1 < priDims; ++d1) {
		for (int d2=0; d2 < priDims; ++d2) {
			priCovData[d1 * priDims + d2] = omxMatrixElement(state->latentCovOut, d1, d2);
		}
	}
	if (forcePositiveSemiDefinite) {
		if (priDims == 1) {
			if (priCovData[0] < BA81_MIN_VARIANCE) priCovData[0] = BA81_MIN_VARIANCE;
		} else {
			Matrix mat(priCovData.data(), priDims, priDims);
			InplaceForcePosSemiDef(mat, NULL, NULL);
		}
	}

	const double Largest = state->LargestDouble;
	std::vector<double> priQarea;
	priQarea.reserve(totalPrimaryPoints);
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, state->quadGridSize, quad);
		double where[priDims];
		pointToWhere(state, quad, where, priDims);
		double den = exp(dmvnorm(priDims, where, state->latentMeanOut->data, priCovData.data()));
		priQarea.push_back(den);
	}

	std::vector<int> priOrder;
	priOrder.reserve(totalPrimaryPoints);
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		priOrder.push_back(qx);
	}
	sortAreaHelper priCmp(priQarea);
	std::sort(priOrder.begin(), priOrder.end(), priCmp);

	state->priQarea.clear();
	state->priQarea.reserve(totalPrimaryPoints);

	double totalArea = 0;
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		double den = priQarea[priOrder[qx]];
		state->priQarea.push_back(den);
		//double prevTotalArea = totalArea;
		totalArea += den;
		// if (totalArea == prevTotalArea) {
		// 	mxLog("%.4g / %.4g = %.4g", den, totalArea, den / totalArea);
		// }
	}

	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		state->priQarea[qx] *= Largest;
		state->priQarea[qx] /= totalArea;
		//mxLog("%.5g,", state->priQarea[qx]);
	}

	for (int sgroup=0; sgroup < numSpecific; sgroup++) {
		totalArea = 0;
		int covCell = (priDims + sgroup) * state->maxAbilities + priDims + sgroup;
		double mean = state->latentMeanOut->data[priDims + sgroup];
		double var = state->latentCovOut->data[covCell];
		if (forcePositiveSemiDefinite && var < BA81_MIN_VARIANCE) var = BA81_MIN_VARIANCE;
		//mxLog("setup[%d] %.2f %.2f", sx, mean, var);
		for (int qx=0; qx < gridsize; qx++) {
			double den = dnorm(state->Qpoint[qx], mean, sqrt(var), FALSE);
			state->speQarea[sIndex(state, sgroup, qx)] = den;
			totalArea += den;
		}
		for (int qx=0; qx < gridsize; qx++) {
			state->speQarea[sIndex(state, sgroup, qx)] *= Largest;
			state->speQarea[sIndex(state, sgroup, qx)] /= totalArea;
		}
		//pda(state->speQarea.data() + sIndex(state, sgroup, 0), 1, gridsize);
	}

	state->wherePrep.clear();
	state->wherePrep.reserve(totalQuadPoints * maxDims);

	if (numSpecific == 0) {
		for (int qx=0; qx < totalPrimaryPoints; qx++) {
			int sortq = priOrder[qx] * maxDims;
			for (int dx=0; dx < maxDims; ++dx) {
				state->wherePrep.push_back(wherePrep[sortq + dx]);
			}
		}
	} else {
		for (int qx=0; qx < totalPrimaryPoints; ++qx) {
			int sortq = priOrder[qx] * gridsize;
			for (int sx=0; sx < gridsize; ++sx) {
				int base = (sortq + sx) * maxDims;
				for (int dx=0; dx < maxDims; ++dx) {
					state->wherePrep.push_back(wherePrep[base + dx]);
				}
			}
		}
	}

	// recompute whereGram because the order might have changed
	std::vector<double> &whereGram = state->whereGram;
	whereGram.resize(totalQuadPoints * triangleLoc1(maxDims));

	for (int qx=0; qx < totalQuadPoints; qx++) {
		double *wh = state->wherePrep.data() + qx * maxDims;
		gramProduct(wh, maxDims, whereGram.data() + qx * triangleLoc1(maxDims));
	}

	//pda(wherePrep.data(), maxDims, totalQuadPoints);

	// The idea here is to avoid denormalized values if they are
	// enabled (5e-324 vs 2e-308).  It would be bad if results
	// changed depending on the denormalization setting.
	// Moreover, we don't lose too much even if denormalized
	// values are disabled.

	state->SmallestPatternLik = 1e16 * std::numeric_limits<double>::min();

	state->expected = Realloc(state->expected, state->totalOutcomes * totalQuadPoints, double);
	state->latentParamVersion = getLatentVersion(state);
}

static void
ba81compute(omxExpectation *oo, const char *what, const char *how)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;

	if (what) {
		if (strcmp(what, "latentDistribution")==0 && how && strcmp(how, "copy")==0) {
			omxCopyMatrix(state->latentMeanOut, state->estLatentMean);
			omxCopyMatrix(state->latentCovOut, state->estLatentCov);
			return;
		}

		if (strcmp(what, "scores")==0) {
			state->type = EXPECTATION_AUGMENTED;
		} else if (strcmp(what, "nothing")==0) {
			state->type = EXPECTATION_OBSERVED;
		} else {
			omxRaiseErrorf("%s: don't know how to predict '%s'",
				       oo->name, what);
		}
		return;
	}

	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	bool itemClean = state->itemParamVersion == omxGetMatrixVersion(state->itemParam) && latentClean;

	if (state->verbose >= 1) {
		mxLog("%s: Qinit %d itemClean %d latentClean %d (1=clean)",
		      oo->name, state->Qpoint.size() != 0, itemClean, latentClean);
	}

	if (!latentClean) ba81SetupQuadrature(oo);

	if (!itemClean) {
		ba81OutcomeProb(state, TRUE, FALSE);
		if (state->numSpecific == 0) {
			if (oo->dynamicDataSource) {
				BA81Engine<BA81Dense, BA81LatentSummary, BA81Estep> engine;
				engine.ba81Estep1(state);
			} else {
				BA81Engine<BA81Dense, BA81LatentFixed, BA81Estep> engine;
				engine.ba81Estep1(state);
			}
		} else {
			if (oo->dynamicDataSource) {
				BA81Engine<BA81TwoTier, BA81LatentSummary, BA81Estep> engine;
				engine.ba81Estep1(state);
			} else {
				BA81Engine<BA81TwoTier, BA81LatentFixed, BA81Estep> engine;
				engine.ba81Estep1(state);
			}
		}
	}

	state->itemParamVersion = omxGetMatrixVersion(state->itemParam);
}

/**
 * MAP is not affected by the number of items. EAP is. Likelihood can
 * get concentrated in a single quadrature ordinate. For 3PL, response
 * patterns can have a bimodal likelihood. This will confuse MAP and
 * is a key advantage of EAP (Thissen & Orlando, 2001, p. 136).
 *
 * Thissen, D. & Orlando, M. (2001). IRT for items scored in two
 * categories. In D. Thissen & H. Wainer (Eds.), \emph{Test scoring}
 * (pp 73-140). Lawrence Erlbaum Associates, Inc.
 */
static void
ba81PopulateAttributes(omxExpectation *oo, SEXP robj)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	int maxAbilities = state->maxAbilities;
	const int numUnique = (int) state->rowMap.size();

	Rf_setAttrib(robj, Rf_install("numStats"), Rf_ScalarReal(numUnique - 1)); // missingness? latent params? TODO

	if (state->debugInternal) {
		const double LogLargest = state->LogLargestDouble;
		int totalOutcomes = state->totalOutcomes;
		SEXP Rlik;
		SEXP Rexpected;

		Rf_protect(Rlik = Rf_allocVector(REALSXP, numUnique));
		memcpy(REAL(Rlik), state->patternLik, sizeof(double) * numUnique);
		double *lik_out = REAL(Rlik);
		for (int px=0; px < numUnique; ++px) {
			// Must return value in log units because it may not be representable otherwise
			lik_out[px] = log(lik_out[px]) - LogLargest;
		}

		Rf_protect(Rexpected = Rf_allocVector(REALSXP, state->totalQuadPoints * totalOutcomes));
		memcpy(REAL(Rexpected), state->expected, sizeof(double) * totalOutcomes * state->totalQuadPoints);

		MxRList dbg;
		dbg.add("patternLikelihood", Rlik);
		dbg.add("em.expected", Rexpected);

		SEXP Rmean, Rcov;
		Rf_protect(Rmean = Rf_allocVector(REALSXP, maxAbilities));
		memcpy(REAL(Rmean), state->estLatentMean->data, maxAbilities * sizeof(double));

		Rf_protect(Rcov = Rf_allocMatrix(REALSXP, maxAbilities, maxAbilities));
		memcpy(REAL(Rcov), state->estLatentCov->data, maxAbilities * maxAbilities * sizeof(double));

		dbg.add("mean", Rmean);
		dbg.add("cov", Rcov);

		Rf_setAttrib(robj, Rf_install("debug"), dbg.asR());
	}

	if (state->scores == SCORES_OMIT) return;

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	/*
	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		Rf_warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 0);
	ba81Estep1(oo);
	*/

	omxData *data = state->data;
	int rows = data->rows;
	int cols = 2 * maxAbilities + triangleLoc1(maxAbilities);
	state->scoresOut.clear();
	state->scoresOut.reserve(cols);
	SEXP Rscores;
	Rf_protect(Rscores = Rf_allocVector(VECSXP, cols));
	for (int cx=0; cx < cols; ++cx) {
		SEXP vec = Rf_allocVector(REALSXP, rows);
		SET_VECTOR_ELT(Rscores, cx, vec);
		state->scoresOut.push_back(REAL(vec));
	}

	const int SMALLBUF = 10;
	char buf[SMALLBUF];
	SEXP names;
	Rf_protect(names = Rf_allocVector(STRSXP, cols));
	for (int nx=0; nx < maxAbilities; ++nx) {
		SET_STRING_ELT(names, nx, Rf_mkChar(state->latentCovOut->rownames[nx]));
		snprintf(buf, SMALLBUF, "se%d", nx+1);
		SET_STRING_ELT(names, maxAbilities + nx, Rf_mkChar(buf));
	}
	for (int nx=0; nx < triangleLoc1(maxAbilities); ++nx) {
		snprintf(buf, SMALLBUF, "cov%d", nx+1);
		SET_STRING_ELT(names, maxAbilities*2 + nx, Rf_mkChar(buf));
	}
	Rf_setAttrib(Rscores, R_NamesSymbol, names);

	SEXP classes;
	Rf_protect(classes = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
	Rf_setAttrib(Rscores, R_ClassSymbol, classes);

	Rf_setAttrib(Rscores, R_RowNamesSymbol, data->getRowNames());

	if (state->numSpecific == 0) {
		BA81Engine<BA81Dense, BA81LatentScores, BA81OmitEstep> engine;
		engine.ba81Estep1(state);
	} else {
		BA81Engine<BA81TwoTier, BA81LatentScores, BA81OmitEstep> engine;
		engine.ba81Estep1(state);
	}

	MxRList out;
	out.add("scores", Rscores);
	Rf_setAttrib(robj, Rf_install("output"), out.asR());
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		mxLog("Freeing %s function.", oo->name);
	}
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	omxFreeMatrix(state->estLatentMean);
	omxFreeMatrix(state->estLatentCov);
	omxFreeMatrix(state->numObsMat);
	Free(state->patternLik);
	Free(state->Sgroup);
	Free(state->expected);
	Free(state->outcomeProb);
	if (state->ownWeights) Free(state->rowWeight);
	delete state;
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    Rf_protect(matrixDims = Rf_getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    Rf_unprotect(1);
}

static void ignoreSetVarGroup(omxExpectation*, FreeVarGroup *)
{}

static omxMatrix *getComponent(omxExpectation *oo, omxFitFunction*, const char *what)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;

	if (strcmp(what, "covariance")==0) {
		return state->estLatentCov;
	} else if (strcmp(what, "mean")==0) {
		return state->estLatentMean;
	} else if (strcmp(what, "numObs")==0) {
		return state->numObsMat;
	} else {
		return NULL;
	}
}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		mxLog("Initializing %s.", oo->name);
	}
	if (!rpf_model) {
		if (0) {
			const int wantVersion = 3;
			int version;
			get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model_GPL");
			(*get_librpf)(&version, &rpf_numModels, &rpf_model);
			if (version < wantVersion) Rf_error("librpf binary API %d installed, at least %d is required",
							 version, wantVersion);
		} else {
			rpf_numModels = librpf_numModels;
			rpf_model = librpf_model;
		}
	}
	
	BA81Expect *state = new BA81Expect;

	// These two constants should be as identical as possible
	state->name = oo->name;
	state->LogLargestDouble = log(std::numeric_limits<double>::max()) - 1;
	state->LargestDouble = exp(state->LogLargestDouble);
	state->OneOverLargestDouble = 1/state->LargestDouble;

	state->numObsMat = NULL;
	state->estLatentMean = NULL;
	state->estLatentCov = NULL;
	state->numSpecific = 0;
	state->excludedPatterns = 0;
	state->patternLik = NULL;
	state->outcomeProb = NULL;
	state->expected = NULL;
	state->type = EXPECTATION_OBSERVED;
	state->scores = SCORES_OMIT;
	state->itemParam = NULL;
	state->EitemParam = NULL;
	state->itemParamVersion = 0;
	state->latentParamVersion = 0;
	state->quadGridSize = 0;
	oo->argStruct = (void*) state;

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf("%s unable to handle data type %s", oo->name, omxDataType(state->data));
		return;
	}

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("ItemSpec")));
	for (int sx=0; sx < Rf_length(tmp); ++sx) {
		SEXP model = VECTOR_ELT(tmp, sx);
		if (!OBJECT(model)) {
			Rf_error("Item models must inherit rpf.base");
		}
		SEXP spec;
		Rf_protect(spec = R_do_slot(model, Rf_install("spec")));
		state->itemSpec.push_back(REAL(spec));
	}

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("design")));
	if (!Rf_isNull(tmp)) {
		int rows, cols;
		getMatrixDims(tmp, &rows, &cols);
		state->design.resize(rows, cols);
		memcpy(state->design.data(), INTEGER(tmp), sizeof(int) * rows * cols);
	}

	state->latentMeanOut = omxNewMatrixFromSlot(rObj, currentState, "mean");
	if (!state->latentMeanOut) Rf_error("Failed to retrieve mean matrix");

	state->latentCovOut  = omxNewMatrixFromSlot(rObj, currentState, "cov");
	if (!state->latentCovOut) Rf_error("Failed to retrieve cov matrix");

	state->itemParam = omxNewMatrixFromSlot(rObj, globalState, "ItemParam");

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("EItemParam")));
	if (!Rf_isNull(tmp)) {
		int rows, cols;
		getMatrixDims(tmp, &rows, &cols);
		if (rows != state->itemParam->rows || cols != state->itemParam->cols) {
			Rf_error("EItemParam must have same dimensions as ItemParam");
		}
		state->EitemParam = REAL(tmp);
	}

	oo->computeFun = ba81compute;
	oo->setVarGroup = ignoreSetVarGroup;
	oo->destructFun = ba81Destroy;
	oo->populateAttrFun = ba81PopulateAttributes;
	oo->componentFun = getComponent;
	oo->canDuplicate = false;
	
	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we get better cache behavior.
	// The following summary stats would be cheaper to calculate too.

	omxData *data = state->data;

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("weightColumn")));
	int weightCol = INTEGER(tmp)[0];
	state->ownWeights = weightCol == NA_INTEGER;
	if (state->ownWeights) {
		// Should rowMap be part of omxData? This is essentially a
		// generic compression step that shouldn't be specific to IFA models.
		state->rowWeight = Realloc(NULL, data->rows, double);
		state->rowMap.resize(data->rows);
		int numUnique = 0;
		for (int rx=0; rx < data->rows; ) {
			int rw = omxDataNumIdenticalRows(state->data, rx);
			state->rowWeight[numUnique] = rw;
			state->rowMap[numUnique] = rx;
			rx += rw;
			++numUnique;
		}
		state->rowMap.resize(numUnique);
		state->weightSum = state->data->rows;
	}
	else {
		if (omxDataColumnIsFactor(data, weightCol)) {
			omxRaiseErrorf("%s: weightColumn %d is a factor", oo->name, 1 + weightCol);
			return;
		}
		state->rowWeight = omxDoubleDataColumn(data, weightCol);
		state->weightSum = 0;
		for (int rx=0; rx < data->rows; ++rx) { state->weightSum += state->rowWeight[rx]; }
		state->rowMap.resize(data->rows);
		for (size_t rx=0; rx < state->rowMap.size(); ++rx) {
			state->rowMap[rx] = rx;
		}
	}
	// complain about non-integral rowWeights (EAP can't work) TODO

	std::vector<int> &rowMap = state->rowMap;

	const int numItems = state->itemParam->cols;
	if (state->itemSpec.size() == 1) {
		for (int ix=1; ix < numItems; ++ix) {
			state->itemSpec.push_back(state->itemSpec[0]);
		}
	}

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("dataColumns")));
	if (Rf_length(tmp) != numItems) Rf_error("dataColumns must be length %d", numItems);
	state->colMap = INTEGER(tmp);
	const int *colMap = state->colMap;

	int maxSpec = 0;
	int maxParam = 0;
	int maxItemDims = 0;

	std::vector<int> &itemOutcomes = state->itemOutcomes;
	std::vector<int> &cumItemOutcomes = state->cumItemOutcomes;
	itemOutcomes.resize(numItems);
	cumItemOutcomes.resize(numItems);
	int totalOutcomes = 0;
	for (int cx = 0; cx < numItems; cx++) {
		if (!omxDataColumnIsFactor(data, colMap[cx])) {
			omxRaiseErrorf("%s: column %d is not a factor", oo->name, 1 + colMap[cx]);
			return;
		}

		const double *spec = state->itemSpec[cx];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		if (maxItemDims < dims)
			maxItemDims = dims;

		int no = spec[RPF_ISpecOutcomes];
		itemOutcomes[cx] = no;
		cumItemOutcomes[cx] = totalOutcomes;
		totalOutcomes += no;

		// TODO this summary stat should be available from omxData
		int dataMax=0;
		for (int rx=0; rx < data->rows; rx++) {
			int pick = omxIntDataElementUnsafe(data, rx, colMap[cx]);
			if (dataMax < pick)
				dataMax = pick;
		}
		if (dataMax > no) {
			Rf_error("Data for item %d has %d outcomes, not %d", cx+1, dataMax, no);
		}

		int numSpec = (*rpf_model[id].numSpec)(spec);
		if (maxSpec < numSpec)
			maxSpec = numSpec;

		int numParam = (*rpf_model[id].numParam)(spec);
		if (maxParam < numParam)
			maxParam = numParam;
	}

	state->totalOutcomes = totalOutcomes;

	if (int(state->itemSpec.size()) != numItems) {
		omxRaiseErrorf("ItemSpec must contain %d item model specifications",
			       data->cols);
		return;
	}

	if (state->design.rows() == 0) {
		state->maxDims = maxItemDims;
		state->maxAbilities = maxItemDims;
		state->design.resize(state->maxDims, numItems);
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			for (int dx=0; dx < state->maxDims; dx++) {
				state->design(dx, ix) = dx < dims? dx+1 : NA_INTEGER;
			}
		}
	} else {
		Eigen::MatrixXi &design = state->design;
		if (design.cols() != numItems) {
			omxRaiseErrorf("Design matrix should have %d columns", numItems);
			return;
		}

		state->maxAbilities = design.maxCoeff();
		maxItemDims = 0;
		for (int ix=0; ix < design.cols(); ix++) {
			int ddim = 0;
			for (int rx=0; rx < design.rows(); rx++) {
				if (design(rx,ix) != NA_INTEGER) ddim += 1;
			}
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			if (ddim > dims) Rf_error("Item %d has %d dims but design assigns %d dims", ix, dims, ddim);
			if (maxItemDims < ddim) {
				maxItemDims = ddim;
			}
		}
		state->maxDims = maxItemDims;
	}
	if (state->maxAbilities <= state->maxDims) {
		state->Sgroup = Calloc(numItems, int);
	} else {
		// Not sure if this is correct, revisit TODO
		int Sgroup0 = -1;
		state->Sgroup = Realloc(NULL, numItems, int);
		for (int dx=0; dx < state->maxDims; dx++) {
			for (int ix=0; ix < numItems; ix++) {
				int ability = state->design(dx, ix);
				if (dx < state->maxDims - 1) {
					if (Sgroup0 <= ability)
						Sgroup0 = ability+1;
					continue;
				}
				int ss=-1;
				if (ability >= Sgroup0) {
					if (ss == -1) {
						ss = ability;
					} else {
						omxRaiseErrorf("Item %d cannot belong to more than "
							       "1 specific dimension (both %d and %d)",
							       ix, ss, ability);
						return;
					}
				}
				if (ss == -1) ss = Sgroup0;
				state->Sgroup[ix] = ss - Sgroup0;
			}
		}
		state->numSpecific = state->maxAbilities - state->maxDims + 1;
	}

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("naAction")));
	bool naFail = strEQ(CHAR(Rf_asChar(tmp)), "fail");

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("minItemsPerScore")));
	int minItemsPerScore = Rf_asInteger(tmp);
	if (minItemsPerScore > numItems) {
		omxRaiseErrorf("%s: minItemsPerScore (=%d) cannot be larger than the number of items (=%d)",
			       oo->name, minItemsPerScore, numItems);
		return;
	}

	state->rowSkip.assign(rowMap.size(), false);

	// Rows with no information about an ability will obtain the
	// prior distribution as an ability estimate. This will
	// throw off multigroup latent distribution estimates.
	for (size_t rx=0; rx < rowMap.size(); rx++) {
		std::vector<int> contribution(state->maxAbilities);
		for (int ix=0; ix < numItems; ix++) {
			int pick = omxIntDataElementUnsafe(data, rowMap[rx], colMap[ix]);
			if (pick == NA_INTEGER) continue;
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			int dr = 0;
			for (int dx=0; dx < dims; dx++) {
				int ability = state->design(dr + dx, ix);
				while (ability == NA_INTEGER) {
					++dr;
					ability = state->design(dr + dx, ix);
				}
				// assume factor loadings are the first item parameters
				if (omxMatrixElement(state->itemParam, dx, ix) == 0) continue;
				contribution[ability - 1] += 1;
			}
		}
		for (int ax=0; ax < state->maxAbilities; ++ax) {
			if (contribution[ax] < minItemsPerScore) {
				if (naFail) {
					int dest = omxDataIndex(data, state->rowMap[rx]);
					omxRaiseErrorf("Data row %d has no information about ability %d", 1+dest, 1+ax);
				}
				// We could compute the other scores, but estimation of the
				// latent distribution is in the hot code path. We can reconsider
				// this choice when we try generating scores instead of the
				// score distribution.
				state->rowSkip[rx] = true;
			}
		}
	}

	if (isErrorRaised()) return;

	if (state->latentMeanOut->rows * state->latentMeanOut->cols != state->maxAbilities) {
		Rf_error("The mean matrix '%s' must be 1x%d or %dx1", state->latentMeanOut->name,
		      state->maxAbilities, state->maxAbilities);
	}
	if (state->latentCovOut->rows != state->maxAbilities ||
	    state->latentCovOut->cols != state->maxAbilities) {
		Rf_error("The cov matrix '%s' must be %dx%d",
		      state->latentCovOut->name, state->maxAbilities, state->maxAbilities);
	}

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("verbose")));
	state->verbose = Rf_asInteger(tmp);

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("debugInternal")));
	state->debugInternal = Rf_asLogical(tmp);

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("qpoints")));
	state->targetQpoints = Rf_asReal(tmp);

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("qwidth")));
	state->Qwidth = Rf_asReal(tmp);

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("scores")));
	const char *score_option = CHAR(Rf_asChar(tmp));
	if (strEQ(score_option, "omit")) state->scores = SCORES_OMIT;
	if (strEQ(score_option, "full")) state->scores = SCORES_FULL;

	state->ElatentVersion = 0;
	state->estLatentMean = omxInitMatrix(state->maxAbilities, 1, TRUE, currentState);
	state->estLatentCov = omxInitMatrix(state->maxAbilities, state->maxAbilities, TRUE, currentState);
	omxCopyMatrix(state->estLatentMean, state->latentMeanOut); // rename matrices TODO
	omxCopyMatrix(state->estLatentCov, state->latentCovOut);
	state->numObsMat = omxInitMatrix(1, 1, TRUE, currentState);
	omxSetVectorElement(state->numObsMat, 0, data->rows);
}
