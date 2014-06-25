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
	const int totalQuadPoints = state->quad.totalQuadPoints;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	const size_t numItems = state->itemSpec.size();
	std::vector<int> &rowMap = state->rowMap;
	double *oProb = state->outcomeProb;
	std::vector<double> &priQarea = state->quad.priQarea;
	std::vector<int*> &dataColumns = state->dataColumns;

	for (int qx=0; qx < totalQuadPoints; ++qx) {
		out[qx] = priQarea[qx];
	}

	const int row = rowMap[px];
	for (size_t ix=0; ix < numItems; ix++) {
		int pick = dataColumns[ix][row];
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		for (int qx=0; qx < totalQuadPoints; ++qx) {
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
	const int totalQuadPoints = state->quad.totalQuadPoints;
	const int totalPrimaryPoints = state->quad.totalPrimaryPoints;
	const int specificPoints = state->quad.quadGridSize;
	const size_t numItems = state->itemSpec.size();
	const double OneOverLargest = state->OneOverLargestDouble;
	std::vector<int> &rowMap = state->rowMap;
	std::vector<double> &speQarea = state->quad.speQarea;
	std::vector<double> &priQarea = state->quad.priQarea;
	std::vector<int*> &dataColumns = state->dataColumns;

	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sx=0; sx < specificPoints * numSpecific; sx++) {
			lxk[qloc] = speQarea[sx];
			++qloc;
		}
	}

	const int row = rowMap[px];
	for (size_t ix=0; ix < numItems; ix++) {
		int pick = dataColumns[ix][row];
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;
		int Sgroup = state->Sgroup[ix];
		double *out1 = lxk;
		for (int qx=0; qx < state->quad.totalQuadPoints; qx++) {
			out1[Sgroup] *= oProb[pick];
			oProb += itemOutcomes[ix];
			out1 += numSpecific;
		}
	}

	for (int qx=0; qx < totalPrimaryPoints * numSpecific; ++qx) Eis[qx] = 0;
	for (int qx=0; qx < totalPrimaryPoints; ++qx) Ei[qx] = priQarea[qx];

	int eisloc = 0;
	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sx=0; sx < specificPoints; sx++) {
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
}

// Depends on item parameters, but not latent distribution
void ba81OutcomeProb(BA81Expect *state, double *param, bool wantLog)
{
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	std::vector<int> &cumItemOutcomes = state->cumItemOutcomes;
	int maxParam = state->maxParam;
	const int maxDims = state->maxDims;
	const size_t numItems = state->itemSpec.size();
	state->outcomeProb = Realloc(state->outcomeProb, state->totalOutcomes * state->quad.totalQuadPoints, double);

#pragma omp parallel for num_threads(Global->numThreads)
	for (size_t ix=0; ix < numItems; ix++) {
		double *qProb = state->outcomeProb + cumItemOutcomes[ix] * state->quad.totalQuadPoints;
		const double *spec = state->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = param + ix * maxParam;
		rpf_prob_t prob_fn = wantLog? rpf_model[id].logprob : rpf_model[id].prob;

		for (int qx=0; qx < state->quad.totalQuadPoints; qx++) {
			double *where = state->quad.wherePrep.data() + qx * maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, maxDims-1)];
			}

			(*prob_fn)(spec, iparam, ptheta.data(), qProb);
			qProb += itemOutcomes[ix];
		}
	}
}

void BA81LatentFixed::normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double patternLik1, int thrId)
{
	double weight = state->rowWeight[px] / patternLik1;
	for (int qx=0; qx < state->ptsPerThread; ++qx) {
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
	omxData *data = state->data;

	// NOTE: Qweight remains unnormalized

	double *scorePad = thrScore.data() + numLatents * thrId;
	OMXZERO(scorePad, numLatents);

	state->quad.EAP(Qweight, 1/patternLik1, scorePad);

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
	for (int qx=0; qx < state->ptsPerThread; ++qx) {
		double tmp = Qweight[qx] * weight;
		Dweight[qx] += tmp;
		Qweight[qx] = tmp;
	}
}

void BA81LatentSummary::end(struct BA81Expect *state)
{
	for (int tx=1; tx < Global->numThreads; ++tx) {
		double *Dweight = thrDweight.data() + state->ptsPerThread * tx;
		double *dest = thrDweight.data();
		for (int qx=0; qx < state->ptsPerThread; ++qx) {
			dest[qx] += Dweight[qx];
		}
	}

	state->quad.EAP(thrDweight.data(), 1/state->weightSum, latentDist.data());

	omxMatrix *meanOut = state->estLatentMean;
	omxMatrix *covOut = state->estLatentCov;
	const int maxAbilities = state->maxAbilities;
	const int primaryDims = state->primaryDims;

	double *latentDist1 = latentDist.data();
	for (int d1=0; d1 < maxAbilities; d1++) {
		omxSetVectorElement(meanOut, d1, latentDist1[d1]);
	}

	for (int d1=0; d1 < primaryDims; d1++) {
		int cx = maxAbilities + triangleLoc1(d1);
		for (int d2=0; d2 <= d1; d2++) {
			double cov = latentDist1[cx];
			omxSetMatrixElement(covOut, d1, d2, cov);
			if (d1 != d2) omxSetMatrixElement(covOut, d2, d1, cov);
			++cx;
		}
	}
	for (int d1=primaryDims; d1 < maxAbilities; d1++) {
		int loc = maxAbilities + triangleLoc0(d1);
		omxSetMatrixElement(covOut, d1, d1, latentDist1[loc]);
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
	return state->quad.totalQuadPoints;
}

template <>
int BA81Config<BA81TwoTier>::getPrimaryPoints(struct BA81Expect *state)
{
	return state->quad.totalPrimaryPoints;
}

template <typename CovType>
void BA81Estep<CovType>::begin(BA81Expect *state)
{
	totalQuadPoints = state->quad.totalQuadPoints;
	numItems = int(state->itemSpec.size());
	colMap = state->colMap;
	data = state->data;
	thrExpected.assign(state->totalOutcomes * state->quad.totalQuadPoints * Global->numThreads, 0.0);
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

		for (int qx=0; qx < totalQuadPoints; ++qx) {
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
		for (int qx=0; qx < totalQuadPoints; ++qx) {
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
	const int expectedSize = state->quad.totalQuadPoints * state->totalOutcomes;
	double *e1 = thrExpected.data();

	state->expected = Realloc(state->expected, state->totalOutcomes * state->quad.totalQuadPoints, double);
	memcpy(state->expected, e1, sizeof(double) * expectedSize);
	e1 += expectedSize;

	for (int tx=1; tx < numThreads; ++tx) {
		for (int ex=0; ex < expectedSize; ++ex) {
			state->expected[ex] += *e1;
			++e1;
		}
	}
}

template <typename CovType>
void BA81OmitEstep<CovType>::recordTable(struct BA81Expect *state)
{
	Free(state->expected);
}

template <
  typename CovTypePar,
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
void BA81EngineBase<CovTypePar, LatentPolicy, EstepPolicy>::engineDone(struct BA81Expect *state)
{
	state->expectedUsed = false;

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
	state->ptsPerThread = state->quad.totalQuadPoints;
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

	if (EstepPolicy<CovType>::hasEnd() && LatentPolicy::hasEnd()) {
#pragma omp parallel sections
		{
		{ EstepPolicy<CovType>::recordTable(state); }
#pragma omp section
		{ LatentPolicy::end(state); }
		}
	} else {
		EstepPolicy<CovType>::recordTable(state);
		LatentPolicy::end(state);
	}

	BA81Engine<CovType, LatentPolicy, EstepPolicy>::engineDone(state);
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
	state->ptsPerThread = state->quad.totalQuadPoints * numSpecific;
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

	const int totalPrimaryPoints = state->quad.totalPrimaryPoints;
	const int specificPoints = state->quad.quadGridSize;
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

		// Can omit rest if we only want patternLik TODO

		for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
			for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
				Eis[qloc] = Ei[qx] / Eis[qloc];
				++qloc;
			}
		}

		for (int qloc=0, eisloc=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
			for (int sx=0; sx < specificPoints; sx++) {
				for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
					Qweight[qloc] *= Eis[eisloc + Sgroup];
					++qloc;
				}
			}
		}

		LatentPolicy::normalizeWeights(state, px, Qweight, patternLik1, thrId);
		EstepPolicy<CovType>::addRow(state, px, Qweight, thrId);
	}

	if (EstepPolicy<CovType>::hasEnd() && LatentPolicy::hasEnd()) {
#pragma omp parallel sections
		{
		{ EstepPolicy<CovType>::recordTable(state); }
#pragma omp section
		{ LatentPolicy::end(state); }
		}
	} else {
		EstepPolicy<CovType>::recordTable(state);
		LatentPolicy::end(state);
	}

	BA81Engine<CovType, LatentPolicy, EstepPolicy>::engineDone(state);
}

static int getLatentVersion(BA81Expect *state)
{
	return omxGetMatrixVersion(state->latentMeanOut) + omxGetMatrixVersion(state->latentCovOut);
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
void ba81SetupQuadrature(omxExpectation* oo)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	if (state->quad.Qpoint.size() == 0 && latentClean) return;

	if (state->verbose >= 1) {
		mxLog("%s: quadrature(%d)", oo->name, getLatentVersion(state));
		if (state->verbose >= 2) {
			pda(state->latentMeanOut->data, 1, state->maxAbilities);
			pda(state->latentCovOut->data, state->maxAbilities, state->maxAbilities);
		}
	}

	const int maxDims = state->maxDims;
	if (maxDims == 0) {
		state->quad.setup0();
		state->latentParamVersion = getLatentVersion(state);
		return;
	}

	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);

	omxMatrix *inputCov = state->latentCovOut;
	Eigen::MatrixXd cov;
	cov.resize(priDims, priDims);

	for (int d1=0; d1 < priDims; ++d1) {
		for (int d2=0; d2 < priDims; ++d2) {
			cov(d1, d2) = omxMatrixElement(inputCov, d1, d2);
		}
	}

	// This is required because the EM acceleration can push the
	// covariance matrix to be slightly non-pd when predictors
	// are highly correlated.

	if (priDims == 1) {
		if (cov(0,0) < BA81_MIN_VARIANCE) cov(0,0) = BA81_MIN_VARIANCE;
	} else {
		Matrix mat(cov.data(), priDims, priDims);
		InplaceForcePosSemiDef(mat, NULL, NULL);
	}

	//pda(inputCov->data, inputCov->rows, inputCov->cols);

	Eigen::VectorXd sVar;
	sVar.resize(numSpecific);
	for (int sx=0; sx < numSpecific; ++sx) {
		int loc = priDims + sx;
		double tmp = inputCov->data[loc * inputCov->rows + loc];
		if (tmp < BA81_MIN_VARIANCE) tmp = BA81_MIN_VARIANCE;
		sVar(sx) = tmp;
	}

	state->quad.setup(state->Qwidth, state->targetQpoints, state->latentMeanOut->data, cov, sVar);

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
			state->expectedUsed = true;
		} else if (strcmp(what, "nothing")==0) {
			state->type = EXPECTATION_OBSERVED;
		} else {
			omxRaiseErrorf("%s: don't know how to predict '%s'",
				       oo->name, what);
		}

		if (state->verbose >= 1) {
			mxLog("%s: predict %s", oo->name, what);
		}
		return;
	}

	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	bool itemClean = state->itemParamVersion == omxGetMatrixVersion(state->itemParam) && latentClean;

	if (state->verbose >= 1) {
		mxLog("%s: Qinit %d itemClean %d latentClean %d (1=clean) expectedUsed=%d",
		      oo->name, state->quad.Qpoint.size() != 0, itemClean, latentClean, state->expectedUsed);
	}

	if (!latentClean) ba81SetupQuadrature(oo);

	if (!itemClean) {
		double *param = state->EitemParam? state->EitemParam : state->itemParam->data;
		ba81OutcomeProb(state, param, FALSE);

		if (state->expectedUsed) {
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
		} else {
			if (state->numSpecific == 0) {
				BA81Engine<BA81Dense, BA81LatentFixed, BA81OmitEstep> engine;
				engine.ba81Estep1(state);
			} else {
				BA81Engine<BA81TwoTier, BA81LatentFixed, BA81OmitEstep> engine;
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

		MxRList dbg;
		dbg.add("patternLikelihood", Rlik);

		if (state->expected) {
			Rf_protect(Rexpected = Rf_allocVector(REALSXP, state->quad.totalQuadPoints * totalOutcomes));
			memcpy(REAL(Rexpected), state->expected, sizeof(double) * totalOutcomes * state->quad.totalQuadPoints);
			dbg.add("em.expected", Rexpected);
		}

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
	state->quad.setOne(state->LargestDouble);
	state->OneOverLargestDouble = 1/state->LargestDouble;

	// The idea here is to avoid denormalized values if they are
	// enabled (5e-324 vs 2e-308).  It would be bad if results
	// changed depending on the denormalization setting.
	// Moreover, we don't lose too much even if denormalized
	// values are disabled. This mainly affects models with
	// more than a thousand items.
	state->SmallestPatternLik = 1e16 * std::numeric_limits<double>::min();
	state->expectedUsed = true;

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

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("verbose")));
	state->verbose = Rf_asInteger(tmp);

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
	// This is probably not worth it
	if (0 && state->itemSpec.size() == 1) {
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

	state->dataColumns.reserve(numItems);

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
		int no = spec[RPF_ISpecOutcomes];
		itemOutcomes[cx] = no;
		cumItemOutcomes[cx] = totalOutcomes;
		totalOutcomes += no;

		int *col = omxIntDataColumnUnsafe(data, colMap[cx]);
		state->dataColumns.push_back(col);

		// TODO this summary stat should be available from omxData
		int dataMax=0;
		for (int rx=0; rx < data->rows; rx++) {
			int pick = col[rx];
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

	if (state->itemParam->rows != maxParam) {
		omxRaiseErrorf("ItemParam matrix must have %d rows", maxParam);
		return;
	}

	state->maxParam = maxParam;
	state->totalOutcomes = totalOutcomes;

	if (int(state->itemSpec.size()) != numItems) {
		omxRaiseErrorf("ItemSpec must contain %d item model specifications",
			       data->cols);
		return;
	}

	state->maxAbilities = state->latentMeanOut->rows * state->latentMeanOut->cols;

	if (state->latentCovOut->rows != state->maxAbilities ||
	    state->latentCovOut->cols != state->maxAbilities) {
		Rf_error("The cov matrix '%s' must be %dx%d",
		      state->latentCovOut->name, state->maxAbilities, state->maxAbilities);
	}

	for (int cx = 0; cx < numItems; cx++) {
		const double *spec = state->itemSpec[cx];
		int dims = spec[RPF_ISpecDims];
		int loadings = 0;
		for (int dx=0; dx < dims; ++dx) {
			if (omxMatrixElement(state->itemParam, dx, cx) != 0) loadings += 1;
		}
		if (loadings > state->maxAbilities) {
			omxRaiseErrorf("Item %d has more factor loadings (%d) than there are factors (%d)",
				       1+cx, loadings, state->maxAbilities);
		}
	}

	// detect two-tier covariance structure
	std::vector<int> orthogonal;
	if (state->maxAbilities >= 3) {
		int mlen = state->maxAbilities;
		Eigen::Map<Eigen::MatrixXd> Ecov(state->latentCovOut->data, mlen, mlen);
		Eigen::Matrix<Eigen::DenseIndex, Eigen::Dynamic, 1> numCov((Ecov.array() != 0.0).matrix().colwise().count());
		std::vector<int> candidate;
		for (int fx=0; fx < numCov.rows(); ++fx) {
			if (numCov(fx) == 1) candidate.push_back(fx);
		}
		if (candidate.size() > 1) {
			std::vector<bool> mask(numItems);
			for (int cx=candidate.size() - 1; cx >= 0; --cx) {
				std::vector<bool> loading(numItems);
				for (int ix=0; ix < numItems; ++ix) {
					loading[ix] = state->itemParam->data[ix * maxParam + candidate[cx]] != 0;
				}
				std::vector<bool> overlap(loading.size());
				std::transform(loading.begin(), loading.end(),
					       mask.begin(), overlap.begin(),
					       std::logical_and<bool>());
				if (std::find(overlap.begin(), overlap.end(), true) == overlap.end()) {
					std::transform(loading.begin(), loading.end(),
						       mask.begin(), mask.begin(),
						       std::logical_or<bool>());
					orthogonal.push_back(candidate[cx]);
				}
			}
		}
		std::reverse(orthogonal.begin(), orthogonal.end());
	}
	if (orthogonal.size() && orthogonal[0] != state->maxAbilities - int(orthogonal.size())) {
		Rf_error("%s: Independent factors must be given after dense factors in %s",
			 oo->name, state->itemParam->name);
	}
	if (orthogonal.size() == 0) {
		state->maxDims = state->maxAbilities;
	} else {
		state->numSpecific = orthogonal.size();
		state->maxDims = state->maxAbilities - orthogonal.size() + 1;
		if (state->verbose >= 1) mxLog("%s: Two-tier structure detected; "
					       "%d abilities reduced to %d dimensions",
					       oo->name, state->maxAbilities, state->maxDims);
		state->Sgroup.assign(numItems, 0);
		for (int ix=0; ix < numItems; ix++) {
			for (int dx=orthogonal[0]; dx < state->maxAbilities; ++dx) {
				if (state->itemParam->data[ix * maxParam + dx] != 0) {
					state->Sgroup[ix] = dx - orthogonal[0];
					continue;
				}
			}
		}
	}

	// TODO: Items with zero loadings can be replaced with equivalent items
	// with fewer factors. This would speed up calculation of derivatives.

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

	if (state->maxAbilities) {
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
				for (int dx=0; dx < dims; dx++) {
					// assume factor loadings are the first item parameters
					if (omxMatrixElement(state->itemParam, dx, ix) == 0) continue;
					contribution[dx] += 1;
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
	}

	if (isErrorRaised()) return;

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("debugInternal")));
	state->debugInternal = Rf_asLogical(tmp);

	Rf_protect(tmp = R_do_slot(rObj, Rf_install("qpoints")));
	state->targetQpoints = Rf_asInteger(tmp);

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
