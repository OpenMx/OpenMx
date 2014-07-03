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

#ifndef _BA81QUAD_H_
#define _BA81QUAD_H_

#include "glue.h"
#include "Eigen/Core"
#include "libifa-rpf.h"

class ba81NormalQuad {
 private:
	inline void pointToWhere(const int *quad, double *where, int upto);
	inline void decodeLocation(int qx, const int dims, int *quad);
	double One, ReciprocalOfOne;

	inline int sIndex(int sx, int qx) {
		//if (sx < 0 || sx >= state->numSpecific) Rf_error("Out of domain");
		//if (qx < 0 || qx >= state->quadGridSize) Rf_error("Out of domain");
		return qx * numSpecific + sx;
	};

	inline void mapDenseSpace(double piece, const double *where,
				  const double *whereGram, double *latentDist);
	inline void mapSpecificSpace(int sgroup, double piece, const double *where,
				     const double *whereGram, double *latentDist);

 public:
	int quadGridSize;                     // rename to gridSize TODO
	int maxDims;
	int primaryDims;
	int numSpecific;
	int maxAbilities;
	std::vector<double> Qpoint;           // quadGridSize
	int totalQuadPoints;                  // quadGridSize ^ maxDims
	int totalPrimaryPoints;               // totalQuadPoints except for specific dim
	int weightTableSize;                  // dense: totalQuadPoints; 2tier: totalQuadPoints * numSpecific
	std::vector<double> priQarea;         // totalPrimaryPoints
	std::vector<double> speQarea;         // quadGridSize * numSpecific
	std::vector<double> wherePrep;        // totalQuadPoints * maxDims
	Eigen::MatrixXd whereGram;            // triangleLoc1(maxDims) x totalQuadPoints

	ba81NormalQuad();
	void setOne(double one) { One = one; ReciprocalOfOne = 1/one; }
	void setup0();
	void setup(double Qwidth, int Qpoints, double *means,
		   Eigen::MatrixXd &priCov, Eigen::VectorXd &sVar);
	inline double getReciprocalOfOne() const { return ReciprocalOfOne; };

	// For dense cov, Dweight is size totalQuadPoints
	// For two-tier, Dweight is numSpecific x totalQuadPoints
	inline void EAP(double *thrDweight, double scalingFactor, double *scorePad);
};

void ba81NormalQuad::mapDenseSpace(double piece, const double *where,
				   const double *whereGram, double *latentDist)
{
	const int pmax = primaryDims;
	int gx = 0;
	int cx = maxAbilities;
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

void ba81NormalQuad::mapSpecificSpace(int sgroup, double piece, const double *where,
				      const double *whereGram, double *latentDist)
{
	const int pmax = primaryDims;

	int sdim = pmax + sgroup;
	double piece_w1 = piece * where[pmax];
	latentDist[sdim] += piece_w1;

	double piece_var = piece * whereGram[triangleLoc0(pmax)];
	int to = maxAbilities + triangleLoc0(sdim);
	latentDist[to] += piece_var;
}

void ba81NormalQuad::EAP(double *thrDweight, double scalingFactor, double *scorePad)
{
	if (numSpecific == 0) { // use template to handle this branch at compile time? TODO
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			mapDenseSpace(thrDweight[qx], &wherePrep[qx * maxDims],
				      &whereGram.coeffRef(0, qx), scorePad);
		}
	} else {
		int qloc=0;
		for (int qx=0; qx < totalQuadPoints; qx++) {
			const double *whPrep = &wherePrep[qx * maxDims];
			const double *whGram = &whereGram.coeffRef(0, qx);
			mapDenseSpace(thrDweight[qloc], whPrep, whGram, scorePad);
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapSpecificSpace(Sgroup, thrDweight[qloc], whPrep, whGram, scorePad);
				++qloc;
			}
		}
	}

	const int padSize = maxAbilities + triangleLoc1(maxAbilities);
	for (int d1=0; d1 < padSize; d1++) {
		scorePad[d1] *= scalingFactor;
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
	for (int sx=0; sx < numSpecific; sx++) {
		int sdim = primaryDims + sx;
		double ma1 = scorePad[sdim];
		scorePad[maxAbilities + triangleLoc0(sdim)] -= ma1 * ma1;
	}
}

class ifaGroup {
 private:
	SEXP Rdata;
	void verifyFactorNames(SEXP mat, const char *matName);
 public:
	const int numThreads;

	// item description related
	std::vector<const double*> spec;
	int maxItemDims;
	int numItems() const { return (int) spec.size(); }
	int paramRows;
	double *param;  // itemParam->data
	std::vector<const char*> itemNames;
	std::vector<int> itemOutcomes;
	std::vector<int> cumItemOutcomes;
	int totalOutcomes;
	std::vector<int> Sgroup;       // item's specific group 0..numSpecific-1

	// latent distribution
	double qwidth;
	int qpoints;
	ba81NormalQuad quad;
	ba81NormalQuad &getQuad() { return quad; };
	bool twotier;  // rename to detectTwoTier TODO
	int maxAbilities;
	int numSpecific;
	double *mean;
	double *cov;
	std::vector<const char*> factorNames;

	// data related
	SEXP dataRowNames;
	std::vector<const int*> dataColumns;
	std::vector<int> rowMap;       // row index into MxData
	int getNumUnique() const { return (int) rowMap.size(); }
	const char *weightColumnName;
	double *rowWeight;
 private:
	int minItemsPerScore;
 public:
	void setMinItemsPerScore(int mips);
	std::vector<bool> rowSkip;     // whether to treat the row as NA

	// workspace
	double *outcomeProb;                  // totalOutcomes * totalQuadPoints
	static const double SmallestPatternLik;
	int excludedPatterns;
	Eigen::ArrayXd patternLik;            // numUnique

	inline static bool validPatternLik(double pl)
	{ return std::isfinite(pl) && pl > SmallestPatternLik; }

	// TODO:
	// scores

	ifaGroup(int cores, bool _twotier);
	~ifaGroup();
	void setGridFineness(double width, int points);
	void import(SEXP Rlist);
	void importSpec(SEXP slotValue);
	void learnMaxAbilities();
	void setLatentDistribution(int dims, double *mean, double *cov);
	inline double *getItemParam(int ix) { return param + paramRows * ix; }
	inline const int *dataColumn(int col) { return dataColumns[col]; };
	void detectTwoTier();
	void buildRowSkip();
	void sanityCheck();
	inline void ba81OutcomeProb(double *param, bool wantLog);
	inline void ba81LikelihoodSlow2(const int px, double *out);
	inline void cai2010EiEis(const int px, double *lxk, double *Eis, double *Ei);
	inline void cai2010part2(double *Qweight, double *Eis, double *Ei);
};

// Depends on item parameters, but not latent distribution
void ifaGroup::ba81OutcomeProb(double *param, bool wantLog)
{
	const int maxDims = quad.maxDims;
	outcomeProb = Realloc(outcomeProb, totalOutcomes * quad.totalQuadPoints, double);

#pragma omp parallel for num_threads(numThreads)
	for (int ix=0; ix < numItems(); ix++) {
		double *qProb = outcomeProb + cumItemOutcomes[ix] * quad.totalQuadPoints;
		const double *ispec = spec[ix];
		int id = ispec[RPF_ISpecID];
		int dims = ispec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = param + paramRows * ix;
		rpf_prob_t prob_fn = wantLog? librpf_model[id].logprob : librpf_model[id].prob;

		for (int qx=0; qx < quad.totalQuadPoints; qx++) {
			double *where = quad.wherePrep.data() + qx * maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, maxDims-1)];
			}

			(*prob_fn)(ispec, iparam, ptheta.data(), qProb);
			qProb += itemOutcomes[ix];
		}
	}
}

void ifaGroup::ba81LikelihoodSlow2(const int px, double *out)
{
	const int totalQuadPoints = quad.totalQuadPoints;
	double *oProb = outcomeProb;
	std::vector<double> &priQarea = quad.priQarea;

	for (int qx=0; qx < totalQuadPoints; ++qx) {
		out[qx] = priQarea[qx];
	}

	const int row = rowMap[px];
	for (int ix=0; ix < numItems(); ix++) {
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

void ifaGroup::cai2010EiEis(const int px, double *lxk, double *Eis, double *Ei)
{
	double *oProb = outcomeProb;
	const int totalQuadPoints = quad.totalQuadPoints;
	const int totalPrimaryPoints = quad.totalPrimaryPoints;
	const int specificPoints = quad.quadGridSize;
	std::vector<double> &speQarea = quad.speQarea;
	std::vector<double> &priQarea = quad.priQarea;

	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sx=0; sx < specificPoints * numSpecific; sx++) {
			lxk[qloc] = speQarea[sx];
			++qloc;
		}
	}

	const int row = rowMap[px];
	for (int ix=0; ix < numItems(); ix++) {
		int pick = dataColumns[ix][row];
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;
		int Sgroup1 = Sgroup[ix];
		double *out1 = lxk;
		for (int qx=0; qx < quad.totalQuadPoints; qx++) {
			out1[Sgroup1] *= oProb[pick];
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
			Ei[qx] *= Eis[eisloc + sgroup] * quad.getReciprocalOfOne();
		}
		eisloc += numSpecific;
	}
}

void ifaGroup::cai2010part2(double *Qweight, double *Eis, double *Ei)
{
	const int totalPrimaryPoints = quad.totalPrimaryPoints;
	const int specificPoints = quad.quadGridSize;

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
}

struct BA81Dense {};
struct BA81TwoTier {};

struct BA81EngineBase {
	inline int getPrimaryPoints(class ifaGroup *state) { return state->quad.totalPrimaryPoints; };
	inline double getPatLik(class ifaGroup *state, int px, double *lxk);
};


double BA81EngineBase::getPatLik(class ifaGroup *state, int px, double *lxk)
{
	const int pts = getPrimaryPoints(state);
	Eigen::ArrayXd &patternLik = state->patternLik;
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
	if (!ifaGroup::validPatternLik(patternLik1)) {
#pragma omp atomic
		state->excludedPatterns += 1;
		patternLik[px] = 0;
		return 0;
	}

	patternLik[px] = patternLik1;
	return patternLik1;
}

template <typename T, typename CovType>
struct BA81OmitEstep {
	void begin(class ifaGroup *state, T extraData) {};
	void addRow(class ifaGroup *state, T extraData, int px, double *Qweight, int thrId) {};
	void recordTable(class ifaGroup *state, T extraData) {};
	bool hasEnd() { return false; }
};

template <
  typename T,
  typename CovTypePar,
  template <typename> class LatentPolicy,
  template <typename, typename> class EstepPolicy
>
struct BA81Engine : LatentPolicy<T>, EstepPolicy<T, CovTypePar>, BA81EngineBase {
	void ba81Estep1(class ifaGroup *state, T extraData);
};


template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename, typename> class EstepPolicy
>
struct BA81Engine<T, BA81Dense, LatentPolicy, EstepPolicy> :
	LatentPolicy<T>, EstepPolicy<T, BA81Dense>, BA81EngineBase {
	typedef BA81Dense CovType;
	void ba81Estep1(class ifaGroup *state, T extraData);
};

template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename, typename> class EstepPolicy
>
void BA81Engine<T, BA81Dense, LatentPolicy, EstepPolicy>::ba81Estep1(class ifaGroup *state, T extraData)
{
	ba81NormalQuad &quad = state->getQuad();
	const int numUnique = state->getNumUnique();
	const int numThreads = state->numThreads;
	Eigen::VectorXd thrQweight;
	thrQweight.resize(quad.weightTableSize * numThreads);
	state->excludedPatterns = 0;
	state->patternLik.resize(numUnique);
	Eigen::ArrayXd &patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<T, CovType>::begin(state, extraData);
	LatentPolicy<T>::begin(state, extraData);

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omp_get_thread_num();
		double *Qweight = thrQweight.data() + quad.weightTableSize * thrId;
		state->ba81LikelihoodSlow2(px, Qweight);

		double patternLik1 = getPatLik(state, px, Qweight);
		if (patternLik1 == 0) continue;

		LatentPolicy<T>::normalizeWeights(state, extraData, px, Qweight, patternLik1, thrId);
		EstepPolicy<T, CovType>::addRow(state, extraData, px, Qweight, thrId);
	}

	if (EstepPolicy<T, CovType>::hasEnd() && LatentPolicy<T>::hasEnd()) {
#pragma omp parallel sections
		{
			{ EstepPolicy<T, CovType>::recordTable(state, extraData); }
#pragma omp section
			{ LatentPolicy<T>::end(state, extraData); }
		}
	} else {
		EstepPolicy<T, CovType>::recordTable(state, extraData);
		LatentPolicy<T>::end(state, extraData);
	}
}

template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename, typename> class EstepPolicy
>
struct BA81Engine<T, BA81TwoTier, LatentPolicy, EstepPolicy> :
	LatentPolicy<T>, EstepPolicy<T, BA81TwoTier>, BA81EngineBase {
	typedef BA81TwoTier CovType;
	void ba81Estep1(class ifaGroup *state, T extraData);
};

template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename, typename> class EstepPolicy
>
void BA81Engine<T, BA81TwoTier, LatentPolicy, EstepPolicy>::ba81Estep1(class ifaGroup *state, T extraData)
{
	ba81NormalQuad &quad = state->getQuad();
	const int numSpecific = quad.numSpecific;
	const int numUnique = state->getNumUnique();
	const int numThreads = state->numThreads;
	Eigen::VectorXd thrQweight;
	thrQweight.resize(quad.weightTableSize * numThreads);
	state->excludedPatterns = 0;
	state->patternLik.resize(numUnique);
	Eigen::ArrayXd &patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<T, CovType>::begin(state, extraData);
	LatentPolicy<T>::begin(state, extraData);

	const int totalPrimaryPoints = quad.totalPrimaryPoints;
	Eigen::ArrayXXd thrEi(totalPrimaryPoints, numThreads);
	Eigen::ArrayXXd thrEis(totalPrimaryPoints * numSpecific, numThreads);

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omp_get_thread_num();
		double *Qweight = thrQweight.data() + quad.weightTableSize * thrId;
		double *Ei = &thrEi.coeffRef(0, thrId);
		double *Eis = &thrEis.coeffRef(0, thrId);
		state->cai2010EiEis(px, Qweight, Eis, Ei);

		double patternLik1 = getPatLik(state, px, Ei);
		if (patternLik1 == 0) continue;

		if (!EstepPolicy<T, CovType>::hasEnd() && !LatentPolicy<T>::hasEnd()) continue;

		state->cai2010part2(Qweight, Eis, Ei);

		LatentPolicy<T>::normalizeWeights(state, extraData, px, Qweight, patternLik1, thrId);
		EstepPolicy<T, CovType>::addRow(state, extraData, px, Qweight, thrId);
	}

	if (EstepPolicy<T, CovType>::hasEnd() && LatentPolicy<T>::hasEnd()) {
#pragma omp parallel sections
		{
			{ EstepPolicy<T, CovType>::recordTable(state, extraData); }
#pragma omp section
			{ LatentPolicy<T>::end(state, extraData); }
		}
	} else {
		EstepPolicy<T, CovType>::recordTable(state, extraData);
		LatentPolicy<T>::end(state, extraData);
	}
}

#endif
