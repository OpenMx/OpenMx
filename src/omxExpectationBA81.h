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

#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

#include "omxExpectation.h"
#include "omxOpenmpWrap.h"

enum score_option {
	SCORES_OMIT,
	SCORES_FULL
};

enum expectation_type {
	EXPECTATION_AUGMENTED, // E-M
	EXPECTATION_OBSERVED,  // regular
};

struct BA81Dense {};
struct BA81TwoTier {};

template <typename CovType>
struct BA81Config {
	int getPrimaryPoints(struct BA81Expect *state);
};

template <typename CovType>
struct BA81PatLik : BA81Config<CovType> {
	double getPatLik(struct BA81Expect *state, int px, double *lxk);
};

template <typename CovType>
struct BA81RefreshPatLik : BA81Config<CovType> {
	void begin(struct BA81Expect *state);
	bool skipRow(struct BA81Expect *state, int px);
	double getPatLik(struct BA81Expect *state, int px, double *lxk);
};

template <typename CovType>
struct BA81Estep {
	std::vector<double> thrExpected;
	int numItems;
	int totalQuadPoints;
	const int *colMap;
	omxData *data;

	void begin(struct BA81Expect *state);
	void addRow(struct BA81Expect *state, int px, double *Qweight, int thrId);
	void recordTable(struct BA81Expect *state);
};

template <typename CovType>
struct BA81OmitEstep {
	void begin(struct BA81Expect *state) {};
	void addRow(struct BA81Expect *state, int px, double *Qweight, int thrId) {};
	void recordTable(struct BA81Expect *state) {};
};

struct BA81LatentFixed {
	void begin(struct BA81Expect *state) {}
	void normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double weight, int thrid);
	void end(struct BA81Expect *state) {};
};

struct BA81LatentEstimate {
	void mapDenseSpace(struct BA81Expect *state, double piece, const double *where,
			   const double *whereGram, double *latentDist);
	void mapSpecificSpace(struct BA81Expect *state, int sgroup, double piece, const double *where,
			      const double *whereGram, double *latentDist);
	void mapSpace(struct BA81Expect *state, double *thrDweight, double *latentDist);
};

struct BA81LatentSummary : BA81LatentEstimate {
	int numLatents;
	std::vector<double> thrDweight;
	std::vector<double> latentDist;

	void begin(struct BA81Expect *state);
	void normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double weight, int thrId);
	void end(struct BA81Expect *state);
};

struct BA81LatentScores : BA81LatentEstimate {
	int numLatents;
	Eigen::VectorXd thrScore;

	void begin(struct BA81Expect *state);
	void normalizeWeights(struct BA81Expect *state, int px, double *Qweight, double weight, int thrId);
	void end(struct BA81Expect *state);
};

template <
  typename CovTypePar,
  typename LatentPolicy,
  template <typename> class EstepPolicy
>
struct BA81Engine : LatentPolicy, EstepPolicy<CovTypePar>, BA81PatLik<CovTypePar> {
	double getPatLik(struct BA81Expect *state, int px, double *lxk);
	void ba81Estep1(struct BA81Expect *state);
};

struct BA81Expect {
	double LogLargestDouble;       // should be const but need constexpr
	double LargestDouble;          // should be const but need constexpr
	double OneOverLargestDouble;   // should be const but need constexpr

	// data characteristics
	omxData *data;
	const int *colMap;             // item column to data column mapping
	std::vector<int> rowMap;       // row index into MxData
	double *rowWeight;
	double weightSum;              // sum of rowWeight
	bool ownWeights;
	std::vector<bool> rowSkip;     // whether to treat the row as NA

	// item description related
	std::vector<const double*> itemSpec;
	std::vector<int> itemOutcomes;
	std::vector<int> cumItemOutcomes;
	int maxDims;
	int maxAbilities;
	int numSpecific;
	int *Sgroup;                   // item's specific group 0..numSpecific-1
	Eigen::MatrixXi design;        // items * maxDims

	// quadrature related
	double Qwidth;
	double targetQpoints;
	long quadGridSize;  // change to int type TODO
	long totalQuadPoints;                 // quadGridSize ^ maxDims
	long totalPrimaryPoints;              // totalQuadPoints except for specific dim TODO
	std::vector<double> wherePrep;        // totalQuadPoints * maxDims
	std::vector<double> whereGram;        // totalQuadPoints * triangleLoc1(maxDims)
	std::vector<double> Qpoint;           // quadGridSize
	std::vector<double> priQarea;         // totalPrimaryPoints
	std::vector<double> speQarea;         // quadGridSize * numSpecific

	// estimation related
	omxMatrix *itemParam;
	double *EitemParam;
	double *patternLik;                   // numUnique
	double SmallestPatternLik;
	int excludedPatterns;
	int totalOutcomes;
	double *outcomeProb;                  // totalOutcomes * totalQuadPoints
	double *expected;                     // totalOutcomes * totalQuadPoints
	int ElatentVersion;
	omxMatrix *latentMeanOut;
	omxMatrix *latentCovOut;
	omxMatrix *estLatentMean;
	omxMatrix *estLatentCov;
	omxMatrix *numObsMat; // this is dumb

	int itemParamVersion;
	int latentParamVersion;
	enum expectation_type type;
	enum score_option scores;
	int verbose;
	bool debugInternal;
	struct omxFitFunction *fit;  // weak pointer

	// workspace
	int ptsPerThread;
	int primaryDims;
	std::vector<double*> scoresOut;
};

extern const struct rpf *rpf_model;
extern int rpf_numModels;

void ba81OutcomeProb(BA81Expect *state, bool estep, bool wantLog);

// state->speQarea[sIndex(state, sgroup, sx)]
OMXINLINE static
int sIndex(BA81Expect *state, int sx, int qx)
{
	//if (sx < 0 || sx >= state->numSpecific) Rf_error("Out of domain");
	//if (qx < 0 || qx >= state->quadGridSize) Rf_error("Out of domain");
	return qx * state->numSpecific + sx;
}

OMXINLINE static void
gramProduct(double *vec, size_t len, double *out)
{
	int cell = 0;
	for (size_t v1=0; v1 < len; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			out[cell] = vec[v1] * vec[v2];
			++cell;
		}
	}
}

OMXINLINE static bool
validPatternLik(BA81Expect *state, double pl)
{
	return std::isfinite(pl) && pl > state->SmallestPatternLik;
}

void ba81SetupQuadrature(omxExpectation* oo);
void ba81LikelihoodSlow2(BA81Expect *state, int px, double *out);
void cai2010EiEis(BA81Expect *state, int px, double *lxk, double *Eis, double *Ei);
static const double BA81_MIN_VARIANCE = .01;

#endif
