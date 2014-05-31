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

// http://en.wikipedia.org/wiki/Curiously_Recurring_Template_Pattern

struct BA81Engine {
	struct BA81Expect *conf;

	BA81Engine(struct BA81Expect *_conf) : conf(_conf) {};
	virtual void ba81Estep1() = 0;
	virtual ~BA81Engine() {};
};

template <typename Derived>
struct BA81EngineBase : BA81Engine {
	typedef BA81Engine super;
	int numSpecific;
	int ptsPerThread;
	Eigen::VectorXd thrQweight;
	const int numThreads;
	const int totalQuadPoints;
	const int maxDims;
	const int maxAbilities;
	const double *wherePrep;
	const double *whereGram;
	const int whereGramSize;
	const int primaryDims;

	BA81EngineBase(BA81Expect *_conf);
	virtual void ba81Estep1();
	void startEstep() { static_cast<Derived*>(this)->startEstep(); };
	double *getThrQweight(int thrId);
	void normalizeWeights(int thrId, double weight) {
		static_cast<Derived*>(this)->normalizeWeights(thrId, weight);
	};
	void weightsToLatentDistribution() {
		static_cast<Derived*>(this)->weightsToLatentDistribution();
	};
	void recordLatentDistribution() {
		static_cast<Derived*>(this)->recordLatentDistribution();
	};
};

struct BA81EngineLatentFixed : BA81EngineBase<BA81EngineLatentFixed> {
	typedef BA81EngineBase<BA81EngineLatentFixed> super;
 	BA81EngineLatentFixed(BA81Expect *_conf) : super(_conf) {};

	void startEstep();
	void normalizeWeights(int thrId, double weight);
	void weightsToLatentDistribution() {};
	void recordLatentDistribution() {};
};

struct BA81EngineLatentFree : BA81EngineBase<BA81EngineLatentFree> {
	typedef BA81EngineBase<BA81EngineLatentFree> super;

	const int numLatents;
	std::vector<double> thrDweight;
	std::vector<double> latentDist;

	BA81EngineLatentFree(BA81Expect *_conf) : super(_conf),
		numLatents(maxAbilities + triangleLoc1(maxAbilities)) {};
	void startEstep();
	void normalizeWeights(int thrId, double weight);
	void weightsToLatentDistribution();
	void mapLatentSpace(int sgroup, double piece, const double *where,
			    const double *whereGram, double *latentDist);
	void recordLatentDistribution();
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
	int *Sgroup;              // item's specific group 0..numSpecific-1
	Eigen::MatrixXi design;        // items * maxDims

	// quadrature related
	double Qwidth;
	double targetQpoints;
	long quadGridSize;
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
