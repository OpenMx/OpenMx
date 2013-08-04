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

#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

#include "omxExpectation.h"
#include "omxOpenmpWrap.h"

enum score_option {
	SCORES_OMIT,
	SCORES_UNIQUE,
	SCORES_FULL
};

enum expectation_type {
	EXPECTATION_UNINITIALIZED,
	EXPECTATION_AUGMENTED, // E-M
	EXPECTATION_OBSERVED,  // regular
};

typedef struct {

	// data characteristics
	omxData *data;
	int numUnique;
	int *numIdentical;        // length numUnique
	int *rowMap;              // length numUnique, index of first instance of pattern

	// item description related
	std::vector<const double*> itemSpec;
	std::vector<int> itemOutcomes;
	int maxDims;
	int maxAbilities;
	int numSpecific;
	int *Sgroup;              // item's specific group 0..numSpecific-1
	omxMatrix *design;        // items * maxDims

	// quadrature related
	double Qwidth;
	double targetQpoints;
	long quadGridSize;
	long totalQuadPoints;     // quadGridSize ^ maxDims
	long totalPrimaryPoints;  // totalQuadPoints except for specific dim TODO
	std::vector<double> Qpoint;           // quadGridSize
	std::vector<double> priQarea;         // totalPrimaryPoints
	std::vector<double> speQarea;         // quadGridSize * numSpecific

	// estimation related
	omxMatrix *customPrior;
	omxMatrix *itemParam;
	omxMatrix *EitemParam;    // E step version
	int cacheLXK;
	bool LXKcached;
	double *lxk;              // wo/cache, numUnique * thread
	double *allElxk;          // numUnique * thread
	double *Eslxk;            // numUnique * #specific dimensions * thread
	double *patternLik;       // numUnique
	int excludedPatterns;
	int totalOutcomes;
	double *outcomeProb;      // totalOutcomes * totalQuadPoints
	double *expected;         // totalOutcomes * totalQuadPoints
	std::vector<double> ElatentMean;      // maxAbilities
	std::vector<double> ElatentCov;       // maxAbilities * maxAbilities
	omxMatrix *latentMeanOut;
	omxMatrix *latentCovOut;
	double *EiCache;          // totalPrimaryPoints * numUnique (two-tier only)

	int itemParamVersion;
	int latentParamVersion;
	enum expectation_type type;
	enum score_option scores;
	bool verbose;
} BA81Expect;

extern const struct rpf *rpf_model;
extern int rpf_numModels;

void computeRPF(BA81Expect *state, omxMatrix *itemParam, const int *quad,
		const bool wantlog, double *outcomeProb);
void cai2010(omxExpectation* oo, const int thrId, int recompute, const long qx);
double *ba81LikelihoodFast1(omxExpectation *oo, const int thrId, int specific, const long qx);

OMXINLINE static int
triangleLoc1(int diag)
{
	//if (diag < 1) error("Out of domain");
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

OMXINLINE static int
triangleLoc0(int diag)
{
	//if (diag < 0) error("Out of domain");
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

OMXINLINE static double *
eBase(BA81Expect *state, int thr)
{
	return state->allElxk + thr * state->numUnique;
}

OMXINLINE static double *
esBase(BA81Expect *state, int thr)
{
	return state->Eslxk + thr * state->numUnique * state->numSpecific;
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

OMXINLINE static double
areaProduct(BA81Expect *state, long qx, int sx, const int sg)
{
	if (state->numSpecific == 0) {
		return state->priQarea[qx];
	} else {
		if (sx == -1) {
			sx = qx % state->quadGridSize;
			qx /= state->quadGridSize;
		}
		return (state->priQarea[qx] * state->speQarea[sg * state->quadGridSize + sx]);
	}
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
validPatternLik(double pl)
{
	return isfinite(pl) && pl > 1e-300;
}

// debug tools
void pda(const double *ar, int rows, int cols);
void pia(const int *ar, int rows, int cols);

#endif
