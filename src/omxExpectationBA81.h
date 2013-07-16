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
	double *logNumIdentical;  // length numUnique
	int *rowMap;              // length numUnique, index of first instance of pattern

	// item description related
	std::vector<const double*> itemSpec;
	int maxOutcomes;
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
	std::vector<double> priLogQarea;      // totalPrimaryPoints
	std::vector<double> speLogQarea;      // quadGridSize * numSpecific

	// estimation related
	omxMatrix *customPrior;
	omxMatrix *itemParam;
	omxMatrix *EitemParam;    // E step version
	int cacheLXK;
	bool LXKcached;
	double *lxk;              // wo/cache, numUnique * thread
	std::vector<double> allElxk;          // numUnique * thread
	std::vector<double> Eslxk;            // numUnique * #specific dimensions * thread
	double *patternLik;       // numUnique
	double *_logPatternLik;   // numUnique
	int totalOutcomes;
	double *expected;         // totalOutcomes * totalQuadPoints
	double *ElatentMean;      // maxAbilities * numUnique
	double *ElatentCov;       // maxAbilities * maxAbilities * numUnique ; only lower triangle is used
	omxMatrix *latentMeanOut;
	omxMatrix *latentCovOut;

	enum expectation_type type;
	enum score_option scores;
} BA81Expect;

extern const struct rpf *rpf_model;
extern int rpf_numModels;

void ba81buildLXKcache(omxExpectation *oo);
double *computeRPF(BA81Expect *state, omxMatrix *itemParam, const int *quad);
void ba81SetupQuadrature(omxExpectation* oo, int gridsize, int flat);
void ba81Estep1(omxExpectation *oo);
void cai2010(omxExpectation* oo, int recompute, const int *primaryQuad);
double *ba81LikelihoodFast(omxExpectation *oo, int specific, const int *quad);
double *getLogPatternLik(omxExpectation* oo);

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

// state->allElxk[eIndex(state, px)]
OMXINLINE static int
eIndex(BA81Expect *state, int px)
{
	return omx_absolute_thread_num() * state->numUnique + px;
}

// state->Eslxk[esIndex(state, sx, px)]
OMXINLINE static int
esIndex(BA81Expect *state, int sx, int px)
{
	return (omx_absolute_thread_num() * state->numUnique * state->numSpecific +
		state->numUnique * sx + px);
}

OMXINLINE static void
pointToWhere(BA81Expect *state, const int *quad, double *where, int upto)
{
	for (int dx=0; dx < upto; dx++) {
		where[dx] = state->Qpoint[quad[dx]];
	}
}

OMXINLINE static long
encodeLocation(const int dims, const long grid, const int *quad)
{
	long qx = 0;
	for (int dx=dims-1; dx >= 0; dx--) {
		qx = qx * grid;
		qx += quad[dx];
	}
	return qx;
}

OMXINLINE static void
decodeLocation(long qx, const int dims, const long grid, int *quad)
{
	for (int dx=0; dx < dims; dx++) {
		quad[dx] = qx % grid;
		qx = qx / grid;
	}
}

OMXINLINE static double
logAreaProduct(BA81Expect *state, const int *quad, const int sg)
{
	int maxDims = state->maxDims;
	if (state->numSpecific == 0) {
		long qloc = encodeLocation(maxDims, state->quadGridSize, quad);
		return state->priLogQarea[qloc];
	} else {
		long priloc = encodeLocation(maxDims-1, state->quadGridSize, quad);
		return (state->priLogQarea[priloc] +
			state->speLogQarea[sg * state->quadGridSize + quad[maxDims - 1]]);
	}
}

// debug tools
void pda(const double *ar, int rows, int cols);
void pia(const int *ar, int rows, int cols);

#endif
