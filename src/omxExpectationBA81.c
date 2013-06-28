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

// Consider replacing log() with log2() in some places? Not worth it?

#include <Rmath.h>
#include "omxExpectation.h"
#include "omxOpenmpWrap.h"
#include "npsolWrap.h"
#include "libirt-rpf.h"
#include "omxOptimizer.h"  // remove TODO
#include "dmvnorm.h"

static const char *NAME = "ExpectationBA81";

static const struct rpf *rpf_model = NULL;
static int rpf_numModels;
static const double MIN_PATTERNLIK = 1e-100;

typedef struct {

	// data characteristics
	omxData *data;
	int numUnique;
	int *numIdentical;        // length numUnique
	double *logNumIdentical;  // length numUnique
	int *rowMap;              // length numUnique

	// item description related
	omxMatrix *itemSpec;
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
	double *Qpoint;           // quadGridSize
	double *priLogQarea;      // totalPrimaryPoints
	double *speLogQarea;      // quadGridSize * numSpecific

	// estimation related
	omxMatrix *EitemParam;    // E step version
	int userEitemParam;
	int rescale;
	omxMatrix *itemParam;     // M step version
	omxMatrix *customPrior;   // TODO remove?
	int derivPadSize;         // maxParam + maxParam*(1+maxParam)/2
	double *thrDeriv;         // itemParam->cols * derivPadSize * thread
	int *paramMap;            // itemParam->cols * derivPadSize -> index of free parameter
	int cacheLXK;
	double *lxk;              // wo/cache, numUnique * thread
	double *allSlxk;          // numUnique * thread
	double *Slxk;             // numUnique * #specific dimensions * thread
	double *patternLik;       // numUnique
	int totalOutcomes;
	double *expected;         // totalOutcomes * totalQuadPoints
	double *latentMean;       // maxAbilities * numUnique
	double *latentCov;        // maxAbilities * maxAbilities * numUnique ; only lower triangle is used
	double *latentMean1;      // maxDims
	double *latentCov1;       // maxDims * maxDims ; only lower triangle is used
	double *latentMeanOut;    // maxAbilities
	double *latentCovOut;     // maxAbilities * maxAbilities ; only lower triangle is used
	// also need corresponding label matrices for equating
	int *freeMean;
	int *freeCov;
	int choleskyError;

	int gradientCount;
	int fitCount;
	double lastEMLL;
} omxBA81State;

static void
pda(const double *ar, int rows, int cols) {   // column major order
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g, ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}
}

static void
pia(const int *ar, int rows, int cols) {   // column major order
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%d, ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}
}

static int
getNumThreads(omxExpectation* oo)
{
	int numThreads = oo->currentState->numThreads;
	if (numThreads < 1) numThreads = 1;
	return numThreads;
}

static void buildParamMap(omxExpectation* oo)
{
	omxState* currentState = oo->currentState;
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	omxMatrix *itemParam = state->itemParam;
	int size = itemParam->cols * state->derivPadSize;

	state->paramMap = Realloc(NULL, size, int);
	for (int px=0; px < size; px++) {
		state->paramMap[px] = -1;
	}

	int numFreeParams = currentState->numFreeParams;
	int *pRow = Realloc(NULL, numFreeParams, int);
	int *pCol = Realloc(NULL, numFreeParams, int);

	for (int px=0; px < numFreeParams; px++) {
		omxFreeVar *fv = currentState->freeVarList + px;
		for (int lx=0; lx < fv->numLocations; lx++) {
			if (~fv->matrices[lx] == itemParam->matrixNumber) {
				pRow[px] = fv->row[lx];
				pCol[px] = fv->col[lx];
				int at = pCol[px] * state->derivPadSize + pRow[px];
				state->paramMap[at] = px;
			}
		}
	}

	for (int p1=0; p1 < numFreeParams; p1++) {
		for (int p2=p1; p2 < numFreeParams; p2++) {
			if (pCol[p1] != pCol[p2]) continue;
			const double *spec = omxMatrixColumn(state->itemSpec, pCol[p1]);
			int id = spec[RPF_ISpecID];
			int numParam = (*rpf_model[id].numParam)(spec);
			int r1 = pRow[p1];
			int r2 = pRow[p2];
			if (r1 > r2) { int tmp=r1; r1=r2; r2=tmp; }
			int rowOffset = 0;
			for (int rx=1; rx <= r2; rx++) rowOffset += rx;
			int at = pCol[p1] * state->derivPadSize + numParam + rowOffset + r1;
			state->paramMap[at] = numFreeParams + p1 * numFreeParams + p2;
		}
	}

	Free(pRow);
	Free(pCol);

	state->thrDeriv = Realloc(NULL, itemParam->cols * state->derivPadSize * getNumThreads(oo), double);
}

OMXINLINE static void
pointToWhere(omxBA81State *state, const int *quad, double *where, int upto)
{
	for (int dx=0; dx < upto; dx++) {
		where[dx] = state->Qpoint[quad[dx]];
	}
}

OMXINLINE static void
assignDims(omxMatrix *itemSpec, omxMatrix *design, int dims, int maxDims, int ix,
	   const double *restrict theta, double *restrict ptheta)
{
	for (int dx=0; dx < dims; dx++) {
		int ability = (int)omxMatrixElement(design, dx, ix) - 1;
		if (ability >= maxDims) ability = maxDims-1;
		ptheta[dx] = theta[ability];
	}
}

/**
 * \param theta Vector of ability parameters, one per ability
 * \returns A numItems by maxOutcomes colMajor vector of doubles. Caller must Free it.
 */
static double *
computeRPF(omxExpectation *oo, omxMatrix *itemParam, const int *quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	int numItems = itemSpec->cols;
	omxMatrix *design = state->design;
	int maxDims = state->maxDims;

	double theta[maxDims];
	pointToWhere(state, quad, theta, maxDims);

	double *outcomeProb = Realloc(NULL, numItems * state->maxOutcomes, double);
	//double *outcomeProb = Calloc(numItems * state->maxOutcomes, double);

	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * state->maxOutcomes;
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double ptheta[dims];
		assignDims(itemSpec, design, dims, maxDims, ix, theta, ptheta);
		(*rpf_model[id].logprob)(spec, iparam, ptheta, out);
#if 0
		for (int ox=0; ox < spec[RPF_ISpecOutcomes]; ox++) {
			if (!isfinite(out[ox]) || out[ox] > 0) {
				Rprintf("spec\n");
				pda(spec, itemSpec->rows, 1);
				Rprintf("item param\n");
				pda(iparam, itemParam->rows, 1);
				Rprintf("where\n");
				pda(ptheta, dims, 1);
				error("RPF returned %20.20f", out[ox]);
			}
		}
#endif
	}

	return outcomeProb;
}

OMXINLINE static long
encodeLocation(const int dims, const long grid, const int *restrict quad)
{
	long qx = 0;
	for (int dx=dims-1; dx >= 0; dx--) {
		qx = qx * grid;
		qx += quad[dx];
	}
	return qx;
}

OMXINLINE static double *
getLXKcache(omxBA81State *state, const int *quad, const int specific)
{
	long ordinate;
	if (state->numSpecific == 0) {
		ordinate = encodeLocation(state->maxDims, state->quadGridSize, quad);
	} else {
		ordinate = (specific * state->totalQuadPoints +
			    encodeLocation(state->maxDims, state->quadGridSize, quad));
	}
	return state->lxk + state->numUnique * ordinate;
}

OMXINLINE static double *
ba81Likelihood(omxExpectation *oo, int specific, const int *restrict quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int numUnique = state->numUnique;
	int maxOutcomes = state->maxOutcomes;
	omxData *data = state->data;
	int numItems = state->itemSpec->cols;
	int *restrict Sgroup = state->Sgroup;
	double *restrict lxk;

	if (!state->cacheLXK) {
		lxk = state->lxk + numUnique * omx_absolute_thread_num();
	} else {
		lxk = getLXKcache(state, quad, specific);
	}

	const double *outcomeProb = computeRPF(oo, state->EitemParam, quad);
	if (!outcomeProb) {
		OMXZERO(lxk, numUnique);
		return lxk;
	}

	const int *rowMap = state->rowMap;
	for (int px=0; px < numUnique; px++) {
		double lxk1 = 0;
		for (int ix=0; ix < numItems; ix++) {
			if (specific != Sgroup[ix]) continue;
			int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
			if (pick == NA_INTEGER) continue;
			double piece = outcomeProb[ix * maxOutcomes + pick-1];  // move -1 elsewhere TODO
			lxk1 += piece;
		}
#if 0
#pragma omp critical(ba81LikelihoodDebug1)
		if (!isfinite(lxk1) || lxk1 > numItems) {
			Rprintf("where\n");
			double where[state->maxDims];
			pointToWhere(state, quad, where, state->maxDims);
			pda(where, state->maxDims, 1);
			Rprintf("prob\n");
			pda(outcomeProb, numItems, maxOutcomes);
			error("Likelihood of row %d is %f", rowMap[px], lxk1);
		}
#endif
		lxk[px] = lxk1;
	}

	Free(outcomeProb);

	return lxk;
}

OMXINLINE static double *
ba81LikelihoodFast(omxExpectation *oo, int specific, const int *restrict quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	if (!state->cacheLXK) {
		return ba81Likelihood(oo, specific, quad);
	} else {
		return getLXKcache(state, quad, specific);
	}

}

OMXINLINE static void
mapLatentSpace(omxBA81State *state, int px, int sgroup, double piece, const double *where)
{
	double *latentMean = state->latentMean;
	double *latentCov = state->latentCov;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int pmax = maxDims;
	if (state->numSpecific) pmax -= 1;

	if (sgroup == 0) {
		for (int d1=0; d1 < pmax; d1++) {
			double piece_w1 = piece * where[d1];
			int mloc = px * maxAbilities + d1;
#pragma omp atomic
			latentMean[mloc] += piece_w1;
			for (int d2=0; d2 <= d1; d2++) {
				int loc = px * maxAbilities * maxAbilities + d2 * maxAbilities + d1;
				double piece_cov = piece_w1 * where[d2];
#pragma omp atomic
				latentCov[loc] += piece_cov;
			}
		}
	}

	if (state->numSpecific) {
		int sdim = maxDims + sgroup - 1;

		double piece_w1 = piece * where[maxDims-1];
		int mloc = px * maxAbilities + sdim;
#pragma omp atomic
		latentMean[mloc] += piece_w1;

		int loc = px * maxAbilities * maxAbilities + sdim * maxAbilities + sdim;
		double piece_var = piece_w1 * where[maxDims-1];
#pragma omp atomic
		latentCov[loc] += piece_var;
	}
}

OMXINLINE static double
logAreaProduct(omxBA81State *state, const int *restrict quad, const int sg)
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

#define CALC_ALLSLXK(state, numUnique) \
	(state->allSlxk + omx_absolute_thread_num() * (numUnique))

#define CALC_SLXK(state, numUnique, numSpecific) \
	(state->Slxk + omx_absolute_thread_num() * (numUnique) * (numSpecific))

OMXINLINE static void
cai2010(omxExpectation* oo, int recompute, const int *restrict primaryQuad,
	double *restrict allSlxk, double *restrict Slxk)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int sDim = maxDims-1;

	int quad[maxDims];
	memcpy(quad, primaryQuad, sizeof(int)*sDim);

	OMXZERO(Slxk, numUnique * numSpecific);
	OMXZERO(allSlxk, numUnique);

	for (int sx=0; sx < numSpecific; sx++) {
		double *eis = Slxk + numUnique * sx;
		int quadGridSize = state->quadGridSize;

		for (int qx=0; qx < quadGridSize; qx++) {
			quad[sDim] = qx;
			double where[maxDims];
			pointToWhere(state, quad, where, maxDims);

			double *lxk;
			if (recompute) {
				lxk = ba81Likelihood(oo, sx, quad);
			} else {
				lxk = getLXKcache(state, quad, sx);
			}

			for (int ix=0; ix < numUnique; ix++) {
				eis[ix] += exp(lxk[ix] + state->priLogQarea[qx]);
			}
		}

		for (int px=0; px < numUnique; px++) {
			eis[px] = log(eis[px]);
			allSlxk[px] += eis[px];
		}
	}
}

OMXINLINE static void
decodeLocation(long qx, const int dims, const long grid,
	       int *restrict quad)
{
	for (int dx=0; dx < dims; dx++) {
		quad[dx] = qx % grid;
		qx = qx / grid;
	}
}

static void
ba81Estep1(omxExpectation *oo) {
	if(OMX_DEBUG_MML) {Rprintf("Beginning %s Computation.\n", NAME);}

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	double *patternLik = state->patternLik;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	double *latentMean = state->latentMean;
	double *latentCov = state->latentCov;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int primaryDims = maxDims;

	OMXZERO(patternLik, numUnique);
	OMXZERO(latentMean, numUnique * maxAbilities);
	OMXZERO(latentCov, numUnique * maxAbilities * maxAbilities);

	// E-step, marginalize person ability
	//
	// Note: In the notation of Bock & Aitkin (1981) and
	// Cai~(2010), these loops are reversed.  That is, the inner
	// loop is over quadrature points and the outer loop is over
	// all response patterns.
	//
	if (numSpecific == 0) {
#pragma omp parallel for num_threads(oo->currentState->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(state, quad, where, maxDims);

			double *lxk = ba81Likelihood(oo, 0, quad);

			double logArea = state->priLogQarea[qx];
#pragma omp critical(EstepUpdate)
			for (int px=0; px < numUnique; px++) {
				double tmp = exp(lxk[px] + logArea);
#if 0
				if (!isfinite(tmp)) {
					Rprintf("where\n");
					pda(where, maxDims, 1);
					error("Row %d lxk %f logArea %f tmp %f",
					      state->rowMap[px], lxk[px], logArea, tmp);
				}
#endif
				patternLik[px] += tmp;
				mapLatentSpace(state, px, 0, tmp, where);
			}
		}
	} else {
		primaryDims -= 1;
		int sDim = primaryDims;
		long specificPoints = state->quadGridSize;

#pragma omp parallel for num_threads(oo->currentState->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, primaryDims, state->quadGridSize, quad);

			double *allSlxk = CALC_ALLSLXK(state, numUnique);
			double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
			cai2010(oo, TRUE, quad, allSlxk, Slxk);

			for (int sgroup=0; sgroup < numSpecific; sgroup++) {
				double *eis = Slxk + numUnique * sgroup;
				for (long sx=0; sx < specificPoints; sx++) {
					quad[sDim] = sx;
					double where[maxDims];
					pointToWhere(state, quad, where, maxDims);
					double logArea = logAreaProduct(state, quad, sgroup);
					double *lxk = ba81LikelihoodFast(oo, sgroup, quad);
					for (int px=0; px < numUnique; px++) {
						double tmp = exp((allSlxk[px] - eis[px]) + lxk[px] + logArea);
						mapLatentSpace(state, px, sgroup, tmp, where);
					}
				}
			}

			double priLogArea = state->priLogQarea[qx];
#pragma omp critical(EstepUpdate)
			for (int px=0; px < numUnique; px++) {
				double tmp = exp(allSlxk[px] + priLogArea);
				patternLik[px] += tmp;  // is it faster to make this line atomic? TODO
			}
		}
	}

	int *numIdentical = state->numIdentical;

	if(0) {
		Rprintf("weight\n");
		for (int px=0; px < numUnique; px++) {
			double weight = numIdentical[px] / patternLik[px];
			Rprintf("%20.20f\n", weight);
		}

		Rprintf("per item mean\n");
		for (int px=0; px < numUnique; px++) {
			Rprintf("[%d] %20.20f\n", px, latentMean[px * maxAbilities]);
		}
	}

	for (int px=0; px < numUnique; px++) {
		if (patternLik[px] < MIN_PATTERNLIK) {
			patternLik[px] = MIN_PATTERNLIK;
			warning("Likelihood of pattern %d is 0, forcing to %.3g",
				px, MIN_PATTERNLIK);
		}

		double weight = numIdentical[px] / patternLik[px];
		for (int d1=0; d1 < primaryDims; d1++) {
			latentMean[px * maxAbilities + d1] *= weight;
			for (int d2=0; d2 <= d1; d2++) {
				int loc = px * maxAbilities * maxAbilities + d2 * maxAbilities + d1;
				latentCov[loc] *= weight;
			}
		}
		for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
			latentMean[px * maxAbilities + sdim] *= weight;
			int loc = px * maxAbilities * maxAbilities + sdim * maxAbilities + sdim;
			latentCov[loc] *= weight;
		}
#if 0
		if (!isfinite(patternLik[px])) {
			error("Likelihood of row %d is %f", state->rowMap[px], patternLik[px]);
		}
#endif
		patternLik[px] = log(patternLik[px]);
	}

	for (int px=1; px < numUnique; px++) {
		for (int d1=0; d1 < primaryDims; d1++) {
			latentMean[d1] += latentMean[px * maxAbilities + d1];
			for (int d2=0; d2 <= d1; d2++) {
				int cell = d2 * maxAbilities + d1;
				int loc = px * maxAbilities * maxAbilities + cell;
				latentCov[cell] += latentCov[loc];
			}
		}
		for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
			latentMean[sdim] += latentMean[px * maxAbilities + sdim];
			int cell = sdim * maxAbilities + sdim;
			int loc = px * maxAbilities * maxAbilities + cell;
			latentCov[cell] += latentCov[loc];
		}
	}

	//pda(latentMean, state->maxAbilities, 1);
	//pda(latentCov, state->maxAbilities, state->maxAbilities);

	omxData *data = state->data;
	for (int d1=0; d1 < maxAbilities; d1++) {
		latentMean[d1] /= data->rows;
	}

	for (int d1=0; d1 < primaryDims; d1++) {
		for (int d2=0; d2 <= d1; d2++) {
			int cell = d2 * maxAbilities + d1;
			int tcell = d1 * maxAbilities + d2;
			latentCov[tcell] = latentCov[cell] =
				latentCov[cell] / data->rows - latentMean[d1] * latentMean[d2];
		}
	}
	for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
		int cell = sdim * maxAbilities + sdim;
		latentCov[cell] = latentCov[cell] / data->rows - latentMean[sdim] * latentMean[sdim];
	}

	//Rprintf("E-step\n");
	//pda(latentMean, state->maxAbilities, 1);
	//pda(latentCov, state->maxAbilities, state->maxAbilities);
}

static void
schilling_bock_2005_rescale(omxExpectation *oo)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *design = state->design;
	double *latentMean = state->latentMean;
	double *latentCov = state->latentCov;
	double *latentMean1 = state->latentMean1;
	double *latentCov1 = state->latentCov1;
	int maxAbilities = state->maxAbilities;
	int maxDims = state->maxDims;

	//Rprintf("schilling bock\n");
	//pda(latentMean, maxAbilities, 1);
	//pda(latentCov, maxAbilities, maxAbilities);
	//omxPrint(design, "design");

	const char triangle = 'L';
	F77_CALL(dpotrf)(&triangle, &maxAbilities, latentCov, &maxAbilities, &state->choleskyError);
	if (state->choleskyError != 0) {
		warning("Cholesky failed with %d; rescaling disabled", state->choleskyError); // make error TODO?
		return;
	}

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
				latentMean1[d1] = 0;
			} else {
				latentMean1[d1] = latentMean[idesign[d1]];
			}
			for (int d2=0; d2 <= d1; d2++) {
				int cell = idesign[d2] * maxAbilities + idesign[d1];
				if (idesign[d1] == -1 || idesign[d2] == -1) {
					latentCov1[d2 * maxDims + d1] = d1==d2? 1 : 0;
				} else {
					latentCov1[d2 * maxDims + d1] = latentCov[cell];
				}
			}
		}
		if (1) {  // ease debugging, make optional TODO
			for (int d1=idx; d1 < maxDims; d1++) latentMean1[d1] = nan("");
			for (int d1=0; d1 < maxDims; d1++) {
				for (int d2=0; d2 < maxDims; d2++) {
					if (d1 < idx && d2 < idx) continue;
					latentCov1[d2 * maxDims + d1] = nan("");
				}
			}
		}
		double *iparam = omxMatrixColumn(itemParam, ix);
		int *mask = state->paramMap + state->derivPadSize * ix;
		rpf_model[id].rescale(spec, iparam, mask, latentMean1, latentCov1);
	}

	// This is an egregious hack :-)
	omxState* currentState = oo->currentState;
	int numFreeParams = currentState->numFreeParams;
	double param[numFreeParams];
	for (int rx=0; rx < itemParam->rows; rx++) {
		for (int cx=0; cx < itemParam->cols; cx++) {
			int vx = state->paramMap[cx * state->derivPadSize + rx];
			if (vx >= 0 && vx < numFreeParams) {
				// Look for gradient entries to determine what to copy.
				// I said it was a hack, no?
				param[vx] = omxMatrixElement(itemParam, rx, cx);
			}
		}
	}
	handleFreeVarList(currentState, param, numFreeParams);
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
static void
ba81SetupQuadrature(omxExpectation* oo, int gridsize, int flat)
{
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	int numUnique = state->numUnique;
	int numThreads = getNumThreads(oo);
	int maxDims = state->maxDims;
	int Qwidth = state->Qwidth;
	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);

	// try starting small and increasing to the cap TODO
	state->quadGridSize = gridsize;

	state->totalQuadPoints = 1;
	for (int dx=0; dx < maxDims; dx++) {
		state->totalQuadPoints *= state->quadGridSize;
	}

	state->totalPrimaryPoints = state->totalQuadPoints;

	if (numSpecific) {
		state->totalPrimaryPoints /= state->quadGridSize;
		state->speLogQarea = Realloc(state->speLogQarea, state->quadGridSize * gridsize, double);
	}

	state->Qpoint = Realloc(state->Qpoint, state->quadGridSize, double);
	state->priLogQarea = Realloc(state->priLogQarea, state->totalPrimaryPoints, double);

	for (int px=0; px < state->quadGridSize; px ++) {
		state->Qpoint[px] = Qwidth - px * 2 * Qwidth / (state->quadGridSize-1);
	}

	if (flat) {
		double flatd = log(1) - log(state->totalPrimaryPoints);
		for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
			state->priLogQarea[qx] = flatd;
		}
		flatd = log(1) - log(state->quadGridSize);
		for (int sx=0; sx < numSpecific; sx++) {
			for (int qx=0; qx < state->quadGridSize; qx++) {
				state->speLogQarea[ sx * state->quadGridSize + qx] = flatd;
			}
		}
	} else {
		double totalArea = 0;
		for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[priDims];
			decodeLocation(qx, priDims, state->quadGridSize, quad);
			double where[priDims];
			pointToWhere(state, quad, where, priDims);
			state->priLogQarea[qx] = dmvnorm(priDims, where, state->latentMeanOut, state->latentCovOut);
			totalArea += exp(state->priLogQarea[qx]);
		}
		totalArea = log(totalArea);
		for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
			state->priLogQarea[qx] -= totalArea;
			//Rprintf("%.5g,", state->priLogQarea[qx]);
		}
		//Rprintf("\n");

		for (int sx=0; sx < numSpecific; sx++) {
			totalArea = 0;
			for (int qx=0; qx < state->quadGridSize; qx++) {
				int covCell = (priDims + sx) * state->maxAbilities + priDims + sx;
				double den = dnorm(state->Qpoint[qx],
						   state->latentMeanOut[priDims + sx],
						   state->latentCovOut[covCell], TRUE);
				state->speLogQarea[sx * state->quadGridSize + qx] = den;
				totalArea += exp(den);
			}
			totalArea = log(totalArea);
			for (int qx=0; qx < state->quadGridSize; qx++) {
				state->speLogQarea[sx * state->quadGridSize + qx] -= totalArea;
			}
		}
	}

	if (!state->cacheLXK) {
		state->lxk = Realloc(state->lxk, numUnique * numThreads, double);
	} else {
		int ns = state->numSpecific;
		if (ns == 0) ns = 1;
		long numOrdinate = ns * state->totalQuadPoints;
		state->lxk = Realloc(state->lxk, numUnique * numOrdinate, double);
	}

	state->expected = Realloc(state->expected, state->totalOutcomes * state->totalQuadPoints, double);
}

OMXINLINE static void
updateLatentParam(omxExpectation* oo)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int maxAbilities = state->maxAbilities;

	for (int d1=0; d1 < maxAbilities; d1++) {
		if (state->freeMean[d1]) {
			state->latentMeanOut[d1] = state->latentMean[d1];
		}

		for (int d2=0; d2 <= d1; d2++) {
			int cell = d2 * maxAbilities + d1;
			if (state->freeCov[cell]) {
				state->latentCovOut[cell] = state->latentCov[cell];
			}
		}
	}
	//Rprintf("updateLatentParam\n");
	//pda(state->latentMeanOut, maxAbilities, 1);
	//pda(state->latentCovOut, maxAbilities, maxAbilities);

	ba81SetupQuadrature(oo, state->targetQpoints, 0);
}

OMXINLINE static void
expectedUpdate(omxData *restrict data, const int *rowMap, const int px, const int item,
	       const double observed, const int outcomes, double *out)
{
	int pick = omxIntDataElementUnsafe(data, rowMap[px], item);
	if (pick == NA_INTEGER) {
		double slice = exp(observed - log(outcomes));
		for (int ox=0; ox < outcomes; ox++) {
			out[ox] += slice;
		}
	} else {
		out[pick-1] += exp(observed);
	}
}

OMXINLINE static void
ba81Expected(omxExpectation* oo)
{
	omxState* currentState = oo->currentState;  // only used in #pragma omp
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxData *data = state->data;
	int numSpecific = state->numSpecific;
	const int *rowMap = state->rowMap;
	double *patternLik = state->patternLik;
	double *logNumIdentical = state->logNumIdentical;
	int numUnique = state->numUnique;
	int maxDims = state->maxDims;
	int numItems = state->itemParam->cols;
	omxMatrix *itemSpec = state->itemSpec;
	int totalOutcomes = state->totalOutcomes;

	OMXZERO(state->expected, totalOutcomes * state->totalQuadPoints);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(currentState->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double *lxk = ba81LikelihoodFast(oo, 0, quad);
			for (int px=0; px < numUnique; px++) {
				double *out = state->expected + qx * totalOutcomes;
				double observed = logNumIdentical[px] + lxk[px] - patternLik[px];
				for (int ix=0; ix < numItems; ix++) {
					const double *spec = omxMatrixColumn(itemSpec, ix);
					int outcomes = spec[RPF_ISpecOutcomes];
					expectedUpdate(data, rowMap, px, ix, observed, outcomes, out);
					out += outcomes;
				}
			}
		}
	} else {
		int sDim = state->maxDims-1;
		long specificPoints = state->quadGridSize;

#pragma omp parallel for num_threads(currentState->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);

			// allSlxk, Slxk only depend on the ordinate of the primary dimensions
			double *allSlxk = CALC_ALLSLXK(state, numUnique);
			double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
			cai2010(oo, !state->cacheLXK, quad, allSlxk, Slxk);

			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				long qloc = encodeLocation(state->maxDims, state->quadGridSize, quad);

				for (int sgroup=0; sgroup < numSpecific; sgroup++) {
					double *eis = Slxk + numUnique * sgroup;
					double *lxk = ba81LikelihoodFast(oo, sgroup, quad);

					for (int px=0; px < numUnique; px++) {
						double *out = state->expected + totalOutcomes * qloc;

						for (int ix=0; ix < numItems; ix++) {
							const double *spec = omxMatrixColumn(itemSpec, ix);
							int outcomes = spec[RPF_ISpecOutcomes];
							if (state->Sgroup[ix] == sgroup) {
								double observed = logNumIdentical[px] + (allSlxk[px] - eis[px]) +
									(lxk[px] - patternLik[px]);
								expectedUpdate(data, rowMap, px, ix, observed, outcomes, out);
							}
							out += outcomes;
						}
					}
				}
			}
		}
	}
	//pda(state->expected, state->totalOutcomes, state->totalQuadPoints);
}

static void
ba81Estep(omxExpectation *oo, enum ComputeExpectationContext ctx) {
	if (ctx == COMPUTE_EXPECT_INITIALIZE) return;

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	if (state->userEitemParam) {
		state->userEitemParam = FALSE;
	} else {
		omxCopyMatrix(state->EitemParam, state->itemParam);
	}
	ba81Estep1(oo);
	if (ctx == COMPUTE_EXPECT_PREFIT) {
		if (state->rescale) schilling_bock_2005_rescale(oo);
	} else {
		updateLatentParam(oo);
	}
	ba81Expected(oo);
}

OMXINLINE static double
ba81Fit1Ordinate(omxExpectation* oo, const int *quad, const double *weight, int want)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	omxMatrix *itemParam = state->itemParam;
	int numItems = itemParam->cols;
	int maxOutcomes = state->maxOutcomes;
	int maxDims = state->maxDims;
	double *myDeriv = state->thrDeriv + itemParam->cols * state->derivPadSize * omx_absolute_thread_num();
	int do_deriv = want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN);

	double where[maxDims];
	pointToWhere(state, quad, where, maxDims);

	double *outcomeProb = computeRPF(oo, itemParam, quad); // avoid malloc/free? TODO
	if (!outcomeProb) return 0;

	double thr_ll = 0;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		int id = spec[RPF_ISpecID];
		int iOutcomes = spec[RPF_ISpecOutcomes];

		double area = exp(logAreaProduct(state, quad, state->Sgroup[ix]));   // avoid exp() here? TODO
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
ba81ComputeMFit1(omxExpectation* oo, int want, double *gradient, double *hessian)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxState* currentState = oo->currentState;  // only used in #pragma omp
	omxMatrix *customPrior = state->customPrior;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *itemSpec = state->itemSpec;
	int maxDims = state->maxDims;
	const int totalOutcomes = state->totalOutcomes;

	double ll = 0;
	if (customPrior) {
		omxRecompute(customPrior);
		ll = customPrior->data[0];
	}

	if (!isfinite(ll)) {
		omxPrint(itemParam, "item param");
		error("Bayesian prior returned %g; do you need to add a lbound/ubound?", ll);
	}

#pragma omp parallel for num_threads(currentState->numThreads)
	for (long qx=0; qx < state->totalQuadPoints; qx++) {
		//double area = exp(state->priLogQarea[qx]);  // avoid exp() here? TODO
		int quad[maxDims];
		decodeLocation(qx, maxDims, state->quadGridSize, quad);
		double *weight = state->expected + qx * totalOutcomes;
		double thr_ll = ba81Fit1Ordinate(oo, quad, weight, want);
		
#pragma omp atomic
		ll += thr_ll;
	}

	if (!customPrior && gradient) {
		double *deriv0 = state->thrDeriv;

		int perThread = itemParam->cols * state->derivPadSize;
		for (int th=1; th < getNumThreads(oo); th++) {
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

		int numFreeParams = currentState->numFreeParams;
		int numParams = itemParam->cols * state->derivPadSize;
		for (int ox=0; ox < numParams; ox++) {
			int to = state->paramMap[ox];
			if (to == -1) continue;

			// Need to check because this can happen if
			// lbounds/ubounds are not set appropriately.
			if (0 && !isfinite(deriv0[ox])) {
				int item = ox / itemParam->rows;
				Rprintf("item parameters:\n");
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
				if (hessian) hessian[to - numFreeParams] -= deriv0[ox];
			}
		}
	}

	return -ll;
}

double
ba81ComputeFit(omxExpectation* oo, int want, double *gradient, double *hessian)
{
	if (!want) return 0;

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxState* currentState = oo->currentState;

	if (want & FF_COMPUTE_FIT) {
		++state->fitCount;
	}

	if (want & (FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN)) {
		// M-step
		++state->gradientCount;

		int numFreeParams = currentState->numFreeParams;
		OMXZERO(gradient, numFreeParams);
		if (hessian) OMXZERO(hessian, numFreeParams * numFreeParams);

		omxMatrix *itemParam = state->itemParam;
		OMXZERO(state->thrDeriv, state->derivPadSize * itemParam->cols * getNumThreads(oo));

		double got = ba81ComputeMFit1(oo, want, gradient, hessian);
		state->lastEMLL = got;
		return got;
	} else {
		// Major EM iteration, note completely different LL calculation

		double *patternLik = state->patternLik;
		int *numIdentical = state->numIdentical;
		int numUnique = state->numUnique;
		double got = 0;
		for (int ux=0; ux < numUnique; ux++) {
			got += numIdentical[ux] * patternLik[ux];
		}
		return -2 * got;
	}
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
omxRListElement *
ba81EAP(omxExpectation *oo, int *numReturns)
{
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	int maxDims = state->maxDims;
	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);
	int maxAbilities = state->maxAbilities;

	*numReturns = 4; // + (maxDims > 1) + (numSpecific > 1);
	omxRListElement *out = (omxRListElement*) R_alloc(*numReturns, sizeof(omxRListElement));
	int ox=0;

	out[ox].numValues = -1;
	strcpy(out[ox].label, "latent.cov");
	out[ox].rows = maxAbilities;
	out[ox].cols = maxAbilities;
	out[ox].values = state->latentCovOut;
	++ox;

	out[ox].numValues = maxAbilities;
	strcpy(out[ox].label, "latent.mean");
	out[ox].values = state->latentMeanOut;
	++ox;

	out[ox].numValues = 1;
	strcpy(out[ox].label, "EM.LL");
	out[ox].values = &state->lastEMLL;
	++ox;

	omxData *data = state->data;
	int numUnique = state->numUnique;

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 1);
	ba81Estep1(oo);   // recalc patternLik with a flat prior

	/*
	double *cov = NULL;
	if (maxDims > 1) {
		strcpy(out[2].label, "ability.cov");
		out[2].numValues = -1;
		out[2].rows = maxDims;
		out[2].cols = maxDims;
		out[2].values = (double*) R_alloc(out[2].rows * out[2].cols, sizeof(double));
		cov = out[2].values;
		OMXZERO(cov, out[2].rows * out[2].cols);
	}
	*/

	// Need a separate work space because the destination needs
	// to be in unsorted order with duplicated rows.
	double *ability = Calloc(numUnique * maxAbilities * 2, double);

	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, state->quadGridSize, quad);
		double where[priDims];
		pointToWhere(state, quad, where, priDims);
		double logArea = state->priLogQarea[qx];

		double *lxk;
		if (numSpecific == 0) {
			lxk = ba81LikelihoodFast(oo, 0, quad);
		} else {
			double *allSlxk = CALC_ALLSLXK(state, numUnique);
			double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
			cai2010(oo, FALSE, quad, allSlxk, Slxk);
			lxk = allSlxk;
		}

		double *row = ability;
		for (int px=0; px < numUnique; px++) {
			double plik = exp(logArea + lxk[px]);
			for (int dx=0; dx < priDims; dx++) {
				double piece = where[dx] * plik;
				row[dx*2] += piece;
				row[dx*2 + 1] += where[dx] * piece;
				// ignore cov, for now
			}
			row += 2 * maxAbilities;
		}
	}

	double *ris = Realloc(NULL, numUnique, double);
	for (int sx=0; sx < numSpecific; sx++) {
		for (int sqx=0; sqx < state->quadGridSize; sqx++) {
			double area = exp(state->speLogQarea[sx * state->quadGridSize + sqx]);
			double ptArea = area * state->Qpoint[sqx];
			OMXZERO(ris, numUnique);
			for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
				int quad[maxDims];
				decodeLocation(qx, priDims, state->quadGridSize, quad);
				quad[priDims] = sqx;

				double *allSlxk = CALC_ALLSLXK(state, numUnique);
				double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
				cai2010(oo, FALSE, quad, allSlxk, Slxk);

				double *eis = Slxk + numUnique * sx;
				double *lxk = ba81LikelihoodFast(oo, sx, quad);

				double logArea = state->priLogQarea[qx];
				for (int px=0; px < numUnique; px++) {
					ris[px] += exp(logArea + lxk[px] + allSlxk[px] - eis[px]);
				}
			}
			double *row = ability;
			for (int px=0; px < numUnique; px++) {
				double piece = ris[px] * ptArea;
			        row[(priDims + sx) * 2] += piece;
			        row[(priDims + sx) * 2 + 1] += piece * state->Qpoint[sqx];
				row += 2 * maxAbilities;
		        }
		}
	}
	Free(ris);

	double *patternLik = state->patternLik;
	double *row = ability;
	for (int px=0; px < numUnique; px++) {
		double denom = exp(patternLik[px]);
		for (int ax=0; ax < maxAbilities; ax++) {
			row[ax * 2] /= denom;
			row[ax * 2 + 1] /= denom;
			row[ax * 2 + 1] -= row[ax * 2] * row[ax * 2];
		}
		row += 2 * maxAbilities;
        }

	/*
	// make symmetric
	for (int d1=0; d1 < maxDims; d1++) {
		for (int d2=0; d2 < d1; d2++) {
			cov[d2 * maxDims + d1] = cov[d1 * maxDims + d2];
		}
	}
	*/

	for (int px=0; px < numUnique; px++) {
		double *arow = ability + px * 2 * maxAbilities;
		for (int dx=0; dx < maxAbilities; dx++) {
			arow[dx*2+1] = sqrt(arow[dx*2+1]);
		}
	}

	strcpy(out[ox].label, "ability");
	out[ox].numValues = -1;
	out[ox].rows = 2 * maxAbilities;
	out[ox].cols = data->rows;
	out[ox].values = (double*) R_alloc(out[ox].rows * out[ox].cols, sizeof(double));

	for (int rx=0; rx < numUnique; rx++) {
		double *pa = ability + rx * 2 * maxAbilities;

		int dups = omxDataNumIdenticalRows(state->data, state->rowMap[rx]);
		for (int dup=0; dup < dups; dup++) {
			int dest = omxDataIndex(data, state->rowMap[rx]+dup);
			memcpy(out[ox].values + dest * out[ox].rows, pa, sizeof(double) * 2 * maxAbilities);
		}
	}
	Free(ability);
	++ox;

	for (int ix=0; ix < state->itemParam->cols; ix++) {
		double *spec = omxMatrixColumn(state->itemSpec, ix);
		int id = spec[RPF_ISpecID];
		double *param = omxMatrixColumn(state->itemParam, ix);
		rpf_model[id].postfit(spec, param);
	}

	return out;
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		Rprintf("Freeing %s function.\n", NAME);
	}
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	//Rprintf("fit %d gradient %d\n", state->fitCount, state->gradientCount);
	omxFreeAllMatrixData(state->itemSpec);
	omxFreeAllMatrixData(state->itemParam);
	omxFreeAllMatrixData(state->EitemParam);
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->customPrior);
	Free(state->logNumIdentical);
	Free(state->numIdentical);
	Free(state->Qpoint);
	Free(state->priLogQarea);
	Free(state->rowMap);
	Free(state->patternLik);
	Free(state->lxk);
	Free(state->Slxk);
	Free(state->allSlxk);
	Free(state->Sgroup);
	Free(state->expected);
	Free(state->paramMap);
	Free(state->thrDeriv);
	Free(state->latentMean);
	Free(state->latentCov);
	Free(state->latentMean1);
	Free(state->latentCov1);
	Free(state->latentMeanOut);
	Free(state->latentCovOut);
	Free(state);
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    PROTECT(matrixDims = getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    UNPROTECT(1);
}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		Rprintf("Initializing %s.\n", NAME);
	}
	if (!rpf_model) {
		const int wantVersion = 3;
		int version;
		get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model_GPL");
		(*get_librpf)(&version, &rpf_numModels, &rpf_model);
		if (version < wantVersion) error("librpf binary API %d installed, at least %d is required",
						 version, wantVersion);
	}
	
	omxBA81State *state = Calloc(1, omxBA81State);
	oo->argStruct = (void*) state;

	PROTECT(tmp = GET_SLOT(rObj, install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf(currentState, "%s unable to handle data type %s", NAME, omxDataType(state->data));
		return;
	}

	state->itemSpec =
		omxNewMatrixFromSlot(rObj, currentState, "ItemSpec");
	state->design =
		omxNewMatrixFromSlot(rObj, currentState, "Design");
	state->itemParam =
		omxNewMatrixFromSlot(rObj, currentState, "ItemParam");
	state->EitemParam =
		omxNewMatrixFromSlot(rObj, currentState, "EItemParam");
	if (state->EitemParam) {
		state->userEitemParam = TRUE;
	} else {
		state->EitemParam =
			omxInitTemporaryMatrix(NULL, state->itemParam->rows, state->itemParam->cols,
					       TRUE, currentState);
	}

	state->customPrior =
		omxNewMatrixFromSlot(rObj, currentState, "CustomPrior");
	
	oo->computeFun = ba81Estep;
	oo->destructFun = ba81Destroy;
	
	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we don't waste RAM.
	// The following summary stats would be cheaper to calculate too.

	int numUnique = 0;
	omxData *data = state->data;
	if (omxDataNumFactor(data) != data->cols) {
		// verify they are ordered factors TODO
		omxRaiseErrorf(currentState, "%s: all columns must be factors", NAME);
		return;
	}

	for (int rx=0; rx < data->rows;) {
		rx += omxDataNumIdenticalRows(state->data, rx);
		++numUnique;
	}
	state->numUnique = numUnique;

	state->rowMap = Realloc(NULL, numUnique, int);
	state->numIdentical = Realloc(NULL, numUnique, int);
	state->logNumIdentical = Realloc(NULL, numUnique, double);

	int numItems = state->itemParam->cols;
	if (data->cols != numItems) {
		error("Data has %d columns for %d items", data->cols, numItems);
	}

	for (int rx=0, ux=0; rx < data->rows; ux++) {
		if (rx == 0) {
			// all NA rows will sort to the top
			int na=0;
			for (int ix=0; ix < numItems; ix++) {
				if (omxIntDataElement(data, 0, ix) == NA_INTEGER) { ++na; }
			}
			if (na == numItems) {
				omxRaiseErrorf(currentState, "Remove rows with all NAs");
				return;
			}
		}
		int dups = omxDataNumIdenticalRows(state->data, rx);
		state->numIdentical[ux] = dups;
		state->logNumIdentical[ux] = log(dups);
		state->rowMap[ux] = rx;
		rx += dups;
	}

	state->patternLik = Realloc(NULL, numUnique, double);

	int numThreads = getNumThreads(oo);

	int maxSpec = 0;
	int maxParam = 0;
	state->maxDims = 0;
	state->maxOutcomes = 0;

	for (int ix=0; ix < numItems; ix++) {
		double *spec = omxMatrixColumn(state->itemSpec, ix);
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= rpf_numModels) {
			omxRaiseErrorf(currentState, "ItemSpec column %d has unknown item model %d", ix, id);
			return;
		}

		double *param = omxMatrixColumn(state->itemParam, ix);
		rpf_model[id].prefit(spec, param);
	}

	int totalOutcomes = 0;
	for (int cx = 0; cx < data->cols; cx++) {
		const double *spec = omxMatrixColumn(state->itemSpec, cx);
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		if (state->maxDims < dims)
			state->maxDims = dims;

		int no = spec[RPF_ISpecOutcomes];
		totalOutcomes += no;
		if (state->maxOutcomes < no)
			state->maxOutcomes = no;

		// TODO this summary stat should be available from omxData
		int dataMax=0;
		for (int rx=0; rx < data->rows; rx++) {
			int pick = omxIntDataElementUnsafe(data, rx, cx);
			if (dataMax < pick)
				dataMax = pick;
		}
		if (dataMax > no) {
			error("Data for item %d has %d outcomes, not %d", cx+1, dataMax, no);
		} else if (dataMax < no) {
			warning("Data for item %d has only %d outcomes, not %d", cx+1, dataMax, no);
			// promote to error?
			// should complain if an outcome is not represented in the data TODO
		}

		int numSpec = (*rpf_model[id].numSpec)(spec);
		if (maxSpec < numSpec)
			maxSpec = numSpec;

		int numParam = (*rpf_model[id].numParam)(spec);
		if (maxParam < numParam)
			maxParam = numParam;
	}

	state->totalOutcomes = totalOutcomes;

	if (state->itemSpec->cols != data->cols || state->itemSpec->rows != maxSpec) {
		omxRaiseErrorf(currentState, "ItemSpec must have %d item columns and %d rows",
			       data->cols, maxSpec);
		return;
	}
	if (state->itemParam->rows != maxParam) {
		omxRaiseErrorf(currentState, "ItemParam should have %d rows", maxParam);
		return;
	}

	state->derivPadSize = maxParam + maxParam*(1+maxParam)/2;

	if (state->design == NULL) {
		state->maxAbilities = state->maxDims;
		state->design = omxInitTemporaryMatrix(NULL, state->maxDims, numItems,
				       TRUE, currentState);
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = omxMatrixColumn(state->itemSpec, ix);
			int dims = spec[RPF_ISpecDims];
			for (int dx=0; dx < state->maxDims; dx++) {
				omxSetMatrixElement(state->design, dx, ix, dx < dims? (double)dx+1 : nan(""));
			}
		}
	} else {
		omxMatrix *design = state->design;
		if (design->cols != numItems ||
		    design->rows != state->maxDims) {
			omxRaiseErrorf(currentState, "Design matrix should have %d rows and %d columns",
				       state->maxDims, numItems);
			return;
		}

		state->maxAbilities = 0;
		for (int ix=0; ix < design->rows * design->cols; ix++) {
			double got = design->data[ix];
			if (!R_FINITE(got)) continue;
			if (round(got) != (int)got) error("Design matrix can only contain integers"); // TODO better way?
			if (state->maxAbilities < got)
				state->maxAbilities = got;
		}
		for (int ix=0; ix < design->cols; ix++) {
			const double *idesign = omxMatrixColumn(design, ix);
			int ddim = 0;
			for (int rx=0; rx < design->rows; rx++) {
				if (isfinite(idesign[rx])) ddim += 1;
			}
			const double *spec = omxMatrixColumn(state->itemSpec, ix);
			int dims = spec[RPF_ISpecDims];
			if (ddim > dims) error("Item %d has %d dims but design assigns %d", ix, dims, ddim);
		}
	}
	if (state->maxAbilities <= state->maxDims) {
		state->Sgroup = Calloc(numItems, int);
	} else {
		// Not sure if this is correct, revisit TODO
		int Sgroup0 = -1;
		state->Sgroup = Realloc(NULL, numItems, int);
		for (int dx=0; dx < state->maxDims; dx++) {
			for (int ix=0; ix < numItems; ix++) {
				int ability = omxMatrixElement(state->design, dx, ix);
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
						omxRaiseErrorf(currentState, "Item %d cannot belong to more than "
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
		state->allSlxk = Realloc(NULL, numUnique * numThreads, double);
		state->Slxk = Realloc(NULL, numUnique * state->numSpecific * numThreads, double);
	}

	PROTECT(tmp = GET_SLOT(rObj, install("free.mean")));
	state->freeMean = Realloc(NULL, state->maxAbilities, int);
	if (isNull(tmp)) {
		OMXZERO(state->freeMean, state->maxAbilities);
	} else {
		if (length(tmp) != state->maxAbilities) {
			error("Specify for %d means whether to freely estimate", state->maxAbilities);
		}
		memcpy(state->freeMean, LOGICAL(tmp), sizeof(int)*state->maxAbilities);
	}

	PROTECT(tmp = GET_SLOT(rObj, install("free.cov")));
	state->freeCov = Realloc(NULL, state->maxAbilities * state->maxAbilities, int);
	if (isNull(tmp)) {
		OMXZERO(state->freeCov, state->maxAbilities * state->maxAbilities);
	} else {
		int rows;
		int cols;
		getMatrixDims(tmp, &rows, &cols);
		if (rows != state->maxAbilities && cols != state->maxAbilities) {
			error("free.cov must be a %d by %d matrix of logicals",
			      state->maxAbilities, state->maxAbilities);
		}
		memcpy(state->freeCov, LOGICAL(tmp),
		       sizeof(int)*state->maxAbilities * state->maxAbilities);
	}

	PROTECT(tmp = GET_SLOT(rObj, install("cache")));
	state->cacheLXK = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("rescale")));
	state->rescale = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qpoints")));
	state->targetQpoints = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qwidth")));
	state->Qwidth = asReal(tmp);

	state->latentMean = Realloc(NULL, state->maxAbilities * numUnique, double);
	state->latentCov = Realloc(NULL, state->maxAbilities * state->maxAbilities * numUnique, double);
	state->latentMean1 = Realloc(NULL, state->maxDims, double);
	state->latentCov1 = Realloc(NULL, state->maxDims * state->maxDims, double);
	state->latentMeanOut = Calloc(state->maxAbilities, double);
	state->latentCovOut = Realloc(NULL, state->maxAbilities * state->maxAbilities, double);

	for (int d1=0; d1 < state->maxAbilities; d1++) {
		for (int d2=0; d2 < state->maxAbilities; d2++) {
			state->latentCovOut[d1 * state->maxAbilities + d2] = d1==d2? 1 : 0;
		}
	}

	ba81SetupQuadrature(oo, state->targetQpoints, 0);

	buildParamMap(oo);

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData for column summary stats?
}
