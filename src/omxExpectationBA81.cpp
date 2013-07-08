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
#include "omxExpectationBA81.h"
#include "omxOpenmpWrap.h"
#include "npsolWrap.h"
#include "libifa-rpf.h"
#include "dmvnorm.h"

static const char *NAME = "ExpectationBA81";

const struct rpf *rpf_model = NULL;
int rpf_numModels;
static const double MIN_PATTERNLIK = 1e-100;

void pda(const double *ar, int rows, int cols)
{
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
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%d, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

OMXINLINE static int
triangleLoc0(int diag)
{ return (diag+1) * (diag+2) / 2 - 1; }

OMXINLINE static int
triangleLoc1(int diag)
{ return (diag) * (diag+1) / 2; }

OMXINLINE static void
assignDims(omxMatrix *itemSpec, omxMatrix *design, int dims, int maxDims, int ix,
	   const double *theta, double *ptheta)
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
double *
computeRPF(omxMatrix *itemSpec, omxMatrix *design, omxMatrix *itemParam,
	   int maxDims, int maxOutcomes, const int *quad, const double *Qpoint)
{
	int numItems = itemSpec->cols;

	double theta[maxDims];
	pointToWhere(Qpoint, quad, theta, maxDims);

	double *outcomeProb = Realloc(NULL, numItems * maxOutcomes, double);
	//double *outcomeProb = Calloc(numItems * maxOutcomes, double);

	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * maxOutcomes;
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double ptheta[dims];
		assignDims(itemSpec, design, dims, maxDims, ix, theta, ptheta);
		(*rpf_model[id].logprob)(spec, iparam, ptheta, out);
#if 0
		for (int ox=0; ox < spec[RPF_ISpecOutcomes]; ox++) {
			if (!isfinite(out[ox]) || out[ox] > 0) {
				mxLog("spec");
				pda(spec, itemSpec->rows, 1);
				mxLog("item param");
				pda(iparam, itemParam->rows, 1);
				mxLog("where");
				pda(ptheta, dims, 1);
				error("RPF returned %20.20f", out[ox]);
			}
		}
#endif
	}

	return outcomeProb;
}

OMXINLINE static double *
getLXKcache(BA81Expect *state, const int *quad, const int specific)
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
ba81Likelihood(omxExpectation *oo, int specific, const int *quad)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	int numUnique = state->numUnique;
	int maxOutcomes = state->maxOutcomes;
	omxData *data = state->data;
	int numItems = state->itemSpec->cols;
	int *Sgroup = state->Sgroup;
	double *lxk;

	if (!state->cacheLXK) {
		lxk = state->lxk + numUnique * omx_absolute_thread_num();
	} else {
		lxk = getLXKcache(state, quad, specific);
	}

	const double *outcomeProb = computeRPF(state->itemSpec, state->design, state->EitemParam,
					       state->maxDims, state->maxOutcomes, quad, state->Qpoint);
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
			mxLog("where");
			double where[state->maxDims];
			pointToWhere(state->Qpoint, quad, where, state->maxDims);
			pda(where, state->maxDims, 1);
			mxLog("prob");
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
ba81LikelihoodFast(omxExpectation *oo, int specific, const int *quad)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (!state->cacheLXK) {
		return ba81Likelihood(oo, specific, quad);
	} else {
		return getLXKcache(state, quad, specific);
	}

}

OMXINLINE static void
mapLatentSpace(BA81Expect *state, int px, int sgroup, double piece, const double *where)
{
	double *ElatentMean = state->ElatentMean;
	double *ElatentCov = state->ElatentCov;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int pmax = maxDims;
	if (state->numSpecific) pmax -= 1;

	if (sgroup == 0) {
		for (int d1=0; d1 < pmax; d1++) {
			double piece_w1 = piece * where[d1];
			int mloc = px * maxAbilities + d1;
#pragma omp atomic
			ElatentMean[mloc] += piece_w1;
			for (int d2=0; d2 <= d1; d2++) {
				int loc = px * maxAbilities * maxAbilities + d2 * maxAbilities + d1;
				double piece_cov = piece_w1 * where[d2];
#pragma omp atomic
				ElatentCov[loc] += piece_cov;
			}
		}
	}

	if (state->numSpecific) {
		int sdim = maxDims + sgroup - 1;

		double piece_w1 = piece * where[maxDims-1];
		int mloc = px * maxAbilities + sdim;
#pragma omp atomic
		ElatentMean[mloc] += piece_w1;

		int loc = px * maxAbilities * maxAbilities + sdim * maxAbilities + sdim;
		double piece_var = piece_w1 * where[maxDims-1];
#pragma omp atomic
		ElatentCov[loc] += piece_var;
	}
}

#define CALC_ALLSLXK(state, numUnique) \
	(state->allSlxk + omx_absolute_thread_num() * (numUnique))

#define CALC_SLXK(state, numUnique, numSpecific) \
	(state->Slxk + omx_absolute_thread_num() * (numUnique) * (numSpecific))

OMXINLINE static void
cai2010(omxExpectation* oo, int recompute, const int *primaryQuad,
	double *allSlxk, double *Slxk)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int sDim = maxDims-1;
	int quadGridSize = state->quadGridSize;
	int quad[maxDims];
	memcpy(quad, primaryQuad, sizeof(int)*sDim);

	OMXZERO(Slxk, numUnique * numSpecific);
	OMXZERO(allSlxk, numUnique);

	for (int sx=0; sx < numSpecific; sx++) {
		double *eis = Slxk + numUnique * sx;

		for (int qx=0; qx < quadGridSize; qx++) {
			quad[sDim] = qx;

			double *lxk;     // a.k.a. "L_is"
			if (recompute) {
				lxk = ba81Likelihood(oo, sx, quad);
			} else {
				lxk = ba81LikelihoodFast(oo, sx, quad);
			}

			for (int ix=0; ix < numUnique; ix++) {
				eis[ix] += exp(lxk[ix] + state->speLogQarea[sx * quadGridSize + qx]);
			}
		}

		for (int px=0; px < numUnique; px++) {
			eis[px] = log(eis[px]);
			allSlxk[px] += eis[px];  // allSlxk a.k.a. "E_i"
		}
	}
}

static void
ba81Estep1(omxExpectation *oo) {
	if(OMX_DEBUG) {mxLog("Beginning %s Computation.", NAME);}

	BA81Expect *state = (BA81Expect*) oo->argStruct;
	double *patternLik = state->patternLik;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	double *ElatentMean = state->ElatentMean;
	double *ElatentCov = state->ElatentCov;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int primaryDims = maxDims;

	OMXZERO(patternLik, numUnique);
	OMXZERO(ElatentMean, numUnique * maxAbilities);
	OMXZERO(ElatentCov, numUnique * maxAbilities * maxAbilities);

	// E-step, marginalize person ability
	//
	// Note: In the notation of Bock & Aitkin (1981) and
	// Cai~(2010), these loops are reversed.  That is, the inner
	// loop is over quadrature points and the outer loop is over
	// all response patterns.
	//
	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(state->Qpoint, quad, where, maxDims);

			double *lxk = ba81Likelihood(oo, 0, quad);

			double logArea = state->priLogQarea[qx];
#pragma omp critical(EstepUpdate)
			for (int px=0; px < numUnique; px++) {
				double tmp = exp(lxk[px] + logArea);
#if 0
				if (!isfinite(tmp)) {
					mxLog("where");
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

#pragma omp parallel for num_threads(Global->numThreads)
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
					pointToWhere(state->Qpoint, quad, where, maxDims);
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

	if (0) {
		mxLog("weight");
		for (int px=0; px < numUnique; px++) {
			double weight = numIdentical[px] / patternLik[px];
			mxLog("%20.20f", weight);
		}

		mxLog("per item mean");
		pda(ElatentMean, maxAbilities, numUnique);
	}

	for (int px=0; px < numUnique; px++) {
		if (patternLik[px] < MIN_PATTERNLIK) {
			patternLik[px] = MIN_PATTERNLIK;
			warning("Likelihood of pattern %d is 0, forcing to %.3g",
				px, MIN_PATTERNLIK);
		}

		double weight = numIdentical[px] / patternLik[px];
		for (int d1=0; d1 < primaryDims; d1++) {
			ElatentMean[px * maxAbilities + d1] *= weight;
			for (int d2=0; d2 <= d1; d2++) {
				int loc = px * maxAbilities * maxAbilities + d2 * maxAbilities + d1;
				ElatentCov[loc] *= weight;
			}
		}
		for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
			ElatentMean[px * maxAbilities + sdim] *= weight;
			int loc = px * maxAbilities * maxAbilities + sdim * maxAbilities + sdim;
			ElatentCov[loc] *= weight;
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
			ElatentMean[d1] += ElatentMean[px * maxAbilities + d1];
			for (int d2=0; d2 <= d1; d2++) {
				int cell = d2 * maxAbilities + d1;
				int loc = px * maxAbilities * maxAbilities + cell;
				ElatentCov[cell] += ElatentCov[loc];
			}
		}
		for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
			ElatentMean[sdim] += ElatentMean[px * maxAbilities + sdim];
			int cell = sdim * maxAbilities + sdim;
			int loc = px * maxAbilities * maxAbilities + cell;
			ElatentCov[cell] += ElatentCov[loc];
		}
	}

	//pda(ElatentMean, 1, state->maxAbilities);
	//pda(ElatentCov, state->maxAbilities, state->maxAbilities);

	omxData *data = state->data;
	for (int d1=0; d1 < maxAbilities; d1++) {
		ElatentMean[d1] /= data->rows;
	}

	for (int d1=0; d1 < primaryDims; d1++) {
		for (int d2=0; d2 <= d1; d2++) {
			int cell = d2 * maxAbilities + d1;
			int tcell = d1 * maxAbilities + d2;
			ElatentCov[tcell] = ElatentCov[cell] =
				ElatentCov[cell] / data->rows - ElatentMean[d1] * ElatentMean[d2];
		}
	}
	for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
		int cell = sdim * maxAbilities + sdim;
		ElatentCov[cell] = ElatentCov[cell] / data->rows - ElatentMean[sdim] * ElatentMean[sdim];
	}

	//mxLog("E-step");
	//pda(ElatentMean, 1, state->maxAbilities);
	//pda(ElatentCov, state->maxAbilities, state->maxAbilities);
	state->validExpectation = TRUE;
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
static void
ba81SetupQuadrature(omxExpectation* oo, int gridsize, int flat)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	int numUnique = state->numUnique;
	int numThreads = Global->numThreads;
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

	double qgs = state->quadGridSize-1;
	for (int px=0; px < state->quadGridSize; px ++) {
		state->Qpoint[px] = Qwidth - px * 2 * Qwidth / qgs;
	}

	if (flat) {
		// not sure why this is useful, remove? TODO
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
		//pda(state->latentMeanOut->data, 1, state->maxAbilities);
		//pda(state->latentCovOut->data, state->maxAbilities, state->maxAbilities);

		double totalArea = 0;
		for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[priDims];
			decodeLocation(qx, priDims, state->quadGridSize, quad);
			double where[priDims];
			pointToWhere(state->Qpoint, quad, where, priDims);
			state->priLogQarea[qx] = dmvnorm(priDims, where,
							 state->latentMeanOut->data,
							 state->latentCovOut->data);
			totalArea += exp(state->priLogQarea[qx]);
		}
		totalArea = log(totalArea);
		for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
			state->priLogQarea[qx] -= totalArea;
			//mxLog("%.5g,", state->priLogQarea[qx]);
		}
		//mxLog("");

		for (int sx=0; sx < numSpecific; sx++) {
			totalArea = 0;
			for (int qx=0; qx < state->quadGridSize; qx++) {
				int covCell = (priDims + sx) * state->maxAbilities + priDims + sx;
				double den = dnorm(state->Qpoint[qx],
						   state->latentMeanOut->data[priDims + sx],
						   state->latentCovOut->data[covCell], TRUE);
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
expectedUpdate(omxData *data, const int *rowMap, const int px, const int item,
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
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	omxData *data = state->data;
	int numSpecific = state->numSpecific;
	const int *rowMap = state->rowMap;
	double *patternLik = state->patternLik;
	double *logNumIdentical = state->logNumIdentical;
	int numUnique = state->numUnique;
	int maxDims = state->maxDims;
	int numItems = state->EitemParam->cols;
	omxMatrix *itemSpec = state->itemSpec;
	int totalOutcomes = state->totalOutcomes;

	OMXZERO(state->expected, totalOutcomes * state->totalQuadPoints);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
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

#pragma omp parallel for num_threads(Global->numThreads)
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
EAPinternal(omxExpectation *oo, std::vector<double> *mean, std::vector<double> *cov)
{
	// add openmp parallelization stuff TODO
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	int numUnique = state->numUnique;
	int maxAbilities = state->maxAbilities;
	int covEntries = triangleLoc1(maxAbilities);
	int maxDims = state->maxDims;
	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);

	mean->assign(numUnique * maxAbilities, 0);
	cov->assign(numUnique * covEntries, 0);

	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, state->quadGridSize, quad);
		double where[priDims];
		pointToWhere(state->Qpoint, quad, where, priDims);
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

		for (int px=0; px < numUnique; px++) {
			double plik = exp(logArea + lxk[px]);
			int cx=0;
			for (int d1=0; d1 < priDims; d1++) {
				double piece = where[d1] * plik;
				(*mean)[px * maxAbilities + d1] += piece;
				for (int d2=0; d2 <= d1; d2++) {
					(*cov)[px * covEntries + cx] += where[d2] * piece;
					++cx;
				}
			}
		}
	}

	if (numSpecific) {
		std::vector<double> ris(numUnique);
		for (int sx=0; sx < numSpecific; sx++) {
			for (int sqx=0; sqx < state->quadGridSize; sqx++) {
				double area = exp(state->speLogQarea[sx * state->quadGridSize + sqx]);
				double ptArea = area * state->Qpoint[sqx];
				ris.assign(numUnique, 0);
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
				for (int px=0; px < numUnique; px++) {
					double piece = ris[px] * ptArea;
					int dim = priDims + sx;
					(*mean)[px * maxAbilities + dim] += piece;
					(*cov)[px * covEntries + triangleLoc0(dim)] += piece * state->Qpoint[sqx];
				}
			}
		}
	}
	//mxLog("step1");
	//pda(mean->data(), maxAbilities, 2);
	//pda(cov->data(), covEntries, numUnique);

	double *patternLik = state->patternLik;
	for (int px=0; px < numUnique; px++) {
		double denom = exp(patternLik[px]);
		for (int ax=0; ax < maxAbilities; ax++) {
			(*mean)[px * maxAbilities + ax] /= denom;
		}
		for (int cx=0; cx < triangleLoc1(priDims); ++cx) {
			(*cov)[px * covEntries + cx] /= denom;
		}
		for (int sx=0; sx < numSpecific; sx++) {
			(*cov)[px * covEntries + triangleLoc0(priDims + sx)] /= denom;
		}

		int cx=0;
		for (int a1=0; a1 < priDims; ++a1) {
			for (int a2=0; a2 <= a1; ++a2) {
				double ma1 = (*mean)[px * maxAbilities + a1];
				double ma2 = (*mean)[px * maxAbilities + a2];
				(*cov)[px * covEntries + cx] -= ma1 * ma2;
				++cx;
			}
		}
		for (int sx=0; sx < numSpecific; sx++) {
			int sdim = priDims + sx;
			double ma1 = (*mean)[px * maxAbilities + sdim];
			(*cov)[px * covEntries + triangleLoc0(sdim)] -= ma1 * ma1;
		}
        }
	//mxLog("step2");
	//pda(mean->data(), maxAbilities, 2);
	//pda(cov->data(), covEntries, numUnique);
}

static void
ba81Estep(omxExpectation *oo, const char *context) {
	if (!context) return;

	BA81Expect *state = (BA81Expect *) oo->argStruct;
	omxRecompute(state->EitemParam);
	omxRecompute(state->latentMeanOut);
	omxRecompute(state->latentCovOut);

	ba81Estep1(oo);
	if (strcmp(context, "E")==0) {
		// for E-M LL
		ba81Expected(oo);
	} else if (strcmp(context, "M")==0) {
		// for regular LL
		BA81Expect *state = (BA81Expect *) oo->argStruct;
		ba81SetupQuadrature(oo, state->targetQpoints, 0);
	} else {
		omxRaiseErrorf(globalState, "Unknown context '%s'", context);
	}
}

static void
copyScore(int rows, int maxAbilities, std::vector<double> &mean,
	  std::vector<double> &cov, const int rx, double *scores, const int dest)
{
	for (int ax=0; ax < maxAbilities; ++ax) {
		scores[rows * ax + dest] = mean[maxAbilities * rx + ax];
	}
	for (int ax=0; ax < maxAbilities; ++ax) {
		scores[rows * (maxAbilities + ax) + dest] =
			sqrt(cov[triangleLoc1(maxAbilities) * rx + triangleLoc0(ax)]);
	}
	for (int ax=0; ax < triangleLoc1(maxAbilities); ++ax) {
		scores[rows * (2*maxAbilities + ax) + dest] =
			cov[triangleLoc1(maxAbilities) * rx + ax];
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
static void
ba81PopulateAttributes(omxExpectation *oo, SEXP robj)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;

	if (state->scores == SCORES_OMIT || !state->validExpectation) return;

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	/*
	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 0);
	ba81Estep1(oo);
	*/

	std::vector<double> mean;
	std::vector<double> cov;
	EAPinternal(oo, &mean, &cov);

	int numUnique = state->numUnique;
	omxData *data = state->data;
	int maxAbilities = state->maxAbilities;
	int rows = state->scores == SCORES_FULL? data->rows : numUnique;
	int cols = 2 * maxAbilities + triangleLoc1(maxAbilities);
	SEXP Rscores;
	PROTECT(Rscores = allocMatrix(REALSXP, rows, cols));
	double *scores = REAL(Rscores);

	const int SMALLBUF = 10;
	char buf[SMALLBUF];
	SEXP names;
	PROTECT(names = allocVector(STRSXP, cols));
	for (int nx=0; nx < maxAbilities; ++nx) {
		snprintf(buf, SMALLBUF, "s%d", nx+1);
		SET_STRING_ELT(names, nx, mkChar(buf));
		snprintf(buf, SMALLBUF, "se%d", nx+1);
		SET_STRING_ELT(names, maxAbilities + nx, mkChar(buf));
	}
	for (int nx=0; nx < triangleLoc1(maxAbilities); ++nx) {
		snprintf(buf, SMALLBUF, "cov%d", nx+1);
		SET_STRING_ELT(names, maxAbilities*2 + nx, mkChar(buf));
	}
	SEXP dimnames;
	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 1, names);
	setAttrib(Rscores, R_DimNamesSymbol, dimnames);

	if (state->scores == SCORES_FULL) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (int rx=0; rx < numUnique; rx++) {
			int dups = omxDataNumIdenticalRows(state->data, state->rowMap[rx]);
			for (int dup=0; dup < dups; dup++) {
				int dest = omxDataIndex(data, state->rowMap[rx]+dup);
				copyScore(rows, maxAbilities, mean, cov, rx, scores, dest);
			}
		}
	} else {
#pragma omp parallel for num_threads(Global->numThreads)
		for (int rx=0; rx < numUnique; rx++) {
			copyScore(rows, maxAbilities, mean, cov, rx, scores, rx);
		}
	}

	setAttrib(robj, install("scores.out"), Rscores);
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		mxLog("Freeing %s function.", NAME);
	}
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	omxFreeAllMatrixData(state->itemSpec);
	omxFreeAllMatrixData(state->EitemParam);
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->latentMeanOut);
	omxFreeAllMatrixData(state->latentCovOut);
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
	Free(state->ElatentMean);
	Free(state->ElatentCov);
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

static void ignoreSetVarGroup(omxExpectation*, FreeVarGroup *)
{}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		mxLog("Initializing %s.", NAME);
	}
	if (!rpf_model) {
		if (0) {
			const int wantVersion = 3;
			int version;
			get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model_GPL");
			(*get_librpf)(&version, &rpf_numModels, &rpf_model);
			if (version < wantVersion) error("librpf binary API %d installed, at least %d is required",
							 version, wantVersion);
		} else {
			rpf_numModels = librpf_numModels;
			rpf_model = librpf_model;
		}
	}
	
	BA81Expect *state = Calloc(1, BA81Expect);
	oo->argStruct = (void*) state;

	PROTECT(tmp = GET_SLOT(rObj, install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf(currentState, "%s unable to handle data type %s", NAME, omxDataType(state->data));
		return;
	}

	// change to regular matrices instead of MxMatrix TODO
	state->itemSpec =
		omxNewMatrixFromSlot(rObj, currentState, "ItemSpec");
	state->design =
		omxNewMatrixFromSlot(rObj, currentState, "Design");

	state->latentMeanOut = omxNewMatrixFromSlot(rObj, currentState, "mean"); // move to FitFunction? TODO
	if (!state->latentMeanOut) error("Failed to retrieve mean matrix");
	state->latentCovOut  = omxNewMatrixFromSlot(rObj, currentState, "cov");
	if (!state->latentCovOut) error("Failed to retrieve cov matrix");

	state->EitemParam =
		omxNewMatrixFromSlot(rObj, currentState, "EItemParam");
	if (!state->EitemParam) error("Must supply EItemParam");

	oo->computeFun = ba81Estep;
	oo->setVarGroup = ignoreSetVarGroup;
	oo->destructFun = ba81Destroy;
	oo->populateAttrFun = ba81PopulateAttributes;
	
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

	int numItems = state->EitemParam->cols;
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

	int numThreads = Global->numThreads;

	int maxSpec = 0;
	int maxParam = 0;
	state->maxDims = 0;
	state->maxOutcomes = 0;

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
	if (state->EitemParam->rows != maxParam) {
		omxRaiseErrorf(currentState, "ItemParam should have %d rows", maxParam);
		return;
	}

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

	if (state->latentMeanOut->rows * state->latentMeanOut->cols != state->maxAbilities) {
		error("The mean matrix '%s' must be 1x%d or %dx1", state->latentMeanOut->name,
		      state->maxAbilities, state->maxAbilities);
	}
	if (state->latentCovOut->rows != state->maxAbilities ||
	    state->latentCovOut->cols != state->maxAbilities) {
		error("The cov matrix '%s' must be %dx%d",
		      state->latentCovOut->name, state->maxAbilities, state->maxAbilities);
	}

	PROTECT(tmp = GET_SLOT(rObj, install("cache")));
	state->cacheLXK = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qpoints")));
	state->targetQpoints = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qwidth")));
	state->Qwidth = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("scores")));
	const char *score_option = CHAR(asChar(tmp));
	if (strcmp(score_option, "omit")==0) state->scores = SCORES_OMIT;
	if (strcmp(score_option, "unique")==0) state->scores = SCORES_UNIQUE;
	if (strcmp(score_option, "full")==0) state->scores = SCORES_FULL;

	state->ElatentMean = Realloc(NULL, state->maxAbilities * numUnique, double);
	state->ElatentCov = Realloc(NULL, state->maxAbilities * state->maxAbilities * numUnique, double);

	ba81SetupQuadrature(oo, state->targetQpoints, 0);

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData for column summary stats?
}
