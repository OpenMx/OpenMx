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

#include <Rmath.h>
#include "omxExpectationBA81.h"
#include "glue.h"
#include "libifa-rpf.h"
#include "dmvnorm.h"

const struct rpf *rpf_model = NULL;
int rpf_numModels;

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

// state->speQarea[sIndex(state, sx, qx)]
OMXINLINE static
int sIndex(BA81Expect *state, int sx, int qx)
{
	//if (sx < 0 || sx >= state->numSpecific) error("Out of domain");
	//if (qx < 0 || qx >= state->quadGridSize) error("Out of domain");
	return sx * state->quadGridSize + qx;
}

/**
 * \param theta Vector of ability parameters, one per ability
 * \returns A numItems by maxOutcomes colMajor vector of doubles. Caller must Free it.
 */
double *computeRPF(BA81Expect *state, omxMatrix *itemParam, const int *quad, const bool wantlog)
{
	omxMatrix *design = state->design;
	int maxDims = state->maxDims;
	int maxOutcomes = state->maxOutcomes;
	size_t numItems = state->itemSpec.size();

	double theta[maxDims];
	pointToWhere(state, quad, theta, maxDims);

	double *outcomeProb = Realloc(NULL, numItems * maxOutcomes, double);
	//double *outcomeProb = Calloc(numItems * maxOutcomes, double);

	for (size_t ix=0; ix < numItems; ix++) {
		const double *spec = state->itemSpec[ix];
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * maxOutcomes;
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double ptheta[dims];

		for (int dx=0; dx < dims; dx++) {
			int ability = (int)omxMatrixElement(design, dx, ix) - 1;
			if (ability >= maxDims) ability = maxDims-1;
			ptheta[dx] = theta[ability];
		}

		if (wantlog) {
			(*rpf_model[id].logprob)(spec, iparam, ptheta, out);
		} else {
			(*rpf_model[id].prob)(spec, iparam, ptheta, out);
		}
#if 0
		for (int ox=0; ox < spec[RPF_ISpecOutcomes]; ox++) {
			if (!isfinite(out[ox]) || out[ox] > 0) {
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
ba81Likelihood(omxExpectation *oo, const int thrId, int specific, const int *quad, const double *outcomeProb)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	int numUnique = state->numUnique;
	int maxOutcomes = state->maxOutcomes;
	omxData *data = state->data;
	size_t numItems = state->itemSpec.size();
	int *Sgroup = state->Sgroup;
	double *lxk;

	if (!state->cacheLXK) {
		lxk = state->lxk + numUnique * thrId;
	} else {
		lxk = getLXKcache(state, quad, specific);
	}

	const int *rowMap = state->rowMap;
	for (int px=0; px < numUnique; px++) {
		double lxk1 = 1;
		for (size_t ix=0; ix < numItems; ix++) {
			if (specific != Sgroup[ix]) continue;
			int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
			if (pick == NA_INTEGER) continue;
			double piece = outcomeProb[ix * maxOutcomes + pick-1];  // move -1 elsewhere TODO
			lxk1 *= piece;
			//mxLog("%d pick %d piece %.3f", ix, pick, piece);
		}
#if 0
#pragma omp critical(ba81LikelihoodDebug1)
		if (!isfinite(lxk1) || lxk1 > numItems) {
			mxLog("where");
			double where[state->maxDims];
			pointToWhere(state, quad, where, state->maxDims);
			pda(where, state->maxDims, 1);
			mxLog("prob");
			pda(outcomeProb, numItems, maxOutcomes);
			error("Likelihood of row %d is %f", rowMap[px], lxk1);
		}
#endif
		lxk[px] = lxk1;
	}

	//mxLog("L.is(%d) at (%d %d) %.2f", specific, quad[0], quad[1], lxk[0]);
	return lxk;
}

double *ba81LikelihoodFast(omxExpectation *oo, const int thrId, int specific, const int *quad)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (!state->cacheLXK) {
		const double *outcomeProb = computeRPF(state, state->EitemParam, quad, FALSE);
		double *ret = ba81Likelihood(oo, thrId, specific, quad, outcomeProb);
		Free(outcomeProb);
		return ret;
	} else {
		return getLXKcache(state, quad, specific);
	}

}

OMXINLINE static void
mapLatentSpace(BA81Expect *state, int sgroup, double piece, const double *where,
	       const std::vector<double> &whereGram, double *latentDist)
{
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int pmax = maxDims;
	if (state->numSpecific) pmax -= 1;

	if (sgroup == 0) {
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

	if (state->numSpecific) {
		int sdim = pmax + sgroup;
		double piece_w1 = piece * where[pmax];
		latentDist[sdim] += piece_w1;

		double piece_var = piece * whereGram[triangleLoc0(pmax)];
		int to = maxAbilities + triangleLoc0(sdim);
		latentDist[to] += piece_var;
	}
}

// Eslxk, allElxk (Ei, Eis) depend on the ordinate of the primary dimensions
void cai2010(omxExpectation* oo, const int thrId, int recompute, const int *primaryQuad)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int sDim = maxDims-1;
	int quadGridSize = state->quadGridSize;
	int quad[maxDims];
	memcpy(quad, primaryQuad, sizeof(int)*sDim);
	double *allElxk = eBase(state, thrId);
	double *Eslxk = esBase(state, thrId);

	for (int px=0; px < numUnique; px++) {
		allElxk[px] = 1;
		for (int sx=0; sx < numSpecific; sx++) {
			Eslxk[sx * numUnique + px] = 0;
		}
	}

	if (!state->cacheLXK) recompute = TRUE;

	for (int qx=0; qx < quadGridSize; qx++) {
		quad[sDim] = qx;
		double *outcomeProb = NULL;
		if (recompute) {
			outcomeProb = computeRPF(state, state->EitemParam, quad, FALSE);
		}
		for (int sx=0; sx < numSpecific; sx++) {
			double *myEslxk = Eslxk + sx * numUnique;
			double *lxk;     // a.k.a. "L_is"
			if (recompute) {
				lxk = ba81Likelihood(oo, thrId, sx, quad, outcomeProb);
			} else {
				lxk = getLXKcache(state, quad, sx);
			}

			for (int ix=0; ix < numUnique; ix++) {
				double area = state->speQarea[sIndex(state, sx, qx)];
				double piece = lxk[ix] * area;
				//mxLog("E.is(%d) at (%d, %d) %.2f + %.2f = %.2f",
				//  sx, primaryQuad[0], qx, lxk[ix], area, piece);
				myEslxk[ix] += piece;
			}
		}
		Free(outcomeProb);
	}

	for (int sx=0; sx < numSpecific; sx++) {
		for (int px=0; px < numUnique; px++) {
			//mxLog("E.is(%d) at (%d) %.2f", sx, primaryQuad[0], state->Eslxk[esIndex(state, sx, 0)]);
			allElxk[px] *= Eslxk[sx * numUnique + px];  // allSlxk a.k.a. "E_i"
		}
	}
}

static void ba81Estep1(omxExpectation *oo)
{
	if(OMX_DEBUG) {mxLog("Beginning %s Computation.", oo->name);}

	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (state->verbose) {
		mxLog("%s: lxk(%d) patternLik ElatentMean ElatentCov",
		      oo->name, omxGetMatrixVersion(state->EitemParam));
	}

	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int primaryDims = maxDims;

	state->patternLik = Realloc(state->patternLik, numUnique, double);
	double *patternLik = state->patternLik;
	OMXZERO(patternLik, numUnique);

	int numLatents = maxAbilities + triangleLoc1(maxAbilities);
	int numLatentsPerThread = numUnique * numLatents;
	double *latentDist = Calloc(numUnique * numLatents * Global->numThreads, double);

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
			int thrId = omx_absolute_thread_num();
			double *thrLatentDist = latentDist + thrId * numLatentsPerThread;
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(state, quad, where, maxDims);
			std::vector<double> whereGram(triangleLoc1(maxDims));
			gramProduct(where, maxDims, whereGram.data());

			const double *outcomeProb = computeRPF(state, state->EitemParam, quad, FALSE);
			double *lxk = ba81Likelihood(oo, thrId, 0, quad, outcomeProb);
			Free(outcomeProb);

			double area = state->priQarea[qx];
			for (int px=0; px < numUnique; px++) {
				double tmp = lxk[px] * area;
#if 0
				if (!isfinite(tmp)) {
					mxLog("where");
					pda(where, maxDims, 1);
					error("Row %d lxk %f logArea %f tmp %f",
					      state->rowMap[px], lxk[px], logArea, tmp);
				}
#endif
#pragma omp atomic
				patternLik[px] += tmp;
				mapLatentSpace(state, 0, tmp, where, whereGram,
					       thrLatentDist + px * numLatents);
			}
		}
	} else {
		primaryDims -= 1;
		int sDim = primaryDims;
		long specificPoints = state->quadGridSize;

#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int thrId = omx_absolute_thread_num();
			double *thrLatentDist = latentDist + thrId * numLatentsPerThread;
			int quad[maxDims];
			decodeLocation(qx, primaryDims, state->quadGridSize, quad);

			cai2010(oo, thrId, TRUE, quad);
			double *allElxk = eBase(state, thrId);
			double *Eslxk = esBase(state, thrId);

			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				double where[maxDims];
				pointToWhere(state, quad, where, maxDims);
				std::vector<double> whereGram(triangleLoc1(maxDims));
				gramProduct(where, maxDims, whereGram.data());

				for (int sgroup=0; sgroup < numSpecific; sgroup++) {
					double area = areaProduct(state, quad, sgroup);
					double *lxk = ba81LikelihoodFast(oo, thrId, sgroup, quad);
					for (int px=0; px < numUnique; px++) {
						double Ei = allElxk[px];
						double Eis = Eslxk[sgroup * numUnique + px];
						double tmp = ((Ei / Eis) * lxk[px] * area);
						mapLatentSpace(state, sgroup, tmp, where, whereGram,
							       thrLatentDist + px * numLatents);
					}
				}
			}

			double priArea = state->priQarea[qx];
			for (int px=0; px < numUnique; px++) {
				double Ei = allElxk[px];
				double tmp = (Ei * priArea);
#pragma omp atomic
				patternLik[px] += tmp;
			}
		}
	}

	int *numIdentical = state->numIdentical;

	//mxLog("raw latent");
	//pda(latentDist, numLatents, numUnique);

#pragma omp parallel for num_threads(Global->numThreads) schedule(dynamic)
	for (int lx=0; lx < maxAbilities + triangleLoc1(primaryDims); ++lx) {
		for (int tx=1; tx < Global->numThreads; ++tx) {
			double *thrLatentDist = latentDist + tx * numLatentsPerThread;
			for (int px=0; px < numUnique; px++) {
				int loc = px * numLatents + lx;
				latentDist[loc] += thrLatentDist[loc];
			}
		}
	}

#pragma omp parallel for num_threads(Global->numThreads)
	for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
		for (int tx=1; tx < Global->numThreads; ++tx) {
			double *thrLatentDist = latentDist + tx * numLatentsPerThread;
			for (int px=0; px < numUnique; px++) {
				int loc = px * numLatents + maxAbilities + triangleLoc0(sdim);
				latentDist[loc] += thrLatentDist[loc];
			}
		}
	}

#pragma omp parallel for num_threads(Global->numThreads)
	for (int px=0; px < numUnique; px++) {
		if (!isfinite(patternLik[px])) {
			omxRaiseErrorf(globalState, "Likelihood of pattern %d is %.3g",
				       px, patternLik[px]);
		}

		double *latentDist1 = latentDist + px * numLatents;
		double weight = numIdentical[px] / patternLik[px];
		int cx = maxAbilities;
		for (int d1=0; d1 < primaryDims; d1++) {
			latentDist1[d1] *= weight;
			for (int d2=0; d2 <= d1; d2++) {
				latentDist1[cx] *= weight;
				++cx;
			}
		}
		for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
			latentDist1[sdim] *= weight;
			int loc = maxAbilities + triangleLoc0(sdim);
			latentDist1[loc] *= weight;
		}
#if 0
		if (!isfinite(patternLik[px])) {
			error("Likelihood of row %d is %f", state->rowMap[px], patternLik[px]);
		}
#endif
	}

	//mxLog("raw latent after weighting");
	//pda(latentDist, numLatents, numUnique);

	std::vector<double> &ElatentMean = state->ElatentMean;
	std::vector<double> &ElatentCov = state->ElatentCov;
	
	ElatentMean.assign(ElatentMean.size(), 0.0);
	ElatentCov.assign(ElatentCov.size(), 0.0);

#pragma omp parallel for num_threads(Global->numThreads)
	for (int d1=0; d1 < maxAbilities; d1++) {
		for (int px=0; px < numUnique; px++) {
			double *latentDist1 = latentDist + px * numLatents;
			int cx = maxAbilities + triangleLoc1(d1);
			if (d1 < primaryDims) {
				ElatentMean[d1] += latentDist1[d1];
				for (int d2=0; d2 <= d1; d2++) {
					int cell = d2 * maxAbilities + d1;
					ElatentCov[cell] += latentDist1[cx];
					++cx;
				}
			} else {
				ElatentMean[d1] += latentDist1[d1];
				int cell = d1 * maxAbilities + d1;
				int loc = maxAbilities + triangleLoc0(d1);
				ElatentCov[cell] += latentDist1[loc];
			}
		}
	}

	//pda(ElatentMean.data(), 1, state->maxAbilities);
	//pda(ElatentCov.data(), state->maxAbilities, state->maxAbilities);

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

	if (state->cacheLXK) state->LXKcached = TRUE;

	Free(latentDist);

	//mxLog("E-step");
	//pda(ElatentMean.data(), 1, state->maxAbilities);
	//pda(ElatentCov.data(), state->maxAbilities, state->maxAbilities);
}

static int getLatentVersion(BA81Expect *state)
{
	return omxGetMatrixVersion(state->latentMeanOut) + omxGetMatrixVersion(state->latentCovOut);
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
static void ba81SetupQuadrature(omxExpectation* oo, int gridsize)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	if (state->verbose) {
		mxLog("%s: quadrature(%d)", oo->name, getLatentVersion(state));
	}
	int numUnique = state->numUnique;
	int numThreads = Global->numThreads;
	int maxDims = state->maxDims;
	double Qwidth = state->Qwidth;
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
		state->speQarea.resize(gridsize * numSpecific);
	}

	state->Qpoint.resize(gridsize);
	state->priQarea.resize(state->totalPrimaryPoints);

	double qgs = state->quadGridSize-1;
	for (int px=0; px < state->quadGridSize; px ++) {
		state->Qpoint[px] = Qwidth - px * 2 * Qwidth / qgs;
	}

	//pda(state->latentMeanOut->data, 1, state->maxAbilities);
	//pda(state->latentCovOut->data, state->maxAbilities, state->maxAbilities);

	double totalArea = 0;
	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, state->quadGridSize, quad);
		double where[priDims];
		pointToWhere(state, quad, where, priDims);
		state->priQarea[qx] = exp(dmvnorm(priDims, where,
						  state->latentMeanOut->data,
						  state->latentCovOut->data));
		totalArea += state->priQarea[qx];
	}
	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		state->priQarea[qx] /= totalArea;
		//mxLog("%.5g,", state->priQarea[qx]);
	}

	for (int sx=0; sx < numSpecific; sx++) {
		totalArea = 0;
		int covCell = (priDims + sx) * state->maxAbilities + priDims + sx;
		double mean = state->latentMeanOut->data[priDims + sx];
		double var = state->latentCovOut->data[covCell];
		//mxLog("setup[%d] %.2f %.2f", sx, mean, var);
		for (int qx=0; qx < state->quadGridSize; qx++) {
			double den = dnorm(state->Qpoint[qx], mean, sqrt(var), FALSE);
			state->speQarea[sIndex(state, sx, qx)] = den;
			totalArea += den;
		}
		for (int qx=0; qx < state->quadGridSize; qx++) {
			state->speQarea[sIndex(state, sx, qx)] /= totalArea;
		}
		//pda(state->speQarea.data() + sIndex(state, sx, 0), 1, state->quadGridSize);
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

static void ba81buildLXKcache(omxExpectation *oo)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	if (!state->cacheLXK || state->LXKcached) return;
	
	ba81Estep1(oo);
}

OMXINLINE static void
expectedUpdate(omxData *data, const int *rowMap, const int px, const int item,
	       const double observed, double *out)
{
	int pick = omxIntDataElementUnsafe(data, rowMap[px], item);
	if (pick != NA_INTEGER) {
		out[pick-1] += observed;
	}
}

OMXINLINE static void
ba81Expected(omxExpectation* oo)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (state->verbose) mxLog("%s: EM.expected", oo->name);

	omxData *data = state->data;
	int numSpecific = state->numSpecific;
	const int *rowMap = state->rowMap;
	double *patternLik = state->patternLik;
	int *numIdentical = state->numIdentical;
	int numUnique = state->numUnique;
	int maxDims = state->maxDims;
	int numItems = state->EitemParam->cols;
	int totalOutcomes = state->totalOutcomes;
	std::vector<int> &itemOutcomes = state->itemOutcomes;

	OMXZERO(state->expected, totalOutcomes * state->totalQuadPoints);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double *lxk = ba81LikelihoodFast(oo, thrId, 0, quad);
			for (int px=0; px < numUnique; px++) {
				double *out = state->expected + qx * totalOutcomes;
				double observed = numIdentical[px] * lxk[px] / patternLik[px];
				for (int ix=0; ix < numItems; ix++) {
					const int outcomes = itemOutcomes[ix];
					expectedUpdate(data, rowMap, px, ix, observed, out);
					out += outcomes;
				}
			}
		}
	} else {
		int sDim = state->maxDims-1;
		long specificPoints = state->quadGridSize;

#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);

			cai2010(oo, thrId, FALSE, quad);
			double *allElxk = eBase(state, thrId);
			double *Eslxk = esBase(state, thrId);

			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				long qloc = encodeLocation(state->maxDims, state->quadGridSize, quad);

				for (int sgroup=0; sgroup < numSpecific; sgroup++) {
					double *lxk = ba81LikelihoodFast(oo, thrId, sgroup, quad);
					double *myEslxk = Eslxk + sgroup * numUnique;

					for (int px=0; px < numUnique; px++) {
						double *out = state->expected + totalOutcomes * qloc;

						for (int ix=0; ix < numItems; ix++) {
							const int outcomes = itemOutcomes[ix];
							if (state->Sgroup[ix] == sgroup) {
								double Ei = allElxk[px];
								double Eis = myEslxk[px];
								double observed = (numIdentical[px] * (Ei / Eis) *
										   (lxk[px] / patternLik[px]));
								expectedUpdate(data, rowMap, px, ix, observed, out);
							}
							out += outcomes;
						}
					}
				}
			}
		}
	}

	if (!state->checkedBadData) {
		std::vector<double> byOutcome(totalOutcomes, 0);
		for (int ox=0; ox < totalOutcomes; ++ox) {
			for (long qx=0; qx < state->totalQuadPoints; qx++) {
				byOutcome[ox] += state->expected[totalOutcomes * qx + ox];
			}
			if (byOutcome[ox] == 0) {
				int uptoItem = 0;
				for (size_t cx = 0; cx < itemOutcomes.size(); cx++) {
					if (ox < uptoItem + itemOutcomes[cx]) {
						int bad = ox - uptoItem;
						omxRaiseErrorf(globalState, "Item %lu outcome %d is never endorsed.\n"
							       "You must collapse categories or omit this item to estimate item parameters.", 1+cx, 1+bad);
						break;
					}
					uptoItem += itemOutcomes[cx];
				}
			}
		}
		state->checkedBadData = TRUE;
	}
	//pda(state->expected, state->totalOutcomes, state->totalQuadPoints);
}

OMXINLINE static void
accumulateScores(BA81Expect *state, int px, int sgroup, double piece, const double *where,
		 int primaryDims, int covEntries, std::vector<double> *mean, std::vector<double> *cov)
{
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;

	if (sgroup == 0) {
		int cx=0;
		for (int d1=0; d1 < primaryDims; d1++) {
			double piece_w1 = piece * where[d1];
			double &dest1 = (*mean)[px * maxAbilities + d1];
#pragma omp atomic
			dest1 += piece_w1;
			for (int d2=0; d2 <= d1; d2++) {
				double &dest2 = (*cov)[px * covEntries + cx];
#pragma omp atomic
				dest2 += where[d2] * piece_w1;
				++cx;
			}
		}
	}

	if (state->numSpecific) {
		int sdim = maxDims + sgroup - 1;
		double piece_w1 = piece * where[primaryDims];
		double &dest3 = (*mean)[px * maxAbilities + sdim];
#pragma omp atomic
		dest3 += piece_w1;

		double &dest4 = (*cov)[px * covEntries + triangleLoc0(sdim)];
#pragma omp atomic
		dest4 += piece_w1 * where[primaryDims];
	}
}

static void
EAPinternalFast(omxExpectation *oo, std::vector<double> *mean, std::vector<double> *cov)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (state->verbose) mxLog("%s: EAP", oo->name);

	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int primaryDims = maxDims;
	int covEntries = triangleLoc1(maxAbilities);

	mean->assign(numUnique * maxAbilities, 0);
	cov->assign(numUnique * covEntries, 0);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(state, quad, where, maxDims);

			double *lxk = ba81LikelihoodFast(oo, thrId, 0, quad);

			double area = state->priQarea[qx];
			for (int px=0; px < numUnique; px++) {
				double tmp = lxk[px] * area;
				accumulateScores(state, px, 0, tmp, where, primaryDims, covEntries, mean, cov);
			}
		}
	} else {
		primaryDims -= 1;
		int sDim = primaryDims;
		long specificPoints = state->quadGridSize;

#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, primaryDims, state->quadGridSize, quad);

			cai2010(oo, thrId, FALSE, quad);
			double *allElxk = eBase(state, thrId);
			double *Eslxk = esBase(state, thrId);

			for (int sgroup=0; sgroup < numSpecific; sgroup++) {
				for (long sx=0; sx < specificPoints; sx++) {
					quad[sDim] = sx;
					double where[maxDims];
					pointToWhere(state, quad, where, maxDims);
					double area = areaProduct(state, quad, sgroup);
					double *lxk = ba81LikelihoodFast(oo, thrId, sgroup, quad);
					for (int px=0; px < numUnique; px++) {
						double Ei = allElxk[px];
						double Eis = Eslxk[sgroup * numUnique + px];
						double tmp = ((Ei / Eis) * lxk[px] * area);
						accumulateScores(state, px, sgroup, tmp, where, primaryDims,
								 covEntries, mean, cov);
					}
				}
			}
		}
	}

	double *patternLik = state->patternLik;
	for (int px=0; px < numUnique; px++) {
		double denom = patternLik[px];
		for (int ax=0; ax < maxAbilities; ax++) {
			(*mean)[px * maxAbilities + ax] /= denom;
		}
		for (int cx=0; cx < triangleLoc1(primaryDims); ++cx) {
			(*cov)[px * covEntries + cx] /= denom;
		}
		for (int sx=0; sx < numSpecific; sx++) {
			(*cov)[px * covEntries + triangleLoc0(primaryDims + sx)] /= denom;
		}
		int cx=0;
		for (int a1=0; a1 < primaryDims; ++a1) {
			for (int a2=0; a2 <= a1; ++a2) {
				double ma1 = (*mean)[px * maxAbilities + a1];
				double ma2 = (*mean)[px * maxAbilities + a2];
				(*cov)[px * covEntries + cx] -= ma1 * ma2;
				++cx;
			}
		}
		for (int sx=0; sx < numSpecific; sx++) {
			int sdim = primaryDims + sx;
			double ma1 = (*mean)[px * maxAbilities + sdim];
			(*cov)[px * covEntries + triangleLoc0(sdim)] -= ma1 * ma1;
		}
        }
}

static void recomputePatternLik(omxExpectation *oo)
{
	BA81Expect *estate = (BA81Expect*) oo->argStruct;
	if (estate->verbose) mxLog("%s: patternLik", oo->name);

	int numUnique = estate->numUnique;
	int numSpecific = estate->numSpecific;
	int maxDims = estate->maxDims;
	int primaryDims = maxDims;
	double *patternLik = estate->patternLik;
	OMXZERO(patternLik, numUnique);

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < estate->totalQuadPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, maxDims, estate->quadGridSize, quad);
			double where[maxDims];
			pointToWhere(estate, quad, where, maxDims);
			double area = estate->priQarea[qx];
			double *lxk = ba81LikelihoodFast(oo, thrId, 0, quad);

			for (int px=0; px < numUnique; px++) {
				double tmp = (lxk[px] * area);
#pragma omp atomic
				patternLik[px] += tmp;
			}
		}
	} else {
		primaryDims -= 1;

#pragma omp parallel for num_threads(Global->numThreads)
		for (long qx=0; qx < estate->totalPrimaryPoints; qx++) {
			const int thrId = omx_absolute_thread_num();
			int quad[maxDims];
			decodeLocation(qx, primaryDims, estate->quadGridSize, quad);

			cai2010(oo, thrId, FALSE, quad);
			double *allElxk = eBase(estate, thrId);

			double priArea = estate->priQarea[qx];
			for (int px=0; px < numUnique; px++) {
				double Ei = allElxk[px];
				double tmp = (Ei * priArea);
#pragma omp atomic
				patternLik[px] += tmp;
			}
		}
	}
}

static void
ba81compute(omxExpectation *oo, const char *context)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;

	if (context) {
		if (strcmp(context, "EM")==0) {
			state->type = EXPECTATION_AUGMENTED;
		} else if (context[0] == 0) {
			state->type = EXPECTATION_OBSERVED;
		} else {
			omxRaiseErrorf(globalState, "Unknown context '%s'", context);
			return;
		}
	}

	omxRecompute(state->EitemParam);

	bool itemClean = state->itemParamVersion == omxGetMatrixVersion(state->EitemParam);
	bool latentClean = state->latentParamVersion == getLatentVersion(state);

	if (state->verbose) {
		mxLog("%s: Qinit %d itemClean %d latentClean %d (1=clean)",
		      oo->name, state->Qpoint.size() != 0, itemClean, latentClean);
	}

	if (state->Qpoint.size() == 0 || !latentClean) {
		ba81SetupQuadrature(oo, state->targetQpoints);
	}
	if (itemClean) {
		ba81buildLXKcache(oo);
		if (!latentClean) recomputePatternLik(oo);
	} else {
		ba81Estep1(oo);
	}

	if (state->type == EXPECTATION_AUGMENTED) {
		ba81Expected(oo);
	}

	state->itemParamVersion = omxGetMatrixVersion(state->EitemParam);
	state->latentParamVersion = getLatentVersion(state);
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
	int maxAbilities = state->maxAbilities;

	SEXP Rmean, Rcov;
	PROTECT(Rmean = allocVector(REALSXP, maxAbilities));
	memcpy(REAL(Rmean), state->ElatentMean.data(), maxAbilities * sizeof(double));

	PROTECT(Rcov = allocMatrix(REALSXP, maxAbilities, maxAbilities));
	memcpy(REAL(Rcov), state->ElatentCov.data(), maxAbilities * maxAbilities * sizeof(double));

	setAttrib(robj, install("empirical.mean"), Rmean);
	setAttrib(robj, install("empirical.cov"), Rcov);

	if (state->type == EXPECTATION_AUGMENTED) {
		int numUnique = state->numUnique;
		int totalOutcomes = state->totalOutcomes;
		SEXP Rlik;
		SEXP Rexpected;

		PROTECT(Rlik = allocVector(REALSXP, numUnique));
		memcpy(REAL(Rlik), state->patternLik, sizeof(double) * numUnique);

		PROTECT(Rexpected = allocMatrix(REALSXP, totalOutcomes, state->totalQuadPoints));
		memcpy(REAL(Rexpected), state->expected, sizeof(double) * totalOutcomes * state->totalQuadPoints);

		setAttrib(robj, install("patternLikelihood"), Rlik);
		setAttrib(robj, install("em.expected"), Rexpected);
	}

	if (state->scores == SCORES_OMIT || state->type == EXPECTATION_UNINITIALIZED) return;

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
	EAPinternalFast(oo, &mean, &cov);

	int numUnique = state->numUnique;
	omxData *data = state->data;
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
		mxLog("Freeing %s function.", oo->name);
	}
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	omxFreeAllMatrixData(state->EitemParam);
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->latentMeanOut);
	omxFreeAllMatrixData(state->latentCovOut);
	omxFreeAllMatrixData(state->customPrior);
	omxFreeAllMatrixData(state->itemParam);
	Free(state->numIdentical);
	Free(state->rowMap);
	Free(state->patternLik);
	Free(state->lxk);
	Free(state->Eslxk);
	Free(state->allElxk);
	Free(state->Sgroup);
	Free(state->expected);
	delete state;
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
		mxLog("Initializing %s.", oo->name);
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
	
	BA81Expect *state = new BA81Expect;
	state->checkedBadData = FALSE;
	state->numSpecific = 0;
	state->numIdentical = NULL;
	state->rowMap = NULL;
	state->design = NULL;
	state->lxk = NULL;
	state->patternLik = NULL;
	state->Eslxk = NULL;
	state->allElxk = NULL;
	state->expected = NULL;
	state->type = EXPECTATION_UNINITIALIZED;
	state->scores = SCORES_OMIT;
	state->itemParam = NULL;
	state->customPrior = NULL;
	state->itemParamVersion = 0;
	state->latentParamVersion = 0;
	oo->argStruct = (void*) state;

	PROTECT(tmp = GET_SLOT(rObj, install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf(currentState, "%s unable to handle data type %s", oo->name, omxDataType(state->data));
		return;
	}

	PROTECT(tmp = GET_SLOT(rObj, install("ItemSpec")));
	for (int sx=0; sx < length(tmp); ++sx) {
		SEXP model = VECTOR_ELT(tmp, sx);
		if (!OBJECT(model)) {
			error("Item models must inherit rpf.base");
		}
		SEXP spec;
		PROTECT(spec = GET_SLOT(model, install("spec")));
		state->itemSpec.push_back(REAL(spec));
	}

	PROTECT(tmp = GET_SLOT(rObj, install("design")));
	if (!isNull(tmp)) {
		// better to demand integers and not coerce to real TODO
		state->design = omxNewMatrixFromRPrimitive(tmp, globalState, FALSE, 0);
	}

	state->latentMeanOut = omxNewMatrixFromSlot(rObj, currentState, "mean");
	if (!state->latentMeanOut) error("Failed to retrieve mean matrix");
	state->latentCovOut  = omxNewMatrixFromSlot(rObj, currentState, "cov");
	if (!state->latentCovOut) error("Failed to retrieve cov matrix");

	state->EitemParam =
		omxNewMatrixFromSlot(rObj, currentState, "EItemParam");
	if (!state->EitemParam) error("Must supply EItemParam");

	state->itemParam =
		omxNewMatrixFromSlot(rObj, globalState, "ItemParam");

	if (state->EitemParam->rows != state->itemParam->rows ||
	    state->EitemParam->cols != state->itemParam->cols) {
		error("ItemParam and EItemParam must be of the same dimension");
	}

	oo->computeFun = ba81compute;
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
		omxRaiseErrorf(currentState, "%s: all columns must be factors", oo->name);
		return;
	}

	for (int rx=0; rx < data->rows;) {
		rx += omxDataNumIdenticalRows(state->data, rx);
		++numUnique;
	}
	state->numUnique = numUnique;

	state->rowMap = Realloc(NULL, numUnique, int);
	state->numIdentical = Realloc(NULL, numUnique, int);

	state->customPrior =
		omxNewMatrixFromSlot(rObj, globalState, "CustomPrior");
	
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
		state->rowMap[ux] = rx;
		rx += dups;
	}

	int numThreads = Global->numThreads;

	int maxSpec = 0;
	int maxParam = 0;
	state->maxDims = 0;
	state->maxOutcomes = 0;

	std::vector<int> &itemOutcomes = state->itemOutcomes;
	itemOutcomes.resize(numItems);
	int totalOutcomes = 0;
	for (int cx = 0; cx < data->cols; cx++) {
		const double *spec = state->itemSpec[cx];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		if (state->maxDims < dims)
			state->maxDims = dims;

		int no = spec[RPF_ISpecOutcomes];
		itemOutcomes[cx] = no;
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

	if (int(state->itemSpec.size()) != data->cols) {
		omxRaiseErrorf(currentState, "ItemSpec must contain %d item model specifications",
			       data->cols);
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
			const double *spec = state->itemSpec[ix];
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
			const double *spec = state->itemSpec[ix];
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
		state->allElxk = Realloc(NULL, numUnique * numThreads, double);
		state->Eslxk = Realloc(NULL, numUnique * state->numSpecific * numThreads, double);
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

	PROTECT(tmp = GET_SLOT(rObj, install("verbose")));
	state->verbose = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("cache")));
	state->cacheLXK = asLogical(tmp);
	state->LXKcached = FALSE;

	PROTECT(tmp = GET_SLOT(rObj, install("qpoints")));
	state->targetQpoints = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qwidth")));
	state->Qwidth = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("scores")));
	const char *score_option = CHAR(asChar(tmp));
	if (strcmp(score_option, "omit")==0) state->scores = SCORES_OMIT;
	if (strcmp(score_option, "unique")==0) state->scores = SCORES_UNIQUE;
	if (strcmp(score_option, "full")==0) state->scores = SCORES_FULL;

	state->ElatentMean.resize(state->maxAbilities);
	state->ElatentCov.resize(state->maxAbilities * state->maxAbilities);

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData for column summary stats?
}
