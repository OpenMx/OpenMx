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

#include "omxExpectation.h"
#include "omxOpenmpWrap.h"
#include "npsolWrap.h"
#include "libirt-rpf.h"
#include "merge.h"

static const char *NAME = "ExpectationBA81";

static const struct rpf *rpf_model = NULL;
static int rpf_numModels;

typedef double *(*rpf_fn_t)(omxExpectation *oo, omxMatrix *itemParam, const int *quad);

typedef struct {

	// data characteristics
	omxData *data;
	int numUnique;
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
	int numQpoints;
	double *Qpoint;
	double *Qarea;
	double *logQarea;
	long *quadGridSize;       // maxDims
	long totalPrimaryPoints;  // product of quadGridSize except specific dim
	long totalQuadPoints;     // product of quadGridSize

	// estimation related
	omxMatrix *EitemParam;    // E step version
	omxMatrix *itemParam;     // M step version
	SEXP rpf;
	rpf_fn_t computeRPF;
	omxMatrix *customPrior;
	int *paramMap;
	double *thrGradient;      // thread * numFreeParams
	double *thrGradient1;     // thread * numFreeParams
	int cacheLXK;		  // w/cache,  numUnique * #specific quad points * totalQuadPoints
	double *lxk;              // wo/cache, numUnique * thread
	double *allSlxk;          // numUnique * thread
	double *Slxk;             // numUnique * #specific dimensions * thread
	double *patternLik;       // length numUnique
	double ll;                // the most recent finite ll

	int gradientCount;
	int fitCount;
} omxBA81State;

static void
pda(const double *ar, int rows, int cols) {
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}

}

static int
findFreeVarLocation(omxMatrix *itemParam, const omxFreeVar *fv)
{
	for (int lx=0; lx < fv->numLocations; lx++) {
		if (~fv->matrices[lx] == itemParam->matrixNumber) return lx;
	}
	return -1;
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

	state->paramMap = Realloc(NULL, itemParam->rows * itemParam->cols, int);
	for (int px=0; px < itemParam->rows * itemParam->cols; px++) {
		state->paramMap[px] = -1;
	}

	int numFreeParams = currentState->numFreeParams;
	for (int px=0; px < numFreeParams; px++) {
		omxFreeVar *fv = currentState->freeVarList + px;
		int loc = findFreeVarLocation(itemParam, fv);
		if (loc == -1) continue;
		state->paramMap[fv->col[loc] * itemParam->rows + fv->row[loc]] = px;
	}

	state->thrGradient = Realloc(NULL, numFreeParams * getNumThreads(oo), double);
	state->thrGradient1 = Realloc(NULL, numFreeParams * getNumThreads(oo), double);
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
 * This is the main function needed to generate simulated data from
 * the model. It could be argued that the rest of the estimation
 * machinery belongs in the fitfunction.
 *
 * \param theta Vector of ability parameters, one per ability
 * \returns A numItems by maxOutcomes colMajor vector of doubles. Caller must Free it.
 */
static double *
standardComputeRPF(omxExpectation *oo, omxMatrix *itemParam, const int *quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	int numItems = itemSpec->cols;
	omxMatrix *design = state->design;
	int maxDims = state->maxDims;

	double theta[maxDims];
	pointToWhere(state, quad, theta, maxDims);

	double *outcomeProb = Realloc(NULL, numItems * state->maxOutcomes, double);

	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * state->maxOutcomes;
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double ptheta[dims];
		assignDims(itemSpec, design, dims, maxDims, ix, theta, ptheta);
		(*rpf_model[id].logprob)(spec, iparam, ptheta, out);
	}

	return outcomeProb;
}

static double *
RComputeRPF1(omxExpectation *oo, omxMatrix *itemParam, const int *quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int maxOutcomes = state->maxOutcomes;
	omxMatrix *design = state->design;
	omxMatrix *itemSpec = state->itemSpec;
	int maxDims = state->maxDims;

	double theta[maxDims];
	pointToWhere(state, quad, theta, maxDims);

	SEXP invoke;
	PROTECT(invoke = allocVector(LANGSXP, 4));
	SETCAR(invoke, state->rpf);
	SETCADR(invoke, omxExportMatrix(itemParam));
	SETCADDR(invoke, omxExportMatrix(itemSpec));

	SEXP where;
	PROTECT(where = allocMatrix(REALSXP, maxDims, itemParam->cols));
	double *ptheta = REAL(where);
	for (int ix=0; ix < itemParam->cols; ix++) {
		int dims = omxMatrixElement(itemSpec, RPF_ISpecDims, ix);
		assignDims(itemSpec, design, dims, maxDims, ix, theta, ptheta + ix*maxDims);
		for (int dx=dims; dx < maxDims; dx++) {
			ptheta[ix*maxDims + dx] = NA_REAL;
		}
	}
	SETCADDDR(invoke, where);

	SEXP matrix;
	PROTECT(matrix = eval(invoke, R_GlobalEnv));

	if (!isMatrix(matrix)) {
		omxRaiseError(oo->currentState, -1,
			      "RPF must return an item by outcome matrix");
		return NULL;
	}

	SEXP matrixDims;
	PROTECT(matrixDims = getAttrib(matrix, R_DimSymbol));
	int *dimList = INTEGER(matrixDims);
	int numItems = state->itemSpec->cols;
	if (dimList[0] != maxOutcomes || dimList[1] != numItems) {
		const int errlen = 200;
		char errstr[errlen];
		snprintf(errstr, errlen, "RPF must return a %d outcomes by %d items matrix",
			 maxOutcomes, numItems);
		omxRaiseError(oo->currentState, -1, errstr);
		return NULL;
	}

	// Unlikely to be of type INTSXP, but just to be safe
	PROTECT(matrix = coerceVector(matrix, REALSXP));
	double *restrict got = REAL(matrix);

	// Need to copy because threads cannot share SEXP
	double *restrict outcomeProb = Realloc(NULL, numItems * maxOutcomes, double);

	// Double check there aren't NAs in the wrong place
	for (int ix=0; ix < numItems; ix++) {
		int numOutcomes = omxMatrixElement(state->itemSpec, RPF_ISpecOutcomes, ix);
		for (int ox=0; ox < numOutcomes; ox++) {
			int vx = ix * maxOutcomes + ox;
			if (isnan(got[vx])) {
				const int errlen = 200;
				char errstr[errlen];
				snprintf(errstr, errlen, "RPF returned NA in [%d,%d]", ox,ix);
				omxRaiseError(oo->currentState, -1, errstr);
			}
			outcomeProb[vx] = got[vx];
		}
	}

	return outcomeProb;
}

static double *
RComputeRPF(omxExpectation *oo, omxMatrix *itemParam, const int *quad)
{
	omx_omp_set_lock(&GlobalRLock);
	PROTECT_INDEX pi = omxProtectSave();
	double *ret = RComputeRPF1(oo, itemParam, quad);
	omxProtectRestore(pi);
	omx_omp_unset_lock(&GlobalRLock);  // hope there was no exception!
	return ret;
}

OMXINLINE static long
encodeLocation(const int dims, const long *restrict grid, const int *restrict quad)
{
	long qx = 0;
	for (int dx=dims-1; dx >= 0; dx--) {
		qx = qx * grid[dx];
		qx += quad[dx];
	}
	return qx;
}

#define CALC_LXK_CACHED(state, numUnique, quad, tqp, specific) \
	((state)->lxk + \
	 (numUnique) * encodeLocation((state)->maxDims, (state)->quadGridSize, quad) + \
	 (numUnique) * (tqp) * (specific))

OMXINLINE static double *
ba81Likelihood(omxExpectation *oo, int specific, const int *restrict quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int numUnique = state->numUnique;
	int maxOutcomes = state->maxOutcomes;
	omxData *data = state->data;
	int numItems = state->itemSpec->cols;
	rpf_fn_t rpf_fn = state->computeRPF;
	int *restrict Sgroup = state->Sgroup;
	double *restrict lxk;

	if (!state->cacheLXK) {
		lxk = state->lxk + numUnique * omx_absolute_thread_num();
	} else {
		lxk = CALC_LXK_CACHED(state, numUnique, quad, state->totalQuadPoints, specific);
	}

	const double *outcomeProb = (*rpf_fn)(oo, state->EitemParam, quad);
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
			lxk1 += outcomeProb[ix * maxOutcomes + pick-1];
		}
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
		return ba81LikelihoodFast(oo, specific, quad);
	} else {
		return CALC_LXK_CACHED(state, state->numUnique, quad, state->totalQuadPoints, specific);
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
		int quadGridSize = state->quadGridSize[sDim];

		for (int qx=0; qx < quadGridSize; qx++) {
			quad[sDim] = qx;
			double *lxk;
			if (recompute) {
				lxk = ba81Likelihood(oo, sx, quad);
			} else {
				lxk = CALC_LXK_CACHED(state, numUnique, quad, state->totalQuadPoints, sx);
			}

			for (int ix=0; ix < numUnique; ix++) {
				eis[ix] += exp(lxk[ix] + state->logQarea[qx]);
			}
		}

		for (int px=0; px < numUnique; px++) {
			eis[px] = log(eis[px]);
			allSlxk[px] += eis[px];
		}
	}
}

OMXINLINE static double
logAreaProduct(omxBA81State *state, const int *restrict quad, const int upto)
{
	double logArea = 0;
	for (int dx=0; dx < upto; dx++) {
		logArea += state->logQarea[quad[dx]];
	}
	return logArea;
}

// The idea of this API is to allow passing in a number larger than 1.
OMXINLINE static void
areaProduct(omxBA81State *state, const int *restrict quad, const int upto, double *restrict out)
{
	for (int dx=0; dx < upto; dx++) {
		*out *= state->Qarea[quad[dx]];
	}
}

OMXINLINE static void
decodeLocation(long qx, const int dims, const long *restrict grid,
	       int *restrict quad)
{
	for (int dx=0; dx < dims; dx++) {
		quad[dx] = qx % grid[dx];
		qx = qx / grid[dx];
	}
}

static void
ba81Estep(omxExpectation *oo) {
	if(OMX_DEBUG_MML) {Rprintf("Beginning %s Computation.\n", NAME);}

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	double *patternLik = state->patternLik;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;

	omxCopyMatrix(state->EitemParam, state->itemParam);

	OMXZERO(patternLik, numUnique);

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
			int quad[state->maxDims];
			decodeLocation(qx, state->maxDims, state->quadGridSize, quad);

			double *lxk = ba81Likelihood(oo, 0, quad);

			double logArea = logAreaProduct(state, quad, state->maxDims);
#pragma omp critical(EstepUpdate)
			for (int px=0; px < numUnique; px++) {
				double tmp = exp(lxk[px] + logArea);
				patternLik[px] += tmp;
			}
		}
	} else {
		int sDim = state->maxDims-1;

#pragma omp parallel for num_threads(oo->currentState->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[state->maxDims];
			decodeLocation(qx, sDim, state->quadGridSize, quad);

			double *allSlxk = CALC_ALLSLXK(state, numUnique);
			double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
			cai2010(oo, TRUE, quad, allSlxk, Slxk);

			double logArea = logAreaProduct(state, quad, sDim);
#pragma omp critical(EstepUpdate)
			for (int px=0; px < numUnique; px++) {
				double tmp = exp(allSlxk[px] + logArea);
				patternLik[px] += tmp;
			}
		}
	}

	for (int px=0; px < numUnique; px++) {
		patternLik[px] = log(patternLik[px]);
	}
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

/** 
 * \param quad a vector that indexes into a multidimensional quadrature
 * \param out points to an array numOutcomes wide
 */
OMXINLINE static void
ba81Weight(omxExpectation* oo, const int item, const int *quad, int outcomes, double *out)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxData *data = state->data;
	const int *rowMap = state->rowMap;
	int specific = state->Sgroup[item];
	double *patternLik = state->patternLik;
	double *logNumIdentical = state->logNumIdentical;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int sDim = state->maxDims-1;

	OMXZERO(out, outcomes);

	if (numSpecific == 0) {
		double *lxk = ba81LikelihoodFast(oo, specific, quad);
		for (int px=0; px < numUnique; px++) {
			double observed = logNumIdentical[px] + lxk[px] - patternLik[px];
			expectedUpdate(data, rowMap, px, item, observed, outcomes, out);
		}
	} else {
		double *allSlxk = CALC_ALLSLXK(state, numUnique);
		double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
		if (quad[sDim] == 0) {
			// allSlxk, Slxk only depend on the ordinate of the primary dimensions
			cai2010(oo, !state->cacheLXK, quad, allSlxk, Slxk);
		}
		double *eis = Slxk + numUnique * specific;

		// Avoid recalc with modest buffer? TODO
		double *lxk = ba81LikelihoodFast(oo, specific, quad);

		for (int px=0; px < numUnique; px++) {
			double observed = logNumIdentical[px] + (allSlxk[px] - eis[px]) +
				(lxk[px] - patternLik[px]);
			expectedUpdate(data, rowMap, px, item, observed, outcomes, out);
		}
	}
}

OMXINLINE static double
ba81Fit1Ordinate(omxExpectation* oo, const int *quad, int want)
{
	omxState* currentState = oo->currentState;
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	omxMatrix *itemParam = state->itemParam;
	int numItems = itemParam->cols;
	rpf_fn_t rpf_fn = state->computeRPF;
	int maxOutcomes = state->maxOutcomes;
	int maxDims = state->maxDims;
	int numFreeParams = currentState->numFreeParams;
	double *gradient = state->thrGradient1 + numFreeParams * omx_absolute_thread_num();
	int do_gradient = want & FF_COMPUTE_GRADIENT;

	double where[maxDims];
	pointToWhere(state, quad, where, maxDims);

	double *outcomeProb = (*rpf_fn)(oo, itemParam, quad); // avoid malloc/free? TODO
	if (!outcomeProb) return 0;

	double thr_ll = 0;
	for (int ix=0; ix < numItems; ix++) {
		const double *spec = omxMatrixColumn(itemSpec, ix);
		int id = spec[RPF_ISpecID];
		int outcomes = spec[RPF_ISpecOutcomes];

		double weight[outcomes];
		ba81Weight(oo, ix, quad, outcomes, weight);
		for (int ox=0; ox < outcomes; ox++) {
#if 0
			if (!isnormal(outcomeProb[ix * maxOutcomes + ox])) {
				pda(itemParam->data, itemParam->rows, itemParam->cols);
				pda(outcomeProb, outcomes, numItems);
				error("RPF produced NAs");
			}
#endif
			double got = weight[ox] * outcomeProb[ix * maxOutcomes + ox];
			areaProduct(state, quad, maxDims, &got);
			thr_ll += got;
		}

		if (do_gradient) {
			int *mask = state->paramMap + itemParam->rows * ix;
			double *iparam = omxMatrixColumn(itemParam, ix);
			(*rpf_model[id].gradient)(spec, iparam, mask, where, weight, gradient);
		}
	}

	Free(outcomeProb);

	if (do_gradient) {
		double *thrG = state->thrGradient + numFreeParams * omx_absolute_thread_num();

		for (int ox=0; ox < numFreeParams; ox++) {
#if 0
			if (!isnormal(gradient[ox])) {
				Rprintf("item spec:\n");
				pda(spec, (*rpf_model[id].numSpec)(spec), 1);
				Rprintf("item parameters:\n");
				pda(iparam, numParam, 1);
				Rprintf("where:\n");
				pda(where, maxDims, 1);
				Rprintf("weight:\n");
				pda(weight, outcomes, 1);
				error("Gradient for item %d param %d is %f; are you missing a lbound/ubound?",
				      item, ox, gradient[ox]);
			}
#endif

			areaProduct(state, quad, maxDims, gradient+ox);
			thrG[ox] += gradient[ox];
		}
	}

	return thr_ll;
}

static double
ba81ComputeFit1(omxExpectation* oo, int want, double *gradient)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxState* currentState = oo->currentState;
	omxMatrix *customPrior = state->customPrior;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *itemSpec = state->itemSpec;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;

	double ll = 0;
	if (customPrior) {
		omxRecompute(customPrior);
		ll = customPrior->data[0];
	} else {
		int numItems = itemSpec->cols;
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = omxMatrixColumn(itemSpec, ix);
			int id = spec[RPF_ISpecID];
			double *iparam = omxMatrixColumn(itemParam, ix);
			ll += (*rpf_model[id].prior)(spec, iparam);
		}
	}

	if (!isnormal(ll)) {
		omxPrint(itemParam, "item param");
		error("Bayesian prior returned NA; do you need to add a lbound/ubound?");
	}

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(currentState->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);
			double thr_ll = ba81Fit1Ordinate(oo, quad, want);

#pragma omp atomic
			ll += thr_ll;
		}
	} else {
		int sDim = state->maxDims-1;
		long *quadGridSize = state->quadGridSize;

#pragma omp parallel for num_threads(currentState->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, quadGridSize, quad);

			double thr_ll = 0;
			long specificPoints = quadGridSize[sDim];
			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				thr_ll += ba81Fit1Ordinate(oo, quad, want);
			}
#pragma omp atomic
			ll += thr_ll;
		}
	}

	if (!customPrior && gradient) {
		int numFreeParams = currentState->numFreeParams;

		for (int th=0; th < getNumThreads(oo); th++) {
			double *thrG = state->thrGradient + numFreeParams * th;

			for (int ox=0; ox < numFreeParams; ox++) {
				gradient[ox] += thrG[ox];
			}
		}

		int numItems = itemParam->cols;
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = omxMatrixColumn(itemSpec, ix);
			int id = spec[RPF_ISpecID];
			int *mask = state->paramMap + itemParam->rows * ix;
			double *iparam = omxMatrixColumn(itemParam, ix);
			(*rpf_model[id].gradient)(spec, iparam, mask, NULL, NULL, gradient);
		}

		for (int ox=0; ox < numFreeParams; ox++) {
			// Need to check because this can happen if
			// lbounds/ubounds are not set appropriately.
			if (!isnormal(gradient[ox])) {
				omxFreeVar *fv = currentState->freeVarList + ox;
				int loc = findFreeVarLocation(itemParam, fv);
				if (loc == -1) continue;

				Rprintf("item parameters:\n");
				const double *spec = omxMatrixColumn(itemSpec, fv->col[loc]);
				int id = spec[RPF_ISpecID];
				int numParam = (*rpf_model[id].numParam)(spec);
				double *iparam = omxMatrixColumn(itemParam, fv->col[loc]);
				pda(iparam, numParam, 1);
				error("Gradient for item %d is %f; are you missing a lbound/ubound?",
				      fv->col[loc], gradient[ox]);
			}

			gradient[ox] = -2 * gradient[ox];
		}
	}

	if (isinf(ll)) {
		// This is a hack to avoid the need to specify
		// ubound/lbound on parameters. Bounds are necessary
		// mainly for debugging derivatives.
		// Perhaps bounds can be pull in from librpf? TODO
		return 2*state->ll;
	} else {
		ll = -2 * ll;
		state->ll = ll;
		return ll;
	}
}

double
ba81ComputeFit(omxExpectation* oo, int want, double *gradient)
{
	if (!want) return 0;

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxState* currentState = oo->currentState;

	if (!state->paramMap) buildParamMap(oo);

	if (want & FF_COMPUTE_FIT) {
		++state->fitCount;
	}

	if (want & FF_COMPUTE_GRADIENT) {
		++state->gradientCount;

		int numFreeParams = currentState->numFreeParams;
		OMXZERO(gradient, numFreeParams);
		OMXZERO(state->thrGradient, numFreeParams * getNumThreads(oo));
	}

	double got = ba81ComputeFit1(oo, want, gradient);
	return got;
}

static void
ba81SetupQuadrature(omxExpectation* oo, int numPoints, double *points, double *area)
{
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	int numUnique = state->numUnique;
	int numThreads = getNumThreads(oo);

	state->numQpoints = numPoints;

	Free(state->Qpoint);
	Free(state->Qarea);
	state->Qpoint = Realloc(NULL, numPoints, double);
	state->Qarea = Realloc(NULL, numPoints, double);
	memcpy(state->Qpoint, points, sizeof(double)*numPoints);
	memcpy(state->Qarea, area, sizeof(double)*numPoints);

	Free(state->logQarea);

	state->logQarea = Realloc(NULL, state->numQpoints, double);
	for (int px=0; px < state->numQpoints; px++) {
		state->logQarea[px] = log(state->Qarea[px]);
	}

	state->totalQuadPoints = 1;
	state->totalPrimaryPoints = 1;
	state->quadGridSize = (long*) R_alloc(state->maxDims, sizeof(long));
	for (int dx=0; dx < state->maxDims; dx++) {
		state->quadGridSize[dx] = state->numQpoints;
		state->totalQuadPoints *= state->quadGridSize[dx];
		if (dx < state->maxDims-1) {
			state->totalPrimaryPoints *= state->quadGridSize[dx];
		}
	}

	Free(state->lxk);

	if (!state->cacheLXK) {
		state->lxk = Realloc(NULL, numUnique * numThreads, double);
	} else {
		int ns = state->numSpecific;
		if (ns == 0) ns = 1;
		state->lxk = Realloc(NULL, numUnique * state->totalQuadPoints * ns, double);
	}
}

static void
ba81EAP1(omxExpectation *oo, double *workspace, long qx, int maxDims, int numUnique,
	 double *ability, double *cov, double *spstats)
{
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	double *patternLik = state->patternLik;
	int quad[maxDims];
	decodeLocation(qx, maxDims, state->quadGridSize, quad);
	double where[maxDims];
	pointToWhere(state, quad, where, maxDims);
	double logArea = logAreaProduct(state, quad, maxDims);
	double *lxk = ba81LikelihoodFast(oo, 0, quad);
	double *myspace = workspace + 2 * maxDims * numUnique * omx_absolute_thread_num();

	for (int px=0; px < numUnique; px++) {
		double *piece = myspace + px * 2 * maxDims;
		double plik = exp(lxk[px] + logArea - patternLik[px]);
		for (int dx=0; dx < maxDims; dx++) {
			piece[dx] = where[dx] * plik;
		}
		/*
		for (int d1=0; d1 < maxDims; d1++) {
			for (int d2=0; d2 <= d1; d2++) {
				covPiece[d1 * maxDims + d2] += piece[d1] * piece[d2];
			}
		}
		*/
	}
#pragma omp critical(EAP1Update)
	for (int px=0; px < numUnique; px++) {
		double *piece = myspace + px * 2 * maxDims;
		double *arow = ability + px * 2 * maxDims;
		for (int dx=0; dx < maxDims; dx++) {
			arow[dx*2] += piece[dx];
		}
		/*
		for (int d1=0; d1 < maxDims; d1++) {
			for (int d2=0; d2 <= d1; d2++) {
				int loc = d1 * maxDims + d2;
				cov[loc] += covPiece[loc];
			}
		}
		*/
	}
}

static void
ba81EAP2(omxExpectation *oo, double *workspace, long qx, int maxDims, int numUnique,
	 double *ability, double *spstats)
{
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	double *patternLik = state->patternLik;
	int quad[maxDims];
	decodeLocation(qx, maxDims, state->quadGridSize, quad);
	double where[maxDims];
	pointToWhere(state, quad, where, maxDims);
	double logArea = logAreaProduct(state, quad, maxDims);
	double *lxk = ba81LikelihoodFast(oo, 0, quad);

	for (int px=0; px < numUnique; px++) {
		double psd[maxDims];
		double *arow = ability + px * 2 * maxDims;
		for (int dx=0; dx < maxDims; dx++) {
			// is this just sqrt(variance) and redundant with the covariance matrix? TODO
			double ldiff = log(fabs(where[dx] - arow[dx*2]));
			psd[dx] = exp(2 * ldiff + lxk[px] + logArea - patternLik[px]);
		}
#pragma omp critical(EAP1Update)
		for (int dx=0; dx < maxDims; dx++) {
			arow[dx*2+1] += psd[dx];
		}
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
	//int numSpecific = state->numSpecific;

	*numReturns = 2; // + (maxDims > 1) + (numSpecific > 1);
	omxRListElement *out = (omxRListElement*) R_alloc(*numReturns, sizeof(omxRListElement));

	out[0].numValues = 1;
	out[0].values = (double*) R_alloc(1, sizeof(double));
	strcpy(out[0].label, "Minus2LogLikelihood");
	out[0].values[0] = state->ll;

	omxData *data = state->data;
	int numUnique = state->numUnique;

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	int numQpoints = state->numQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	double Qpoint[numQpoints];
	double Qarea[numQpoints];
	const double Qwidth = 4;
	for (int qx=0; qx < numQpoints; qx++) {
		Qpoint[qx] = Qwidth - qx * Qwidth*2 / (numQpoints-1);
		Qarea[qx] = 1.0/numQpoints;
	}
	ba81SetupQuadrature(oo, numQpoints, Qpoint, Qarea);
	ba81Estep(oo);   // recalc patternLik with a flat prior

	double *cov = NULL;
	/*
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
	double *spstats = NULL;
	/*
	if (numSpecific) {
		strcpy(out[3].label, "specific");
		out[3].numValues = -1;
		out[3].rows = numSpecific;
		out[3].cols = 2;
		out[3].values = (double*) R_alloc(out[3].rows * out[3].cols, sizeof(double));
		spstats = out[3].values;
	}
	*/

	// allocation of workspace could be optional
	int numThreads = getNumThreads(oo);
	double *workspace = Realloc(NULL, numUnique * maxDims * 2 * numThreads, double);

	// Need a separate work space because the destination needs
	// to be in unsorted order with duplicated rows.
	double *ability = Calloc(numUnique * maxDims * 2, double);

#pragma omp parallel for num_threads(oo->currentState->numThreads)
	for (long qx=0; qx < state->totalQuadPoints; qx++) {
		ba81EAP1(oo, workspace, qx, maxDims, numUnique, ability, cov, spstats);
	}

	/*
	// make symmetric
	for (int d1=0; d1 < maxDims; d1++) {
		for (int d2=0; d2 < d1; d2++) {
			cov[d2 * maxDims + d1] = cov[d1 * maxDims + d2];
		}
	}
	*/

#pragma omp parallel for num_threads(oo->currentState->numThreads)
	for (long qx=0; qx < state->totalQuadPoints; qx++) {
		ba81EAP2(oo, workspace, qx, maxDims, numUnique, ability, spstats);
	}

	for (int px=0; px < numUnique; px++) {
		double *arow = ability + px * 2 * maxDims;
		for (int dx=0; dx < maxDims; dx++) {
			arow[dx*2+1] = sqrt(arow[dx*2+1]);
		}
	}

	strcpy(out[1].label, "ability");
	out[1].numValues = -1;
	out[1].rows = data->rows;
	out[1].cols = 2 * maxDims;
	out[1].values = (double*) R_alloc(out[1].rows * out[1].cols, sizeof(double));

	for (int rx=0; rx < numUnique; rx++) {
		double *pa = ability + rx * 2 * maxDims;

		int dups = omxDataNumIdenticalRows(state->data, state->rowMap[rx]);
		for (int dup=0; dup < dups; dup++) {
			int dest = omxDataIndex(data, state->rowMap[rx]+dup);
			int col=0;
			for (int dx=0; dx < maxDims; dx++) {
				out[1].values[col * out[1].rows + dest] = pa[col]; ++col;
				out[1].values[col * out[1].rows + dest] = pa[col]; ++col;
			}
		}
	}
	Free(ability);
	Free(workspace);
	return out;
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		Rprintf("Freeing %s function.\n", NAME);
	}
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	Rprintf("fit %d gradient %d\n", state->fitCount, state->gradientCount);
	omxFreeAllMatrixData(state->itemSpec);
	omxFreeAllMatrixData(state->itemParam);
	omxFreeAllMatrixData(state->EitemParam);
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->customPrior);
	Free(state->logNumIdentical);
	Free(state->Qpoint);
	Free(state->Qarea);
	Free(state->logQarea);
	Free(state->rowMap);
	Free(state->patternLik);
	Free(state->lxk);
	Free(state->Slxk);
	Free(state->allSlxk);
	Free(state->Sgroup);
	Free(state->paramMap);
	Free(state->thrGradient);
	Free(state->thrGradient1);
	Free(state);
}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		Rprintf("Initializing %s.\n", NAME);
	}
	if (!rpf_model) {
		const int wantVersion = 2;
		int version;
		get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model");
		(*get_librpf)(&version, &rpf_numModels, &rpf_model);
		if (version < wantVersion) error("librpf binary API %d installed, at least %d is required",
						 version, wantVersion);
	}
	
	omxBA81State *state = Calloc(1, omxBA81State);
	oo->argStruct = (void*) state;

	state->ll = 1e20;   // finite but big
	
	PROTECT(tmp = GET_SLOT(rObj, install("data")));
	state->data = omxNewDataFromMxDataPtr(tmp, currentState);
        UNPROTECT(1);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf(currentState, "%s unable to handle data type %s", NAME, omxDataType(state->data));
		return;
	}

	PROTECT(state->rpf = GET_SLOT(rObj, install("RPF")));
	if (state->rpf == R_NilValue) {
		state->computeRPF = standardComputeRPF;
	} else {
		state->computeRPF = RComputeRPF;
	}

	state->itemSpec =
		omxNewMatrixFromIndexSlot(rObj, currentState, "ItemSpec");
	state->design =
		omxNewMatrixFromIndexSlot(rObj, currentState, "Design");
	state->itemParam =
		omxNewMatrixFromIndexSlot(rObj, currentState, "ItemParam");
	state->EitemParam =
		omxInitTemporaryMatrix(NULL, state->itemParam->rows, state->itemParam->cols,
				       TRUE, currentState);
	state->customPrior =
		omxNewMatrixFromIndexSlot(rObj, currentState, "CustomPrior");
	
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
	state->logNumIdentical = Realloc(NULL, numUnique, double);

	int numItems = state->itemParam->cols;

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

	for (int cx = 0; cx < data->cols; cx++) {
		const double *spec = omxMatrixColumn(state->itemSpec, cx);
		int id = spec[RPF_ISpecID];
		if (id < 0 || id >= rpf_numModels) {
			omxRaiseErrorf(currentState, "ItemSpec column %d has unknown item model %d", cx, id);
			return;
		}

		int dims = spec[RPF_ISpecDims];
		if (state->maxDims < dims)
			state->maxDims = dims;

		// TODO verify that item model can have requested number of outcomes
		int no = spec[RPF_ISpecOutcomes];
		if (state->maxOutcomes < no)
			state->maxOutcomes = no;

		int numSpec = (*rpf_model[id].numSpec)(spec);
		if (maxSpec < numSpec)
			maxSpec = numSpec;

		int numParam = (*rpf_model[id].numParam)(spec);
		if (maxParam < numParam)
			maxParam = numParam;
	}

	if (state->itemSpec->cols != data->cols || state->itemSpec->rows != maxSpec) {
		omxRaiseErrorf(currentState, "ItemSpec must have %d item columns and %d rows",
			       data->cols, maxSpec);
		return;
	}
	if (state->itemParam->rows != maxParam) {
		omxRaiseErrorf(currentState, "ItemParam should have %d rows", maxParam);
		return;
	}

	if (state->design == NULL) {
		state->maxAbilities = state->maxDims;
		state->design = omxInitTemporaryMatrix(NULL, state->maxDims, numItems,
				       TRUE, currentState);
		for (int ix=0; ix < numItems; ix++) {
			for (int dx=0; dx < state->maxDims; dx++) {
				omxSetMatrixElement(state->design, dx, ix, (double)dx+1);
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
			if (round(got) != got) error("Design matrix can only contain integers"); // TODO better way?
			if (state->maxAbilities < got)
				state->maxAbilities = got;
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

	PROTECT(tmp = GET_SLOT(rObj, install("cache")));
	state->cacheLXK = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("GHpoints")));
	double *qpoints = REAL(tmp);
	int numQPoints = length(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("GHarea")));
	double *qarea = REAL(tmp);
	if (numQPoints != length(tmp)) error("length(GHpoints) != length(GHarea)");

	ba81SetupQuadrature(oo, numQPoints, qpoints, qarea);

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData for column summary stats?
}
