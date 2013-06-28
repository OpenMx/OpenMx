/*
  Copyright 2012 Joshua Nathaniel Pritikin and contributors

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

#include "omxExpectation.h"
#include "omxOpenmpWrap.h"
#include "npsolWrap.h"
#include "libirt-rpf.h"

static const char *NAME = "ExpectationBA81";

typedef double *(*rpf_fn_t)(omxExpectation *oo, omxMatrix *itemParam, const int *quad);

typedef int (*rpf_numParam_t)(const int numDims, const int numOutcomes);
typedef void (*rpf_logprob_t)(const int numDims, const double *restrict param,
			      const double *restrict th,
			      const int numOutcomes, double *restrict out);
struct rpf {
	const char name[8];
	rpf_numParam_t numParam;
	rpf_logprob_t logprob;
};

static const struct rpf rpf_table[] = {
	{ "drm1",  irt_rpf_1dim_drm_numParam,  irt_rpf_1dim_drm_logprob },
	{ "drm",   irt_rpf_mdim_drm_numParam,  irt_rpf_mdim_drm_logprob },
	{ "gpcm1", irt_rpf_1dim_gpcm_numParam, irt_rpf_1dim_gpcm_logprob }
};
static const int numStandardRPF = (sizeof(rpf_table) / sizeof(struct rpf));

typedef struct {

	omxData *data;
	int numUnique;
	omxMatrix *itemSpec;
	int *Sgroup;              // item's specific group 0..numSpecific-1
	int maxOutcomes;
	int maxDims;
	int numGHpoints;
	double *GHpoint;
	double *GHarea;
	long *quadGridSize;       // maxDims
	long totalPrimaryPoints;  // product of quadGridSize except specific dim
	long totalQuadPoints;     // product of quadGridSize
	int maxAbilities;
	int numSpecific;
	omxMatrix *design;        // items * maxDims
	omxMatrix *itemPrior;
	omxMatrix *itemParam;     // M step version
	omxMatrix *EitemParam;    // E step version
	SEXP rpf;
	rpf_fn_t computeRPF;

	int cacheLXK;		  // w/cache,  numUnique * #specific quad points * totalQuadPoints
	double *lxk;              // wo/cache, numUnique * thread
	double *allSlxk;          // numUnique * thread
	double *Slxk;             // numUnique * #specific dimensions * thread

	double *patternLik;       // length numUnique
	double *logNumIdentical;  // length numUnique
	double ll;                // the most recent finite ll

} omxBA81State;

enum ISpecRow {
	ISpecID,
	ISpecOutcomes,
	ISpecDims,
	ISpecRowCount
};

/*
static void
pda(const double *ar, int rows, int cols) {
	for (int rx=0; rx < rows; rx++) {
		for (int cx=0; cx < cols; cx++) {
			Rprintf("%.6g ", ar[cx * rows + rx]);
		}
		Rprintf("\n");
	}

}
*/

OMXINLINE static void
pointToWhere(omxBA81State *state, const int *quad, double *where, int upto)
{
	for (int dx=0; dx < upto; dx++) {
		where[dx] = state->GHpoint[quad[dx]];
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
		int outcomes = omxMatrixElement(itemSpec, ISpecOutcomes, ix);
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * state->maxOutcomes;
		int id = omxMatrixElement(itemSpec, ISpecID, ix);
		int dims = omxMatrixElement(itemSpec, ISpecDims, ix);
		double ptheta[dims];
		assignDims(itemSpec, design, dims, maxDims, ix, theta, ptheta);
		(*rpf_table[id].logprob)(dims, iparam, ptheta, outcomes, out);
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
		int dims = omxMatrixElement(itemSpec, ISpecDims, ix);
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
		int numOutcomes = omxMatrixElement(state->itemSpec, ISpecOutcomes, ix);
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

	for (int px=0, row=0; px < numUnique; px++) {
		double lxk1 = 0;
		for (int ix=0; ix < numItems; ix++) {
			if (specific != Sgroup[ix]) continue;
			int pick = omxIntDataElementUnsafe(data, row, ix);
			if (pick == NA_INTEGER) continue;
			lxk1 += outcomeProb[ix * maxOutcomes + pick-1];
		}
		lxk[px] = lxk1;
		row += omxDataNumIdenticalRows(data, row);
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
				eis[ix] += exp(lxk[ix] + state->GHarea[qx]);
			}
		}

		for (int px=0; px < numUnique; px++) {
			eis[px] = log(eis[px]);
			allSlxk[px] += eis[px];
		}
	}
}

OMXINLINE static double
areaProduct(omxBA81State *state, const int *restrict quad, const int upto)
{
	double logArea = 0;
	for (int dx=0; dx < upto; dx++) {
		logArea += state->GHarea[quad[dx]];
	}
	return logArea;
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

			double logArea = areaProduct(state, quad, state->maxDims);
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

			double logArea = areaProduct(state, quad, sDim);
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
expectedUpdate(omxData *restrict data, int *restrict row, const int item,
	       const double observed, const int outcomes, double *out)
{
	int pick = omxIntDataElementUnsafe(data, *row, item);
	if (pick == NA_INTEGER) {
		double slice = exp(observed - log(outcomes));
		for (int ox=0; ox < outcomes; ox++) {
			out[ox] += slice;
		}
	} else {
		out[pick-1] += exp(observed);
	}
	*row += omxDataNumIdenticalRows(data, *row);
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
	int specific = state->Sgroup[item];
	double *patternLik = state->patternLik;
	double *logNumIdentical = state->logNumIdentical;
	int numUnique = state->numUnique;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;
	int sDim = state->maxDims-1;

	OMXZERO(out, outcomes);

	if (numSpecific == 0) {
		double *lxk = ba81LikelihoodFast(oo, specific, quad);
		for (int px=0, row=0; px < numUnique; px++) {
			double observed = logNumIdentical[px] + lxk[px] - patternLik[px];
			expectedUpdate(data, &row, item, observed, outcomes, out);
		}
	} else {
		double *allSlxk = CALC_ALLSLXK(state, numUnique);
		double *Slxk = CALC_SLXK(state, numUnique, numSpecific);
		if (quad[sDim] == 0) {
			// allSlxk, Slxk only depend on the ordinate of the primary dimensions
			cai2010(oo, !state->cacheLXK, quad, allSlxk, Slxk);
		}
		double *eis = Slxk + numUnique * specific;
		double *lxk = ba81LikelihoodFast(oo, specific, quad);

		for (int px=0, row=0; px < numUnique; px++) {
			double observed = logNumIdentical[px] + (allSlxk[px] - eis[px]) +
				(lxk[px] - patternLik[px]);
			expectedUpdate(data, &row, item, observed, outcomes, out);
		}
	}

	double logArea = areaProduct(state, quad, maxDims);

	for (int ox=0; ox < outcomes; ox++) {
		out[ox] = log(out[ox]) + logArea;
	}
}

OMXINLINE static double
ba81Fit1Ordinate(omxExpectation* oo, const int *quad)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemParam = state->itemParam;
	int numItems = itemParam->cols;
	rpf_fn_t rpf_fn = state->computeRPF;
	int maxOutcomes = state->maxOutcomes;

	double *outcomeProb = (*rpf_fn)(oo, itemParam, quad);
	if (!outcomeProb) return 0;

	double thr_ll = 0;
	for (int ix=0; ix < numItems; ix++) {
		int outcomes = omxMatrixElement(state->itemSpec, ISpecOutcomes, ix);
		double out[outcomes];
		ba81Weight(oo, ix, quad, outcomes, out);
		for (int ox=0; ox < outcomes; ox++) {
			double got = exp(out[ox]) * outcomeProb[ix * maxOutcomes + ox];
			thr_ll += got;
		}
	}

	Free(outcomeProb);
	return thr_ll;
}

static double
ba81ComputeFit1(omxExpectation* oo)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemPrior = state->itemPrior;
	int numSpecific = state->numSpecific;
	int maxDims = state->maxDims;

	omxRecompute(itemPrior);
	double ll = itemPrior->data[0];

	if (numSpecific == 0) {
#pragma omp parallel for num_threads(oo->currentState->numThreads)
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);

			double thr_ll = ba81Fit1Ordinate(oo, quad);

#pragma omp atomic
			ll += thr_ll;
		}
	} else {
		int sDim = state->maxDims-1;
		long *quadGridSize = state->quadGridSize;

#pragma omp parallel for num_threads(oo->currentState->numThreads)
		for (long qx=0; qx < state->totalPrimaryPoints; qx++) {
			int quad[maxDims];
			decodeLocation(qx, maxDims, state->quadGridSize, quad);

			double thr_ll = 0;
			long specificPoints = quadGridSize[sDim];
			for (long sx=0; sx < specificPoints; sx++) {
				quad[sDim] = sx;
				thr_ll += ba81Fit1Ordinate(oo, quad);
			}
#pragma omp atomic
			ll += thr_ll;
		}
	}

	if (isinf(ll)) {
		return 2*state->ll;
	} else {
		// TODO need to *2 also?
		state->ll = -ll;
		return -ll;
	}
}

double
ba81ComputeFit(omxExpectation* oo)
{
	double got = ba81ComputeFit1(oo);
	return got;
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		Rprintf("Freeing %s function.\n", NAME);
	}
	omxBA81State *state = (omxBA81State *) oo->argStruct;
	omxFreeAllMatrixData(state->itemSpec);
	omxFreeAllMatrixData(state->itemParam);
	omxFreeAllMatrixData(state->EitemParam);
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->itemPrior);
	Free(state->logNumIdentical);
	Free(state->patternLik);
	Free(state->lxk);
	Free(state->Slxk);
	Free(state->allSlxk);
	Free(state->Sgroup);
	Free(state);
}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		Rprintf("Initializing %s.\n", NAME);
	}
	
	omxBA81State *state = Calloc(1, omxBA81State);
	state->ll = 10^9;   // finite but big
	
	PROTECT(tmp = GET_SLOT(rObj, install("GHpoints")));
	state->numGHpoints = length(tmp);
	state->GHpoint = REAL(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("GHarea")));
	if (state->numGHpoints != length(tmp)) error("length(GHpoints) != length(GHarea)");
	state->GHarea = REAL(tmp);

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
		// and analytic gradient TODO
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
	state->itemPrior =
		omxNewMatrixFromIndexSlot(rObj, currentState, "ItemPrior");
	
	oo->computeFun = ba81Estep;
	//	oo->gradientFun = ba81Gradient;
	oo->destructFun = ba81Destroy;
	
	oo->argStruct = (void*) state;

	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we don't waste RAM.
	// The following summary stats would be cheaper to calculate too.

	int numUnique = 0;
	omxData *data = state->data;
	if (omxDataNumFactor(data) != data->cols) {
		omxRaiseErrorf(currentState, "%s: all columns must be factors", NAME);
		return;
	}

	for (int rx=0; rx < data->rows;) {
		rx += omxDataNumIdenticalRows(state->data, rx);
		++numUnique;
	}
	state->numUnique = numUnique;

	state->logNumIdentical = Realloc(NULL, numUnique, double);

	int numItems = state->itemParam->cols;

	for (int rx=0, ux=0; rx < data->rows;) {
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
		state->logNumIdentical[ux++] = log(dups);
		rx += dups;
	}

	state->patternLik = Realloc(NULL, numUnique, double);

	int numThreads = oo->currentState->numThreads;
	if (numThreads < 1) numThreads = 1;

	if (state->itemSpec->cols != data->cols || state->itemSpec->rows != ISpecRowCount) {
		omxRaiseErrorf(currentState, "ItemSpec must have %d item columns and %d rows",
			       data->cols, ISpecRowCount);
		return;
	}

	int maxParam = 0;
	state->maxDims = 0;
	state->maxOutcomes = 0;

	for (int cx = 0; cx < data->cols; cx++) {
		int id = omxMatrixElement(state->itemSpec, ISpecID, cx);
		if (id < 0 || id >= numStandardRPF) {
			omxRaiseErrorf(currentState, "ItemSpec column %d has unknown item model %d", cx, id);
			return;
		}

		int dims = omxMatrixElement(state->itemSpec, ISpecDims, cx);
		if (state->maxDims < dims)
			state->maxDims = dims;

		// TODO verify that item model can have requested number of outcomes
		int no = omxMatrixElement(state->itemSpec, ISpecOutcomes, cx);
		if (state->maxOutcomes < no)
			state->maxOutcomes = no;

		int numParam = (*rpf_table[id].numParam)(dims, no);
		if (maxParam < numParam)
			maxParam = numParam;
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
		int Sgroup0 = state->maxDims;
		state->Sgroup = Realloc(NULL, numItems, int);
		for (int ix=0; ix < numItems; ix++) {
			int ss=-1;
			for (int dx=0; dx < state->maxDims; dx++) {
				int ability = omxMatrixElement(state->design, dx, ix);
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
			}
			if (ss == -1) ss = 0;
			state->Sgroup[ix] = ss - Sgroup0;
		}
		state->numSpecific = state->maxAbilities - state->maxDims + 1;
		state->allSlxk = Realloc(NULL, numUnique * numThreads, double);
		state->Slxk = Realloc(NULL, numUnique * state->numSpecific * numThreads, double);
	}

	state->totalQuadPoints = 1;
	state->totalPrimaryPoints = 1;
	state->quadGridSize = (long*) R_alloc(state->maxDims, sizeof(long));
	for (int dx=0; dx < state->maxDims; dx++) {
		state->quadGridSize[dx] = state->numGHpoints;
		state->totalQuadPoints *= state->quadGridSize[dx];
		if (dx < state->maxDims-1) {
			state->totalPrimaryPoints *= state->quadGridSize[dx];
		}
	}

	PROTECT(tmp = GET_SLOT(rObj, install("cache")));
	state->cacheLXK = asLogical(tmp);

	if (!state->cacheLXK) {
		state->lxk = Realloc(NULL, numUnique * numThreads, double);
	} else {
		int ns = state->numSpecific;
		if (ns == 0) ns = 1;
		state->lxk = Realloc(NULL, numUnique * state->totalQuadPoints * ns, double);
	}

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData?
}

SEXP omx_get_rpf_names()
{
	SEXP outsxp;
	PROTECT(outsxp = allocVector(STRSXP, numStandardRPF));
	for (int sx=0; sx < numStandardRPF; sx++) {
		SET_STRING_ELT(outsxp, sx, mkChar(rpf_table[sx].name));
	}
	UNPROTECT(1);
	return outsxp;
}
