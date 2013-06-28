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

typedef double *(*rpf_fn_t)(omxExpectation *oo, omxMatrix *itemParam, const double *theta);

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
	const struct rpf *stdrpf;
	int maxOutcomes;
	int maxDims;
	int *numQuadPoints;
	long totalQuadPoints;  // product of numQuadPoints
	int maxAbilities;
	omxMatrix *design;
	omxMatrix *itemPrior;
	omxMatrix *itemParam;
	omxMatrix *EitemParam;
	SEXP rpf;
	rpf_fn_t computeRPF;

	double **lxk;  // per thread vector of numUnique doubles

	double *patternLik; // array(dim=c(o@numPatterns))
	double *logNumIdentical; // array(dim=c(o@numPatterns))
	double ll;

} omxBA81State;

enum ISpecRow {
	ISpecID,
	ISpecOutcomes,
	ISpecDims,
	ISpecRowCount
};

// This acts like a prior distribution for difficulty.
// PC-BILOG does not impose any additional prior
// on difficulty estimates (Baker & Kim, 2004, p. 196).
//
// quadrature points and weights via BILOG
static const int numQuadratures = 10;
static const double quadPoint[] = {
	-4.000, -3.111,-2.222,-1.333,-0.4444,
	0.4444, 1.333, 2.222, 3.111, 4.000
};
static const double quadArea[] = {
	0.000119, 0.002805, 0.03002, 0.1458, 0.3213,
	0.3213, 0.1458, 0.03002, 0.002805, 0.000119
};
static const double quadLogArea[] = {
	-9.036387, -5.876352, -3.505891, -1.925519, -1.135380,
	-1.135380, -1.925519, -3.505891, -5.876352, -9.036387
};

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		Rprintf("Freeing %s function.\n", NAME);
	}
	omxBA81State *mml = (omxBA81State *) oo->argStruct;
	if (mml->data) omxFreeData(mml->data);
	if (mml->itemSpec) omxFreeMatrixData(mml->itemSpec);
	if (mml->itemParam) omxFreeMatrixData(mml->itemParam);
	if (mml->EitemParam) omxFreeMatrixData(mml->EitemParam);
	if (mml->design) omxFreeMatrixData(mml->design);
	if (mml->itemPrior) omxFreeMatrixData(mml->itemPrior);
}

OMXINLINE static double *
ba81Likelihood(omxExpectation *oo, const int *quad)   // TODO optimize more
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	double *restrict lxk = state->lxk[omx_absolute_thread_num()];
	int numUnique = state->numUnique;
	int maxOutcomes = state->maxOutcomes;
	omxData *data = state->data;
	int numItems = state->itemSpec->cols;
	rpf_fn_t rpf_fn = state->computeRPF;

	double where[state->maxDims];
	for (int dx=0; dx < state->maxDims; dx++) {
		where[dx] = quadPoint[quad[dx]];
	}

	const double *outcomeProb = (*rpf_fn)(oo, state->EitemParam, where);
	if (!outcomeProb) {
		OMXZERO(lxk, numUnique);
		return lxk;
	}

	for (int px=0, row=0; px < numUnique; px++) {
		double lxk1 = 0;
		for (int ix=0; ix < numItems; ix++) {
			int pick = omxIntDataElement(data, row, ix);
			if (pick == NA_INTEGER) continue;
			lxk1 += outcomeProb[ix * maxOutcomes + pick-1];
		}
		lxk[px] = lxk1;
		row += omxDataNumIdenticalRows(data, row);
	}

	Free(outcomeProb);

	return lxk;
}

/** 
 * This is the weight used in the likelihood function.
 *
 * \param quad a vector that indexes into a multidimensional quadrature
 * \param out points to an array numOutcomes wide
 */
static void
ba81LogExpected(omxExpectation* oo, const int item, const int *quad, double *out)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxData *data = state->data;
	int numItems = state->itemSpec->cols;
	if (item < 0 || item >= numItems) error("ba81LogExpected: item %d out of range", item);

	// If lxk was cached then it would save 1 evaluation per
	// quadrature per fit.  However, this is only a good tradeoff
	// if there is enough RAM. This is nice potential
	// optimization, but larger estimation problems would not
	// benefit.

	double *lxk = ba81Likelihood(oo, quad);

	double *patternLik = state->patternLik;
	double *logNumIdentical = state->logNumIdentical;
	int numUnique = state->numUnique;
	int outcomes = omxMatrixElement(state->itemSpec, ISpecOutcomes, item);

	OMXZERO(out, outcomes);

	for (int px=0, row=0; px < numUnique; px++) {
		double nt = logNumIdentical[px] + lxk[px] - patternLik[px];
		int pick = omxIntDataElement(data, row, item);
		if (pick == NA_INTEGER) {
			double slice = exp(nt - log(outcomes));
			for (int ox=0; ox < outcomes; ox++) {
				out[ox] += slice;
			}
		} else {
			out[pick-1] += exp(nt);
		}
		row += omxDataNumIdenticalRows(data, row);
	}

	double logArea = 0;
	for (int dx=0; dx < state->maxDims; dx++) {
		logArea += quadLogArea[quad[dx]];
	}

	for (int ox=0; ox < outcomes; ox++) {
		out[ox] = log(out[ox]) + logArea;
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
standardComputeRPF(omxExpectation *oo, omxMatrix *itemParam, const double *theta)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemSpec = state->itemSpec;
	int numItems = itemSpec->cols;
	omxMatrix *design = state->design;

	double *outcomeProb = Realloc(NULL, numItems * state->maxOutcomes, double);

	for (int ix=0; ix < numItems; ix++) {
		int outcomes = omxMatrixElement(itemSpec, ISpecOutcomes, ix);
		double *iparam = omxMatrixColumn(itemParam, ix);
		double *out = outcomeProb + ix * state->maxOutcomes;
		int id = omxMatrixElement(itemSpec, ISpecID, ix);
		int dims = omxMatrixElement(itemSpec, ISpecDims, ix);
		double ptheta[dims];
		for (int dx=0; dx < dims; dx++) {
			ptheta[dx] = theta[((int)omxMatrixElement(design, dx, ix)-1) % dims];
		}
		(*rpf_table[id].logprob)(dims, iparam, ptheta, outcomes, out);
	}

	return outcomeProb;
}

static double *
RComputeRPF1(omxExpectation *oo, omxMatrix *itemParam, const double *theta)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	int maxOutcomes = state->maxOutcomes;
	omxMatrix *design = state->design;
	omxMatrix *itemSpec = state->itemSpec;

	SEXP invoke;
	PROTECT(invoke = allocVector(LANGSXP, 4));
	SETCAR(invoke, state->rpf);
	SETCADR(invoke, omxExportMatrix(itemParam));
	SETCADDR(invoke, omxExportMatrix(itemSpec));

	int dims = state->maxDims;
	SEXP where;
	PROTECT(where = allocMatrix(REALSXP, dims, itemParam->cols));
	double *ptheta = REAL(where);
	for (int ix=0; ix < itemParam->cols; ix++) {
		int idims = omxMatrixElement(itemSpec, ISpecDims, ix);
		for (int dx=0; dx < dims; dx++) {
			double out;
			if (dx < idims) {
				out = theta[((int)omxMatrixElement(design, dx, ix)-1) % dims];
			} else {
				out = NA_REAL;
			}
			ptheta[ix * dims + dx] = out;
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
RComputeRPF(omxExpectation *oo, omxMatrix *itemParam, const double *theta)
{
	omx_omp_set_lock(&GlobalRLock);
	PROTECT_INDEX pi = omxProtectSave();
	double *ret = RComputeRPF1(oo, itemParam, theta);
	omxProtectRestore(pi);
	omx_omp_unset_lock(&GlobalRLock);  // hope there was no exception!
	return ret;
}

OMXINLINE static void
decodeLocation(long qx, const int dims, const int *restrict grid,
	       int *restrict quad)
{
	for (int dx=0; dx < dims; dx++) {
		quad[dx] = qx % grid[dx];
		qx = qx / grid[dx];
	}
}

static void ba81Estep(omxExpectation *oo) {
	if(OMX_DEBUG_MML) {Rprintf("Beginning %s Computation.\n", NAME);}

	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxData *data = state->data;
	double *patternLik = state->patternLik;
	int numUnique = state->numUnique;

	OMXZERO(patternLik, numUnique);
	omxCopyMatrix(state->EitemParam, state->itemParam);

	// E-step, marginalize person ability
	//
	// Note: In Cai (2010) and Bock & Aitkin (1981), these loops
	// are reversed.  That is, the inner loop is over quadrature
	// points and the outer loop is over all response patterns.
	//
#pragma omp parallel for num_threads(oo->currentState->numThreads)
	for (long qx=0; qx < state->totalQuadPoints; qx++) {
		int quad[state->maxDims];
		decodeLocation(qx, state->maxDims, state->numQuadPoints, quad);
		double *lxk = ba81Likelihood(oo, quad);

		// with more threads, it might be better to make this whole loop atomic?
		for (int px=0, row=0; px < numUnique; px++) {
			double logArea = 0;
			for (int dx=0; dx < state->maxDims; dx++) {
				logArea += quadLogArea[quad[dx]];
			}
			double tmp = exp(lxk[px] + logArea);

#pragma omp atomic
			patternLik[px] += tmp;

			row += omxDataNumIdenticalRows(data, row);
		}
	}

	for (int px=0; px < numUnique; px++) {
		patternLik[px] = log(patternLik[px]);
	}
}

static double
ba81ComputeFit1(omxExpectation* oo)
{
	omxBA81State *state = (omxBA81State*) oo->argStruct;
	omxMatrix *itemPrior = state->itemPrior;
	omxMatrix *itemParam = state->itemParam;
	int numItems = itemParam->cols;

	omxRecompute(itemPrior);
	double ll = itemPrior->data[0];
	int maxOutcomes = state->maxOutcomes;
	rpf_fn_t rpf_fn = state->computeRPF;

#pragma omp parallel for num_threads(oo->currentState->numThreads)
	for (long qx=0; qx < state->totalQuadPoints; qx++) {
		int quad[state->maxDims];
		decodeLocation(qx, state->maxDims, state->numQuadPoints, quad);

		double where[state->maxDims];
		for (int dx=0; dx < state->maxDims; dx++) {
			where[dx] = quadPoint[quad[dx]];
		}

		double *outcomeProb = (*rpf_fn)(oo, itemParam, where);
		if (!outcomeProb) continue;

		double thr_ll = 0;

		for (int ix=0; ix < numItems; ix++) {
			int outcomes = omxMatrixElement(state->itemSpec, ISpecOutcomes, ix);
			double out[outcomes];
			ba81LogExpected(oo, ix, quad, out);
			for (int ox=0; ox < outcomes; ox++) {
				double got = exp(out[ox]) * outcomeProb[ix * maxOutcomes + ox];
				thr_ll += got;
			}
		}

		Free(outcomeProb);

		#pragma omp atomic
		ll += thr_ll;
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

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP nextMatrix;
	SEXP rObj = oo->rObj;
	
	if(OMX_DEBUG) {
		Rprintf("Initializing %s.\n", NAME);
	}
	
	omxBA81State *state = (omxBA81State*) R_alloc(1, sizeof(*state));
	state->ll = 1;
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));
	state->data =
		omxNewDataFromMxDataPtr(nextMatrix, currentState);
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
	oo->destructFun = ba81Destroy;
	
	oo->argStruct = (void*) state;

	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we don't waste RAM.
	// The following summary stats would be cheaper to calculate too.

	int numUnique = 0;
	omxData *data = state->data;
	for (int rx=0; rx < data->rows;) {
		rx += omxDataNumIdenticalRows(state->data, rx);
		++numUnique;
	}
	state->numUnique = numUnique;

	state->logNumIdentical =
		(double*) R_alloc(numUnique, sizeof(double));

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

	state->patternLik =
		(double*) R_alloc(numUnique, sizeof(double));

	int numThreads = oo->currentState->numThreads;
	if (numThreads < 1) numThreads = 1;

	state->lxk = (double**) R_alloc(numThreads, sizeof(double*));
	for (int th=0; th < numThreads; th++) {
		state->lxk[th] = (double*) R_alloc(numUnique, sizeof(double));
	}

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
	if (state->maxAbilities > state->maxDims)
		error("Bi-factor style models are not supported yet");

	state->totalQuadPoints = 1;
	state->numQuadPoints = (int*) R_alloc(state->maxDims, sizeof(int));
	for (int dx=0; dx < state->maxDims; dx++) {
		state->numQuadPoints[dx] = numQuadratures;   // make configurable TODO
		state->totalQuadPoints *= numQuadratures;
	}

	// verify data bounded between 1 and numOutcomes TODO
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
