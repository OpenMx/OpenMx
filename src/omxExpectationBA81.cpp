/*
 * Copyright 2012-2017 Joshua Nathaniel Pritikin and contributors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <limits>
#include <Rmath.h>

#include "omxExpectationBA81.h"
#include "glue.h"
#include <libifa-rpf.h>
#include "dmvnorm.h"
#include "omxBuffer.h"
#include "matrix.h"
#include "EnableWarnings.h"

#define USE_EXTERNAL_LIBRPF 1

const struct rpf *Glibrpf_model = NULL;
int Glibrpf_numModels;

void pda(const double *ar, int rows, int cols)
{
	if (rows == 0 || cols == 0) return;
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
	if (rows == 0 || cols == 0) return;
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%d, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

template <typename T>
void BA81LatentFixed<T>::normalizeWeights(class ifaGroup *grp, T extraData,
					  int px, double patternLik1, int thrId)
{
	double weight = grp->rowWeight[px] / patternLik1;
	grp->quad.weightBy(thrId, weight);
}

template <typename T>
void BA81LatentSummary<T>::normalizeWeights(class ifaGroup *grp, T extraData,
					    int px, double patternLik1, int thrId)
{
	double weight = grp->rowWeight[px] / patternLik1;
	grp->quad.weightByAndSummarize(thrId, weight);
}

static void exportLatentDistToOMX(ba81NormalQuad &quad, double *latentDist1, omxMatrix *meanOut, omxMatrix *covOut)
{
	const int maxAbilities = quad.abilities();

	if (meanOut) {
		for (int d1=0; d1 < maxAbilities; d1++) {
			omxSetVectorElement(meanOut, d1, latentDist1[d1]);
		}
	}

	if (covOut) {
		for (int d1=0; d1 < maxAbilities; d1++) {
			int cx = maxAbilities + triangleLoc1(d1);
			for (int d2=0; d2 <= d1; d2++) {
				double cov = latentDist1[cx];
				omxSetMatrixElement(covOut, d1, d2, cov);
				if (d1 != d2) omxSetMatrixElement(covOut, d2, d1, cov);
				++cx;
			}
		}
	}
}

template <typename T>
void BA81LatentSummary<T>::end(class ifaGroup *grp, T extraData)
{
	ba81NormalQuad &quad = grp->quad;
	int dim = quad.abilities();
	int numLatents = dim + triangleLoc1(dim);
	Eigen::ArrayXd latentDist(numLatents);
	quad.EAP(extraData->weightSum, latentDist);
	for (int d1=quad.abilities(); d1 < numLatents; d1++) {
		latentDist[d1] *= extraData->weightSum / (extraData->weightSum - 1.0);
	}
	exportLatentDistToOMX(quad, latentDist.data(), extraData->estLatentMean, extraData->estLatentCov);

	++extraData->ElatentVersion;
}

void ba81AggregateDistributions(std::vector<struct omxExpectation *> &expectation,
				int *version, omxMatrix *meanMat, omxMatrix *covMat)
{
	int allVer = 0;
	for (size_t ex=0; ex < expectation.size(); ++ex) {
		BA81Expect *ba81 = (BA81Expect *) expectation[ex];
		allVer += ba81->ElatentVersion;
	}
	if (*version == allVer) return;
	*version = allVer;

	BA81Expect *exemplar = (BA81Expect *) expectation[0];
	ba81NormalQuad &quad = exemplar->getQuad();
	ba81NormalQuad combined(quad);

	int got = 0;
	for (size_t ex=0; ex < expectation.size(); ++ex) {
		BA81Expect *ba81 = (BA81Expect *) expectation[ex];
		// double weight = 1/ba81->weightSum; ?
		combined.addSummary(ba81->grp.quad);
		++got;
	}
	if (got == 0) return;

	int dim = quad.abilities();
	int numLatents = dim + triangleLoc1(dim);
	Eigen::ArrayXd latentDist(numLatents);
	combined.EAP(got, latentDist);
	for (int d1=quad.abilities(); d1 < numLatents; d1++) {
		latentDist[d1] *= got / (got - 1.0);
	}
	exportLatentDistToOMX(quad, latentDist.data(), meanMat, covMat);
}

template <typename T>
void BA81Estep<T>::begin(ifaGroup *state)
{
	state->quad.allocEstep(Global->numThreads);
}

template <typename T>
void BA81Estep<T>::addRow(class ifaGroup *state, int mpx, int thrId)
{
	state->quad.addToExpected(thrId, mpx);
}

template <typename T>
void BA81Estep<T>::recordTable(class ifaGroup *state)
{
	state->quad.prepExpectedTable();
}

static unsigned getLatentVersion(BA81Expect *state)
{
	unsigned vv = 1;  // to ensure it doesn't match on the first test
	if (state->_latentMeanOut) vv += omxGetMatrixVersion(state->_latentMeanOut);
	if (state->_latentCovOut) vv += omxGetMatrixVersion(state->_latentCovOut);
	return vv;
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
void ba81RefreshQuadrature(omxExpectation* oo)
{
	BA81Expect *state = (BA81Expect *) oo;
	ba81NormalQuad &quad = state->getQuad();

	Eigen::VectorXd mean;
	Eigen::MatrixXd fullCov;
	state->getLatentDistribution(NULL, mean, fullCov);

	if (state->verbose >= 1) {
		mxLog("%s: refresh quadrature", oo->name);
		if (state->verbose >= 2) {
			int dim = mean.rows();
			pda(mean.data(), 1, dim);
			pda(fullCov.data(), dim, dim);
		}
	}

	quad.refresh(mean, fullCov);
}

void BA81Expect::refreshPatternLikelihood(bool hasFreeLatent)
{
	if (hasFreeLatent) {
		BA81Engine<BA81Expect*, BA81LatentSummary, BA81OmitEstep> engine;
		engine.ba81Estep1(&this->grp, this);
	} else {
		BA81Engine<BA81Expect*, BA81LatentFixed, BA81OmitEstep> engine;
		engine.ba81Estep1(&this->grp, this);
	}
}

void BA81Expect::compute(FitContext *fc, const char *what, const char *how)
{
	omxExpectation *oo = this;
	BA81Expect *state = (BA81Expect *) oo;

	if (what) {
		if (strcmp(what, "latentDistribution")==0 && how && strcmp(how, "copy")==0) {
			omxCopyMatrix(state->_latentMeanOut, state->estLatentMean);
			omxCopyMatrix(state->_latentCovOut, state->estLatentCov);

			double sampleSizeAdj = (state->weightSum - 1.0) / state->weightSum;
			int covSize = state->_latentCovOut->rows * state->_latentCovOut->cols;
			for (int cx=0; cx < covSize; ++cx) {
				state->_latentCovOut->data[cx] *= sampleSizeAdj;
			}
			return;
		}

		if (strcmp(what, "scores")==0) {
			state->expectedUsed = true;
			state->type = EXPECTATION_AUGMENTED;
		} else if (strcmp(what, "nothing")==0) {
			state->type = EXPECTATION_OBSERVED;
		} else {
			omxRaiseErrorf("%s: don't know how to predict '%s'",
				       oo->name, what);
		}

		if (state->verbose >= 1) {
			mxLog("%s: predict %s", oo->name, what);
		}
		return;
	}

	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	bool itemClean = state->itemParamVersion == omxGetMatrixVersion(state->itemParam) && latentClean;

	ba81NormalQuad &quad = state->getQuad();

	if (state->verbose >= 1) {
		mxLog("%s: Qinit %d itemClean %d latentClean %d (1=clean) expectedUsed=%d",
		      oo->name, (int)quad.isAllocated(), itemClean, latentClean, state->expectedUsed);
	}

	if (!latentClean) {
		ba81RefreshQuadrature(oo);
		state->latentParamVersion = getLatentVersion(state);
	}

	if (!itemClean) {
		double *param = state->EitemParam? state->EitemParam : state->itemParam->data;
		state->grp.quad.cacheOutcomeProb(param, FALSE);

		bool estep = state->expectedUsed;
		if (estep) {
			if (oo->dynamicDataSource) {
				BA81Engine<BA81Expect*, BA81LatentSummary, BA81Estep> engine;
				engine.ba81Estep1(&state->grp, state);
			} else {
				BA81Engine<BA81Expect*, BA81LatentFixed, BA81Estep> engine;
				engine.ba81Estep1(&state->grp, state);
			}
		} else {
			state->grp.quad.releaseEstep();
			state->refreshPatternLikelihood(oo->dynamicDataSource);
		}
		if (oo->dynamicDataSource && state->verbose >= 2) {
			mxLog("%s: empirical distribution mean and cov:", state->name);
			omxPrint(state->estLatentMean, "mean");
			omxPrint(state->estLatentCov, "cov");
		}
		if (state->verbose >= 1) {
			const int numUnique = state->getNumUnique();
			mxLog("%s: estep<%s, %s> %d/%d rows excluded",
			      oo->name,
			      (estep && oo->dynamicDataSource? "summary":"fixed"),
			      (estep? "estep":"omitEstep"),
			      state->grp.excludedPatterns, numUnique);
		}
	}

	state->itemParamVersion = omxGetMatrixVersion(state->itemParam);
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
void BA81Expect::populateAttr(SEXP robj)
{
	if (!debugInternal) return;

	ba81NormalQuad &quad = getQuad();
	int maxAbilities = quad.abilities();
	const int numUnique = getNumUnique();

	const double LogLargest = LogLargestDouble;
	SEXP Rlik;

	if (grp.patternLik.size() != numUnique) {
		refreshPatternLikelihood(dynamicDataSource);
	}

	Rf_protect(Rlik = Rf_allocVector(REALSXP, numUnique));
	memcpy(REAL(Rlik), grp.patternLik.data(), sizeof(double) * numUnique);
	double *lik_out = REAL(Rlik);
	for (int px=0; px < numUnique; ++px) {
		// Must return value in log units because it may not be representable otherwise
		lik_out[px] = log(lik_out[px]) - LogLargest;
	}

	MxRList dbg;
	dbg.add("patternLikelihood", Rlik);

	if (quad.getEstepTableSize(0)) {
		SEXP Rexpected;
		Rf_protect(Rexpected = Rf_allocVector(REALSXP, quad.getEstepTableSize(0)));
		Eigen::Map< Eigen::ArrayXd > box(REAL(Rexpected), quad.getEstepTableSize(0));
		quad.exportEstepTable(0, box);
		dbg.add("em.expected", Rexpected);
	}

	SEXP Rmean, Rcov;
	if (estLatentMean) {
		Rf_protect(Rmean = Rf_allocVector(REALSXP, maxAbilities));
		memcpy(REAL(Rmean), estLatentMean->data, maxAbilities * sizeof(double));
		dbg.add("mean", Rmean);
	}
	if (estLatentCov) {
		Rf_protect(Rcov = Rf_allocMatrix(REALSXP, maxAbilities, maxAbilities));
		memcpy(REAL(Rcov), estLatentCov->data, maxAbilities * maxAbilities * sizeof(double));
		dbg.add("cov", Rcov);
	}

	Rf_setAttrib(robj, Rf_install("debug"), dbg.asR());
}

BA81Expect::~BA81Expect()
{
	if(OMX_DEBUG) {
		mxLog("Freeing %s function.", name);
	}
	omxFreeMatrix(estLatentMean);
	omxFreeMatrix(estLatentCov);
}

omxMatrix *BA81Expect::getComponent(const char *what)
{
	if (strcmp(what, "covariance")==0) {
		return estLatentCov;
	} else if (strcmp(what, "mean")==0) {
		return estLatentMean;
	} else {
		return NULL;
	}
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    ScopedProtect p1(matrixDims, Rf_getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
}

omxExpectation *omxInitExpectationBA81() { return new BA81Expect; }

void BA81Expect::init() {
	SEXP tmp;
	
	if(OMX_DEBUG) {
		mxLog("Initializing %s.", name);
	}
	if (!Glibrpf_model) {
#if USE_EXTERNAL_LIBRPF
		get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model_GPL");
		(*get_librpf)(LIBIFA_RPF_API_VERSION, &Glibrpf_numModels, &Glibrpf_model);
#else
		// if linking against included source code
		Glibrpf_numModels = librpf_numModels;
		Glibrpf_model = librpf_model;
#endif
	}
	
	BA81Expect *state = this;

	// These two constants should be as identical as possible
	state->name = name;
	if (0) {
		state->LogLargestDouble = 0.0;
		state->LargestDouble = 1.0;
	} else {
		state->LogLargestDouble = log(std::numeric_limits<double>::max()) - 1;
		state->LargestDouble = exp(state->LogLargestDouble);
		ba81NormalQuad &quad = state->getQuad();
		quad.setOne(state->LargestDouble);
	}

	state->expectedUsed = false;

	state->estLatentMean = NULL;
	state->estLatentCov = NULL;
	state->type = EXPECTATION_OBSERVED;
	state->itemParam = NULL;
	state->EitemParam = NULL;
	state->itemParamVersion = 0;
	state->latentParamVersion = 0;

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);
	}

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf("%s unable to handle data type %s", name, omxDataType(state->data));
		return;
	}

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("verbose")));
	state->verbose = Rf_asInteger(tmp);
	}

	int targetQpoints;
	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("qpoints")));
		targetQpoints = Rf_asInteger(tmp);
	}

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("qwidth")));
	state->grp.setGridFineness(Rf_asReal(tmp), targetQpoints);
	}

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("ItemSpec")));
	state->grp.importSpec(tmp);
	if (state->verbose >= 2) mxLog("%s: found %d item specs", name, state->numItems());
	}

	state->_latentMeanOut = omxNewMatrixFromSlot(rObj, currentState, "mean");
	state->_latentCovOut  = omxNewMatrixFromSlot(rObj, currentState, "cov");

	state->itemParam = omxNewMatrixFromSlot(rObj, currentState, "item");
	state->grp.param = state->itemParam->data; // algebra not allowed yet TODO

	const int numItems = state->itemParam->cols;
	if (state->numItems() != numItems) {
		omxRaiseErrorf("ItemSpec length %d must match the number of item columns (%d)",
			       state->numItems(), numItems);
		return;
	}
	if (state->itemParam->rows != state->grp.impliedParamRows) {
		omxRaiseErrorf("item matrix must have %d rows", state->grp.impliedParamRows);
		return;
	}
	state->grp.paramRows = state->itemParam->rows;

	// for algebra item param, will need to defer until later?
	state->grp.learnMaxAbilities();

	int maxAbilities = state->grp.itemDims;
	state->grp.setFactorNames(state->itemParam->rownames);

	{
		ProtectedSEXP tmp2(R_do_slot(rObj, Rf_install(".detectIndependence")));
		state->grp.detectIndependence = Rf_asLogical(tmp2);
	}

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("EstepItem")));
	if (!Rf_isNull(tmp)) {
		int rows, cols;
		getMatrixDims(tmp, &rows, &cols);
		if (rows != state->itemParam->rows || cols != state->itemParam->cols) {
			Rf_error("EstepItem must have the same dimensions as the item MxMatrix");
		}
		state->EitemParam = REAL(tmp);
	}
	}

	canDuplicate = false;
	
	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we get better cache behavior.
	// The following summary stats would be cheaper to calculate too.

	if (data->hasDefinitionVariables()) Rf_error("%s: not implemented yet", name);

	std::vector<int> &rowMap = state->grp.rowMap;

	int weightCol;
	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("weightColumn")));
		weightCol = INTEGER(tmp)[0];
	}

	if (weightCol == NA_INTEGER && !data->hasWeight()) {
		state->grp.rowWeight = (double*) R_alloc(data->rows, sizeof(double));
		for (int rx=0; rx < data->rows; ++rx) {
			state->grp.rowWeight[rx] = 1.0;
		}
		state->weightSum = state->data->rows;
	} else if (data->hasWeight()) {
		if (weightCol != NA_INTEGER) {
			Rf_warning("Data '%s' already has a weight column; "
				   "weight column provided to '%s' ignored", data->name, name);
		}
		state->grp.rowWeight = data->getWeightColumn();
		state->weightSum = omxDataNumObs(data);
	} else if (weightCol != NA_INTEGER) {
		if (omxDataColumnIsFactor(data, weightCol)) {
			omxRaiseErrorf("%s: weightColumn %d is a factor", name, 1 + weightCol);
			return;
		}
		state->grp.rowWeight = omxDoubleDataColumn(data, weightCol);
		state->weightSum = 0;
		for (int rx=0; rx < data->rows; ++rx) { state->weightSum += state->grp.rowWeight[rx]; }
	}

	// complain about non-integral rowWeights (EAP can't work) TODO

	rowMap.resize(data->rows);
	for (size_t rx=0; rx < rowMap.size(); ++rx) {
		rowMap[rx] = rx;
	}

	auto colMap = getDataColumns();

	for (int cx = 0; cx < numItems; cx++) {
		int *col = omxIntDataColumnUnsafe(data, colMap[cx]);
		state->grp.dataColumns.push_back(col);
	}

	// sanity check data
	for (int cx = 0; cx < numItems; cx++) {
		if (!omxDataColumnIsFactor(data, colMap[cx])) {
			data->omxPrintData("diagnostic", 3);
			omxRaiseErrorf("%s: column %d is not a factor", name, int(1 + colMap[cx]));
			return;
		}
	}

	// TODO the max outcome should be available from omxData
	for (int rx=0; rx < data->rows; rx++) {
		int cols = 0;
		for (int cx = 0; cx < numItems; cx++) {
			const int *col = state->grp.dataColumns[cx];
			int pick = col[rx];
			if (pick == NA_INTEGER) continue;
			++cols;
			const int no = state->grp.itemOutcomes[cx];
			if (pick > no) {
				Rf_error("Data for item '%s' has at least %d outcomes, not %d",
					 state->itemParam->colnames[cx], pick, no);
			}
		}
		if (cols == 0) {
			Rf_error("Row %d has all NAs", 1+rx);
		}
	}

	if (state->_latentMeanOut && state->_latentMeanOut->rows * state->_latentMeanOut->cols != maxAbilities) {
		Rf_error("The mean matrix '%s' must be a row or column vector of size %d",
			 state->_latentMeanOut->name(), maxAbilities);
	}

	if (state->_latentCovOut && (state->_latentCovOut->rows != maxAbilities ||
				    state->_latentCovOut->cols != maxAbilities)) {
		Rf_error("The cov matrix '%s' must be %dx%d",
			 state->_latentCovOut->name(), maxAbilities, maxAbilities);
	}

	state->grp.setLatentDistribution(state->_latentMeanOut? state->_latentMeanOut->data : NULL,
					 state->_latentCovOut? state->_latentCovOut->data : NULL);

	{
		EigenArrayAdaptor Eparam(state->itemParam);
		Eigen::Map< Eigen::VectorXd > meanVec(state->grp.mean, maxAbilities);
		Eigen::Map< Eigen::MatrixXd > covMat(state->grp.cov, maxAbilities, maxAbilities);
		state->grp.quad.setStructure(state->grp.qwidth, state->grp.qpoints,
					     Eparam, meanVec, covMat);
	}

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("minItemsPerScore")));
	state->grp.setMinItemsPerScore(Rf_asInteger(tmp));
	}

	state->grp.buildRowSkip();

	if (isErrorRaised()) return;

	{ScopedProtect p1(tmp, R_do_slot(rObj, Rf_install("debugInternal")));
	state->debugInternal = Rf_asLogical(tmp);
	}

	state->ElatentVersion = 0;
	if (state->_latentMeanOut) {
		state->estLatentMean = omxInitMatrix(maxAbilities, 1, TRUE, currentState);
		omxCopyMatrix(state->estLatentMean, state->_latentMeanOut); // rename matrices TODO
	}
	if (state->_latentCovOut) {
		state->estLatentCov = omxInitMatrix(maxAbilities, maxAbilities, TRUE, currentState);
		omxCopyMatrix(state->estLatentCov, state->_latentCovOut);
	}
}

void BA81Expect::invalidateCache()
{
	grp.rowWeight = data->getWeightColumn();
}

const char *BA81Expect::getLatentIncompatible(BA81Expect *other)
{
	// NOTE: grp.quad not initialized yet
	// make method of ifaGroup ?
	if (grp.itemOutcomes != other->grp.itemOutcomes) return "items";
	if (grp.itemDims != other->grp.itemDims) return "number of factors";
	if (grp.qpoints != other->grp.qpoints) return "qpoints";
	if (grp.qwidth != other->grp.qwidth) return "qwidth";
	return 0;
}

