/*
 *  Copyright 2013 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxDefines.h"
#include "Compute.h"
#include "omxState.h"
#include "omxExportBackendState.h"
#include "omxRFitFunction.h"
#include "matrix.h"
#include "omxBuffer.h"

void FitContext::init()
{
	size_t numParam = varGroup->vars.size();
	wanted = 0;
	mac = parent? parent->mac : 0;
	fit = parent? parent->fit : 0;
	caution = parent? parent->caution : 0;
	est = new double[numParam];
	flavor = new int[numParam];
	grad = new double[numParam];
	hess = new double[numParam * numParam];
	ihess = new double[numParam * numParam];
	changedEstimates = false;
	inform = INFORM_UNINITIALIZED;
	iterations = 0;
}

FitContext::FitContext(std::vector<double> &startingValues)
{
	parent = NULL;
	varGroup = Global->freeGroup[0];
	init();

	size_t numParam = varGroup->vars.size();
	if (startingValues.size() != numParam) {
		error("Got %d starting values for %d parameters",
		      startingValues.size(), numParam);
	}
	memcpy(est, startingValues.data(), sizeof(double) * numParam);

	for (size_t v1=0; v1 < numParam; v1++) {
		grad[v1] = nan("unset");
		for (size_t v2=0; v2 < numParam; v2++) {
			hess[v1 * numParam + v2] = nan("unset");
		}
	}
}

FitContext::FitContext(FitContext *parent, FreeVarGroup *varGroup)
{
	this->parent = parent;
	this->varGroup = varGroup;
	init();

	FreeVarGroup *src = parent->varGroup;
	FreeVarGroup *dest = varGroup;
	size_t svars = parent->varGroup->vars.size();
	size_t dvars = varGroup->vars.size();
	if (dvars == 0) return;
	mapToParent.resize(dvars);

	size_t d1 = 0;
	for (size_t s1=0; s1 < src->vars.size(); ++s1) {
		if (src->vars[s1] != dest->vars[d1]) continue;
		mapToParent[d1] = s1;
		est[d1] = parent->est[s1];

		if (parent->wanted & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_INFO)) {
			grad[d1] = parent->grad[s1];

			size_t d2 = 0;
			for (size_t s2=0; s2 < src->vars.size(); ++s2) {
				if (src->vars[s2] != dest->vars[d2]) continue;
				hess[d1 * dvars + d2] = parent->hess[s1 * svars + s2];
				if (++d2 == dvars) break;
			}
		}

		// ihess TODO?

		if (++d1 == dvars) break;
	}
	if (d1 != dvars) error("Parent free parameter group is not a superset");

	// pda(parent->est, 1, svars);
	// pda(est, 1, dvars);
	// pda(parent->grad, 1, svars);
	// pda(grad, 1, dvars);
	// pda(parent->hess, svars, svars);
	// pda(hess, dvars, dvars);
}

void FitContext::copyParamToModel(omxMatrix *mat)
{ copyParamToModel(mat->currentState); }

void FitContext::copyParamToModel(omxMatrix *mat, double *at)
{ copyParamToModel(mat->currentState, at); }

void FitContext::updateParent()
{
	FreeVarGroup *src = varGroup;
	FreeVarGroup *dest = parent->varGroup;
	size_t svars = varGroup->vars.size();
	size_t dvars = parent->varGroup->vars.size();

	parent->wanted |= wanted;
	parent->fit = fit;
	parent->mac = mac;
	parent->caution = caution;

	// rewrite using mapToParent TODO

	if (svars > 0) {
		size_t s1 = 0;
		for (size_t d1=0; d1 < dest->vars.size(); ++d1) {
			if (dest->vars[d1] != src->vars[s1]) continue;
			parent->est[d1] = est[s1];

			if (wanted & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_INFO)) {
				parent->grad[d1] = grad[s1];

				size_t s2 = 0;
				for (size_t d2=0; d2 < dest->vars.size(); ++d2) {
					if (dest->vars[d2] != src->vars[s2]) continue;
					parent->hess[d1 * dvars + d2] = hess[s1 * svars + s2];
					if (++s2 == svars) break;
				}
			}

			// ihess TODO?

			if (++s1 == svars) break;
		}
		if (wanted & FF_COMPUTE_PARAMFLAVOR) {
			for (size_t s1=0; s1 < src->vars.size(); ++s1) {
				parent->flavor[mapToParent[s1]] = flavor[s1];
			}
		}
	}
	
	// pda(est, 1, svars);
	// pda(parent->est, 1, dvars);
	// pda(grad, 1, svars);
	// pda(parent->grad, 1, dvars);
	// pda(hess, svars, svars);
	// pda(parent->hess, dvars, dvars);
}

void FitContext::updateParentAndFree()
{
	updateParent();
	delete this;
}

void FitContext::log(const char *where)
{
	log(where, wanted);
}

void FitContext::log(const char *where, int what)
{
	size_t count = varGroup->vars.size();
	std::string buf(where);
	buf += " ---\n";
	if (what & FF_COMPUTE_MAXABSCHANGE) buf += string_snprintf("MAC: %.5f\n", mac);
	if (what & FF_COMPUTE_FIT) buf += string_snprintf("fit: %.5f\n", fit);
	if (what & FF_COMPUTE_ESTIMATE) {
		buf += string_snprintf("est %lu: c(", count);
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.5f", est[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_GRADIENT) {
		buf += string_snprintf("grad %lu: c(", count);
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.5f", grad[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & (FF_COMPUTE_HESSIAN | FF_COMPUTE_INFO)) {
		buf += string_snprintf("hess %lux%lu: c(", count, count);
		for (size_t v1=0; v1 < count; ++v1) {
			for (size_t v2=0; v2 < count; ++v2) {
				buf += string_snprintf("%.5f", hess[v1 * count + v2]);
				if (v1 < count-1 || v2 < count-1) buf += ", ";
			}
			buf += "\n";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_IHESSIAN) {
		buf += string_snprintf("ihess %lux%lu: c(", count, count);
		for (size_t v1=0; v1 < count; ++v1) {
			for (size_t v2=0; v2 < count; ++v2) {
				buf += string_snprintf("%.5f", ihess[v1 * count + v2]);
				if (v1 < count-1 || v2 < count-1) buf += ", ";
			}
			buf += "\n";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_HGPROD) {
		buf += string_snprintf("ihess %%*%% grad %lu: list(", hgProd.size());
		for (size_t px=0; px < hgProd.size(); ++px) {
			buf += string_snprintf("c(%d, %d, %d)", hgProd[px].hentry,
					       hgProd[px].gentry, hgProd[px].dest);
			if (px < hgProd.size() - 1) buf += ", ";
		}
		buf += ")\n";
	}
	mxLogBig(buf);
}

static void _fixSymmetry(const char *name, double *mat, size_t numParam, bool force)
{
	for (size_t h1=1; h1 < numParam; h1++) {
		for (size_t h2=0; h2 < h1; h2++) {
			if (!force && mat[h2 * numParam + h1] != 0) {
				omxRaiseErrorf(globalState, "%s is not upper triangular", name);
				break;
			}
			mat[h2 * numParam + h1] = mat[h1 * numParam + h2];
		}
	}
}

void FitContext::fixHessianSymmetry(int want, bool force)
{
	size_t numParam = varGroup->vars.size();

	if (want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_INFO)) {
		_fixSymmetry("Hessian/information", hess, numParam, force);
	}

	if (want & FF_COMPUTE_IHESSIAN) {
		_fixSymmetry("Inverse Hessian", ihess, numParam, force);
	}
}

static void omxRepopulateRFitFunction(omxFitFunction* oo, double* x, int n)
{
	omxRFitFunction* rFitFunction = (omxRFitFunction*)oo->argStruct;

	SEXP theCall, estimate;

	PROTECT(estimate = allocVector(REALSXP, n));
	double *est = REAL(estimate);
	for(int i = 0; i < n ; i++) {
		est[i] = x[i];
	}

	PROTECT(theCall = allocVector(LANGSXP, 4));

	SETCAR(theCall, install("imxUpdateModelValues"));
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->flatModel);
	SETCADDDR(theCall, estimate);

	REPROTECT(rFitFunction->model = eval(theCall, R_GlobalEnv), rFitFunction->modelIndex);

	UNPROTECT(2); // theCall, estimate
}

void FitContext::copyParamToModel(omxState* os)
{
	copyParamToModel(os, est);
}

void FitContext::maybeCopyParamToModel(omxState* os)
{
	if (changedEstimates) {
		copyParamToModel(os, est);
		changedEstimates = false;
	}
}

void FitContext::copyParamToModel(omxState* os, double *at)
{
	size_t numParam = varGroup->vars.size();
	if(OMX_DEBUG) {
		mxLog("Copying %lu free parameter estimates to model %p", numParam, os);
	}

	if(numParam == 0) return;

	// Confidence Intervals & Hessian Calculation probe the parameter space
	// near the best estimate. If stale, we need to restore the best estimate
	// before returning results to the user.
	os->stale = at != est;

	os->computeCount++;

	if(OMX_VERBOSE) {
		std::string buf;
		buf += string_snprintf("Call: %d.%d (%ld) ", os->majorIteration, os->minorIteration, os->computeCount);
		buf += ("Estimates: [");
		for(size_t k = 0; k < numParam; k++) {
			buf += string_snprintf(" %f", at[k]);
		}
		buf += ("]\n");
		mxLogBig(buf);
	}

	for(size_t k = 0; k < numParam; k++) {
		omxFreeVar* freeVar = varGroup->vars[k];
		for(size_t l = 0; l < freeVar->locations.size(); l++) {
			omxFreeVarLocation *loc = &freeVar->locations[l];
			omxMatrix *matrix = os->matrixList[loc->matrix];
			int row = loc->row;
			int col = loc->col;
			omxSetMatrixElement(matrix, row, col, at[k]);
			if(OMX_DEBUG) {
				mxLog("Setting location (%d, %d) of matrix %d to value %f for var %lu",
					row, col, loc->matrix, at[k], k);
			}
		}
	}

	if (RFitFunction) omxRepopulateRFitFunction(RFitFunction, at, numParam);

	varGroup->markDirty(os);

	if (!os->childList) return;

	for(int i = 0; i < Global->numChildren; i++) {
		copyParamToModel(os->childList[i], at);
	}
}

double *FitContext::take(int want)
{
	if (!(want & (wanted | FF_COMPUTE_ESTIMATE))) {
		error("Attempt to take %d but not available", want);
	}

	double *ret = NULL;
	switch(want) {
	case FF_COMPUTE_ESTIMATE:
		ret = est;
		est = NULL;
		break;
	case FF_COMPUTE_HESSIAN:
	case FF_COMPUTE_INFO:
		ret = hess;
		hess = NULL;
		break;
	case FF_COMPUTE_IHESSIAN:
		ret = ihess;
		ihess = NULL;
		break;
	default:
		error("Taking of %d is not implemented", want);
	}
	if (!ret) error("Attempt to take %d, already taken", want);
	return ret;
}

FitContext::~FitContext()
{
	if (est) delete [] est;
	if (flavor) delete [] flavor;
	if (grad) delete [] grad;
	if (hess) delete [] hess;
	if (ihess) delete [] ihess;
}

omxFitFunction *FitContext::RFitFunction = NULL;

void FitContext::setRFitFunction(omxFitFunction *rff)
{
	if (rff) {
		Global->numThreads = 1;
		if (RFitFunction) {
			error("You can only create 1 MxRFitFunction per independent model");
		}
	}
	RFitFunction = rff;
}

Ramsay1975::Ramsay1975(FitContext *fc, int flavor, double caution, int verbose,
		       double minCaution)
{
	this->fc = fc;
	this->flavor = flavor;
	this->verbose = verbose;
	this->caution = caution;
	this->minCaution = minCaution;
	maxCaution = 0.0;
	highWatermark = std::max(0.5, caution);  // arbitrary guess

	numParam = fc->varGroup->vars.size();
	prevAdj1.assign(numParam, 0);
	prevAdj2.resize(numParam);
	prevEst.resize(numParam);
	memcpy(prevEst.data(), fc->est, sizeof(double) * numParam);
}

void Ramsay1975::recordEstimate(int px, double newEst)
{
	omxFreeVar *fv = fc->varGroup->vars[px];
	bool hitBound=false;
	double param = newEst;
	if (param < fv->lbound) {
		hitBound=true;
		param = prevEst[px] - (prevEst[px] - fv->lbound) / 2;
	}
	if (param > fv->ubound) {
		hitBound=true;
		param = prevEst[px] + (fv->ubound - prevEst[px]) / 2;
	}
	
	prevAdj2[px] = prevAdj1[px];
	prevAdj1[px] = param - prevEst[px];
	
	if (verbose >= 4) {
		std::string buf;
		buf += string_snprintf("~%d~%s: %.4f -> %.4f", px, fv->name, prevEst[px], param);
		if (hitBound) {
			buf += string_snprintf(" wanted %.4f but hit bound", newEst);
		}
		if (prevAdj1[px] * prevAdj2[px] < 0) {
			buf += " *OSC*";
		}
		buf += "\n";
		mxLogBig(buf);
	}

	fc->est[px] = param;
	prevEst[px] = param;
}

void Ramsay1975::apply()
{
	for (size_t px=0; px < numParam; ++px) {
		recordEstimate(px, (1 - caution) * fc->est[px] + caution * prevEst[px]);
	}
}

void Ramsay1975::recalibrate(bool *restart)
{
	double normPrevAdj2 = 0;
	double normAdjDiff = 0;
	std::vector<double> adjDiff(numParam);

	// The choice of norm is also arbitrary. Other norms might work better.
	for (size_t px=0; px < numParam; ++px) {
		if (fc->flavor[px] != flavor) continue;
		adjDiff[px] = prevAdj1[px] - prevAdj2[px];
		normPrevAdj2 += prevAdj2[px] * prevAdj2[px];
	}

	for (size_t px=0; px < numParam; ++px) {
		if (fc->flavor[px] != flavor) continue;
		normAdjDiff += adjDiff[px] * adjDiff[px];
	}
	if (normAdjDiff == 0) {
		return;
		//error("Ramsay: no free variables of flavor %d", flavor);
	}

	double ratio = sqrt(normPrevAdj2 / normAdjDiff);
	//if (verbose >= 3) mxLog("Ramsay[%d]: sqrt(%.5f/%.5f) = %.5f",
	// flavor, normPrevAdj2, normAdjDiff, ratio);

	double newCaution = 1 - (1-caution) * ratio;
	if (newCaution > .95) newCaution = .95;  // arbitrary guess
	if (newCaution < 0) newCaution /= 2;     // don't get overconfident
	if (newCaution < minCaution) newCaution = minCaution;
	if (newCaution < caution) {
		caution = newCaution/3 + 2*caution/3;  // don't speed up too fast, arbitrary ratio
	} else {
		caution = newCaution;
	}
	maxCaution = std::max(maxCaution, caution);
	if (caution < highWatermark || (normPrevAdj2 < 1e-3 && normAdjDiff < 1e-3)) {
		if (verbose >= 3) mxLog("Ramsay[%d]: %.2f caution", flavor, caution);
	} else {
		if (verbose >= 3) mxLog("Ramsay[%d]: caution %.2f > %.2f, extreme oscillation, restart recommended",
					flavor, caution, highWatermark);
		*restart = TRUE;
	}
	highWatermark += .02; // arbitrary guess
}

void Ramsay1975::restart()
{
	memcpy(prevEst.data(), fc->est, sizeof(double) * numParam);
	prevAdj1.assign(numParam, 0);
	prevAdj2.assign(numParam, 0);
	highWatermark = 1 - (1 - highWatermark) * .5; // arbitrary guess
	caution = std::max(caution, highWatermark);   // arbitrary guess
	maxCaution = std::max(maxCaution, caution);
	highWatermark = caution;
	if (verbose >= 3) {
		mxLog("Ramsay[%d]: restart with %.2f caution %.2f highWatermark",
		      flavor, caution, highWatermark);
	}
}

omxCompute::omxCompute()
{
	varGroup = NULL;
}

omxCompute::~omxCompute()
{}

void omxCompute::initFromFrontend(SEXP rObj)
{
	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("id")));
	if (length(slotValue) == 1) {
		int id = INTEGER(slotValue)[0];
		varGroup = Global->findVarGroup(id);
	}

	if (!varGroup) {
		if (!R_has_slot(rObj, install("free.set"))) {
			varGroup = Global->freeGroup[0];
		} else {
			PROTECT(slotValue = GET_SLOT(rObj, install("free.set")));
			if (length(slotValue) != 0) {
				// it's a free.set with no free variables
				varGroup = Global->findVarGroup(-1);
			} else {
				varGroup = Global->freeGroup[0];
			}
		}
	}
}

class omxComputeSequence : public omxCompute {
	typedef omxCompute super;
	std::vector< omxCompute* > clist;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~omxComputeSequence();
};

class omxComputeIterate : public omxCompute {
	typedef omxCompute super;
	std::vector< omxCompute* > clist;
	int maxIter;
	double tolerance;
	int verbose;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~omxComputeIterate();
};

class omxComputeOnce : public omxCompute {
	typedef omxCompute super;
	std::vector< omxMatrix* > algebras;
	std::vector< omxExpectation* > expectations;
	int verbose;
	const char *context;
	bool mac;
	bool fit;
	bool gradient;
	bool hessian;
	bool ihessian;
	bool infoMat;
	bool hgprod;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual omxFitFunction *getFitFunction();
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
};

class ComputeEM : public omxCompute {
	typedef omxCompute super;
	std::vector< omxExpectation* > expectations;
	omxCompute *fit1;
	omxCompute *fit2;
	int maxIter;
	int mstepIter;
	int totalMstepIter;
	double tolerance;
	int verbose;
	bool useRamsay;
	bool information;
	std::vector<Ramsay1975*> ramsay;
	std::vector<double*> estHistory;
	FitContext *recentFC;  //nice if can use std::unique_ptr
	std::vector<double> optimum;
	std::vector<double> stdError;
	double bestFit;
 	static const double MIDDLE_START = 0.21072103131565256273; // -log(.9)*2 constexpr
	static const double MIDDLE_END = 0.0020010006671670687271; // -log(.999)*2

	void setExpectationContext(const char *context);

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~ComputeEM();
};

class ComputeStandardError : public omxCompute {
	typedef omxCompute super;
 public:
        virtual void reportResults(FitContext *fc, MxRList *out);
};

static class omxCompute *newComputeSequence()
{ return new omxComputeSequence(); }

static class omxCompute *newComputeIterate()
{ return new omxComputeIterate(); }

static class omxCompute *newComputeOnce()
{ return new omxComputeOnce(); }

static class omxCompute *newComputeEM()
{ return new ComputeEM(); }

static class omxCompute *newComputeStandardError()
{ return new ComputeStandardError(); }

struct omxComputeTableEntry {
        char name[32];
        omxCompute *(*ctor)();
};

static const struct omxComputeTableEntry omxComputeTable[] = {
        {"MxComputeEstimatedHessian", &newComputeEstimatedHessian},
        {"MxComputeGradientDescent", &newComputeGradientDescent},
	{"MxComputeSequence", &newComputeSequence },
	{"MxComputeIterate", &newComputeIterate },
	{"MxComputeOnce", &newComputeOnce },
        {"MxComputeNewtonRaphson", &newComputeNewtonRaphson},
        {"MxComputeEM", &newComputeEM },
	{"MxComputeStandardError", &newComputeStandardError}
};

omxCompute *omxNewCompute(omxState* os, const char *type)
{
        omxCompute *got = NULL;

        for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxComputeTable); fx++) {
                const struct omxComputeTableEntry *entry = omxComputeTable + fx;
                if(strcmp(type, entry->name) == 0) {
                        got = entry->ctor();
                        break;
                }
        }

        if (!got) error("Compute %s is not implemented", type);

        return got;
}

void omxComputeSequence::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("steps")));

	for (int cx = 0; cx < length(slotValue); cx++) {
		SEXP step = VECTOR_ELT(slotValue, cx);
		SEXP s4class;
		PROTECT(s4class = STRING_ELT(getAttrib(step, install("class")), 0));
		omxCompute *compute = omxNewCompute(globalState, CHAR(s4class));
		compute->initFromFrontend(step);
		if (isErrorRaised(globalState)) break;
		clist.push_back(compute);
	}
}

void omxComputeSequence::compute(FitContext *fc)
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		FitContext *context = fc;
		if (fc->varGroup != clist[cx]->varGroup) {
			context = new FitContext(fc, clist[cx]->varGroup);
		}
		clist[cx]->compute(context);
		if (context != fc) context->updateParentAndFree();
		if (isErrorRaised(globalState)) break;
	}
}

void omxComputeSequence::reportResults(FitContext *fc, MxRList *out)
{
	// put this stuff in a new list?
	// merge with Iterate TODO
	for (size_t cx=0; cx < clist.size(); ++cx) {
		FitContext *context = fc;
		if (fc->varGroup != clist[cx]->varGroup) {
			context = new FitContext(fc, clist[cx]->varGroup);
		}
		clist[cx]->reportResults(context, out);
		if (context != fc) context->updateParentAndFree();
		if (isErrorRaised(globalState)) break;
	}
}

double omxComputeSequence::getOptimizerStatus()
{
	// for backward compatibility, not indended to work generally
	for (size_t cx=0; cx < clist.size(); ++cx) {
		double got = clist[cx]->getOptimizerStatus();
		if (got != NA_REAL) return got;
	}
	return NA_REAL;
}

omxComputeSequence::~omxComputeSequence()
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		delete clist[cx];
	}
}

void omxComputeIterate::initFromFrontend(SEXP rObj)
{
	SEXP slotValue;

	super::initFromFrontend(rObj);

	PROTECT(slotValue = GET_SLOT(rObj, install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	PROTECT(slotValue = GET_SLOT(rObj, install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) error("tolerance must be positive");

	PROTECT(slotValue = GET_SLOT(rObj, install("steps")));

	for (int cx = 0; cx < length(slotValue); cx++) {
		SEXP step = VECTOR_ELT(slotValue, cx);
		SEXP s4class;
		PROTECT(s4class = STRING_ELT(getAttrib(step, install("class")), 0));
		omxCompute *compute = omxNewCompute(globalState, CHAR(s4class));
		compute->initFromFrontend(step);
		if (isErrorRaised(globalState)) break;
		clist.push_back(compute);
	}

	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);
}

void omxComputeIterate::compute(FitContext *fc)
{
	int iter = 0;
	double prevFit = 0;
	double mac = tolerance * 10;
	while (1) {
		for (size_t cx=0; cx < clist.size(); ++cx) {
			FitContext *context = fc;
			if (fc->varGroup != clist[cx]->varGroup) {
				context = new FitContext(fc, clist[cx]->varGroup);
			}
			clist[cx]->compute(context);
			if (context != fc) context->updateParentAndFree();
			if (isErrorRaised(globalState)) break;
		}
		if (fc->wanted & FF_COMPUTE_MAXABSCHANGE) {
			if (fc->mac < 0) {
				warning("MAC estimated at %.4f; something is wrong", fc->mac);
				break;
			} else {
				mac = fc->mac;
				if (verbose) mxLog("ComputeIterate: mac %.9g", mac);
			}
		}
		if (fc->wanted & FF_COMPUTE_FIT) {
			if (fc->fit == 0) {
				warning("Fit estimated at 0; something is wrong");
				break;
			}
			if (prevFit != 0) {
				double change = prevFit - fc->fit;
				if (verbose) mxLog("ComputeIterate: fit %.9g change %.9g", fc->fit, change);
				mac = fabs(change);
			} else {
				if (verbose) mxLog("ComputeIterate: initial fit %.9g", fc->fit);
			}
			prevFit = fc->fit;
		}
		if (!(fc->wanted & (FF_COMPUTE_MAXABSCHANGE | FF_COMPUTE_FIT))) {
			omxRaiseErrorf(globalState, "ComputeIterate: neither MAC nor fit available");
		}
		if (isErrorRaised(globalState) || ++iter > maxIter || mac < tolerance) break;
	}
}

void omxComputeIterate::reportResults(FitContext *fc, MxRList *out)
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		FitContext *context = fc;
		if (fc->varGroup != clist[cx]->varGroup) {
			context = new FitContext(fc, clist[cx]->varGroup);
		}
		clist[cx]->reportResults(context, out);
		if (context != fc) context->updateParentAndFree();
		if (isErrorRaised(globalState)) break;
	}
}

double omxComputeIterate::getOptimizerStatus()
{
	// for backward compatibility, not indended to work generally
	for (size_t cx=0; cx < clist.size(); ++cx) {
		double got = clist[cx]->getOptimizerStatus();
		if (got != NA_REAL) return got;
	}
	return NA_REAL;
}

omxComputeIterate::~omxComputeIterate()
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		delete clist[cx];
	}
}

void ComputeEM::initFromFrontend(SEXP rObj)
{
	recentFC = NULL;

	SEXP slotValue;
	SEXP s4class;

	super::initFromFrontend(rObj);

	PROTECT(slotValue = GET_SLOT(rObj, install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	PROTECT(slotValue = GET_SLOT(rObj, install("information")));
	information = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("ramsay")));
	useRamsay = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) error("tolerance must be positive");

	PROTECT(slotValue = GET_SLOT(rObj, install("what")));
	for (int wx=0; wx < length(slotValue); ++wx) {
		int objNum = INTEGER(slotValue)[wx];
		omxExpectation *expectation = globalState->expectationList[objNum];
		setFreeVarGroup(expectation, varGroup);
		omxCompleteExpectation(expectation);
		expectations.push_back(expectation);
	}

	PROTECT(slotValue = GET_SLOT(rObj, install("mstep.fit")));
	PROTECT(s4class = STRING_ELT(getAttrib(slotValue, install("class")), 0));
	fit1 = omxNewCompute(globalState, CHAR(s4class));
	fit1->initFromFrontend(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("fit")));
	PROTECT(s4class = STRING_ELT(getAttrib(slotValue, install("class")), 0));
	fit2 = omxNewCompute(globalState, CHAR(s4class));
	fit2->initFromFrontend(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);
}

void ComputeEM::setExpectationContext(const char *context)
{
	for (size_t wx=0; wx < expectations.size(); ++wx) {
		omxExpectation *expectation = expectations[wx];
		if (verbose >= 4) mxLog("ComputeEM: expectation[%lu] %s context %s", wx, expectation->name, context);
		omxExpectationCompute(expectation, context);
	}
}

void pda(const double *ar, int rows, int cols);

void ComputeEM::compute(FitContext *fc)
{
	int totalMstepIter = 0;
	int iter = 0;
	double prevFit = 0;
	double mac = tolerance * 10;
	bool converged = false;
	const size_t freeVars = fc->varGroup->vars.size();
	const int freeVarsEM = (int) fit1->varGroup->vars.size();
	bool in_middle = false;

	OMXZERO(fc->flavor, freeVars);

	FitContext *tmp = new FitContext(fc, fit1->varGroup);
	for (int vx=0; vx < freeVarsEM; ++vx) {
		fc->flavor[tmp->mapToParent[vx]] = 1;
	}
	tmp->updateParentAndFree();

	ramsay.push_back(new Ramsay1975(fc, int(ramsay.size()), 0, verbose, -1.25)); // other param
	ramsay.push_back(new Ramsay1975(fc, int(ramsay.size()), 0, verbose, -1));    // EM param

	if (verbose >= 1) mxLog("ComputeEM: Welcome, tolerance=%g ramsay=%d info=%d flavors=%ld",
				tolerance, useRamsay, information, ramsay.size());

	while (1) {
		setExpectationContext("EM");

		if (recentFC) delete recentFC;
		recentFC = new FitContext(fc, fit1->varGroup);
		fit1->compute(recentFC);
		if (recentFC->inform == INFORM_ITERATION_LIMIT) {
			fc->inform = INFORM_ITERATION_LIMIT;
			omxRaiseErrorf(globalState, "ComputeEM: iteration limited reached");
			break;
		}
		mstepIter = recentFC->iterations;
		recentFC->updateParent();

		setExpectationContext("");

		{
			FitContext *context = fc;
			if (fc->varGroup != fit2->varGroup) {
				context = new FitContext(fc, fit2->varGroup);
			}

			// For IFA, PREOPTIMIZE updates latent distribution parameters
			omxFitFunction *ff2 = fit2->getFitFunction();
			if (ff2) omxFitFunctionCompute(ff2, FF_COMPUTE_PREOPTIMIZE, context);

			if (!useRamsay) {
				fc->maybeCopyParamToModel(globalState);
			} else {
				context->updateParent();

				bool wantRestart;
				if (iter > 3 && iter % 3 == 0) {
					for (size_t rx=0; rx < ramsay.size(); ++rx) {
						ramsay[rx]->recalibrate(&wantRestart);
					}
				}
				for (size_t rx=0; rx < ramsay.size(); ++rx) {
					ramsay[rx]->apply();
				}
				fc->copyParamToModel(globalState);
			}

			fit2->compute(context);
			if (context != fc) context->updateParentAndFree();
		}

		totalMstepIter += mstepIter;

		if (!(fc->wanted & FF_COMPUTE_FIT)) {
			omxRaiseErrorf(globalState, "ComputeEM: fit not available");
			break;
		}
		if (fc->fit == 0) {
			omxRaiseErrorf(globalState, "Fit estimated at 0; something is wrong");
			break;
		}
		double change = 0;
		if (prevFit != 0) {
			change = prevFit - fc->fit;
			if (0 < change && change < MIDDLE_START) in_middle = true;
			if (verbose >= 2) mxLog("ComputeEM[%d]: msteps %d fit %.9g change %.9g",
						iter, mstepIter, fc->fit, change);
			mac = fabs(change);
		} else {
			if (verbose >= 2) mxLog("ComputeEM: msteps %d initial fit %.9g",
						mstepIter, fc->fit);
		}

		if (in_middle && change > MIDDLE_END) estHistory.push_back(recentFC->take(FF_COMPUTE_ESTIMATE));
		prevFit = fc->fit;
		converged = mac < tolerance;
		if (isErrorRaised(globalState) || ++iter > maxIter || converged) break;
	}

	bestFit = fc->fit;
	if (verbose >= 1) mxLog("ComputeEM: cycles %d/%d total mstep %d fit %f",
				iter, maxIter,totalMstepIter, bestFit);

	if (!converged || !information) return;

	if (estHistory.size() < 2) {
		if (verbose >= 1) mxLog("ComputeEM: history too short to estimate SEs; try increasing EM tolerance");
		return;
	}

	// what about latent distribution parameters? TODO

	recentFC->fixHessianSymmetry(FF_COMPUTE_IHESSIAN);
	double *ihess = recentFC->take(FF_COMPUTE_IHESSIAN);

	double semTolerance = sqrt(tolerance);
	optimum.resize(freeVars);
	memcpy(optimum.data(), fc->est, sizeof(double) * freeVars);
	omxBuffer<double> rij(freeVarsEM * freeVarsEM);
	setExpectationContext("EM");

	// could parallelize by variable instead of at lower level, not sure which is better TODO
	for (int vx=0; vx < freeVarsEM && !isErrorRaised(globalState); ++vx) {
		int base = vx * freeVarsEM;
		bool semConverged;
		for (size_t cx = 0; cx < estHistory.size(); ++cx) {
			if (verbose >= 2) mxLog("ComputeEM: probing param %d from history %ld/%ld",
						vx, cx, estHistory.size());

			memcpy(fc->est, optimum.data(), sizeof(double) * freeVars);
			FitContext *emfc = new FitContext(fc, fit1->varGroup);
			emfc->est[vx] = estHistory[cx][vx];
			emfc->copyParamToModel(globalState);
			fit1->compute(emfc);

			double denom = estHistory[cx][vx] - optimum[emfc->mapToParent[vx]];
			if (verbose >= 1 && fabs(denom) < 1e-5) {
				mxLog("ComputeEM: param %d history %ld denom=%f < 1e-5 (this is bad)", vx, cx, denom);
			}

			semConverged = true;
			for (int v1=0; v1 < freeVarsEM; ++v1) {
				double got = (emfc->est[v1] - optimum[emfc->mapToParent[v1]]) / denom;
				if (fabs(got - rij[base + v1]) >= semTolerance) semConverged = false;
				rij[base + v1] = got;
			}
			//pda(rij.data() + base, 1, freeVarsEM);
			if (verbose >= 3) mxLog("rij[%d]=%f", vx, rij[base + vx]); // useful?
			delete emfc;
			if ((cx > 0 && semConverged) || isErrorRaised(globalState)) break;
		}
	}

	memcpy(fc->est, optimum.data(), sizeof(double) * freeVars);
	fc->copyParamToModel(globalState);

	//pda(rij.data(), freeVarsEM, freeVarsEM);

	// rij = I-rij
	for (int v1=0; v1 < freeVarsEM; ++v1) {
		for (int v2=0; v2 < freeVarsEM; ++v2) {
			int cell = v1 * freeVarsEM + v2;
			double entry = rij[cell];
			if (v1 == v2) entry = 1 - entry;
			else entry = -entry;
			rij[cell] = entry;
		}
	}
	// make symmetric
	for (int v1=1; v1 < freeVarsEM; ++v1) {
		for (int v2=0; v2 < v1; ++v2) {
			int c1 = v1 * freeVarsEM + v2;
			int c2 = v2 * freeVarsEM + v1;
			double mean = (rij[c1] + rij[c2])/2;
			rij[c1] = mean;
			rij[c2] = mean;
		}
	}

	//mxLog("symm");
	//pda(rij.data(), freeVarsEM, freeVarsEM);

	//pda(ihess, freeVarsEM, freeVarsEM);

	// ihess = ihess %*% rij^{-1}
	omxBuffer<int> ipiv(freeVarsEM);
	int info;
	F77_CALL(dgesv)(&freeVarsEM, &freeVarsEM, rij.data(), &freeVarsEM,  // dsysv? TODO
			ipiv.data(), ihess, &freeVarsEM, &info);
	if (info < 0) error("dgesv %d", info);
	if (info > 0) {
		omxRaiseErrorf(globalState, "EM map is not positive definite %d", info);
		return;
	}

	//pda(ihess, freeVarsEM, freeVarsEM);

	// rewrite in terms of ComputeStandardError TODO
	stdError.resize(freeVars);
	for (int v1=0; v1 < freeVarsEM; ++v1) {
		int cell = v1 * freeVarsEM + v1;
		stdError[ recentFC->mapToParent[v1] ] = sqrt(ihess[cell]);
	}	

	delete [] ihess;
}

void ComputeEM::reportResults(FitContext *fc, MxRList *out)
{
	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(bestFit)));
	out->push_back(std::make_pair(mkChar("Minus2LogLikelihood"), ScalarReal(bestFit)));

	size_t numFree = fc->varGroup->vars.size();
	if (!numFree) return;

	if (optimum.size() == numFree) {
		SEXP Rvec;
		PROTECT(Rvec = allocVector(REALSXP, numFree));
		memcpy(REAL(Rvec), optimum.data(), sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("estimate"), Rvec));
	}

	if (stdError.size() == numFree) {
		// make conditional TODO
		SEXP Rvec;
		PROTECT(Rvec = allocVector(REALSXP, numFree));
		memcpy(REAL(Rvec), stdError.data(), sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("standardErrors"), Rvec));
	}
}

double ComputeEM::getOptimizerStatus()
{
	// for backward compatibility, not indended to work generally
	return NA_REAL;
}

ComputeEM::~ComputeEM()
{
	for (size_t rx=0; rx < ramsay.size(); ++rx) {
		delete ramsay[rx];
	}
	ramsay.clear();

	delete fit1;
	delete fit2;

	for (size_t hx=0; hx < estHistory.size(); ++hx) {
		delete [] estHistory[hx];
	}
	estHistory.clear();
	if (recentFC) delete recentFC;
}

void omxComputeOnce::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("what")));
	for (int wx=0; wx < length(slotValue); ++wx) {
		int objNum = INTEGER(slotValue)[wx];
		if (objNum >= 0) {
			omxMatrix *algebra = globalState->algebraList[objNum];
			if (algebra->fitFunction) {
				setFreeVarGroup(algebra->fitFunction, varGroup);
				omxCompleteFitFunction(algebra);
			}
			algebras.push_back(algebra);
		} else {
			omxExpectation *expectation = globalState->expectationList[~objNum];
			setFreeVarGroup(expectation, varGroup);
			omxCompleteExpectation(expectation);
			expectations.push_back(expectation);
		}
	}

	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);

	context = "";

	PROTECT(slotValue = GET_SLOT(rObj, install("context")));
	if (length(slotValue) == 0) {
		// OK
	} else if (length(slotValue) == 1) {
		SEXP elem;
		PROTECT(elem = STRING_ELT(slotValue, 0));
		context = CHAR(elem);
	}

	PROTECT(slotValue = GET_SLOT(rObj, install("maxAbsChange")));
	mac = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("fit")));
	fit = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("gradient")));
	gradient = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("hessian")));
	hessian = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("information")));
	infoMat = asLogical(slotValue);

	if (hessian && infoMat) error("Cannot compute the Hessian and Fisher Information matrix simultaneously");

	PROTECT(slotValue = GET_SLOT(rObj, install("ihessian")));
	ihessian = asLogical(slotValue);

	PROTECT(slotValue = GET_SLOT(rObj, install("hgprod")));
	hgprod = asLogical(slotValue);

	if (algebras.size() == 1 && algebras[0]->fitFunction) {
		omxFitFunction *ff = algebras[0]->fitFunction;
		if (gradient && !ff->gradientAvailable) {
			error("Gradient requested but not available");
		}
		if ((hessian || ihessian || hgprod) && !ff->hessianAvailable) {
			// add a separate flag for hgprod TODO
			error("Hessian requested but not available");
		}
		// add check for information TODO
	}
}

omxFitFunction *omxComputeOnce::getFitFunction()
{
	if (algebras.size() == 1 && algebras[0]->fitFunction) {
		return algebras[0]->fitFunction;
	} else {
		return NULL;
	}
}

void omxComputeOnce::compute(FitContext *fc)
{
	if (algebras.size()) {
		int want = 0;
		size_t numParam = fc->varGroup->vars.size();
		if (mac) {
			want |= FF_COMPUTE_MAXABSCHANGE;
			fc->mac = 0;
		}
		if (fit) {
			want |= FF_COMPUTE_FIT;
			fc->fit = 0;
		}
		if (gradient) {
			want |= FF_COMPUTE_GRADIENT;
			OMXZERO(fc->grad, numParam);
		}
		if (hessian) {
			want |= FF_COMPUTE_HESSIAN;
			OMXZERO(fc->hess, numParam * numParam);
		}
		if (infoMat) {
			want |= FF_COMPUTE_INFO;
			OMXZERO(fc->hess, numParam * numParam);
		}
		if (ihessian) {
			want |= FF_COMPUTE_IHESSIAN;
			OMXZERO(fc->ihess, numParam * numParam);
		}
		if (hgprod) {
			want |= FF_COMPUTE_HGPROD;
			fc->hgProd.resize(0);
		}
		if (!want) return;

		for (size_t wx=0; wx < algebras.size(); ++wx) {
			omxMatrix *algebra = algebras[wx];
			if (algebra->fitFunction) {
				if (verbose) mxLog("ComputeOnce: fit %p want %d",
						   algebra->fitFunction, want);

				omxFitFunctionCompute(algebra->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
				fc->maybeCopyParamToModel(globalState);

				omxFitFunctionCompute(algebra->fitFunction, want, fc);
				fc->fit = algebra->data[0];
				fc->fixHessianSymmetry(want);
			} else {
				if (verbose) mxLog("ComputeOnce: algebra %p", algebra);
				omxForceCompute(algebra);
			}
		}
	} else if (expectations.size()) {
		for (size_t wx=0; wx < expectations.size(); ++wx) {
			omxExpectation *expectation = expectations[wx];
			if (verbose) mxLog("ComputeOnce: expectation[%lu] %p context %s", wx, expectation, context);
			omxExpectationCompute(expectation, context);
		}
	}
}

void omxComputeOnce::reportResults(FitContext *fc, MxRList *out)
{
	if (algebras.size()==0 || algebras[0]->fitFunction == NULL) return;

	omxMatrix *algebra = algebras[0];

	omxPopulateFitFunction(algebra, out);

	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fc->fit)));
	out->push_back(std::make_pair(mkChar("Minus2LogLikelihood"), ScalarReal(fc->fit)));

	size_t numFree = fc->varGroup->vars.size();
	if (numFree) {
		SEXP estimate;
		PROTECT(estimate = allocVector(REALSXP, numFree));
		memcpy(REAL(estimate), fc->est, sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("estimate"), estimate));

		if (gradient) {
			SEXP Rgradient;
			PROTECT(Rgradient = allocVector(REALSXP, numFree));
			memcpy(REAL(Rgradient), fc->grad, sizeof(double) * numFree);
			out->push_back(std::make_pair(mkChar("gradient"), Rgradient));
		}

		if (hessian) {
			SEXP Rhessian;
			PROTECT(Rhessian = allocMatrix(REALSXP, numFree, numFree));
			memcpy(REAL(Rhessian), fc->hess, sizeof(double) * numFree * numFree);
			out->push_back(std::make_pair(mkChar("hessian"), Rhessian));
		}

		if (infoMat) {
			SEXP Rhessian;
			PROTECT(Rhessian = allocMatrix(REALSXP, numFree, numFree));
			memcpy(REAL(Rhessian), fc->hess, sizeof(double) * numFree * numFree);
			out->push_back(std::make_pair(mkChar("information"), Rhessian));
		}

		if (ihessian) {
			SEXP Rihessian;
			PROTECT(Rihessian = allocMatrix(REALSXP, numFree, numFree));
			memcpy(REAL(Rihessian), fc->ihess, sizeof(double) * numFree * numFree);
			out->push_back(std::make_pair(mkChar("ihessian"), Rihessian));
		}

		if (hgprod) {
			// TODO
		}
	}
}

void ComputeStandardError::reportResults(FitContext *fc, MxRList *out)
{
	if (!(fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_INFO | FF_COMPUTE_IHESSIAN))) {
		error("Hessian not available?");
	}

	int numParams = int(fc->varGroup->vars.size());

	if (!(fc->wanted & FF_COMPUTE_IHESSIAN)) {
		// Populate upper triangle
		for(int i = 0; i < numParams; i++) {
			for(int j = 0; j <= i; j++) {
				fc->ihess[i*numParams+j] = fc->hess[i*numParams+j];
			}
		}

		Matrix wmat(fc->ihess, numParams, numParams);
		InvertSymmetricIndef(wmat, 'U');
		fc->fixHessianSymmetry(FF_COMPUTE_IHESSIAN, true);
	}

	// This function calculates the standard errors from the Hessian matrix
	// sqrt(diag(solve(hessian)))

	// We report the fit in -2LL units instead of -LL so we need to adjust here.
	const double scale = sqrt(2); // constexpr

	SEXP stdErrors;
	PROTECT(stdErrors = allocMatrix(REALSXP, numParams, 1));
	double* stdErr = REAL(stdErrors);
	for(int i = 0; i < numParams; i++) {
		stdErr[i] = scale * sqrt(fc->ihess[i * numParams + i]);
	}

	out->push_back(std::make_pair(mkChar("standardErrors"), stdErrors));
}
