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

void FitContext::init()
{
	size_t numParam = varGroup->vars.size();
	wanted = 0;
	mac = parent? parent->mac : 0;
	fit = parent? parent->fit : 0;
	caution = parent? parent->caution : 0;
	est = new double[numParam];
	flavor = new int[numParam];
	forwardDeriv = false;
	grad = new double[numParam];
	hess = new double[numParam * numParam];
	ihess = new double[numParam * numParam];
	changedEstimates = false;
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

	size_t d1 = 0;
	for (size_t s1=0; s1 < src->vars.size(); ++s1) {
		if (src->vars[s1] != dest->vars[d1]) continue;
		est[d1] = parent->est[s1];

		if (forwardDeriv) {
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

void FitContext::updateParentAndFree()
{
	FreeVarGroup *src = varGroup;
	FreeVarGroup *dest = parent->varGroup;
	size_t svars = varGroup->vars.size();
	size_t dvars = parent->varGroup->vars.size();

	parent->wanted = wanted;
	parent->fit = fit;
	parent->mac = mac;
	parent->caution = caution;

	if (svars > 0) {
		size_t s1 = 0;
		for (size_t d1=0; d1 < dest->vars.size(); ++d1) {
			if (dest->vars[d1] != src->vars[s1]) continue;
			parent->est[d1] = est[s1];

			if (forwardDeriv) {
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
	}
	
	// pda(est, 1, svars);
	// pda(parent->est, 1, dvars);
	// pda(grad, 1, svars);
	// pda(parent->grad, 1, dvars);
	// pda(hess, svars, svars);
	// pda(parent->hess, dvars, dvars);

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
	if (what & FF_COMPUTE_HESSIAN) {
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

void FitContext::fixHessianSymmetry(int want)
{
	size_t numParam = varGroup->vars.size();

	if (want & FF_COMPUTE_HESSIAN) {
		for (size_t h1=1; h1 < numParam; h1++) {
			for (size_t h2=0; h2 < h1; h2++) {
				if (hess[h2 * numParam + h1] != 0) {
					omxRaiseErrorf(globalState, "Hessian is not upper triangular");
					break;
				}
				hess[h2 * numParam + h1] = hess[h1 * numParam + h2];
			}
		}
	}

	if (want & FF_COMPUTE_IHESSIAN) {
		for (size_t h1=1; h1 < numParam; h1++) {
			for (size_t h2=0; h2 < h1; h2++) {
				if (ihess[h2 * numParam + h1] != 0) {
					omxRaiseErrorf(globalState, "Inverse Hessian is not upper triangular");
					break;
				}
				ihess[h2 * numParam + h1] = ihess[h1 * numParam + h2];
			}
		}
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

FitContext::~FitContext()
{
	delete [] est;
	delete [] flavor;
	delete [] grad;
	delete [] hess;
	delete [] ihess;
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
	bool hgprod;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
};

static class omxCompute *newComputeSequence()
{ return new omxComputeSequence(); }

static class omxCompute *newComputeIterate()
{ return new omxComputeIterate(); }

static class omxCompute *newComputeOnce()
{ return new omxComputeOnce(); }

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
