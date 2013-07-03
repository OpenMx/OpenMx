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

std::vector<int> FitContext::markMatrices;

void FitContext::init()
{
	size_t numParam = varGroup->vars.size();
	fit = 0;
	est = new double[numParam];
	grad = new double[numParam];
	hess = new double[numParam * numParam];
}

FitContext::FitContext()
{
	parent = NULL;
	varGroup = Global->freeGroup[0];
	init();

	size_t numParam = varGroup->vars.size();
	for (size_t v1=0; v1 < numParam; v1++) {
		est[v1] = Global->freeGroup[0]->vars[v1]->start;
		grad[v1] = nan("unset");
		for (size_t v2=0; v2 < numParam; v2++) {
			hess[v1 * numParam + v2] = nan("unset");
		}
	}
}

FitContext::FitContext(FitContext *parent)
{
	varGroup = parent->varGroup;
	init();
	fit = parent->fit;
	size_t numParam = varGroup->vars.size();
	memcpy(est, parent->est, sizeof(double) * numParam);
	memcpy(grad, parent->grad, sizeof(double) * numParam);
	memcpy(hess, parent->hess, sizeof(double) * numParam * numParam);
}

// arg to control what to copy? usually don't want everything TODO
FitContext::FitContext(FitContext *parent, FreeVarGroup *varGroup)
{
	this->parent = parent;
	this->varGroup = varGroup;
	init();

	FreeVarGroup *src = parent->varGroup;
	FreeVarGroup *dest = varGroup;
	size_t svars = parent->varGroup->vars.size();
	size_t dvars = varGroup->vars.size();

	size_t d1 = 0;
	for (size_t s1=0; s1 < src->vars.size(); ++s1) {
		if (src->vars[s1] != dest->vars[d1]) continue;
		est[d1] = parent->est[s1];
		grad[d1] = parent->grad[s1];

		size_t d2 = 0;
		for (size_t s2=0; s2 < src->vars.size(); ++s2) {
			if (src->vars[s2] != dest->vars[d2]) continue;
			hess[d1 * dvars + d2] = parent->hess[s1 * svars + s2];
			if (++d2 == dvars) break;
		}

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

	size_t s1 = 0;
	for (size_t d1=0; d1 < dest->vars.size(); ++d1) {
		if (dest->vars[d1] != src->vars[s1]) continue;
		parent->est[d1] = est[s1];
		parent->grad[d1] = grad[s1];

		size_t s2 = 0;
		for (size_t d2=0; d2 < dest->vars.size(); ++d2) {
			if (dest->vars[d2] != src->vars[s2]) continue;
			parent->hess[d1 * dvars + d2] = hess[s1 * svars + s2];
			if (++s2 == svars) break;
		}

		if (++s1 == svars) break;
	}
	
	// pda(est, 1, svars);
	// pda(parent->est, 1, dvars);
	// pda(grad, 1, svars);
	// pda(parent->grad, 1, dvars);
	// pda(hess, svars, svars);
	// pda(parent->hess, dvars, dvars);

	delete this;
}

void FitContext::log(int what)
{
	size_t count = varGroup->vars.size();
	std::string buf;
	if (what & FF_COMPUTE_FIT) buf += string_snprintf("fit: %.3f\n", fit);
	if (what & FF_COMPUTE_ESTIMATE) {
		buf += "est: c(";
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.3f", est[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_GRADIENT) {
		buf += "grad: c(";
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.3f", grad[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_HESSIAN) {
		buf += "hess: c(";
		for (size_t v1=0; v1 < count; ++v1) {
			for (size_t v2=0; v2 < count; ++v2) {
				buf += string_snprintf("%.3f", hess[v1 * count + v2]);
				if (v1 < count-1 || v2 < count-1) buf += ", ";
			}
			buf += "\n";
		}
		buf += ")\n";
	}
	mxLogBig(buf);
}

void FitContext::cacheFreeVarDependencies()
{
	omxState *os = globalState;
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	markMatrices.clear();
	markMatrices.resize(numMats + numAlgs, 0);

	// More efficient to use the appropriate group instead of the default group. TODO
	FreeVarGroup *varGroup = Global->freeGroup[0];
	for (size_t freeVarIndex = 0; freeVarIndex < varGroup->vars.size(); freeVarIndex++) {
		omxFreeVar* freeVar = varGroup->vars[freeVarIndex];
		int *deps   = freeVar->deps;
		int numDeps = freeVar->numDeps;
		for (int index = 0; index < numDeps; index++) {
			markMatrices[deps[index] + numMats] = 1;
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

void FitContext::copyParamToModel(omxState* os, double *at)
{
	size_t numParam = varGroup->vars.size();
	if(OMX_DEBUG) {
		mxLog("Copying %d free parameter estimates to model %p", numParam, os);
	}

	if(numParam == 0) return;

	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	os->computeCount++;

	if(OMX_VERBOSE) {
		std::string buf;
		buf += string_snprintf("Call: %d.%d (%d) ", os->majorIteration, os->minorIteration, os->computeCount);
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
				mxLog("Setting location (%d, %d) of matrix %d to value %f for var %d",
					row, col, loc->matrix, at[k], k);
			}
		}
	}

	if (RFitFunction) omxRepopulateRFitFunction(RFitFunction, at, numParam);

	for(size_t i = 0; i < numMats; i++) {
		if (markMatrices[i]) {
			int offset = ~(i - numMats);
			omxMarkDirty(os->matrixList[offset]);
		}
	}

	for(size_t i = 0; i < numAlgs; i++) {
		if (markMatrices[i + numMats]) {
			omxMarkDirty(os->algebraList[i]);
		}
	}

	if (!os->childList) return;

	for(int i = 0; i < Global->numChildren; i++) {
		copyParamToModel(os->childList[i]);
	}
}

FitContext::~FitContext()
{
	delete [] est;
	delete [] grad;
	delete [] hess;
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

omxCompute::~omxCompute()
{}

void omxComputeOperation::initFromFrontend(SEXP rObj)
{
	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("free.group")));
	int paramGroup = INTEGER(slotValue)[0];
	varGroup = Global->freeGroup[paramGroup];
}

class omxComputeAssign : public omxComputeOperation {
	typedef omxComputeOperation super;
	std::vector< int > from;
	std::vector< int > to;
	std::vector< int > pmap;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
};

class omxComputeSequence : public omxCompute {
	std::vector< omxCompute* > clist;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~omxComputeSequence();
};

class omxComputeIterate : public omxCompute {
	std::vector< omxCompute* > clist;
	int maxIter;
	double tolerance;
	bool verbose;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~omxComputeIterate();
};

class omxComputeOnce : public omxComputeOperation {
	typedef omxComputeOperation super;
	omxMatrix *fitMatrix;
	omxExpectation *expectation;
	const char *context;

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

static class omxCompute *newComputeAssign()
{ return new omxComputeAssign(); }

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
	{"MxComputeAssign", &newComputeAssign },
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

void omxComputeAssign::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	SEXP slotValue;

	PROTECT(slotValue = GET_SLOT(rObj, install("from")));
	for (int mx=0; mx < length(slotValue); ++mx) {
		from.push_back(INTEGER(slotValue)[mx]);
	}

	PROTECT(slotValue = GET_SLOT(rObj, install("to")));
	for (int mx=0; mx < length(slotValue); ++mx) {
		to.push_back(INTEGER(slotValue)[mx]);
	}

	if (to.size() != from.size()) {
		error("omxComputeAssign: length of from and to must match");
	}

	pmap.assign(varGroup->vars.size(), -1);
	size_t offset=0;
	for (size_t mx=0; mx < from.size(); mx++) {
		int fmat = from[mx];
		int tmat = to[mx];
		for (size_t p1=0; p1 < varGroup->vars.size(); ++p1) {
			omxFreeVarLocation *l1 = varGroup->vars[p1]->getLocation(fmat);
			if (!l1) continue;
			bool found = FALSE;
			for (size_t r1=0; r1 < varGroup->vars.size(); ++r1) {
				size_t p2 = (r1 + offset) % varGroup->vars.size();
				omxFreeVarLocation *l2 = varGroup->vars[p2]->getLocation(tmat);
				if (!l2) continue;
				if (l1->row != l2->row || l1->col != l2->col) continue;
				if (pmap[p1] != -1) {
					error("omxComputeAssign: cannot copy %s to more than 1 place",
					      globalState->matrixList[~fmat]->name);
				}
				if (p1 == p2) error("omxComputeAssign: cannot copy %s to itself",
						    globalState->matrixList[~fmat]->name);
				pmap[p1] = p2;
				offset = p2+1;  // probably a good guess
				found = TRUE;
				break;
			}
			if (!found) error("omxComputeAssign: %s[%d,%d] is fixed, cannot copy %s",
					  globalState->matrixList[~tmat]->name, 1+l1->row, 1+l1->col,
					  globalState->matrixList[~fmat]->name);
		}
	}
}

void omxComputeAssign::compute(FitContext *fc)
{
	for (size_t px=0; px < varGroup->vars.size(); ++px) {
		if (pmap[px] == -1) continue;
		fc->est[ pmap[px] ] = fc->est[px];
	}
	fc->copyParamToModel(globalState);
}

void omxComputeAssign::reportResults(FitContext *fc, MxRList *out)
{}

void omxComputeSequence::initFromFrontend(SEXP rObj)
{
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
	verbose = asLogical(slotValue);
}

void omxComputeIterate::compute(FitContext *fc)
{
	int iter = 0;
	double prevFit = 0;
	double change = tolerance * 10;
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
		if (prevFit != 0) {
			change = prevFit - fc->fit;
			if (verbose) mxLog("fit %.9g change %.9g", fc->fit, change);
		}
		prevFit = fc->fit;
		if (isErrorRaised(globalState) || ++iter > maxIter || fabs(change) < tolerance) break;
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

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	if (fitMatrix) {
		setFreeVarGroup(fitMatrix->fitFunction, varGroup);
		omxCompleteFitFunction(fitMatrix);
	}

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("expectation")));
	if (length(slotValue)) {
		int expNumber = INTEGER(slotValue)[0];	
		expectation = omxExpectationFromIndex(expNumber, globalState);
		setFreeVarGroup(expectation, varGroup);
		omxCompleteExpectation(expectation);
	}

	if (fitMatrix && expectation) {
		error("Cannot evaluate a fitfunction and expectation simultaneously");
	}
	if (!fitMatrix && !expectation) error("No function specified to evaluate");

	PROTECT(slotValue = GET_SLOT(rObj, install("context")));
	if (length(slotValue) == 0) {
		// OK
	} else if (length(slotValue) == 1) {
		SEXP elem;
		PROTECT(elem = STRING_ELT(slotValue, 0));
		context = CHAR(elem);
	}
}

void omxComputeOnce::compute(FitContext *fc)
{
	if (fitMatrix) {
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, fc);
		fc->fit = fitMatrix->data[0];
	} else if (expectation) {
		omxExpectationCompute(expectation, context);
	}
}

void omxComputeOnce::reportResults(FitContext *fc, MxRList *out)
{
	if (!fitMatrix) return;

	omxPopulateFitFunction(fitMatrix, out);

	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fc->fit)));

	size_t numFree = fc->varGroup->vars.size();
	if (numFree) {
		SEXP estimate;
		PROTECT(estimate = allocVector(REALSXP, numFree));
		memcpy(REAL(estimate), fc->est, sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("estimate"), estimate));
	}
}
