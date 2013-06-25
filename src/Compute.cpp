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
	numParam = varGroup->vars.size();
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
	memcpy(est, parent->est, sizeof(double) * numParam);
	memcpy(grad, parent->grad, sizeof(double) * numParam);
	memcpy(hess, parent->hess, sizeof(double) * numParam * numParam);
}

// arg to control what to copy? usually don't want everything TODO
FitContext::FitContext(FitContext *parent, int group)
{
	this->parent = parent;
	varGroup = Global->freeGroup[group];
	init();

	FreeVarGroup *src = parent->varGroup;
	FreeVarGroup *dest = varGroup;

	size_t d1 = 0;
	for (size_t s1=0; s1 < src->vars.size(); ++s1) {
		if (src->vars[s1] != dest->vars[d1]) continue;
		est[d1] = parent->est[s1];
		grad[d1] = parent->grad[s1];
		++d1;

		size_t d2 = 0;
		for (size_t s2=0; s2 < src->vars.size(); ++s2) {
			if (src->vars[s2] != dest->vars[d2]) continue;
			size_t cell = s1 * numParam + s2;
			hess[cell] = parent->hess[cell];
			++d2;
		}
	}
	if (d1 != numParam-1) error("Parent free parameter group is not a superset");
}

void FitContext::copyParamToModel(omxMatrix *mat)
{ copyParamToModel(mat->currentState); }

void FitContext::copyParamToModel(omxMatrix *mat, double *at)
{ copyParamToModel(mat->currentState, at); }

void FitContext::updateParentAndFree()
{
	FreeVarGroup *src = varGroup;
	FreeVarGroup *dest = parent->varGroup;

	size_t s1 = 0;
	for (size_t d1=0; d1 < dest->vars.size(); ++d1) {
		if (dest->vars[d1] != src->vars[s1]) continue;
		parent->est[d1] = est[s1];
		parent->grad[d1] = grad[s1];
		++s1;

		size_t s2 = 0;
		for (size_t d2=0; d2 < dest->vars.size(); ++d2) {
			if (dest->vars[d2] != src->vars[s2]) continue;
			size_t cell = d1 * numParam + d2;
			parent->hess[cell] = hess[cell];
			++s2;
		}
	}
	
	delete this;
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
	paramGroup = INTEGER(slotValue)[0];
}

class omxComputeSequence : public omxCompute {
	std::vector< omxCompute* > clist;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus();
	virtual ~omxComputeSequence();
};

class omxComputeOnce : public omxComputeOperation {
	typedef omxComputeOperation super;
	omxMatrix *fitMatrix;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *out);
};

static class omxCompute *newComputeSequence()
{ return new omxComputeSequence(); }

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
	{"MxComputeOnce", &newComputeOnce },
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
		int cgroup = clist[cx]->paramGroup;
		if (cgroup) context = new FitContext(fc, cgroup);
		clist[cx]->compute(context);
		if (cgroup) context->updateParentAndFree();
		if (isErrorRaised(globalState)) break;
	}
}

void omxComputeSequence::reportResults(FitContext *fc, MxRList *out)
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		FitContext *context = fc;
		int cgroup = clist[cx]->paramGroup;
		if (cgroup) context = new FitContext(fc, cgroup);
		clist[cx]->reportResults(context, out);
		if (cgroup) context->updateParentAndFree();
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

void omxComputeOnce::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, Global->freeGroup[paramGroup]);
}

void omxComputeOnce::compute(FitContext *fc)
{
        for(size_t index = 0; index < globalState->matrixList.size(); index++) {
            omxMarkDirty(globalState->matrixList[index]);
        }
        for(size_t index = 0; index < globalState->algebraList.size(); index++) {
            omxMarkDirty(globalState->algebraList[index]);
        }
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, NULL);
	fc->fit = fitMatrix->data[0];
}

void omxComputeOnce::reportResults(FitContext *fc, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fc->fit)));

	size_t numFree = fc->numParam;
	if (numFree) {
		SEXP estimate;
		PROTECT(estimate = allocVector(REALSXP, numFree));
		memcpy(REAL(estimate), fc->est, sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("estimate"), estimate));
	}
}
