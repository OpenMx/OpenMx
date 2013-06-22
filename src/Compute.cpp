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

class omxComputeSequence : public omxCompute {
	std::vector< omxCompute* > clist;
	double *est;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(double *startVals);
        virtual void reportResults(MxRList *out);
	virtual double getFit() { return 0; }
	virtual double *getEstimate() { return est; }
	virtual double getOptimizerStatus();
	virtual ~omxComputeSequence();
};

class omxComputeOnce : public omxCompute {
	omxMatrix *fitMatrix;
	double fit;
	double *est;

 public:
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(double *startVals);
        virtual void reportResults(MxRList *out);
	virtual double getFit() { return fit; }
	virtual double *getEstimate() { return est; }
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

void omxComputeSequence::compute(double *startVals)
{
	est = startVals;
	for (size_t cx=0; cx < clist.size(); ++cx) {
		clist[cx]->compute(est);
		est = clist[cx]->getEstimate();
		if (isErrorRaised(globalState)) break;
	}
}

void omxComputeSequence::reportResults(MxRList *out)
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		clist[cx]->reportResults(out);
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
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
}

void omxComputeOnce::compute(double *startVals)
{
	est = startVals;
        for(size_t index = 0; index < globalState->matrixList.size(); index++) {
            omxMarkDirty(globalState->matrixList[index]);
        }
        for(size_t index = 0; index < globalState->algebraList.size(); index++) {
            omxMarkDirty(globalState->algebraList[index]);
        }
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, NULL);
	fit = fitMatrix->data[0];
}

void omxComputeOnce::reportResults(MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fit)));

	if (est) {
		int numFree = Global.numFreeParams;
		SEXP estimate;
		PROTECT(estimate = allocVector(REALSXP, numFree));
		memcpy(REAL(estimate), est, sizeof(double)*numFree);
		out->push_back(std::make_pair(mkChar("estimate"), estimate));
	}
}
