/*
 *  Copyright 2007-2013 The OpenMx Project
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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <sys/types.h>

#include "omxDefines.h"
#include "omxState.h"
#include "Compute.h"
#include "omxRFitFunction.h"

// TODO move code to other FitContext method defs

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
