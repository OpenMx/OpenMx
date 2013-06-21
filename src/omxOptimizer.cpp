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

static std::vector<int> markMatrices;   // constant, therefore thread-safe

void cacheFreeVarDependencies()
{
	omxState *os = globalState;
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	markMatrices.clear();
	markMatrices.resize(numMats + numAlgs, 0);

	for(int freeVarIndex = 0; freeVarIndex < Global.numFreeParams; freeVarIndex++) {
		omxFreeVar* freeVar = Global.freeVarList + freeVarIndex;
		int *deps   = freeVar->deps;
		int numDeps = freeVar->numDeps;
		for (int index = 0; index < numDeps; index++) {
			markMatrices[deps[index] + numMats] = 1;
		}
	}

}

void markFreeVarDependenciesHelper(omxState* os, int varNumber) {

	int numDeps = Global.freeVarList[varNumber].numDeps;
	int *deps = Global.freeVarList[varNumber].deps;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(os->matrixList[~value]);
		} else {
			omxMarkDirty(os->algebraList[value]);
		}
	}

}

void markFreeVarDependencies(omxState* os, int varNumber) {

	int numChildren = Global.numChildren;

	markFreeVarDependenciesHelper(os, varNumber);

	if (!os->childList) return;

	for(int i = 0; i < numChildren; i++) {
		markFreeVarDependencies(os->childList[i], varNumber);
	}
}

void handleFreeVarListHelper(omxState* os, double* x, int numVars)
{
	int numChildren = Global.numChildren;

	if(OMX_DEBUG) {
		mxLog("Processing Free Parameter Estimates.");
		mxLog("Number of free parameters is %d.", numVars);
	}

	if(numVars == 0) return;

	omxFreeVar* freeVarList = Global.freeVarList;
	size_t numMats = os->matrixList.size();
	int numAlgs = os->algebraList.size();

	os->computeCount++;

	if(OMX_VERBOSE) {
		std::string buf;
		buf += "--------------------------\n";
		buf += string_snprintf("Call: %d.%d (%d)", os->majorIteration, os->minorIteration, os->computeCount);
		buf += ("Estimates: [");
		for(int k = 0; k < numVars; k++) {
			buf += string_snprintf(" %f", x[k]);
		}
		buf += ("]\n");
		buf += "--------------------------\n";
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < numVars; k++) {
		omxFreeVar* freeVar = freeVarList + k;
		// if(OMX_DEBUG) { mxLog("%d: %f - %d", k,  x[k], freeVarList[k].numLocations); }
		for(size_t l = 0; l < freeVar->locations.size(); l++) {
			omxFreeVarLocation *loc = &freeVar->locations[l];
			omxMatrix *matrix = os->matrixList[loc->matrix];
			int row = loc->row;
			int col = loc->col;
			omxSetMatrixElement(matrix, row, col, x[k]);
			if(OMX_DEBUG) {
				mxLog("Setting location (%d, %d) of matrix %d to value %f for var %d",
					row, col, loc->matrix, x[k], k);
			}
		}
	}

	for(size_t i = 0; i < numMats; i++) {
		if (markMatrices[i]) {
			int offset = ~(i - numMats);
			omxMarkDirty(os->matrixList[offset]);
		}
	}

	for(int i = 0; i < numAlgs; i++) {
		if (markMatrices[i + numMats]) {
			omxMarkDirty(os->algebraList[i]);
		}
	}

	if (!os->childList) return;

	for(int i = 0; i < numChildren; i++) {
		handleFreeVarListHelper(os->childList[i], x, numVars);
	}
}

/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(omxFitFunction* oo, double* x, int numVars)
{
	handleFreeVarListHelper(oo->matrix->currentState, x, numVars);
}
