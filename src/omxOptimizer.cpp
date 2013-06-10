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

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <sys/types.h>

#include "omxDefines.h"
#include "omxState.h"

void cacheFreeVarDependencies(omxState* os)
{
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	os->markMatrices.clear();
	os->markMatrices.resize(numMats + numAlgs, 0);

	for(int freeVarIndex = 0; freeVarIndex < os->numFreeParams; freeVarIndex++) {
		omxFreeVar* freeVar = os->freeVarList + freeVarIndex;
		int *deps   = freeVar->deps;
		int numDeps = freeVar->numDeps;
		for (int index = 0; index < numDeps; index++) {
			os->markMatrices[deps[index] + numMats] = 1;
		}
	}

}

void markFreeVarDependenciesHelper(omxState* os, int varNumber) {

	int numDeps = os->freeVarList[varNumber].numDeps;
	int *deps = os->freeVarList[varNumber].deps;

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

	int numChildren = os->numChildren;

	markFreeVarDependenciesHelper(os, varNumber);

	for(int i = 0; i < numChildren; i++) {
		markFreeVarDependencies(os->childList[i], varNumber);
	}
}

void handleFreeVarListHelper(omxState* os, double* x, int numVars)
{
	int numChildren = os->numChildren;

	if(OMX_DEBUG && os->parentState == NULL) {
		Rprintf("Processing Free Parameter Estimates.\n");
		Rprintf("Number of free parameters is %d.\n", numVars);
	}

	if(numVars == 0) return;

	omxFreeVar* freeVarList = os->freeVarList;
	size_t numMats = os->matrixList.size();
	int numAlgs = os->algebraList.size();

	os->computeCount++;

	if(OMX_VERBOSE && os->parentState == NULL) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d.%d (%d)\n", os->majorIteration, os->minorIteration, os->computeCount);
		Rprintf("Estimates: [");
		for(int k = 0; k < numVars; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
		Rprintf("--------------------------\n");
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < numVars; k++) {
		omxFreeVar* freeVar = freeVarList + k;
		// if(OMX_DEBUG) { Rprintf("%d: %f - %d\n", k,  x[k], freeVarList[k].numLocations); }
		for(size_t l = 0; l < freeVar->locations.size(); l++) {
			omxFreeVarLocation *loc = &freeVar->locations[l];
			omxMatrix *matrix = os->matrixList[loc->matrix];
			int row = loc->row;
			int col = loc->col;
			omxSetMatrixElement(matrix, row, col, x[k]);
			if(OMX_DEBUG && os->parentState == NULL) {
				Rprintf("Setting location (%d, %d) of matrix %d to value %f for var %d\n",
					row, col, loc->matrix, x[k], k);
			}
		}
	}

	for(size_t i = 0; i < numMats; i++) {
		if (os->markMatrices[i]) {
			int offset = ~(i - numMats);
			omxMarkDirty(os->matrixList[offset]);
		}
	}

	for(int i = 0; i < numAlgs; i++) {
		if (os->markMatrices[i + numMats]) {
			omxMarkDirty(os->algebraList[i]);
		}
	}

	for(int i = 0; i < numChildren; i++) {
		handleFreeVarListHelper(os->childList[i], x, numVars);
	}
}

/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(omxFitFunction* oo, double* x, int numVars)
{
	handleFreeVarListHelper(oo->matrix->currentState, x, numVars);
}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str) {
/* Attribution: modified from the code given in Writing R Extensions */
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

SEXP getVar(SEXP str, SEXP env) {
/* Attribution: modified from the code given in Writing R Extensions */
   SEXP ans;
   if(!isString(str) || length(str) != 1)
        error("getVar: variable name is not a single string");
   if(!isEnvironment(env))
	error("getVar: env should be an environment");
   ans = findVar(install(CHAR(STRING_ELT(str, 0))), env);
   return(ans);
}

