/*
 *  Copyright 2007-2012 The OpenMx Project
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

void markFreeVarDependenciesHelper(omxState* os, int varNumber) {

	int numDeps = os->freeVarList[varNumber].numDeps;
	int *deps = os->freeVarList[varNumber].deps;

	omxMatrix** matrixList = os->matrixList;
	omxMatrix** algebraList = os->algebraList;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(matrixList[~value]);
		} else {
			omxMarkDirty(algebraList[value]);
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

void handleFreeVarListHelper(omxState* os, double* x, int numVars, int parent, int *markMatrices) {

	int numChildren = os->numChildren;

	if(OMX_DEBUG && os->parentState == NULL) {
		Rprintf("Processing Free Parameter Estimates.\n");
		Rprintf("Number of free parameters is %d.\n", numVars);
	}

	if(numVars == 0) return;

	omxFreeVar* freeVarList = os->freeVarList;
	omxMatrix** matrixList  = os->matrixList;
	omxMatrix** algebraList = os->algebraList;
	int numMats = os->numMats;
	int numAlgs = os->numAlgs;

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
		// if(OMX_DEBUG) { Rprintf("%d: %f - %d\n", k,  x[k], freeVarList[k].numLocations); }
		for(int l = 0; l < freeVarList[k].numLocations; l++) {
			omxMatrix *matrix = matrixList[freeVarList[k].matrices[l]];
			int row = freeVarList[k].row[l];
			int col = freeVarList[k].col[l];
			omxSetMatrixElement(matrix, row, col, x[k]);
			if(OMX_DEBUG && os->parentState == NULL) {
				Rprintf("Setting location (%d, %d) of matrix %d to value %f for var %d\n",
					row, col, freeVarList[k].matrices[l], x[k], k);
			}
		}
		if (parent) {
			int *deps   = os->freeVarList[k].deps;
			int numDeps = os->freeVarList[k].numDeps;
			for (int index = 0; index < numDeps; index++) {
				markMatrices[deps[index] + numMats] = 1;
			}
		}
	}

	for(int i = 0; i < numMats; i++) {
		if (markMatrices[i]) {
			int offset = ~(i - numMats);
			omxMarkDirty(matrixList[offset]);
		}
	}

	for(int i = 0; i < numAlgs; i++) {
		if (markMatrices[i + numMats]) {
			omxMarkDirty(algebraList[i]);
		}
	}

	// The if-statement is redundant, but the OpenMP
	// specification is ambiguous on the outcome of num_threads(0)
	if (numChildren > 0) {
		#pragma omp parallel for num_threads(numChildren) 
		for(int i = 0; i < numChildren; i++) {
			handleFreeVarListHelper(os->childList[i], x, numVars, 0, markMatrices);
		}
	}
}

/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(omxState* os, double* x, int numVars) {

	memset(os->markMatrices, 0, sizeof(int) * (os->numMats + os->numAlgs));
	handleFreeVarListHelper(os, x, numVars, 1, os->markMatrices);

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

