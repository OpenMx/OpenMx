/*
 *  Copyright 2007-2011 The OpenMx Project
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


/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(omxState* os, double* x, int numVars) {

	int numChildren = os->numChildren;

	if(OMX_DEBUG && os->parentState == NULL) {
		Rprintf("Processing Free Parameter Estimates.\n");
		Rprintf("Number of free parameters is %d.\n", numVars);
	}

	if(numVars == 0) return;

	omxFreeVar* freeVarList = os->freeVarList;
	omxMatrix** matrixList = os->matrixList;

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
			omxMarkDirty(matrix);
		}
	}

	for(int i = 0; i < numChildren; i++) {
		handleFreeVarList(os->childList[i], x, numVars);
	}
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

