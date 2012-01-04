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

#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "npsolWrap.h"
#include "omxImportFrontendState.h"

/* NPSOL-specific globals */
const double NPSOL_BIGBND = 1e20;
const double NEG_INF = -2e20;
const double INF = 2e20;

/* NPSOL-related functions */
extern void F77_SUB(npoptn)(char* string, int length);

void omxSetNPSOLOpts(SEXP options, int *numHessians, int *calculateStdErrors, 
	int *ciMaxIterations, int *disableOptimizer, int *numThreads) {

		char optionCharArray[250] = "";			// For setting options
		int numOptions = length(options);
		SEXP optionNames;
		PROTECT(optionNames = GET_NAMES(options));
		for(int i = 0; i < numOptions; i++) {
			const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
			const char *nextOptionValue = STRING_VALUE(VECTOR_ELT(options, i));
			int lenName = strlen(nextOptionName);
			int lenValue = strlen(nextOptionValue);
			if(matchCaseInsensitive(nextOptionName, lenName, "Calculate Hessian")) {
				if(OMX_DEBUG) { Rprintf("Found hessian option... Value: %s. ", nextOptionValue);};
				if(!matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Enabling explicit hessian calculation.\n");}
					*numHessians = 1;
				}
			} else if(matchCaseInsensitive(nextOptionName, lenName, "Standard Errors")) {
				if(OMX_DEBUG) { Rprintf("Found standard error option...Value: %s. ", nextOptionValue);};
				if(!matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Enabling explicit standard error calculation.\n");}
					*calculateStdErrors = TRUE;
					*numHessians = 1;
				}
			} else if(matchCaseInsensitive(nextOptionName, lenName, "CI Max Iterations")) { 
				int newvalue = atoi(nextOptionValue);
				if (newvalue > 0) *ciMaxIterations = newvalue;
			} else if(matchCaseInsensitive(nextOptionName, lenName, "useOptimizer")) {
				if(OMX_DEBUG) { Rprintf("Found useOptimizer option...");};	
				if(matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Disabling optimization.\n");}
					*disableOptimizer = 1;
				}
			} else if(matchCaseInsensitive(nextOptionName, lenName, "Number of Threads")) {
				*numThreads = atoi(nextOptionValue);
				if(OMX_DEBUG) { Rprintf("Found Number of Threads option (# = %d)...\n", *numThreads);};
			} else {
				sprintf(optionCharArray, "%s %s", nextOptionName, nextOptionValue);
				F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
				if(OMX_DEBUG) { Rprintf("Option %s \n", optionCharArray); }
			}
		}
		UNPROTECT(1); // optionNames
}

