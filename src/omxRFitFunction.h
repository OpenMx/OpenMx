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

#ifndef _OMX_R_FITFUNCTION_
#define _OMX_R_FITFUNCTION_

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

typedef struct {

	SEXP fitfun;
	SEXP model;
	PROTECT_INDEX modelIndex;
	SEXP flatModel;
	SEXP parameters;
	SEXP state;
	PROTECT_INDEX stateIndex;

} omxRFitFunction;

void omxDestroyRFitFunction(omxFitFunction *oo);

void omxRepopulateRFitFunction(omxFitFunction* oo, double* x, int n);

omxRListElement* omxSetFinalReturnsRFitFunction(omxFitFunction *oo, int *numReturns);

void omxInitRFitFunction(omxFitFunction* oo);

#endif /* _OMX_R_FITFUNCTION_ */
