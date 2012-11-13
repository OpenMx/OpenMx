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

#ifndef _OMX_R_OBJECTIVE_
#define _OMX_R_OBJECTIVE_

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

typedef struct {

	SEXP objfun;
	SEXP model;
	PROTECT_INDEX modelIndex;
	SEXP flatModel;
	SEXP parameters;
	SEXP state;
	PROTECT_INDEX stateIndex;

} omxRObjective;

#ifdef _OPENMP

#include "omp.h"
extern omp_lock_t robjective_lock;

#else // (!defined _OPENMP)

extern void* robjective_lock;

#endif // _OPENMP


void omxDestroyRObjective(omxObjective *oo);
void omxCallRObjective(omxObjective *oo);

void omxRepopulateRObjective(omxObjective* oo, double* x, int n);

omxRListElement* omxSetFinalReturnsRObjective(omxObjective *oo, int *numReturns);

void omxInitRObjective(omxObjective* oo, SEXP rObj);

#endif /* _OMX_R_OBJECTIVE_ */
