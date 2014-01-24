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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"
#include "omxRObjective.h"
#include "omxOpenmpWrap.h"

#ifdef _OPENMP

omp_lock_t robjective_lock;

#else

void* robjective_lock = NULL;

#endif

void omxDestroyRObjective(omxObjective *oo) {

	UNPROTECT(5); 			// objfun, model, flatModel, parameters, and state
}

void omxCallRObjective(omxObjective *oo) {
	omx_omp_set_lock(&robjective_lock);

	omxRObjective* rObjective = (omxRObjective*)oo->argStruct;
	SEXP theCall, theReturn;
	PROTECT(theCall = allocVector(LANGSXP, 3));
	SETCAR(theCall, rObjective->objfun);
	SETCADR(theCall, rObjective->model);
	SETCADDR(theCall, rObjective->state);

	PROTECT(theReturn = eval(theCall, R_GlobalEnv));

	if (LENGTH(theReturn) == 1) {
		oo->matrix->data[0] = REAL(AS_NUMERIC(theReturn))[0];
	} else if (LENGTH(theReturn) == 2) {
		oo->matrix->data[0] = REAL(VECTOR_ELT(theReturn, 0))[0];
		REPROTECT(rObjective->state = VECTOR_ELT(theReturn, 1), rObjective->stateIndex);
	} else {
		// throw an error
	}

	UNPROTECT(2); // theCall and theReturn

	omx_omp_unset_lock(&robjective_lock);
}

// I have no idea what I'm supposed to do here...
unsigned short int omxNeedsUpdateRObjective(omxObjective* oo) {
	return(TRUE);
}

void omxRepopulateRObjective(omxObjective* oo, double* x, int n) {
	omx_omp_set_lock(&robjective_lock);

	omxRObjective* rObjective = (omxRObjective*)oo->argStruct;

	SEXP theCall, estimate;

	PROTECT(estimate = allocVector(REALSXP, n));
	double *est = REAL(estimate);
	for(int i = 0; i < n ; i++) {
		est[i] = x[i];
	}

	PROTECT(theCall = allocVector(LANGSXP, 5));

	SETCAR(theCall, install("imxUpdateModelValues"));
	SETCADR(theCall, rObjective->model);
	SETCADDR(theCall, rObjective->flatModel);
	SETCADDDR(theCall, rObjective->parameters);
	SETCAD4R(theCall, estimate);

	REPROTECT(rObjective->model = eval(theCall, R_GlobalEnv), rObjective->modelIndex);

	UNPROTECT(2); // theCall, estimate
	omx_omp_unset_lock(&robjective_lock);

	omxMarkDirty(oo->matrix);
}

omxRListElement* omxSetFinalReturnsRObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 1;
	omxRListElement* retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = oo->matrix->data[0];

	return retVal;
}

void omxInitRObjective(omxObjective* oo, SEXP rObj) {
	if(OMX_DEBUG) { Rprintf("Initializing R objective function.\n"); }
	omxRObjective *newObj = (omxRObjective*) R_alloc(1, sizeof(omxRObjective));
	
	/* Set Objective Calls to R Objective Calls */
	oo->objectiveFun = omxCallRObjective;
	oo->needsUpdateFun = omxNeedsUpdateRObjective;
	oo->setFinalReturns = omxSetFinalReturnsRObjective;
	oo->destructFun = omxDestroyRObjective;
	oo->repopulateFun = omxRepopulateRObjective;
	oo->argStruct = (void*) newObj;
	
	PROTECT(newObj->objfun = GET_SLOT(rObj, install("objfun")));
	PROTECT_WITH_INDEX(newObj->model = GET_SLOT(rObj, install("model")), &(newObj->modelIndex));
	PROTECT(newObj->flatModel = GET_SLOT(rObj, install("flatModel")));
	PROTECT(newObj->parameters = GET_SLOT(rObj, install("parameters")));
	PROTECT_WITH_INDEX(newObj->state = GET_SLOT(rObj, install("state")), &(newObj->stateIndex));

}


