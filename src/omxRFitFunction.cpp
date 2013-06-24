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
#include "omxAlgebraFunctions.h"
#include "omxRFitFunction.h"
#include "omxOpenmpWrap.h"
#include "npsolWrap.h"
#include "Compute.h"

void omxDestroyRFitFunction(omxFitFunction *off) {
	UNPROTECT(4); 			// fitfun, model, flatModel, and state
}

static void omxCallRFitFunction(omxFitFunction *oo, int want, double *gradient) {
	omxState* currentState = oo->matrix->currentState;
	omxRFitFunction* rFitFunction = (omxRFitFunction*)oo->argStruct;

	SEXP theCall, theReturn;
	PROTECT(theCall = allocVector(LANGSXP, 3));
	SETCAR(theCall, rFitFunction->fitfun);
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->state);

	PROTECT(theReturn = eval(theCall, R_GlobalEnv));

	if (LENGTH(theReturn) < 1) {
		// seems impossible, but report it if it happens
		omxRaiseErrorf(currentState, "FitFunction returned nothing");
	} else if (LENGTH(theReturn) == 1) {
		oo->matrix->data[0] = asReal(theReturn);
	} else if (LENGTH(theReturn) == 2) {
		oo->matrix->data[0] = asReal(VECTOR_ELT(theReturn, 0));
		REPROTECT(rFitFunction->state = VECTOR_ELT(theReturn, 1), rFitFunction->stateIndex);
	} else if (LENGTH(theReturn) > 2) {
		omxRaiseErrorf(currentState, "FitFunction returned more than 2 arguments");
	}

	UNPROTECT(2); // theCall and theReturn
}

omxRListElement* omxSetFinalReturnsRFitFunction(omxFitFunction *oo, int *numReturns) {
	*numReturns = 1;
	omxRListElement* retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = oo->matrix->data[0];

	return retVal;
}

void omxInitRFitFunction(omxFitFunction* oo) {
	FitContext::setRFitFunction(oo);

	if(OMX_DEBUG) { mxLog("Initializing R fit function."); }
	omxRFitFunction *newObj = (omxRFitFunction*) R_alloc(1, sizeof(omxRFitFunction));
	
	SEXP rObj = oo->rObj;

	/* Set Fit Function Calls to RFitFunction Calls */
	oo->computeFun = omxCallRFitFunction;
	oo->setFinalReturns = omxSetFinalReturnsRFitFunction;
	oo->destructFun = omxDestroyRFitFunction;
	oo->argStruct = (void*) newObj;
	
	PROTECT(newObj->fitfun = GET_SLOT(rObj, install("fitfun")));
	PROTECT_WITH_INDEX(newObj->model = GET_SLOT(rObj, install("model")), &(newObj->modelIndex));
	PROTECT(newObj->flatModel = GET_SLOT(rObj, install("flatModel")));
	PROTECT_WITH_INDEX(newObj->state = GET_SLOT(rObj, install("state")), &(newObj->stateIndex));

}


