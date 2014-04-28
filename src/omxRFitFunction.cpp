/*
 *  Copyright 2007-2014 The OpenMx Project
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

#include "omxAlgebraFunctions.h"
#include "omxRFitFunction.h"
#include "omxOpenmpWrap.h"
#include "glue.h"
#include "Compute.h"

void omxDestroyRFitFunction(omxFitFunction *off) {
	Rf_unprotect(4); 			// fitfun, model, flatModel, and state
}

static void omxCallRFitFunction(omxFitFunction *oo, int want, FitContext *) {
	if (want & (FF_COMPUTE_PREOPTIMIZE)) return;

	omxRFitFunction* rFitFunction = (omxRFitFunction*)oo->argStruct;

	SEXP theCall, theReturn;
	Rf_protect(theCall = Rf_allocVector(LANGSXP, 3));
	SETCAR(theCall, rFitFunction->fitfun);
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->state);

	Rf_protect(theReturn = Rf_eval(theCall, R_GlobalEnv));

	if (LENGTH(theReturn) < 1) {
		// seems impossible, but report it if it happens
		omxRaiseErrorf("FitFunction returned nothing");
	} else if (LENGTH(theReturn) == 1) {
		oo->matrix->data[0] = Rf_asReal(theReturn);
	} else if (LENGTH(theReturn) == 2) {
		oo->matrix->data[0] = Rf_asReal(VECTOR_ELT(theReturn, 0));
		R_Reprotect(rFitFunction->state = VECTOR_ELT(theReturn, 1), rFitFunction->stateIndex);
	} else if (LENGTH(theReturn) > 2) {
		omxRaiseErrorf("FitFunction returned more than 2 arguments");
	}

	Rf_unprotect(2); // theCall and theReturn
}

void omxInitRFitFunction(omxFitFunction* oo) {
	FitContext::setRFitFunction(oo);

	if(OMX_DEBUG) { mxLog("Initializing R fit function."); }
	omxRFitFunction *newObj = (omxRFitFunction*) R_alloc(1, sizeof(omxRFitFunction));
	
	SEXP rObj = oo->rObj;

	/* Set Fit Function Calls to RFitFunction Calls */
	oo->computeFun = omxCallRFitFunction;
	oo->destructFun = omxDestroyRFitFunction;
	oo->argStruct = (void*) newObj;
	
	Rf_protect(newObj->fitfun = R_do_slot(rObj, Rf_install("fitfun")));
	R_ProtectWithIndex(newObj->model = R_do_slot(rObj, Rf_install("model")), &(newObj->modelIndex));
	Rf_protect(newObj->flatModel = R_do_slot(rObj, Rf_install("flatModel")));
	R_ProtectWithIndex(newObj->state = R_do_slot(rObj, Rf_install("state")), &(newObj->stateIndex));

}


