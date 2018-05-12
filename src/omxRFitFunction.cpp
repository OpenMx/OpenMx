/*
 *  Copyright 2007-2018 by the individuals mentioned in the source code history
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

#include "glue.h"
#include "Compute.h"
#include "omxRFitFunction.h"
#include "EnableWarnings.h"

void omxRFitFunction::compute(int want, FitContext *)
{
	auto *oo = this;
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;

	omxRFitFunction* rFitFunction = this;

	SEXP theCall;
	ScopedProtect p2(theCall, Rf_allocVector(LANGSXP, 3));
	SETCAR(theCall, rFitFunction->fitfun);
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->state);

	ProtectedSEXP theReturn(Rf_eval(theCall, R_GlobalEnv));

	if (LENGTH(theReturn) < 1) {
		// seems impossible, but report it if it happens
		omxRaiseErrorf("FitFunction returned nothing");
	} else if (LENGTH(theReturn) == 1) {
		oo->matrix->data[0] = Rf_asReal(theReturn);
	} else if (LENGTH(theReturn) == 2) {
		oo->matrix->data[0] = Rf_asReal(VECTOR_ELT(theReturn, 0));
		rFitFunction->state = VECTOR_ELT(theReturn, 1);
		Rf_setAttrib(rObj, Rf_install("state"), rFitFunction->state); //protect it
	} else if (LENGTH(theReturn) > 2) {
		omxRaiseErrorf("FitFunction returned more than 2 arguments");
	}
}

omxFitFunction *omxInitRFitFunction()
{ return new omxRFitFunction; }

void omxRFitFunction::init()
{
	FitContext::setRFitFunction(this);

	if(OMX_DEBUG) { mxLog("Initializing R fit function."); }
	omxRFitFunction *newObj = this;
	
	
	ProtectedSEXP Runit(R_do_slot(rObj, Rf_install("units")));
	setUnitsFromName(CHAR(STRING_ELT(Runit, 0)));

	newObj->fitfun = R_do_slot(rObj, Rf_install("fitfun"));
	newObj->model = R_do_slot(rObj, Rf_install("model"));
	newObj->flatModel = R_do_slot(rObj, Rf_install("flatModel"));
	newObj->state = R_do_slot(rObj, Rf_install("state"));
}
