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

#ifndef _OMX_ALGEBRA_FITFUNCTION_
#define _OMX_ALGEBRA_FITFUNCTION_ TRUE

typedef struct {

	omxMatrix *algebra;

} omxAlgebraFitFunction;

void omxDestroyAlgebraFitFunction(omxFitFunction *off) {

}

static void omxCallAlgebraFitFunction(omxFitFunction *off, int want, FitContext *) {
	if(OMX_DEBUG_ALGEBRA) {mxLog("Beginning Algebra Fit Function Computation.");}
	omxMatrix* algebra = ((omxAlgebraFitFunction*)(off->argStruct))->algebra;

	omxRecompute(algebra);
	
	// This should really be checked elsewhere. TODO
	if(algebra->rows != 1 || algebra->cols != 1) {
		Rf_error("MxAlgebraFitFunction's fit function algebra does not evaluate to a 1x1 matrix.");
	}
	
	off->matrix->data[0] = algebra->data[0];
	
	if(OMX_DEBUG) {mxLog("Algebra Fit Function value is %f.", off->matrix->data[0]);}
}

void omxInitAlgebraFitFunction(omxFitFunction* off) {
	
	if(OMX_DEBUG) {
		mxLog("Initializing Algebra fitFunction function.");
	}
	
	SEXP rObj = off->rObj;
	SEXP newptr;
	
	omxAlgebraFitFunction *newObj = (omxAlgebraFitFunction*) R_alloc(1, sizeof(omxAlgebraFitFunction));
	Rf_protect(newptr = R_do_slot(rObj, Rf_install("algebra")));
	newObj->algebra = omxMatrixLookupFromState1(newptr, off->matrix->currentState);
	
	off->computeFun = omxCallAlgebraFitFunction;
	off->destructFun = omxDestroyAlgebraFitFunction;
	
	off->argStruct = (void*) newObj;
}


#endif /* _OMX_ALGEBRA_FITFUNCTION_ */
