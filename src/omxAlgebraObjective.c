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

#include "omxAlgebraFunctions.h"
#include "omxObjectiveTable.h"

#ifndef _OMX_ALGEBRA_OBJECTIVE_
#define _OMX_ALGEBRA_OBJECTIVE_ TRUE

typedef struct {

	omxMatrix *algebra;

} omxAlgebraObjective;

void omxDestroyAlgebraObjective(omxObjective *oo) {

}

void omxCallAlgebraObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
	if(OMX_DEBUG_ALGEBRA) {Rprintf("Beginning Algebra Objective Computation.\n");}
	omxMatrix* algebra = ((omxAlgebraObjective*)(oo->argStruct))->algebra;

	omxRecompute(algebra);
	
	// This should really be checked elsewhere.
	if(algebra->rows != 1 || algebra->cols != 1) {
		error("MxAlgebraObjective's objective algebra does not evaluate to a 1x1 matrix.");
	}
	
	oo->matrix->data[0] = algebra->data[0];
	
	if(OMX_DEBUG) {Rprintf("Algebra Objective value is %f.\n", oo->matrix->data[0]);}
}

omxRListElement* omxSetFinalReturnsAlgebraObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 1;
	omxRListElement* retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);

	return retVal;
}

void omxInitAlgebraObjective(omxObjective* oo, SEXP rObj) {
	
	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Initializing Algebra objective function.\n");
	}
	
	SEXP newptr;
	
	omxAlgebraObjective *newObj = (omxAlgebraObjective*) R_alloc(1, sizeof(omxAlgebraObjective));
	PROTECT(newptr = GET_SLOT(rObj, install("algebra")));
	newObj->algebra = omxNewMatrixFromMxIndex(newptr, oo->matrix->currentState);
	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Algebra Objective Bound to Algebra %d\n", newObj->algebra);
	}
	UNPROTECT(1);
	
	oo->objectiveFun = omxCallAlgebraObjective;
	oo->setFinalReturns = omxSetFinalReturnsAlgebraObjective;
	oo->destructFun = omxDestroyAlgebraObjective;
	oo->repopulateFun = NULL;
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_R_OBJECTIVE_ */
