/*
 *  Copyright 2007-2009 The OpenMx Project
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

	omxMatrix* algebra = ((omxAlgebraObjective*)(oo->argStruct))->algebra;

	omxRecomputeMatrix(algebra);
	
	// This should really be checked elsewhere.
	if(algebra->rows != 1 || algebra->cols != 1) {
		error("MxAlgebraObjective's objective algebra does not evaluate to a 1x1 matrix.");
	}
	
	oo->myMatrix->data[0] = algebra->data[0];
	
}

unsigned short int omxNeedsUpdateAlgebraObjective(omxObjective *oo) {

	if(oo->myMatrix->data[0] != ((omxAlgebraObjective*)oo->argStruct)->algebra->data[0]) return TRUE;
	return omxNeedsUpdate(((omxAlgebraObjective*)oo->argStruct)->algebra);
}

void omxInitAlgebraObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	if(OMX_DEBUG) { Rprintf("Initializing Algebra objective function.\n"); }
	
	SEXP newptr;
	
	omxAlgebraObjective *newObj = (omxAlgebraObjective*) R_alloc(1, sizeof(omxAlgebraObjective));
	PROTECT(newptr = GET_SLOT(rObj, install("algebra")));
	newObj->algebra = omxNewMatrixFromMxIndex(newptr);
	if(OMX_DEBUG) {Rprintf("Algebra Objective Bound to Algebra %d", newObj->algebra);}
	UNPROTECT(1);
	
	oo->objectiveFun = omxCallAlgebraObjective;
	oo->needsUpdateFun = omxNeedsUpdateAlgebraObjective;
	oo->destructFun = omxDestroyAlgebraObjective;
	oo->repopulateFun = NULL;
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_R_OBJECTIVE_ */
