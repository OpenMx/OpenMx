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

/***********************************************************
* 
*  omxObjective.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Objective objects are a subclass of data matrix that evaluates
*   itself anew at each iteration, so that any changes to
*   free parameters can be incorporated into the update.
*
**********************************************************/

#include "omxMatrix.h"
#include "omxObjective.h"
#include "omxRObjective.h"
#include "omxRAMObjective.h"
#include "omxFIMLObjective.h"
#include "omxAlgebraObjective.h"\

/* Need a better way to deal with these. */
extern omxMatrix** algebraList;
extern omxMatrix** matrixList;

void omxFreeObjectiveArgs(omxObjective *oo) {
	/* Completely destroy the objective function tree */

	oo->destructFun(oo);
	oo->myMatrix = NULL;
	
}

void omxObjectiveCompute(omxObjective *oo) {
	if(OMX_DEBUG) { Rprintf("Objective compute: 0x%0x (needed: %s).\n", oo, (oo->myMatrix->isDirty?"Yes":"No"));}

	oo->objectiveFun(oo);

	oo->myMatrix->isDirty = FALSE;

}

unsigned short omxObjectiveNeedsUpdate(omxObjective *oo)
{
	if(OMX_DEBUG) {Rprintf("omxObjectiveNeedsUpdate:");}
	unsigned short needsIt = TRUE;
	if(!(oo->needsUpdateFun == NULL)) {
		needsIt = oo->needsUpdateFun(oo);
	}
	
	if(OMX_DEBUG) {Rprintf("%s\n", (needsIt?"Yes":"No"));}
	
	return needsIt;

}


void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP rObj, SEXP dataList) {

	const char *objType;
	SEXP objectiveClass;
	omxObjective *obj = (omxObjective*) R_alloc(sizeof(omxObjective), 1);

	/* Register Objective and Matrix with each other */
	obj->myMatrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// Objective matrices MUST be 1x1.
	om->objective = obj;
	
	/* Default NeedsUpdate is NULL */
	obj->needsUpdateFun = NULL;
	
	/* Get Objective Type */
	PROTECT(objectiveClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	objType = CHAR(objectiveClass);
	obj->objType[250] = '\0';
	strncpy(obj->objType, objType, 249);
	
	/* Switch based on objective type. */  // Right now, this is hard-wired.  // TODO: Replace with hash function and table lookup.
	if(strncmp(objType, "MxRAMObjective", 21) == 0) { // Covariance-style optimization.
		obj->initFun = omxInitRAMObjective;
		obj->objectiveFun = omxCallRAMObjective;
		obj->destructFun = omxDestroyRAMObjective;
	} else if(strncmp(objType, "MxFIMLObjective", 15) == 0) {
		obj->initFun = omxInitFIMLObjective;
		obj->objectiveFun = omxCallFIMLObjective;
		obj->destructFun = omxDestroyFIMLObjective;
	} else if(strncmp(objType, "MxAlgebraObjective", 18) == 0) {
		obj->initFun = omxInitAlgebraObjective;
		obj->objectiveFun = omxCallAlgebraObjective;
		obj->destructFun = omxDestroyAlgebraObjective;
	} else if(strncmp(objType, "MxRObjective", 12) == 0) {
		obj->initFun = omxInitRObjective;
		obj->objectiveFun = omxCallRObjective;
		obj->destructFun = omxDestroyRObjective;
	} else {
		error("Objective function type %s not implemented on this kernel.", objType);
	}

	obj->initFun(obj, rObj, dataList);
	
	obj->myMatrix->isDirty = TRUE;

	UNPROTECT(1);	/* objectiveClass */

}

void omxObjectiveGradient(omxObjective* oo, double* gradient) {
	if(!(oo->gradientFun == NULL)) { oo->gradientFun(oo, gradient); }
	return;
}

void omxObjectivePrint(omxObjective* oo, char* d) {
	Rprintf("(Objective, type %s) ", oo->objType);
	omxMatrixPrint(oo->myMatrix, d);
}
