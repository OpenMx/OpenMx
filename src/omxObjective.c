/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
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
*   // Question: Should Objective be a ``subtype'' of 
*   // omxAlgebra or a separate beast entirely?
*
**********************************************************/

#include "omxObjective.h"

void omxInitEmptyObjective(omxObjective *oo) {
	/* Sets everything to NULL to avoid bad pointer calls */
	
	oo->initFun = NULL;
	oo->destructFun = NULL;
	oo->repopulateFun = NULL;
	oo->objectiveFun = NULL;
	oo->needsUpdateFun = NULL;
	oo->gradientFun = NULL;
	oo->argStruct = NULL;
	strncpy(oo->objType, "\0", 1);
	oo->matrix = NULL;
	
}

void omxFreeObjectiveArgs(omxObjective *oo) {
	/* Completely destroy the objective function tree */

	oo->destructFun(oo);
	oo->matrix = NULL;
	
}

void omxObjectiveCompute(omxObjective *oo) {
	if(OMX_DEBUG) { Rprintf("Objective compute: 0x%0x (needed: %s).\n", oo, (oo->matrix->isDirty?"Yes":"No"));}

	oo->objectiveFun(oo);

	oo->matrix->isDirty = FALSE;

}

unsigned short omxObjectiveNeedsUpdate(omxObjective *oo)
{
	if(OMX_DEBUG) {Rprintf("omxObjectiveNeedsUpdate:");}
	unsigned short needsIt = TRUE;   		// Defaults to TRUE if unspecified
	if(!(oo->needsUpdateFun == NULL)) {
		needsIt = oo->needsUpdateFun(oo);
	}
	
	if(OMX_DEBUG) {Rprintf("%s\n", (needsIt?"Yes":"No"));}
	
	return needsIt;

}


void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP rObj, SEXP dataList) {

	int i;
	const char *objType;
	SEXP objectiveClass;
	char errorCode[51];
	omxObjective *obj = (omxObjective*) R_alloc(1, sizeof(omxObjective));
	omxInitEmptyObjective(obj);

	/* Register Objective and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// Objective matrices MUST be 1x1.
	om->objective = obj;
	
	/* Get Objective Type */
	PROTECT(objectiveClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	objType = CHAR(objectiveClass);
	obj->objType[250] = '\0';
	strncpy(obj->objType, objType, 249);
	
	/* Switch based on objective type. */ 
	for(i = 0; i < omxObjectiveTableLength; i++) {
		if(strncmp(objType, omxObjectiveSymbolTable[i].name, 250) == 0) {
			obj->initFun = omxObjectiveSymbolTable[i].initFun;
			break;
		}
	}

	if(i == omxObjectiveTableLength) { 
		error("Objective type %s not supported.\n", obj->objType);
	}

	obj->initFun(obj, rObj, dataList);

	if(obj->objectiveFun == NULL) {								// If initialization fails, error code goes in argStruct
  		strncpy(errorCode, "No error code reported.", 25);		// If no error code is reported, we report that.
		if(obj->argStruct != NULL) {
			strncpy(errorCode, (char*)(obj->argStruct), 51);
		}
		error("Could not initialize objective function %s.  Error: %s\n", obj->objType, errorCode);
	}
	
	obj->matrix->isDirty = TRUE;

	UNPROTECT(1);	/* objectiveClass */

}

void omxObjectiveGradient(omxObjective* oo, double* gradient) {
	if(!(oo->gradientFun == NULL)) { oo->gradientFun(oo, gradient); }
	return;
}

void omxObjectivePrint(omxObjective* oo, char* d) {
	Rprintf("(Objective, type %s) ", oo->objType);
	omxPrintMatrix(oo->matrix, d);
}
