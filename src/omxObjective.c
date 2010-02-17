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

void omxCalculateStdErrorFromHessian(int scale, omxObjective *oo) {
	/* This function calculates the standard errors from the hessian matrix */
	// sqrt(diag(solve(hessian)))

	if(oo->hessian == NULL) return;
	
	int numParams = oo->matrix->currentState->numFreeParams;
	
	if(oo->stdError == NULL) {
		oo->stdError = (double*) R_alloc(numParams, sizeof(double));
	}
	
	double* stdErr = oo->stdError;
	
	double* hessian = oo->hessian;
	double* workspace = (double *) Calloc(numParams * numParams, double);
	
	for(int i = 0; i < numParams; i++)
		for(int j = 0; j <= i; j++)
			workspace[i*numParams+j] = hessian[i*numParams+j];		// Populate upper triangle
	
	char u = 'U';
	int ipiv[numParams];
	int lwork = -1;
	double temp;
	int info = 0;
	
	F77_CALL(dsytrf)(&u, &numParams, workspace, &numParams, ipiv, &temp, &lwork, &info);
	
	lwork = (temp > numParams?temp:numParams);
	
	double* work = (double*) Calloc(lwork, double);
	
	F77_CALL(dsytrf)(&u, &numParams, workspace, &numParams, ipiv, work, &lwork, &info);
	
	if(info != 0) {
		
		oo->stdError = NULL;
		
	} else {
		
		F77_CALL(dsytri)(&u, &numParams, workspace, &numParams, ipiv, work, &info);
	
		if(info != 0) {
			oo->stdError = NULL;
		} else {
			for(int i = 0; i < numParams; i++) {
				stdErr[i] = scale * sqrt(workspace[i * numParams + i]);
			}
		}
	}
	
	Free(workspace);
	Free(work);
	
}

void omxInitEmptyObjective(omxObjective *oo) {
	/* Sets everything to NULL to avoid bad pointer calls */
	
	oo->initFun = NULL;
	oo->destructFun = NULL;
	oo->repopulateFun = NULL;
	oo->objectiveFun = NULL;
	oo->needsUpdateFun = NULL;
	oo->getStandardErrorFun = NULL;
	oo->setFinalReturns = NULL;
	oo->gradientFun = NULL;
	oo->argStruct = NULL;
	oo->objType = (char*) calloc(251, sizeof(char*));
	oo->objType[0] = '\0';
	oo->matrix = NULL;
	oo->hessian = NULL;
	
}

void omxFreeObjectiveArgs(omxObjective *oo) {
	/* Completely destroy the objective function tree */
	free(oo->objType);
	oo->destructFun(oo);
	oo->matrix = NULL;
	
}

void omxObjectiveCompute(omxObjective *oo) {
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Objective compute: 0x%0x (needed: %s).\n", oo, (oo->matrix->isDirty?"Yes":"No"));}

	oo->objectiveFun(oo);

	if(oo->matrix != NULL)
		omxMatrixCompute(oo->matrix);
}

unsigned short omxObjectiveNeedsUpdate(omxObjective *oo)
{
	if(OMX_DEBUG_MATRIX) {Rprintf("omxObjectiveNeedsUpdate:");}
	unsigned short needsIt = TRUE;   		// Defaults to TRUE if unspecified
	if(!(oo->needsUpdateFun == NULL)) {
		if(OMX_DEBUG_MATRIX) {Rprintf("Calling update function 0x%x:", oo->needsUpdateFun);}
		needsIt = oo->needsUpdateFun(oo);
	}
	
	if(OMX_DEBUG_MATRIX) {Rprintf("%s\n", (needsIt?"Yes":"No"));}
	
	return needsIt;
}


void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP rObj) {

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
	strncpy(obj->objType, objType, 250);
	
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

	obj->initFun(obj, rObj);

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

omxMatrix* omxNewMatrixFromIndexSlot(SEXP rObj, omxState* currentState, char* const slotName) {
	SEXP slotValue;
	omxMatrix* newMatrix = NULL;
	if(strncmp(slotName, "", 1) == 0) return NULL;
	PROTECT(slotValue = GET_SLOT(rObj, install(slotName)));
	newMatrix = omxNewMatrixFromMxIndex(slotValue, currentState);
	if(newMatrix != NULL) omxRecompute(newMatrix);
	else if(OMX_DEBUG) Rprintf("No M found.\n");
	UNPROTECT(1);
	return newMatrix;
}

omxData* omxNewDataFromDataSlot(SEXP rObj, omxState* currentState, char* const dataSlotName) {
	
	SEXP slotValue;
	
	PROTECT(slotValue = GET_SLOT(rObj, install(dataSlotName)));
	if(OMX_DEBUG) { Rprintf("Data Element %d.\n", AS_INTEGER(slotValue)); }
	omxData* dataElt = omxNewDataFromMxDataPtr(slotValue, currentState);
	UNPROTECT(1); // newMatrix
	
	return dataElt;
	
}
