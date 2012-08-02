/*
 *  Copyright 2007-2012 The OpenMx Project
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

void omxCalculateStdErrorFromHessian(double scale, omxObjective *oo) {
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
				stdErr[i] = sqrt(scale) * sqrt(workspace[i * numParams + i]);
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
	oo->populateAttrFun = NULL;
	oo->setFinalReturns = NULL;
	oo->gradientFun = NULL;
	oo->sharedArgs = NULL;
	oo->argStruct = NULL;
	oo->subObjective = NULL;
	oo->rObj = NULL;
	oo->objType = Calloc(MAX_STRING_LEN, char);
	oo->matrix = NULL;
	oo->stdError = NULL;
	oo->hessian = NULL;
	oo->gradient = NULL;
	oo->isPrepopulated = FALSE;
}

omxObjective* omxCreateSubObjective(omxObjective *oo) {

	if(OMX_DEBUG) {Rprintf("Creating SubObjective Object....\n");}
    if(oo == NULL) {
	 	if(OMX_DEBUG) {Rprintf("Got Null objective.  Returning.");}
		return NULL;
	}
    omxObjective* subObjective = (omxObjective*) Calloc(1, omxObjective);
    omxInitEmptyObjective(subObjective);
	omxCreateDuplicateObjective(subObjective, oo, oo->matrix->currentState);
	oo->subObjective = subObjective;
	
    return subObjective;

}

void omxFreeObjectiveArgs(omxObjective *oo) {
	if(oo==NULL) return;
    
	/* Completely destroy the objective function tree */
	if(OMX_DEBUG) {Rprintf("Freeing objective object at 0x%x with subobjective 0x%x.\n", oo, oo->subObjective);}
	if(oo->matrix != NULL) {
		if(oo->objType != NULL) {
			Free(oo->objType);
		}
		if(oo->subObjective != NULL) {
			omxFreeObjectiveArgs(oo->subObjective);
			Free(oo->subObjective);
		}
		if(oo->destructFun != NULL) {
			if(OMX_DEBUG) {Rprintf("Calling objective destructor for 0x%x.\n", oo);}
			oo->destructFun(oo);
		}
		oo->matrix = NULL;
	}
}

void omxObjectiveCompute(omxObjective *oo) {
	if(OMX_DEBUG_ALGEBRA) { 
	    Rprintf("Objective compute: 0x%0x (needed: %s).\n", oo, (oo->matrix->isDirty?"Yes":"No"));
	}

	oo->objectiveFun(oo);

	omxMarkClean(oo->matrix);

}

void omxDuplicateObjectiveMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;
	if(src->objective == NULL) return;
    
	omxFillMatrixFromMxObjective(tgt, src->objective->rObj, src->hasMatrixNumber, src->matrixNumber);

}

omxObjective* omxCreateDuplicateObjective(omxObjective *tgt, const omxObjective *src, omxState* newState) {

	if(OMX_DEBUG) {Rprintf("Duplicating objective 0x%x into 0x%x.", src, tgt);}

	if(src == NULL) {
		return NULL;
	}
	
	if(tgt == NULL) {
        tgt = (omxObjective*) R_alloc(1, sizeof(omxObjective));
        omxInitEmptyObjective(tgt);
    }

	tgt->initFun 					= src->initFun;
	tgt->destructFun 				= src->destructFun;
	tgt->repopulateFun 				= src->repopulateFun;
	tgt->objectiveFun 				= src->objectiveFun;
	tgt->needsUpdateFun				= src->needsUpdateFun;
	tgt->getStandardErrorFun 		= src->getStandardErrorFun;
	tgt->populateAttrFun 			= src->populateAttrFun;
	tgt->setFinalReturns 			= src->setFinalReturns;
	tgt->gradientFun 				= src->gradientFun;
	tgt->sharedArgs					= src->sharedArgs;
	tgt->matrix						= src->matrix;
	tgt->rObj						= src->rObj;
	tgt->subObjective	 			= src->subObjective;
	tgt->stdError					= src->stdError;
	tgt->hessian 					= src->hessian;
	tgt->gradient 					= src->gradient;

    strncpy(tgt->objType, src->objType, MAX_STRING_LEN);
	return tgt;

}

unsigned short omxObjectiveNeedsUpdate(omxObjective *oo)
{
	if(OMX_DEBUG_MATRIX) { Rprintf("omxObjectiveNeedsUpdate:"); }
	unsigned short needsIt = TRUE;     // Defaults to TRUE if unspecified
	if(!(oo->needsUpdateFun == NULL)) {
		if(OMX_DEBUG_MATRIX) {Rprintf("Calling update function 0x%x:", oo->needsUpdateFun);}
		needsIt = oo->needsUpdateFun(oo);
		if(!needsIt && !(oo->subObjective == NULL)) {
			needsIt = omxObjectiveNeedsUpdate(oo->subObjective);
		}
	} else if(!(oo->subObjective == NULL)) {
		needsIt = omxObjectiveNeedsUpdate(oo->subObjective);
	}
	
	if(OMX_DEBUG_MATRIX) {Rprintf("%s\n", (needsIt?"Yes":"No"));}
	
	return needsIt;
}

void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP rObj,
	unsigned short hasMatrixNumber, int matrixNumber) {

	int i;
	const char *objType;
	SEXP objectiveClass;
	char errorCode[MAX_STRING_LEN];
	omxObjective *obj = (omxObjective*) R_alloc(1, sizeof(omxObjective));
	omxInitEmptyObjective(obj);

	/* Register Objective and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// Objective matrices MUST be 1x1.
	om->objective = obj;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;
	
	/* Get Objective Type */
	PROTECT(objectiveClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	objType = CHAR(objectiveClass);
	obj->objType[MAX_STRING_LEN - 1] = '\0';
	strncpy(obj->objType, objType, MAX_STRING_LEN);
	
	/* Switch based on objective type. */ 
	for(i = 0; i < omxObjectiveTableLength; i++) {
		if(strncmp(objType, omxObjectiveSymbolTable[i].name, MAX_STRING_LEN) == 0) {
			obj->initFun = omxObjectiveSymbolTable[i].initFun;
			break;
		}
	}

	if(i == omxObjectiveTableLength) {
		char newError[MAX_STRING_LEN];
		sprintf(newError, "Objective function %s not implemented.\n", obj->objType);
		omxRaiseError(om->currentState, -1, newError);
	}
	
	obj->rObj = rObj;
	obj->initFun(obj, rObj);

	if(obj->objectiveFun == NULL) {// If initialization fails, error code goes in argStruct
		if(om->currentState->statusCode != 0) {
			strncpy(errorCode, om->currentState->statusMsg, 150); // Report a status error
		} else {
			// If no error code is reported, we report that.
  			strncpy(errorCode, "No error code reported.", 25);
		}
		if(obj->argStruct != NULL) {
			strncpy(errorCode, (char*)(obj->argStruct), 51);
		}
		char newError[MAX_STRING_LEN];
		sprintf(newError, "Could not initialize objective function %s.  Error: %s\n", 
		    obj->objType, errorCode);
		omxRaiseError(om->currentState, -1, newError);
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


/* Helper functions */
omxMatrix* omxNewMatrixFromIndexSlot(SEXP rObj, omxState* currentState, char* const slotName) {
	SEXP slotValue;
	omxMatrix* newMatrix = NULL;
	if(strncmp(slotName, "", 1) == 0) return NULL;
	PROTECT(slotValue = GET_SLOT(rObj, install(slotName)));
	newMatrix = omxNewMatrixFromMxIndex(slotValue, currentState);
	if(newMatrix != NULL) omxRecompute(newMatrix);
	else if(OMX_DEBUG) Rprintf("No slot %s found.\n", slotName);
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
