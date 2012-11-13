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

typedef struct omxObjectiveTableEntry omxObjectiveTableEntry;

struct omxObjectiveTableEntry {

	char name[32];
	void (*initFun)(omxObjective*, SEXP) ;

};

extern void omxInitAlgebraObjective(omxObjective *oo, SEXP rObj);
extern void omxInitFIMLObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRAMObjective(omxObjective *oo, SEXP rObj);
extern void omxInitLISRELObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRowObjective(omxObjective *oo, SEXP rObj);
extern void omxInitMLObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRObjective(omxObjective *oo, SEXP rObj);
extern void omxInitWLSObjective(omxObjective *oo, SEXP rObj);

static const omxObjectiveTableEntry omxObjectiveSymbolTable[] = {
	{"MxAlgebraObjective", 			&omxInitAlgebraObjective},
	{"MxFIMLObjective",				&omxInitFIMLObjective},
	{"MxRAMObjective", 				&omxInitRAMObjective},
	{"MxWLSObjective",				&omxInitWLSObjective},
	{"MxRowObjective", 				&omxInitRowObjective},
	{"MxMLObjective", 				&omxInitMLObjective},
	{"MxRObjective",				&omxInitRObjective},
	{"MxLISRELObjective",			&omxInitLISRELObjective},
	{"", 0}
};

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
	
  memset(oo, 0, sizeof(omxObjective));
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
	} else if (tgt->objType) {
	  omxRaiseError(newState, -1, "omxCreateDuplicateObjective requested to overwrite target");
	}

	memcpy(tgt, src, sizeof(omxObjective));

	return tgt;

}

void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP rObj,
	unsigned short hasMatrixNumber, int matrixNumber) {

	omxObjective *obj = (omxObjective*) R_alloc(1, sizeof(omxObjective));
	omxInitEmptyObjective(obj);

	/* Register Objective and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// Objective matrices MUST be 1x1.
	om->objective = obj;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;
	
	/* Get Objective Type */
	{
	SEXP objectiveClass;
	PROTECT(objectiveClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	const char *objType = CHAR(objectiveClass);
	
	/* Switch based on objective type. */ 
	const omxObjectiveTableEntry *entry = omxObjectiveSymbolTable;
	while (entry->initFun) {
		if(strncmp(objType, entry->name, MAX_STRING_LEN) == 0) {
			obj->objType = entry->name;
			obj->initFun = entry->initFun;
			break;
		}
		entry += 1;
	}

	UNPROTECT(1);	/* objectiveClass */
	}

	if (!obj->initFun) {
		const int MaxErrorLen = 256;
		char newError[MaxErrorLen];
		snprintf(newError, MaxErrorLen, "Objective function %s not implemented.\n", obj->objType);
		omxRaiseError(om->currentState, -1, newError);
		return;
	}
	
	obj->rObj = rObj;
	obj->initFun(obj, rObj);

	if(obj->objectiveFun == NULL) {// If initialization fails, error code goes in argStruct
	        char *errorCode;
		if(om->currentState->statusCode != 0) {
			errorCode = om->currentState->statusMsg; // Report a status error
		} else {
			// If no error code is reported, we report that.
		        errorCode = "No error code reported.";
		}
		if(obj->argStruct != NULL) {
		        errorCode = (char*) obj->argStruct;
		}
		const int MaxErrorLen = 256;
		char newError[MaxErrorLen];
		snprintf(newError, MaxErrorLen, "Could not initialize objective function %s.  Error: %s\n", 
		    obj->objType, errorCode);
		omxRaiseError(om->currentState, -1, newError);
	}
	
	obj->matrix->isDirty = TRUE;


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
