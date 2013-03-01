/*
 *  Copyright 2007-2013 The OpenMx Project
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
*  omxFitFunction.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	FitFunction objects are a subclass of data matrix that evaluates
*   itself anew at each iteration, so that any changes to
*   free parameters can be incorporated into the update.
*   // Question: Should FitFunction be a ``subtype'' of 
*   // omxAlgebra or a separate beast entirely?
*
**********************************************************/

#include "omxFitFunction.h"

typedef struct omxFitFunctionTableEntry omxFitFunctionTableEntry;

struct omxFitFunctionTableEntry {

	char name[32];
	void (*initFun)(omxFitFunction*, SEXP);

};

extern void omxInitAlgebraFitFunction(omxFitFunction *off, SEXP rObj);
extern void omxInitWLSFitFunction(omxFitFunction *off, SEXP rObj);
extern void omxInitRowFitFunction(omxFitFunction *off, SEXP rObj);
extern void omxInitMLFitFunction(omxFitFunction *off, SEXP rObj);
extern void omxInitRFitFunction(omxFitFunction *off, SEXP rObj);

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction},
	{"MxFitFunctionML", 				&omxInitMLFitFunction},
	{"MxFitFunctionR",					&omxInitRFitFunction},
	{"", 0}
};

void omxCalculateStdErrorFromHessian(double scale, omxFitFunction *off) {
	/* This function calculates the standard errors from the hessian matrix */
	// sqrt(diag(solve(hessian)))

	if(off->hessian == NULL) return;
	
	int numParams = off->matrix->currentState->numFreeParams;
	
	if(off->stdError == NULL) {
		off->stdError = (double*) R_alloc(numParams, sizeof(double));
	}
	
	double* stdErr = off->stdError;
	
	double* hessian = off->hessian;
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
		
		off->stdError = NULL;
		
	} else {
		
		F77_CALL(dsytri)(&u, &numParams, workspace, &numParams, ipiv, work, &info);
	
		if(info != 0) {
			off->stdError = NULL;
		} else {
			for(int i = 0; i < numParams; i++) {
				stdErr[i] = sqrt(scale) * sqrt(workspace[i * numParams + i]);
			}
		}
	}
	
	Free(workspace);
	Free(work);
	
}

void omxInitEmptyFitFunction(omxFitFunction *off) {
	/* Sets everything to NULL to avoid bad pointer calls */
	
	memset(off, 0, sizeof(omxFitFunction));
}

void omxFreeFitFunctionArgs(omxFitFunction *off) {
	if(off==NULL) return;
    
	/* Completely destroy the fit function structures */
	if(OMX_DEBUG) {Rprintf("Freeing fit function object at 0x%x.\n", off);}
	if(off->matrix != NULL) {
		if(off->destructFun != NULL) {
			if(OMX_DEBUG) {Rprintf("Calling fit function destructor for 0x%x.\n", off);}
			off->destructFun(off);
		}
		off->matrix = NULL;
	}
}

void omxFitFunctionCreateChildren(omxState *globalState, int numThreads)
{
	if (numThreads <= 1) return;

	omxMatrix *fm = globalState->fitMatrix;
	if (!fm) return;

	omxFitFunction *ff = fm->fitFunction;
	if (!ff->usesChildModels) return;

	globalState->numChildren = numThreads;

	globalState->childList = (omxState**) Calloc(numThreads, omxState*);

	for(int ii = 0; ii < numThreads; ii++) {
		globalState->childList[ii] = (omxState*) R_alloc(1, sizeof(omxState));
		omxInitState(globalState->childList[ii], globalState);
		omxDuplicateState(globalState->childList[ii], globalState);
	}
}

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;
	if(src->fitFunction == NULL) return;
    
	omxFillMatrixFromMxFitFunction(tgt, src->fitFunction->rObj, src->hasMatrixNumber, src->matrixNumber);

}

omxFitFunction* omxCreateDuplicateFitFunction(omxFitFunction *tgt, const omxFitFunction *src, omxState* newState) {

	if(OMX_DEBUG) {Rprintf("Duplicating fit function 0x%x into 0x%x.", src, tgt);}

	if(src == NULL) {
		return NULL;
	}
	
	if(tgt == NULL) {
        tgt = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
        omxInitEmptyFitFunction(tgt);
    } else {
		omxRaiseError(newState, -1,
			"omxCreateDuplicateFitFunction requested to overwrite target");
		return NULL;
	}

	memcpy(tgt, src, sizeof(omxFitFunction));
	return tgt;

}

void omxFitFunctionCompute(omxFitFunction *off, int want, double* gradient) {
	if(OMX_DEBUG_ALGEBRA) { 
	    Rprintf("FitFunction compute: 0x%0x (needed: %s).\n", off, (off->matrix->isDirty?"Yes":"No"));
	}

	off->computeFun(off, want, gradient);

	omxMarkClean(off->matrix);
}

void omxFillMatrixFromMxFitFunction(omxMatrix* om, SEXP rObj,
	unsigned short hasMatrixNumber, int matrixNumber) {

	SEXP slotValue, fitFunctionClass;
	omxFitFunction *obj = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	omxInitEmptyFitFunction(obj);

	/* Register FitFunction and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// FitFunction matrices MUST be 1x1.
	om->fitFunction = obj;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;
	
	/* Get FitFunction Type */
	PROTECT(fitFunctionClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	{
	  const char *fitType = CHAR(fitFunctionClass);
	
	  /* Switch based on fit function type. */ 
	  const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable;
	  while (entry->initFun) {
	    if(strncmp(fitType, entry->name, MAX_STRING_LEN) == 0) {
	      obj->fitType = entry->name;
	      obj->initFun = entry->initFun;
	      break;
	    }
	    entry += 1;
	  }

	  if(obj->initFun == NULL) {
	    const int MaxErrorLen = 256;
	    char newError[MaxErrorLen];
	    snprintf(newError, MaxErrorLen, "Fit function %s not implemented.\n", fitType);
	    omxRaiseError(om->currentState, -1, newError);
	    return;
	  }
	}
	UNPROTECT(1);	/* fitType */

	PROTECT(slotValue = GET_SLOT(rObj, install("expectation")));
	if (LENGTH(slotValue) != 1) {
	    const int MaxErrorLen = 256;
	    char newError[MaxErrorLen];
	    snprintf(newError, MaxErrorLen, "Fit function %s expectation improperly initialized\n", obj->fitType);
	    error(newError);
	}
	int expNumber = INTEGER(slotValue)[0];	
	if(expNumber == NA_INTEGER) {						// Has no expectation associated with it
		obj->expectation = NULL;
	} else {
		obj->expectation = omxNewExpectationFromExpectationIndex(expNumber, om->currentState);
	}
	UNPROTECT(1);	/* slotValue */
	
	if (om->currentState->statusMsg[0]) return;

	obj->rObj = rObj;
	obj->initFun(obj, rObj);

	if(obj->computeFun == NULL) {// If initialization fails, error code goes in argStruct
		char *errorCode;
		if(om->currentState->statusCode != 0) {
			errorCode = om->currentState->statusMsg;
		} else {
			// If no error code is reported, we report that.
  			errorCode = "No error code reported.";
		}
		if(obj->argStruct != NULL) {
			errorCode = (char*)(obj->argStruct);
		}
        const int MaxErrorLen = 256;
        char newError[MaxErrorLen];
        snprintf(newError, MaxErrorLen, "Could not initialize fit function %s.  Error: %s\n",
			obj->fitType, errorCode); 
		omxRaiseError(om->currentState, -1, newError);
	}
	
	obj->matrix->isDirty = TRUE;

}

void omxFitFunctionPrint(omxFitFunction* off, char* d) {
	Rprintf("(FitFunction, type %s) ", off->fitType);
	omxPrintMatrix(off->matrix, d);
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
