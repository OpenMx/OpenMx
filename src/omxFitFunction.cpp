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
#include "omxOptimizer.h"
#include "fitMultigroup.h"

typedef struct omxFitFunctionTableEntry omxFitFunctionTableEntry;

struct omxFitFunctionTableEntry {

	char name[32];
	void (*initFun)(omxFitFunction*);

};

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction},
	{"MxFitFunctionML", 				&omxInitMLFitFunction},
	{"imxFitFunctionFIML", &omxInitFIMLFitFunction},
	{"MxFitFunctionR",					&omxInitRFitFunction},
	{"MxFitFunctionMultigroup", &initFitMultigroup}
};

void omxFreeFitFunctionArgs(omxFitFunction *off) {
	if(off==NULL) return;
    
	/* Completely destroy the fit function structures */
	if(OMX_DEBUG) {mxLog("Freeing fit function object at 0x%x.", off);}
	if(off->matrix != NULL) {
		if(off->destructFun != NULL) {
			if(OMX_DEBUG) {mxLog("Calling fit function destructor for 0x%x.", off);}
			off->destructFun(off);
		}
		off->matrix = NULL;
	}
}

void omxFitFunctionCreateChildren(omxState *globalState)
{
	if (Global->numThreads <= 1) return;

	int numThreads = Global->numChildren = Global->numThreads;

	globalState->childList = (omxState**) Calloc(numThreads, omxState*);

	for(int ii = 0; ii < numThreads; ii++) {
		globalState->childList[ii] = new omxState;
		omxInitState(globalState->childList[ii]);
		omxDuplicateState(globalState->childList[ii], globalState);
	}
}

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;

	omxFitFunction *ff = src->fitFunction;
	if(ff == NULL) return;
    
	omxFillMatrixFromMxFitFunction(tgt, ff->fitType, src->matrixNumber);
	omxCompleteFitFunction(tgt, ff->rObj);
}

void omxFitFunctionCompute(omxFitFunction *off, int want, double* gradient) {
	if (!off->initialized) error("FitFunction not initialized");

	if(OMX_DEBUG_ALGEBRA) { 
	    mxLog("FitFunction compute: 0x%0x (needed: %s).", off, (off->matrix->isDirty?"Yes":"No"));
	}

	off->computeFun(off, want, gradient);

	omxMarkClean(off->matrix);
}

void omxFillMatrixFromMxFitFunction(omxMatrix* om, const char *fitType, int matrixNumber)
{
	omxFitFunction *obj = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	memset(obj, 0, sizeof(omxFitFunction));

	/* Register FitFunction and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1, FALSE);					// FitFunction matrices MUST be 1x1.
	om->fitFunction = obj;
	om->hasMatrixNumber = TRUE;
	om->matrixNumber = matrixNumber;
	
	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if(strcmp(fitType, entry->name) == 0) {
			obj->fitType = entry->name;
			obj->initFun = entry->initFun;
			break;
		}
	}

	if (obj->initFun == NULL) error("Fit function %s not implemented", fitType);
}

void omxCompleteFitFunction(omxMatrix *om, SEXP rObj)
{
	omxFitFunction *obj = om->fitFunction;
	if (obj->initialized) return;
	obj->rObj = rObj;

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("expectation")));
	if (LENGTH(slotValue) == 1) {
		int expNumber = INTEGER(slotValue)[0];	
		if(expNumber != NA_INTEGER) {
			obj->expectation = omxExpectationFromIndex(expNumber, om->currentState);
		}
	}
	UNPROTECT(1);	/* slotValue */
	
	obj->initFun(obj);

	if(obj->computeFun == NULL) error("Failed to initialize fit function %s", obj->fitType); 
	
	obj->matrix->isDirty = TRUE;
	obj->initialized = TRUE;
}

void omxFitFunctionPrint(omxFitFunction* off, const char* d) {
	mxLog("(FitFunction, type %s)", off->fitType);
	omxPrintMatrix(off->matrix, d);
}


/* Helper functions */
omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* currentState, const char* slotName) {
	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install(slotName)));
	omxMatrix* newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
	UNPROTECT(1);
	return newMatrix;
}

