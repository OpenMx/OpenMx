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

extern void omxInitAlgebraFitFunction(omxFitFunction *off);
extern void omxInitWLSFitFunction(omxFitFunction *off);
extern void omxInitRowFitFunction(omxFitFunction *off);
extern void omxInitMLFitFunction(omxFitFunction *off);
extern void omxInitRFitFunction(omxFitFunction *off);

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction},
	{"MxFitFunctionML", 				&omxInitMLFitFunction},
	{"MxFitFunctionR",					&omxInitRFitFunction},
	{"MxFitFunctionMultigroup", &initFitMultigroup}
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
		globalState->childList[ii] = new omxState;
		omxInitState(globalState->childList[ii], globalState);
		omxDuplicateState(globalState->childList[ii], globalState);
	}
}

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;
	if(src->fitFunction == NULL) return;
    
	omxFillMatrixFromMxFitFunction(src->fitFunction->rObj, tgt, src->currentState);

}

omxFitFunction* omxCreateDuplicateFitFunction(omxFitFunction *tgt, const omxFitFunction *src, omxState* newState) {

	if(OMX_DEBUG) {Rprintf("Duplicating fit function 0x%x into 0x%x.", src, tgt);}

	if(src == NULL) {
		return NULL;
	}
	
	if(tgt == NULL) {
        tgt = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	OMXZERO(tgt, 1);
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

omxFitFunction *omxNewInternalFitFunction(omxState* os, const char *fitType,
					  omxExpectation *expect, SEXP rObj, omxMatrix *matrix)
{
	omxFitFunction *obj = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	OMXZERO(obj, 1);

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if(strcmp(fitType, entry->name) == 0) {
			obj->fitType = entry->name;
			obj->initFun = entry->initFun;
			break;
		}
	}

	if(obj->initFun == NULL) error("Fit function %s not implemented", fitType);

	if (!matrix) {
		obj->matrix = omxInitMatrix(NULL, 1, 1, TRUE, os);
		obj->matrix->hasMatrixNumber = TRUE;
		obj->matrix->matrixNumber = ~os->algebraList.size();
		os->algebraList.push_back(obj->matrix);
	} else {
		obj->matrix = matrix;
	}

	obj->matrix->fitFunction = obj;
	
	obj->rObj = rObj;
	obj->expectation = expect;

	return obj;
}

void omxFillMatrixFromMxFitFunction(SEXP rObj, omxMatrix *matrix, omxState *os)
{
	SEXP slotValue, fitFunctionClass;

	PROTECT(fitFunctionClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	const char *fitType = CHAR(fitFunctionClass);

	omxExpectation *expect = NULL;
	PROTECT(slotValue = GET_SLOT(rObj, install("expectation")));
	if (LENGTH(slotValue) == 1) {
		int expNumber = INTEGER(slotValue)[0];	
		if(expNumber != NA_INTEGER) {
			expect = omxExpectationFromIndex(expNumber, os);
		}
	}

	omxNewInternalFitFunction(os, fitType, expect, rObj, matrix);

	UNPROTECT(2);
}

void omxInitializeFitFunction(omxMatrix *om)
{
	omxFitFunction *obj = om->fitFunction;
	if (!obj) error("Matrix 0x%p has no fit function", om);

	if (obj->initialized) return;
	obj->initialized = TRUE;

	obj->initFun(obj);

	if(obj->computeFun == NULL) error("Could not initialize fit function %s", obj->fitType);
	
	obj->matrix->isDirty = TRUE;
}

void omxFitFunctionPrint(omxFitFunction* off, const char* d) {
	Rprintf("(FitFunction, type %s) ", off->fitType);
	omxPrintMatrix(off->matrix, d);
}


/* Helper functions */
omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* currentState, const char* slotName) {
	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install(slotName)));
	omxMatrix* newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
	if (newMatrix) omxRecompute(newMatrix);
	UNPROTECT(1);
	return newMatrix;
}

