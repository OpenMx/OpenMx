/*
 *  Copyright 2007-2014 The OpenMx Project
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
#include "fitMultigroup.h"

typedef struct omxFitFunctionTableEntry omxFitFunctionTableEntry;

struct omxFitFunctionTableEntry {

	char name[32];
	void (*initFun)(omxFitFunction*);
	void (*setVarGroup)(omxFitFunction*, FreeVarGroup *);  // TODO ugh, just convert to C++

};

static void defaultSetFreeVarGroup(omxFitFunction *ff, FreeVarGroup *fvg)
{
	if (ff->freeVarGroup && ff->freeVarGroup != fvg) {
		Rf_warning("%s: setFreeVarGroup called with different group (%d vs %d)",
			ff->matrix->name, ff->freeVarGroup->id[0], fvg->id[0]);
	}
	ff->freeVarGroup = fvg;
}

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction, defaultSetFreeVarGroup},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction, defaultSetFreeVarGroup},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction, defaultSetFreeVarGroup},
	{"MxFitFunctionML", 				&omxInitMLFitFunction, defaultSetFreeVarGroup},
	{"imxFitFunctionFIML", &omxInitFIMLFitFunction, defaultSetFreeVarGroup},
	{"MxFitFunctionR",					&omxInitRFitFunction, defaultSetFreeVarGroup},
	{"MxFitFunctionMultigroup", &initFitMultigroup, mgSetFreeVarGroup},
};

void omxFreeFitFunctionArgs(omxFitFunction *off) {
	if(off==NULL) return;
    
	/* Completely destroy the fit function structures */
	if(off->matrix != NULL) {
		if (off->destructFun) off->destructFun(off);
		off->matrix = NULL;
	}
}

void omxFitFunctionCreateChildren(omxState *globalState)
{
	if (Global->numThreads <= 1) return;

	for(size_t j = 0; j < globalState->expectationList.size(); j++) {
		if (!globalState->expectationList[j]->canDuplicate) return;
	}

	if (globalState->childList.size()) Rf_error("Children already created");

	int numThreads = Global->numThreads;

	globalState->childList.resize(numThreads);

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
	setFreeVarGroup(tgt->fitFunction, src->fitFunction->freeVarGroup);
	tgt->fitFunction->rObj = ff->rObj;
	omxCompleteFitFunction(tgt);
}

void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc)
{
	if (!off->initialized) Rf_error("FitFunction not initialized");

	off->computeFun(off, want, fc);
	if (fc) fc->wanted |= want;

	if (want & FF_COMPUTE_FIT) {
		omxMarkClean(off->matrix);
	}
}

void ComputeFit(omxMatrix *fitMat, int want, FitContext *fc)
{
	bool doFit = want & FF_COMPUTE_FIT;
	// R_CheckUserInterrupt(); add here? TODO

#pragma omp atomic
	++Global->computeCount; // could avoid lock by keeping in FitContext

	if (doFit) Global->checkpointPrefit(fc, fc->est, false);
	omxFitFunction *ff = fitMat->fitFunction;
	if (ff) {
		omxFitFunctionCompute(ff, want, fc);
	} else {
		if (want != FF_COMPUTE_FIT) Rf_error("Only fit is available");
		omxForceCompute(fitMat);
	}
	if (doFit) {
		if (fitMat->rows != 1) {
			if (strEQ(ff->fitType, "MxFitFunctionML") || strEQ(ff->fitType, "imxFitFunctionFIML")) {
				// NOTE: Floating-point addition is not
				// associative. If we compute this in parallel
				// then we introduce non-determinancy.
				double sum = 0;
				for(int i = 0; i < fitMat->rows; i++) {
					sum += log(omxVectorElement(fitMat, i));
				}
				fc->fit = sum * Global->llScale;
				if (!Global->rowLikelihoodsWarning) {
					Rf_warning("%s does not evaluate to a 1x1 matrix. Fixing model by adding "
						   "mxAlgebra(-2*sum(log(%s)), 'm2ll'), mxFitFunctionAlgebra('m2ll')",
						   fitMat->name, fitMat->name);
					Global->rowLikelihoodsWarning = true;
				}
			} else {
				omxRaiseErrorf("%s of type %s returned %d values instead of 1, not sure how to proceed",
					       fitMat->name, ff->fitType, fitMat->rows);
				fc->fit = nan("unknown");
			}
		} else {
			fc->fit = fitMat->data[0];
		}
		if (std::isfinite(fc->fit)) {
			fc->resetIterationError();
		}
		Global->checkpointPostfit(fc);
	}
}

void defaultAddOutput(omxFitFunction* oo, MxRList *out)
{}

void omxFillMatrixFromMxFitFunction(omxMatrix* om, const char *fitType, int matrixNumber)
{
	omxFitFunction *obj = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	memset(obj, 0, sizeof(omxFitFunction));

	/* Register FitFunction and Matrix with each other */
	obj->matrix = om;
	omxResizeMatrix(om, 1, 1);					// FitFunction matrices MUST be 1x1.
	om->fitFunction = obj;
	om->hasMatrixNumber = TRUE;
	om->matrixNumber = matrixNumber;
	
	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if(strcmp(fitType, entry->name) == 0) {
			obj->fitType = entry->name;
			obj->initFun = entry->initFun;

			// We need to set up the FreeVarGroup before calling initFun
			// because older fit functions expect to know the number of
			// free variables during initFun.
			obj->setVarGroup = entry->setVarGroup; // ugh!
			obj->addOutput = defaultAddOutput;
			break;
		}
	}

	if (obj->initFun == NULL) Rf_error("Fit function %s not implemented", fitType);
}

void omxCompleteFitFunction(omxMatrix *om)
{
	omxFitFunction *obj = om->fitFunction;
	if (obj->initialized) return;
	SEXP rObj = obj->rObj;

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("expectation")));
	if (LENGTH(slotValue) == 1) {
		int expNumber = INTEGER(slotValue)[0];	
		if(expNumber != NA_INTEGER) {
			obj->expectation = omxExpectationFromIndex(expNumber, om->currentState);
			setFreeVarGroup(obj->expectation, obj->freeVarGroup);
			omxCompleteExpectation(obj->expectation);
		}
	}
	Rf_unprotect(1);	/* slotValue */
	
	obj->initFun(obj);

	if(obj->computeFun == NULL) Rf_error("Failed to initialize fit function %s", obj->fitType); 
	
	obj->matrix->data[0] = NA_REAL;
	omxMarkDirty(obj->matrix);
	obj->initialized = TRUE;
}

void setFreeVarGroup(omxFitFunction *ff, FreeVarGroup *fvg)
{
	(*ff->setVarGroup)(ff, fvg);
}

void omxFitFunctionPrint(omxFitFunction* off, const char* d) {
	mxLog("(FitFunction, type %s)", off->fitType);
	omxPrintMatrix(off->matrix, d);
}


/* Helper functions */
omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* currentState, const char* slotName) {
	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install(slotName)));
	omxMatrix* newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
	Rf_unprotect(1);
	return newMatrix;
}

