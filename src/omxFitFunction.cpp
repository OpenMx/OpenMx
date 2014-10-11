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
	if (OMX_DEBUG && ff->freeVarGroup && ff->freeVarGroup != fvg) {
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

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;

	omxFitFunction *ff = src->fitFunction;
	if(ff == NULL) return;
    
	omxFillMatrixFromMxFitFunction(tgt, ff->fitType, src->matrixNumber, ff->rObj);
	setFreeVarGroup(tgt->fitFunction, src->fitFunction->freeVarGroup);
}

void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_DIMS | FF_COMPUTE_INITIAL_FIT)) return;

	if (!off->initialized) Rf_error("FitFunction not initialized");

	off->computeFun(off, want, fc);
	if (fc) fc->wanted |= want;
}

void ComputeFit(omxMatrix *fitMat, int want, FitContext *fc)
{
	bool doFit = want & FF_COMPUTE_FIT;
	R_CheckUserInterrupt();

#pragma omp atomic
	++Global->computeCount; // could avoid lock by keeping in FitContext

	// old version of openmp can't do this as part of the atomic instruction
	int evaluation = Global->computeCount;

	if (doFit) {
		if (OMX_DEBUG) {
			mxLog("%s: starting evaluation %d, want %d", fitMat->name, evaluation, want);
		}
		Global->checkpointPrefit(fc, fc->est, false);
	}
	omxFitFunction *ff = fitMat->fitFunction;
	if (ff) {
		omxFitFunctionCompute(ff, want, fc);
	} else {
		if (want != FF_COMPUTE_FIT) Rf_error("Only fit is available");
		omxRecompute(fitMat, want, fc);
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
		if (OMX_DEBUG) {
			mxLog("%s: completed evaluation %d, fit=%f", fitMat->name, evaluation, fc->fit);
		}
	}
}

void defaultAddOutput(omxFitFunction* oo, MxRList *out)
{}

omxFitFunction *omxNewInternalFitFunction(omxState* os, const char *fitType,
					  omxExpectation *expect, omxMatrix *matrix, bool rowLik)
{
	omxFitFunction *obj = (omxFitFunction*) R_alloc(1, sizeof(omxFitFunction));
	OMXZERO(obj, 1);

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

	if(obj->initFun == NULL) Rf_error("Fit function '%s' not implemented", fitType);

	if (!matrix) {
		obj->matrix = omxInitMatrix(1, 1, TRUE, os);
		obj->matrix->hasMatrixNumber = TRUE;
		obj->matrix->matrixNumber = ~os->algebraList.size();
		os->algebraList.push_back(obj->matrix);
	} else {
		obj->matrix = matrix;
	}

	obj->matrix->fitFunction = obj;

	obj->expectation = expect;

	if (rowLik && expect && expect->data) {
		omxData *dat = expect->data;
		omxResizeMatrix(matrix, dat->rows, 1);
	} else {
		omxResizeMatrix(matrix, 1, 1);
	}

	return obj;
}

void omxFillMatrixFromMxFitFunction(omxMatrix* om, const char *fitType, int matrixNumber, SEXP rObj)
{
	om->hasMatrixNumber = TRUE;
	om->matrixNumber = matrixNumber;

	SEXP slotValue;
	omxExpectation *expect = NULL;
	{
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("expectation")));
		if (Rf_length(slotValue) == 1) {
			int expNumber = Rf_asInteger(slotValue);
			if(expNumber != NA_INTEGER) {
				expect = omxExpectationFromIndex(expNumber, om->currentState);
			}
		}
	}

	bool rowLik = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));

	omxFitFunction *ff =
		omxNewInternalFitFunction(om->currentState, fitType, expect, om, rowLik);
	ff->rObj = rObj;
}

void omxChangeFitType(omxFitFunction *oo, const char *fitType)
{
	if (oo->initialized) {
		Rf_error("%s: cannot omxChangeFitType from %s to %s; already initialized",
			 oo->matrix->name, oo->fitType, fitType);
	}

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if (strEQ(fitType, entry->name)) {
			oo->fitType = entry->name;
			oo->initFun = entry->initFun;
			return;
		}
	}

	Rf_error("Cannot find fit type '%s'", fitType);
}

void omxCompleteFitFunction(omxMatrix *om)
{
	omxFitFunction *obj = om->fitFunction;
	if (obj->initialized) return;

	if (obj->expectation) {
		setFreeVarGroup(obj->expectation, obj->freeVarGroup);
		omxCompleteExpectation(obj->expectation);
	}

	obj->initFun(obj);

	if(obj->computeFun == NULL) Rf_error("Failed to initialize fit function %s", obj->fitType); 
	
	obj->matrix->data[0] = NA_REAL;
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
	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install(slotName)));
	omxMatrix* newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
	return newMatrix;
}

