/*
 *  Copyright 2007-2015 The OpenMx Project
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
	{"MxFitFunctionGREML", &omxInitGREMLFitFunction, defaultSetFreeVarGroup},
};

void omxFitFunction::setUnitsFromName(const char *name)
{
	if (strEQ(name, "-2lnL")) {
		units = FIT_UNITS_MINUS2LL;
		ciFun = loglikelihoodCIFun;
	} else {
		Rf_warning("Unknown units '%s' passed to fit function '%s'",
			   name, matrix->name);
		units = FIT_UNITS_UNKNOWN;
	}
}

const char *fitUnitsToName(int units)
{
	switch (units) {
	case FIT_UNITS_UNINITIALIZED: return "";
	case FIT_UNITS_UNKNOWN: return "?";
	case FIT_UNITS_MINUS2LL: return "-2lnL";
	case FIT_UNITS_SQUARED_RESIDUAL: return "r'Wr";
	default: Rf_error("Don't know how to stringify units %d", units);
	}
}

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

void omxFitFunctionComputeAuto(omxFitFunction *off, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_DIMS | FF_COMPUTE_INITIAL_FIT)) return;

	if (!off->initialized) Rf_error("FitFunction not initialized");

	if (!fc->CI) {
		off->computeFun(off, want, fc);
	} else {
		off->ciFun(off, want, fc);
	}

	if (fc) fc->wanted |= want;
}

void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_DIMS | FF_COMPUTE_INITIAL_FIT)) return;

	if (!off->initialized) Rf_error("FitFunction not initialized");

	off->computeFun(off, want, fc);
	if (fc) fc->wanted |= want;
}

void omxFitFunctionComputeCI(omxFitFunction *off, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_DIMS | FF_COMPUTE_INITIAL_FIT)) return;

	if (!off->initialized) Rf_error("FitFunction not initialized");

	off->ciFun(off, want, fc);
	if (fc) fc->wanted |= want;
}

double totalLogLikelihood(omxMatrix *fitMat)
{
	if (fitMat->rows != 1) {
		omxFitFunction *ff = fitMat->fitFunction;
		if (strEQ(ff->fitType, "MxFitFunctionML") || strEQ(ff->fitType, "imxFitFunctionFIML")) {
			// NOTE: Floating-point addition is not
			// associative. If we compute this in parallel
			// then we introduce non-determinancy.
			double sum = 0;
			for(int i = 0; i < fitMat->rows; i++) {
				sum += log(omxVectorElement(fitMat, i));
			}
			if (!Global->rowLikelihoodsWarning) {
				Rf_warning("%s does not evaluate to a 1x1 matrix. Fixing model by adding "
					   "mxAlgebra(-2*sum(log(%s)), 'm2ll'), mxFitFunctionAlgebra('m2ll')",
					   fitMat->name, fitMat->name);
				Global->rowLikelihoodsWarning = true;
			}
			return sum * Global->llScale;
		} else {
			omxRaiseErrorf("%s of type %s returned %d values instead of 1, not sure how to proceed",
				       fitMat->name, ff->fitType, fitMat->rows);
			return nan("unknown");
		}
	} else {
		return fitMat->data[0];
	}
}

void omxFitFunctionPreoptimize(omxFitFunction *off, FitContext *fc)
{
	omxFitFunctionComputeAuto(off, FF_COMPUTE_PREOPTIMIZE, fc);
	fc->fitUnits = off->units;
}

void ComputeFit(const char *callerName, omxMatrix *fitMat, int want, FitContext *fc)
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
		Global->checkpointPrefit(callerName, fc, fc->est, false);
	}
	omxFitFunction *ff = fitMat->fitFunction;
	if (ff) {
		omxFitFunctionComputeAuto(ff, want, fc);
	} else {
		if (want != FF_COMPUTE_FIT) Rf_error("Only fit is available");
		if (fc->CI) Rf_error("CIs cannot be computed for unitless algebra");
		omxRecompute(fitMat, fc);
	}
	if (doFit) {
		fc->fit = totalLogLikelihood(fitMat);
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

static void defaultCIFun(omxFitFunction* oo, int ffcompute, FitContext *fc)
{
	Rf_error("Confidence intervals are not supported for units %d", oo->units);
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
	if(obj->ciFun == NULL) obj->ciFun = defaultCIFun;

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

void loglikelihoodCIFun(omxFitFunction *ff, int want, FitContext *fc)
{
	const omxConfidenceInterval *CI = fc->CI;

	if (want & FF_COMPUTE_PREOPTIMIZE) {
		fc->targetFit = (fc->lowerBound? CI->lbound : CI->ubound) + fc->fit;
		//mxLog("Set target fit to %f (MLE %f)", fc->targetFit, fc->fit);
		return;
	}

	if (!(want & FF_COMPUTE_FIT)) {
		Rf_error("Not implemented yet");
	}

	omxMatrix *fitMat = ff->matrix;

	// We need to compute the fit here because that's the only way to
	// check our soft feasibility constraints. If parameters don't
	// change between here and the constraint evaluation then we
	// should avoid recomputing the fit again in the constraint. TODO

	omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
	const double fit = totalLogLikelihood(fitMat);
	omxRecompute(CI->matrix, fc);
	double CIElement = omxMatrixElement(CI->matrix, CI->row, CI->col);
	omxResizeMatrix(fitMat, 1, 1);

	if (!std::isfinite(fit) || !std::isfinite(CIElement)) {
		fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
		fitMat->data[0] = nan("infeasible");
		return;
	}

	if (want & FF_COMPUTE_FIT) {
		fitMat->data[0] = (fc->lowerBound? CIElement : -CIElement);
		//mxLog("param at %f", fitMat->data[0]);
	}
	if (want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
		// add deriv adjustments here TODO
	}
}

