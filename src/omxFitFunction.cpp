/*
 *  Copyright 2007-2018 The OpenMx Project
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
#include "Compute.h"
#include "EnableWarnings.h"

typedef struct omxFitFunctionTableEntry omxFitFunctionTableEntry;

struct omxFitFunctionTableEntry {

	char name[32];
	omxFitFunction *(*allocate)();

};

omxFitFunction *omxFitFunction::initMorph()
{
	init();
	return this;
}

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction},
	{"MxFitFunctionML", 				&omxInitMLFitFunction},
	{"imxFitFunctionFIML", &omxInitFIMLFitFunction},
	{"MxFitFunctionR",					&omxInitRFitFunction},
	{"MxFitFunctionMultigroup", &initFitMultigroup},
	{"MxFitFunctionGREML", &omxInitGREMLFitFunction},
	{"imxFitFunctionFellner", &InitFellnerFitFunction},
	{"imxFitFunctionBA81", &omxInitFitFunctionBA81},
	{"imxFitFunciontStateSpace", &ssMLFitInit},
	{"imxFitFunciontHiddenMarkov", &InitMarkovFF},
};

void omxFitFunction::setUnitsFromName(const char *name)
{
	if (strEQ(name, "-2lnL")) {
		units = FIT_UNITS_MINUS2LL;
	} else {
		Rf_warning("Unknown units '%s' passed to fit function '%s'",
			   name, matrix->name());
		units = FIT_UNITS_UNKNOWN;
	}
}

bool fitUnitsIsChiSq(FitStatisticUnits units)
{
	return units == FIT_UNITS_MINUS2LL || units == FIT_UNITS_SQUARED_RESIDUAL_CHISQ;
}

const char *fitUnitsToName(FitStatisticUnits units)
{
	switch (units) {
	case FIT_UNITS_UNINITIALIZED: return "";
	case FIT_UNITS_UNKNOWN: return "?";
	case FIT_UNITS_PROBABILITY: return "Pr";
	case FIT_UNITS_MINUS2LL: return "-2lnL";
	case FIT_UNITS_SQUARED_RESIDUAL: return "r'Wr";
	case FIT_UNITS_SQUARED_RESIDUAL_CHISQ: return "r'Wr";
	default: Rf_error("Don't know how to stringify units %d", units);
	}
}

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;

	omxFitFunction *ff = src->fitFunction;
	if(ff == NULL) return;

	omxFillMatrixFromMxFitFunction(tgt, src->matrixNumber, ff->rObj);
}

static void ciFunction(omxFitFunction *ff, int want, FitContext *fc)
{
	if (fitUnitsIsChiSq(ff->units)) {
		fc->ciobj->evalFit(ff, want, fc);
	} else {
		Rf_error("Confidence intervals are not supported for units %s",
			 fitUnitsToName(ff->units));
	}
}

void omxFitFunctionComputeAuto(omxFitFunction *off, int want, FitContext *fc)
{
	if (!off->initialized) return;

	if (!fc->ciobj) {
		off->compute(want, fc);
	} else {
		ciFunction(off, want, fc);
	}

	if (fc) fc->wanted |= want;
}

void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc)
{
	if (!off->initialized) return;

	off->compute(want, fc);
	if (fc) fc->wanted |= want;
}

void omxFitFunctionComputeCI(omxFitFunction *off, int want, FitContext *fc)
{
	if (!off->initialized) return;

	ciFunction(off, want, fc);
	if (fc) fc->wanted |= want;
}

double totalLogLikelihood(omxMatrix *fitMat)
{
	if (fitMat->rows != 1) {
		omxFitFunction *ff = fitMat->fitFunction;
		if (ff->units == FIT_UNITS_PROBABILITY) {
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
					   fitMat->name(), fitMat->name());
				Global->rowLikelihoodsWarning = true;
			}
			return sum * Global->llScale;
		} else {
			omxRaiseErrorf("%s of type %s returned %d values instead of 1, not sure how to proceed",
				       fitMat->name(), ff->fitType, fitMat->rows);
			return nan("unknown");
		}
	} else {
		return fitMat->data[0];
	}
}

void ComputeFit(const char *callerName, omxMatrix *fitMat, int want, FitContext *fc)
{
	fc->incrComputeCount();
	fc->skippedRows = 0;
	omxFitFunction *ff = fitMat->fitFunction;
	if (ff) {
		omxFitFunctionComputeAuto(ff, want, fc);
	} else {
		if (want != FF_COMPUTE_FIT) Rf_error("Only fit is available");
		if (fc->ciobj) Rf_error("CIs cannot be computed for unitless algebra");
		omxRecompute(fitMat, fc);
	}
	if (ff && want & FF_COMPUTE_FIT) {
		fc->fit = totalLogLikelihood(fitMat);
		if (std::isfinite(fc->fit)) {
			fc->resetIterationError();
		}
		Global->checkpointPostfit(callerName, fc, fc->est, false);
		if (OMX_DEBUG) {
			mxLog("%s: completed evaluation, fit=%.12g skippedRows=%d",
			      fitMat->name(), fc->fit, fc->skippedRows);
		}
	}
}

static omxFitFunction *omxNewInternalFitFunction(omxState* os, const char *fitType,
						 omxExpectation *expect, omxMatrix *matrix, bool rowLik)
{
	omxFitFunction *obj = 0;

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if(strcmp(fitType, entry->name) == 0) {
			obj = entry->allocate();
			obj->fitType = entry->name;
			break;
		}
	}
	if (!obj) Rf_error("omxNewInternalFitFunction: cannot find '%s'", fitType);

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

void omxFillMatrixFromMxFitFunction(omxMatrix* om, int matrixNumber, SEXP rObj)
{
	om->hasMatrixNumber = TRUE;
	om->matrixNumber = matrixNumber;

	ProtectedSEXP fitFunctionClass(STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
	const char *fitType = CHAR(fitFunctionClass);

	omxExpectation *expect = NULL;
	ProtectedSEXP slotValue(R_do_slot(rObj, Rf_install("expectation")));
	if (Rf_length(slotValue) == 1) {
		int expNumber = Rf_asInteger(slotValue);
		if(expNumber != NA_INTEGER) {
			expect = omxExpectationFromIndex(expNumber, om->currentState);
		}
	}

	bool rowLik = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));

	omxFitFunction *ff =
		omxNewInternalFitFunction(om->currentState, fitType, expect, om, rowLik);
	ff->rObj = rObj;
}

omxFitFunction *omxChangeFitType(omxFitFunction *oo, const char *fitType)
{
	if (oo->initialized) {
		Rf_error("%s: cannot omxChangeFitType from %s to %s; already initialized",
			 oo->matrix->name(), oo->fitType, fitType);
	}

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if (strEQ(fitType, entry->name)) {
			auto *newObj = entry->allocate();
			newObj->rObj = oo->rObj;
			newObj->expectation = oo->expectation;
			newObj->fitType = entry->name;
			newObj->matrix = oo->matrix;
			newObj->units = oo->units;
			oo->matrix = 0;
			newObj->matrix->fitFunction = newObj;
			delete oo;
			// Need to call initMorph again? Probably never have 2 levels
			// of specialization?
			newObj->init();
			return newObj;
		}
	}

	Rf_error("Cannot find fit type '%s'", fitType);
}

void omxCompleteFitFunction(omxMatrix *om)
{
	omxFitFunction *obj = om->fitFunction;
	if (obj->initialized) return;

	int depth = Global->mpi->getDepth();

	if (obj->expectation) {
		omxCompleteExpectation(obj->expectation);
	}

	obj = obj->initMorph();

	if (Global->mpi->getDepth() != depth) Rf_error("%s improperly used the R protect stack", om->name());

	obj->matrix->data[0] = NA_REAL;
	obj->initialized = TRUE;
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

