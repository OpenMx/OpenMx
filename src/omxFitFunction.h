/*
 *  Copyright 2007-2017 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/***********************************************************
* 
*  omxFitFunction.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxFitFunction class
*   omxFitFunction objects hold necessary information to simplify
* 	fit function calculation.
*
**********************************************************/

#ifndef _OMXFITFUNCTION_H_
#define _OMXFITFUNCTION_H_

#include "omxDefines.h"
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxData.h"
#include "omxState.h"
#include "omxExpectation.h"

struct omxFitFunction {
	SEXP rObj;
	omxExpectation* expectation;

	// This is always a pointer to a static string.
	// We do not need to allocate or free it.
	const char* fitType;

	omxMatrix* matrix;
	bool initialized;
	bool gradientAvailable;
	bool hessianAvailable;
	FitStatisticUnits units;
	bool canDuplicate;
	bool openmpUser; // can decide this in omxAlgebraPreeval

	omxFitFunction() : rObj(0), expectation(0), initialized(false), gradientAvailable(false),
		hessianAvailable(false), units(FIT_UNITS_UNINITIALIZED), canDuplicate(false), openmpUser(false) {};
	virtual ~omxFitFunction() {};
	virtual omxFitFunction *initMorph();
	virtual void init()=0;
	virtual void compute(int ffcompute, FitContext *fc)=0;
	virtual void invalidateCache() {};

	// addOutput should only be used for returning global results
	virtual void addOutput(MxRList *out) {};

	// populateAttr should be used for returning results specific to fit functions or expectations
	virtual void populateAttr(SEXP algebra) {};

	void setUnitsFromName(const char *name);
	const char *name() const { return matrix->name(); }
};

/* Initialize and Destroy */
void omxFillMatrixFromMxFitFunction(omxMatrix* om, int matrixNumber, SEXP rObj);

omxFitFunction *omxChangeFitType(omxFitFunction *oo, const char *fitType);

void omxCompleteFitFunction(omxMatrix *om);

	void omxGetFitFunctionStandardErrors(omxFitFunction *oo);					// Get Standard Errors

/* FitFunction-specific implementations of matrix functions */
void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc);
void omxFitFunctionComputeAuto(omxFitFunction *off, int want, FitContext *fc);
void omxFitFunctionComputeCI(omxFitFunction *off, int want, FitContext *fc);
	void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* targetState);

void omxFitFunctionPrint(omxFitFunction *source, const char* d);
	
omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* state, const char* slotName);

omxFitFunction *omxInitFIMLFitFunction();
omxFitFunction *omxInitAlgebraFitFunction();
omxFitFunction *omxInitWLSFitFunction();
omxFitFunction *omxInitRowFitFunction();
omxFitFunction *omxInitMLFitFunction();
omxFitFunction *omxInitRFitFunction();
omxFitFunction *omxInitFitFunctionBA81();
omxFitFunction *InitMarkovFF();
omxFitFunction *omxInitGREMLFitFunction();
omxFitFunction *InitFellnerFitFunction();
omxFitFunction *ssMLFitInit();

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg);

void ComputeFit(const char *callerName, omxMatrix *fitMat, int want, FitContext *fc);
void loglikelihoodCIFun(omxFitFunction* oo, int ffcompute, FitContext *fc);

double totalLogLikelihood(omxMatrix *fitMat);

const char *fitUnitsToName(FitStatisticUnits units);

#endif /* _OMXFITFUNCTION_H_ */
