/*
 *  Copyright 2007-2012 The OpenMx Project
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

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"

typedef struct omxFitFunction omxFitFunction;


#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxAlgebraFunctions.h"
#include "omxData.h"
#include "omxState.h"
#include "omxExpectation.h"

struct omxFitFunction {					// A fit function

	/* Fields unique to FitFunction Functions */
	void (*initFun)(omxFitFunction *oo, SEXP rObj);								// Wrapper for initialization function (probably not needed)
	void (*destructFun)(omxFitFunction* oo);									// Wrapper for the destructor object
	// ffcompute is somewhat redundent because grad=NULL when gradients are unwanted
	void (*computeFun)(omxFitFunction* oo, int ffcompute, double* grad);
	void (*repopulateFun)(omxFitFunction* oo, double* x, int n);					// To repopulate any data stored in the fit function
	omxRListElement* (*setFinalReturns)(omxFitFunction* oo, int *numVals);		// Sets any R returns.
	double* (*getStandardErrorFun)(omxFitFunction* oo);							// To calculate standard errors
	void (*populateAttrFun)(omxFitFunction* oo, SEXP algebra);					// Add attributes to the result algebra object
	
	SEXP rObj;																	// Original r Object Pointer
	omxExpectation* expectation;												// Data expectation object
	void* sharedArgs;															// Common argument structure
	void* argStruct;															// Arguments needed for fit function
// This is always a pointer to a static string.
// We do not need to allocate or free it.
	const char* fitType;														// Type of FitFunction Function
	double* hessian;															// Hessian details
	double* gradient;															// Gradient details
	double* stdError;															// Standard Error estimates
	unsigned short int isPrepopulated;											// Object has had some values prepopulated to allow object sharing

	omxMatrix* matrix;															// The (1x1) matrix populated by this fit function

};

/* Initialize and Destroy */
	void omxInitEmptyFitFunction(omxFitFunction *oo);
	void omxFillMatrixFromMxFitFunction(omxMatrix* om, SEXP mxobj, 
		unsigned short hasMatrixNumber, int matrixNumber);						// Create an omxFitFunction from an R MxFitFunction object
	void omxFreeFitFunctionArgs(omxFitFunction* fitFunction);						// Frees all args
	void omxGetFitFunctionStandardErrors(omxFitFunction *oo);					// Get Standard Errors

/* FitFunction-specific implementations of matrix functions */
void omxFitFunctionCompute(omxFitFunction *off, int want, double* gradient);
	void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* targetState);
	omxFitFunction* omxCreateDuplicateFitFunction(omxFitFunction *tgt, const omxFitFunction *src, omxState* newState);
	
	void omxFitFunctionPrint(omxFitFunction *source, char* d);					// Pretty-print a (small) matrix
	
	/* Helper functions */
	void omxCalculateStdErrorFromHessian(double scale, omxFitFunction *oo);	// Does what it says
	
	/* Helpers related to fit function initialization */
	omxMatrix* omxNewMatrixFromIndexSlot(SEXP rObj, omxState* state, char* const slotName);	// Gets a matrix from an R SEXP slot

	omxData* omxNewDataFromDataSlot(SEXP rObj, omxState* state, char* const dataSlotName);	// Gets an mxData object from a data slot

enum omxFFCompute {
  FF_COMPUTE_FIT      = 1<<0,
  FF_COMPUTE_GRADIENT = 1<<1,
  FF_COMPUTE_HESSIAN  = 1<<2
};

#endif /* _OMXFITFUNCTION_H_ */
