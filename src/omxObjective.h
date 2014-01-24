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
*  omxObjective.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxObjective class
*   omxObjective objects hold necessary information to simplify
* 	objective function calculation.
*
**********************************************************/

#ifndef _OMXOBJECTIVE_H_
#define _OMXOBJECTIVE_H_

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"

typedef struct omxObjective omxObjective;
typedef struct omxRListElement omxRListElement;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxAlgebraFunctions.h"
#include "omxData.h"
#include "omxState.h"
#include "omxObjectiveTable.h"

struct omxRListElement {
	char label[250];
	double* values;
	int numValues;
	int rows, cols;
};

struct omxObjective {					// An objective

	/* Fields unique to Objective Functions */
	void (*initFun)(omxObjective *oo, SEXP rObj);								// Wrapper for initialization function (probably not needed)
	void (*destructFun)(omxObjective* oo);										// Wrapper for the destructor object
	void (*objectiveFun)(omxObjective* oo);										// Wrapper for the objective function itself
	void (*repopulateFun)(omxObjective* oo, double* x, int n);					// To repopulate any data stored in the objective function
	omxRListElement* (*setFinalReturns)(omxObjective* oo, int *numVals);		// Sets any R returns.
	unsigned short int (*needsUpdateFun)(omxObjective* oo);						// To calculate recomputation
	double* (*getStandardErrorFun)(omxObjective* oo);							// To calculate standard errors
	void (*gradientFun)(omxObjective* oo, double* grad);						// To calculate gradient
	void (*populateAttrFun)(omxObjective* oo, SEXP algebra);					// Add attributes to the result algebra object
	
	omxObjective* subObjective;													// Subobjective structure in case it is needed.
	
	SEXP rObj;																	// Original r Object Pointer
	void* sharedArgs;                                                           // Common argument structure
	void* argStruct;															// Arguments needed for objective function
	char* objType;																// Type of Objective Function
	double* hessian;															// Hessian details
	double* gradient;															// Gradient details
	double* stdError;															// Standard Error estimates
	unsigned short int isPrepopulated;											// Object has had some values prepopulated to allow object sharing

	omxMatrix* matrix;															// The (1x1) matrix populated by this objective function

};

/* Initialize and Destroy */
	void omxInitEmptyObjective(omxObjective *oo);
	omxObjective* omxCreateSubObjective(omxObjective *oo);					// Generate a subobjective object
	void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP mxobj, 
		unsigned short hasMatrixNumber, int matrixNumber);					// Create an objective function from an R MxObjective object
	void omxFreeObjectiveArgs(omxObjective* objective);						// Frees all args
	void omxGetObjectiveStandardErrors(omxObjective *oo);					// Get Standard Errors

/* Objective-specific implementations of matrix functions */
	void omxObjectiveCompute(omxObjective *oo);
	unsigned short int omxObjectiveNeedsUpdate(omxObjective *oo);
	void omxObjectiveGradient(omxObjective* oo, double* gradient);			// For gradient calculation.  If needed.
	void omxDuplicateObjectiveMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* targetState);
	omxObjective* omxCreateDuplicateObjective(omxObjective *tgt, const omxObjective *src, omxState* newState);
	
	void omxObjectivePrint(omxObjective *source, char* d);					// Pretty-print a (small) matrix
	
	/* Helper functions */
	void omxCalculateStdErrorFromHessian(double scale, omxObjective *oo);	// Does what it says
	
	/* Helpers related to objective initialization */
	omxMatrix* omxNewMatrixFromIndexSlot(SEXP rObj, omxState* state, char* const slotName);	// Gets a matrix from an R SEXP slot
	omxData* omxNewDataFromDataSlot(SEXP rObj, omxState* state, char* const dataSlotName);	// Gets an mxData object from a data slot

#endif /* _OMXOBJECTIVE_H_ */
