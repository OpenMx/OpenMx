/*
 *  Copyright 2007-2009 The OpenMx Project
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
*  omxState.h
*
*  Created: Timothy R. Brick 	Date: 2009-05-23
*
*	Contains header information for the omxState structure
*   omxStates keep the current optimization state.
*
**********************************************************/

#ifndef _OMXSTATE_H_
#define _OMXSTATE_H_

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"

/* Forward declarations for later includes */
typedef struct omxState omxState;
typedef struct omxFreeVar omxFreeVar;
typedef struct omxConstraint omxConstraint;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxObjective.h"
#include "omxData.h"
//#include "omxOptimizer.h"											// omxOptimizer objects coming soon

/* Structure definitions for object evaluation */  // Might be cleaner to give these their own files.
struct omxFreeVar {			// Free Variables
	double lbound, ubound;	// Bounds
	int numLocations;
	double** location;		// And where they go.
	int* matrices;			// Matrix numbers for dirtying.
};

struct omxConstraint {		// Free Variable Constraints
	int size;
	int opCode;
	double* lbound;
	double* ubound;
	omxMatrix* result;
};

#define MAX_STRING_LEN 250

struct omxState {													// The Current State of Optimization

/* Model and Optimizer Pointers */

//	omxOptimizer* optimizer;										// Current Optimizer
	int numMats, numAlgs, numData;									// Number of matrices, algebras, and data elements
	omxMatrix** matrixList;											// Model Matrices
	omxMatrix** algebraList;										// Model Algebras
	omxData** dataList;												// Data Objects
	omxMatrix* objectiveMatrix;										// Objective Algebra

	/* May want to farm these out to the omxObjective object. */
	int numConstraints;
	omxConstraint* conList;											// List of constraints

	int numFreeParams;
	omxFreeVar* freeVarList;										// List of Free Variables and where they go.

	unsigned short int allowShortcut;								// See below:
	/* This is primarily included to allow special case objective functions and optimizers to disable it.
		All matrices are initialized at compute count 0.  The optimizer is requested to run the objective
		function once at compute count 1 to trace dependencies.  All matrices that are updated at compute
		count 1 (that is, are non-static) will advance their lastCompute to 1, and those that are computed
		at each row will advance their lastRow to the total number of rows.
		If Shortcut evaluation is enabled, any matrix/algebra with a lastCompute of 0 will not check to see 
		if it needs updating, ever.  Any matrix/algebra with a lastRow of 0 will not check row-by-row to see
		if it need re-evaluation.  This can be disabled in cases where the dependency tree might change during
		optimization.  NOT YET IMPLEMENTED. */

/* Current Optimization State (optimizer-specific) */
	void* optimizerState;											// Optimizer specific state storage

/* Data members for use by Objective Function and Algebra Calculations */
	long int computeCount;											// How many times have things been evaluated so far?
	long int currentRow;											// If we're calculating row-by-row, what row are we on?
	
	int statusCode;													// Status code, if appropriate
	char statusMsg[250];											// Status/Error message to report
	double saturatedModel;											// Saturated model likelihood, where applicable

};

/* Initialize and Destroy */
	void omxInitState(omxState* state);									// Null Constructor
	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList, omxMatrix** algebraList, omxData** dataList, omxMatrix* objective); 
	void omxFreeState(omxState *oo);									// Destructor
	
	void omxRaiseError(omxState *oo, int errorCode, char* errorMsg);	// Raise an Error 
																		// TODO: Move RaiseError to omxOptimizer.
	
/* Advance a step */
	void omxRaiseError(omxState *oo, int errorCode, char* errorStr);
	void omxStateNextRow(omxState *oo);									// Advance Row
	void omxStateNextEvaluation(omxState *oo);							// Advance Evaluation count
	
#endif /* _OMXSTATE_H_ */


