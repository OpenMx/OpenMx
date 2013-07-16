/*
 *  Copyright 2007-2013 The OpenMx Project
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
*  omxExpectation.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxExpectation class
*   omxExpectation objects hold necessary information to simplify
* 	Expectation function calculation.
*
**********************************************************/

#ifndef _OMXEXPECTATION_H_
#define _OMXEXPECTATION_H_

#include <R.h>
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxDefines.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxAlgebraFunctions.h"
#include "omxData.h"
#include "omxState.h"

/* Def Var and Threshold Structures */
struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
	omxData* source;		// Data source
	int numLocations;		// Num locations
	int* rows;				// row positions
	int* cols;				// column positions
	int* matrices;			// matrix numbers
	int  numDeps;           // number of algebra/matrix dependencies
	int* deps;              // indices of algebra/matrix dependencies

};

/* Expectation structure itself */
struct omxExpectation {					// An Expectation

	/* Fields unique to Expectation Functions */
	void (*initFun)(omxExpectation *ox);
	void (*destructFun)(omxExpectation* ox);									// Wrapper for the destructor object
	void (*computeFun)(omxExpectation* ox, const char *context);
	void (*printFun)(omxExpectation* ox);										// Prints the appropriate pieces of the expectation
	void (*populateAttrFun)(omxExpectation* ox, SEXP algebra);					// Add attributes to the result algebra object
	void (*setVarGroup)(omxExpectation*, FreeVarGroup *);  // TODO remove
	
	// componentfun & mutateFun probably take encapsulation a little too seriously.
	// The Fit function should probably just include the structure definition
	// for the expectation and access fields directly or through object methods.
	omxMatrix* (*componentFun)(omxExpectation*, omxFitFunction*, const char*);
	void (*mutateFun)(omxExpectation*, omxFitFunction*, const char*, omxMatrix*);

	SEXP rObj;																	// Original r Object Pointer
	void* argStruct;															// Arguments needed for Expectation function
        const char* expType;   // pointer to a static string, no need to allocate or free

	omxData* data;																// Not sure if this is appropriate, but the expectation passes the actual data object
	omxMatrix* dataColumns;
	omxThresholdColumn* thresholds;
	omxDefinitionVar* defVars;
	int numOrdinal, numDefs;
	
	/* Replication of some of the structures from Matrix */
	unsigned short isComplete;													// Whether or not this expectation has been initialize
	omxState* currentState;
	int expNum;

	omxExpectation *container;
	int numSubmodels;
	omxExpectation **submodels;

	// omxExpectation should not need to know about free variables.
	FreeVarGroup *freeVarGroup; // TODO remove
};

omxExpectation *
omxNewInternalExpectation(const char *expType, omxState* os);

	void omxCompleteExpectation(omxExpectation *ox);
void setFreeVarGroup(omxExpectation *ox, FreeVarGroup *fvg);

	void omxFreeExpectationArgs(omxExpectation* Expectation);					// Frees all args
omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os);
	omxExpectation* omxNewIncompleteExpectation(SEXP mxobj, int expNum, omxState* os);
	

/* Expectation-specific implementations of matrix functions */
void omxExpectationRecompute(omxExpectation *ox);
void omxExpectationCompute(omxExpectation *ox, const char *context);
	omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState);
	
	void omxExpectationPrint(omxExpectation *source, char* d);					// Pretty-print a (small-form) expectation
	
omxMatrix* omxGetExpectationComponent(omxExpectation *ox, omxFitFunction *off, const char* component);
	
void omxSetExpectationComponent(omxExpectation *ox, omxFitFunction *off, const char* component, omxMatrix *om);

void omxInitNormalExpectation(omxExpectation *ox);
void omxInitLISRELExpectation(omxExpectation *ox);
void omxInitStateSpaceExpectation(omxExpectation *ox);
void omxInitRAMExpectation(omxExpectation *ox);
void omxInitExpectationBA81(omxExpectation* oo);

#endif /* _OMXEXPECTATION_H_ */
