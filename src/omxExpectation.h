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

typedef struct omxExpectation omxExpectation;
typedef struct omxRListElement omxRListElement;
typedef struct omxDefinitionVar omxDefinitionVar;
typedef struct omxThresholdColumn omxThresholdColumn;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxAlgebraFunctions.h"
#include "omxData.h"
#include "omxState.h"

struct omxRListElement {
	char label[250];
	double* values;
	int numValues;
	int rows, cols;
};

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

struct omxThresholdColumn {		 	// Threshold

	omxMatrix* matrix;		// Which Matrix/Algebra it comes from
	int column;				// Which column has the thresholds
	int numThresholds;		// And how many thresholds

};

/* Expectation structure itself */
struct omxExpectation {					// An Expectation

	/* Fields unique to Expectation Functions */
	void (*initFun)(omxExpectation *ox, SEXP rObj);								// Wrapper for initialization function (probably not needed)
	void (*destructFun)(omxExpectation* ox);									// Wrapper for the destructor object
	void (*computeFun)(omxExpectation* ox);										// Wrapper for the Expectation function itself
	void (*printFun)(omxExpectation* ox);										// Prints the appropriate pieces of the expectation
	void (*setFitFun)(omxExpectation ox, omxFitFunction* off);					// To handle fit interactions
	void (*repopulateFun)(omxExpectation* ox, double* x, int n);				// To repopulate any data stored in the Expectation function
	void (*populateAttrFun)(omxExpectation* ox, SEXP algebra);					// Add attributes to the result algebra object
	omxMatrix* (*componentFun)(omxExpectation*, omxFitFunction*, const char*);		// Return component locations to expectation
	void (*mutateFun)(omxExpectation*, omxFitFunction*, const char*, omxMatrix*); // Modify/set/mutate components of expectation
	
	SEXP rObj;																	// Original r Object Pointer
	void* sharedArgs;															// Common argument structure
	void* argStruct;															// Arguments needed for Expectation function
        const char* expType;   // pointer to a static string, no need to allocate or free
	unsigned short int isPrepopulated;											// Object has had some values prepopulated to allow object sharing
	omxData* data;																// Not sure if this is appropriate, but the expectation passes the actual data object
	omxMatrix* dataColumns;
	omxThresholdColumn* thresholds;
	omxDefinitionVar* defVars;
	int numOrdinal, numDefs;
	
	/* Replication of some of the structures from Matrix */
	unsigned short isComplete;													// Whether or not this expectation has been initialize
	omxState* currentState;
	int expNum;

};

/* Initialize and Destroy */
	void omxInitEmptyExpectation(omxExpectation *ox);
	void omxCompleteExpectation(omxExpectation *ox);
	void omxFreeExpectationArgs(omxExpectation* Expectation);					// Frees all args
	omxExpectation* omxNewExpectationFromExpectationIndex(int expIndex, omxState* os);
	omxExpectation* omxNewIncompleteExpectation(SEXP mxobj, int expNum, omxState* os);
	

/* Expectation-specific implementations of matrix functions */
	void omxExpectationRecompute(omxExpectation *ox);
	void omxExpectationCompute(omxExpectation *ox);
	omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState);
	
	void omxExpectationPrint(omxExpectation *source, char* d);					// Pretty-print a (small-form) expectation
	
	omxMatrix* omxGetExpectationComponent(omxExpectation *ox, omxFitFunction *off, char* component);	// Get a component
	
	void omxSetExpectationComponent(omxExpectation *ox, omxFitFunction *off, char* component, omxMatrix *om);	// Set a component

#endif /* _OMXEXPECTATION_H_ */
