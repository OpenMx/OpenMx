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
#include <sys/types.h>

#ifdef WIN32

#include <winsock.h>

#else

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#endif

#include <time.h>
#include <unistd.h>
#include "omxDefines.h"

/* Forward declarations for later includes */
typedef struct omxState omxState;
typedef struct omxFreeVar omxFreeVar;
typedef struct omxConstraint omxConstraint;
typedef struct omxCheckpoint omxCheckpoint;
typedef enum omxCheckpointType omxCheckpointType;
typedef struct omxOptimizerState omxOptimizerState;
typedef struct omxConfidenceInterval omxConfidenceInterval;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxObjective.h"
#include "omxData.h"
//#include "omxOptimizer.h"											// omxOptimizer objects coming soon

/* Structure definitions for object evaluation */  // Might be cleaner to give these their own files.
struct omxFreeVar {			// Free Variables
	double lbound, ubound;	// Bounds
	int numLocations;
	int* matrices;			// Matrix numbers.
	int *row, *col;			// Locations for copying.
	const char* name;
};

struct omxConstraint {		// Free Variable Constraints
	int size;
	int opCode;
	double* lbound;
	double* ubound;
	omxMatrix* result;
};

struct omxOptimizerState {			// For hessian or confidence interval computation
	int currentParameter;			// Which parameter is being examined?
	double offset;					// Current offset of optimization
	short int alpha;				// Parameter multiplier
	// Objective should be:  (3.84 - (-2LL))^2 + alpha * parameter
	// Alpha should generally be +1 to minimize parameter -1 to maximize
};

enum omxCheckpointType {
	OMX_FILE_CHECKPOINT = 0,
	OMX_SOCKET_CHECKPOINT = 1,
	OMX_CONNECTION_CHECKPOINT = 2
};

struct omxCheckpoint {
	omxCheckpointType type;
	time_t time;
	int numIterations;
	unsigned long int lastCheckpoint;	// FIXME: Cannot update at sub-second times.
	FILE* file;						// TODO: Maybe make the connection piece a union instead.
	int socket;
	SEXP connection;
	unsigned short int saveHessian;
};

struct omxConfidenceInterval {		// For Confidence interval request
	omxMatrix* matrix;				// The matrix
	int row, col;					// Location of element to calculate
	double ubound;					// Objective-space upper boundary
	double lbound;					// Objective-space lower boundary
	double max;						// Value at upper bound
	double min;						// Value at lower bound
	int lCode;						// Optimizer code at lower bound
	int uCode;						// Optimizer code at upper bound
	unsigned short calcLower;		// Are we currently calculating lbound?
};

#define MAX_STRING_LEN 250

struct omxState {													// The Current State of Optimization

/* Model and Optimizer Pointers */

//	omxOptimizer* optimizer;										// Current Optimizer
	int numMats, numAlgs, numData, numChildren;						// Number of matrices, algebras, and data elements
	int numDynamic; 			                                    // The current number of dynamic matrices
	int maxDynamic;													// The current capacity of dynamic matrices
	omxMatrix** matrixList;											// Model Matrices
	omxMatrix** algebraList;										// Model Algebras
    omxMatrix** dynamicList;                                        // Dynamic matrices and algebras
	omxData** dataList;												// Data Objects
    omxState** childList;                                           // List of child states
    omxState* parentState;                                          // Parent State
    omxMatrix** parentMatrix;                                       // Parent's Matrix List
    omxMatrix** parentAlgebra;                                      // Parent's Algebra List
	omxConstraint* parentConList;									// Parent's Constraint List
                                                                    // TODO: Need a way to deal with unregistered matrices that have free vars
	omxMatrix* objectiveMatrix;										// Objective Algebra

	/* May want to farm these out to the omxObjective object. */
	int numConstraints;
	omxConstraint* conList;											// List of constraints
	int numIntervals;
	int currentInterval;											// The interval currently being calculated
	omxConfidenceInterval* intervalList;							// List of confidence intervals requested

	int numFreeParams;
	omxFreeVar* freeVarList;										// List of Free Variables and where they go.

	/* Saved Optimum State */ // TODO: Rename saved optimum state
	double* optimalValues;											// Values of the free parameters at the optimum value
	double optimum;													// Objective value at last saved optimum
	double* hessian;												// Current hessian storage
	int optimumStatus;												// Optimizer status of last saved optimum (0=converged, 1=green, -1=error, >1=red)
	char optimumMsg[250];											// Status message of last saved optimum
	omxOptimizerState* optimizerState;								// Current optimum parameters for limit computation

/* Current Optimization State (optimizer-specific) */
//	void* optimizerInfo;											// Optimizer specific storage

/* Data members for use by Objective Function and Algebra Calculations */
	long int computeCount;											// How many times have things been evaluated so far?
	long int currentRow;											// If we're calculating row-by-row, what row are we on?

	/* For Checkpointing */
	int majorIteration;												// Major iteration number
	int minorIteration;												// Minor iteration within major iteration
	time_t startTime;												// Time of first computation
	time_t endTime;													// 'Cause we might as well report it
	omxCheckpoint* checkpointList;									// List of checkpoints
	char *chkptText1, *chkptText2;									// Placeholders for checkpointing text
	int numCheckpoints;												// Number of checkpoints

	int statusCode;													// Status code, if appropriate
	char statusMsg[250];											// Status/Error message to report
	double saturatedModel;											// Saturated model likelihood, where applicable

};

/* Initialize and Destroy */
	void omxInitState(omxState* state, omxState *parentState, int numChildren); // Constructor
	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList, omxMatrix** algebraList, omxData** dataList, omxMatrix* objective);
	void omxFreeState(omxState *oo);									// Destructor
	void omxSaveState(omxState *os, double* freeVals, double minimum);	// Saves the current optimization values //TODO: Rename omxSaveState.
	void omxUpdateState(omxState* tgt, omxState* src);					// Updates the tgt state with the contents of src state
    void omxDuplicateState(omxState *tgt, omxState* src, unsigned short fullCopy); 
                                                                        // Duplicates the current state object
	omxState* omxGetState(omxState *os, int stateNum);					// Retrieve a child by number

    void omxSetMajorIteration(omxState *state, int value);				// Recursively set major iteration number
    void omxSetMinorIteration(omxState *state, int value);				// Recursively set minor iteration number

	void omxAddDynamicMatrix(omxState* state, omxMatrix* matrix);
    omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element);
    
	void omxRaiseError(omxState *oo, int errorCode, char* errorMsg);	// Raise an Error
																		// TODO: Move RaiseError to omxOptimizer.

/* Advance a step */
	void omxStateNextRow(omxState *oo);									// Advance Row
	void omxStateNextEvaluation(omxState *oo);							// Advance Evaluation count

	void omxSaveCheckpoint(omxState* os, double* x, double* f);			// Save out checkpoints

#endif /* _OMXSTATE_H_ */


