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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <sys/types.h>

#ifdef WIN32

#else

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#endif

#include <time.h>
#include <unistd.h>
#include <string>

#include "omxDefines.h"

/* Forward declarations for later includes */
typedef struct omxState omxState;
typedef struct omxFreeVar omxFreeVar;
typedef struct omxConstraint omxConstraint;
typedef struct omxCheckpoint omxCheckpoint;
typedef struct omxConfidenceInterval omxConfidenceInterval;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxData.h"

#include <vector>

struct omxFreeVarLocation {
	int matrix;
	int row, col;
};

struct omxFreeVar {			// Free Variables
	double lbound, ubound;	// Bounds
	std::vector<omxFreeVarLocation> locations;
	int numDeps;            // number of algebra/matrix dependencies
	int *deps;              // indices of algebra/matrix dependencies
	const char* name;
};

struct omxConstraint {		// Free Variable Constraints
	int size;
	int opCode;
	double* lbound;
	double* ubound;
	omxMatrix* result;
};

enum omxCheckpointType {
	OMX_FILE_CHECKPOINT,
	OMX_CONNECTION_CHECKPOINT
};
typedef enum omxCheckpointType omxCheckpointType;

struct omxCheckpoint {
	omxCheckpointType type;
	time_t time;
	int numIterations;
	time_t lastCheckpoint;	// FIXME: Cannot update at sub-second times.
	FILE* file;						// TODO: Maybe make the connection piece a union instead.
	SEXP connection;
	unsigned short int saveHessian;
};

struct omxConfidenceInterval {		// For Confidence interval request
	omxMatrix* matrix;				// The matrix
	int row, col;					// Location of element to calculate
	double ubound;					// Fit-space upper boundary
	double lbound;					// Fit-space lower boundary
	double max;						// Value at upper bound
	double min;						// Value at lower bound
	int lCode;						// Optimizer code at lower bound
	int uCode;						// Optimizer code at upper bound
	unsigned short calcLower;		// Are we currently calculating lbound?
};

#define MAX_STRING_LEN 250

// omxGlobal is for state that is read-only during parallel sections.
struct omxGlobal {
	int ciMaxIterations;
	int numThreads;
	int analyticGradients;
	int numChildren;
};
extern struct omxGlobal Global;

// omxState is for stuff that duplication for thread safety.
struct omxState {
	std::vector< omxMatrix* > matrixList;
	std::vector< omxMatrix* > algebraList;
	std::vector< omxExpectation* > expectationList;
	std::vector< omxData* > dataList;
	std::vector< class omxCompute* > computeList;
	omxState** childList;											// List of child states
	omxState* parentState;											// Parent State

	/* May want to farm these out to the omxFitFunction object. */
	int numConstraints;
	int nclin, ncnln;                                               // Number of linear and nonlinear constraints
	omxConstraint* conList;											// List of constraints
	int numIntervals;
	omxConfidenceInterval* intervalList;							// List of confidence intervals requested

	int numFreeParams;
	omxFreeVar* freeVarList;										// List of Free Variables and where they go.

/* Data members for use by Fit Function and Algebra Calculations */
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

	char statusMsg[MAX_STRING_LEN];
};

extern omxState* globalState;

/* Initialize and Destroy */
	void omxInitState(omxState* state, omxState *parentState);
	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList, omxMatrix** algebraList, omxData** dataList, omxMatrix* fitFunction);
void omxFreeChildStates(omxState *state);
void omxFreeState(omxState *state);
	void omxDuplicateState(omxState *tgt, omxState* src); 
                                                                        // Duplicates the current state object
	void omxSetMajorIteration(omxState *state, int value);				// Recursively set major iteration number
	void omxSetMinorIteration(omxState *state, int value);				// Recursively set minor iteration number

	omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element);

	void omxResetStatus(omxState *state);    
inline bool isErrorRaised(omxState *state) { return state->statusMsg[0] != 0; }
void omxRaiseError(omxState *state, int errorCode, const char* errorMsg); // DEPRECATED
void omxRaiseErrorf(omxState *state, const char* errorMsg, ...);
																		// TODO: Move RaiseError to omxOptimizer.

/* Advance a step */
	void omxStateNextRow(omxState *state);								// Advance Row
	void omxStateNextEvaluation(omxState *state);						// Advance Evaluation count

	void omxWriteCheckpointMessage(omxState *os, char *msg);
	void omxSaveCheckpoint(omxState* os, double* x, double* f, int force);	// Save out checkpoints
void omxExamineFitOutput(omxState *state, omxMatrix *fitMatrix, int *mode);

void mxLog(const char* msg, ...);   // thread-safe
void mxLogBig(const std::string str);
std::string string_snprintf(const std::string fmt, ...);

#endif /* _OMXSTATE_H_ */


