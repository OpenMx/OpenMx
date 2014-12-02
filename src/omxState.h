/*
 *  Copyright 2007-2014 The OpenMx Project
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

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
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
#include <stdarg.h>

#include "omxDefines.h"

/* Forward declarations for later includes */
typedef struct omxState omxState;
typedef struct omxFreeVar omxFreeVar;
typedef struct omxConstraint omxConstraint;
typedef struct omxConfidenceInterval omxConfidenceInterval;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxData.h"

struct omxFreeVarLocation {
	int matrix;
	int row, col;
};

struct omxFreeVar {
	int id;
	double lbound, ubound;
	std::vector<omxFreeVarLocation> locations;
	int numDeps;            // number of algebra/matrix dependencies
	int *deps;              // indices of algebra/matrix dependencies
	const char* name;
	
	omxFreeVarLocation *getLocation(int matrix);
};

#define FREEVARGROUP_ALL      0
#define FREEVARGROUP_NONE    -1
#define FREEVARGROUP_INVALID -2

struct FreeVarGroup {
	std::vector<int> id;
	std::vector< omxFreeVar* > vars;

	// see cacheDependencies
	std::vector<bool> dependencies;
	std::vector<bool> locations;

	int lookupVar(const char *name);  // index or -1 if not found
	void cacheDependencies(omxState *os);
	void markDirty(omxState *os);
	void log(omxState *os);
	bool hasSameVars(FreeVarGroup *g2);
	bool isDisjoint(FreeVarGroup *other);
};

#define NEG_INF -2e20
#define INF 2e20

struct omxConstraint {		// Free Variable Constraints
	const char *name;
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

class omxCheckpoint {
	bool wroteHeader;
	time_t lastCheckpoint;	// FIXME: Cannot update at sub-second times.
	int lastIterations;
	int lastEvaluation;
	bool fitPending;

	void omxWriteCheckpointHeader();
	void _prefit(FitContext *fc, double *est, bool force, const char *context);

 public:
	omxCheckpointType type;
	time_t timePerCheckpoint;
	int iterPerCheckpoint;
	int evalsPerCheckpoint;
	FILE* file;

	omxCheckpoint();
	void message(FitContext *fc, double *est, const char *msg);
	void prefit(FitContext *fc, double *est, bool force);
	void postfit(FitContext *fc);
	~omxCheckpoint();
};

struct omxConfidenceInterval {		// For Confidence interval request
	const char *name;
	omxMatrix* matrix;				// The matrix
	int row, col;					// Location of element to calculate
	double ubound;					// Fit-space upper boundary
	double lbound;					// Fit-space lower boundary
	double max;						// Value at upper bound
	double min;						// Value at lower bound
	int lCode;						// Optimizer code at lower bound
	int uCode;						// Optimizer code at upper bound
	unsigned short calcLower;		// Are we currently calculating lbound?
	bool isWholeAlgebra() const { return row == -1 && col == -1; }
};

// omxGlobal is for state that is read-only during parallel sections.
class omxGlobal {
	bool unpackedConfidenceIntervals;
	std::vector< FreeVarGroup* > freeGroup;

 public:
	int ciMaxIterations;
	int numThreads;
	int analyticGradients;
	double llScale;
	int anonAlgebra;
	bool rowLikelihoodsWarning;

	double maxptsa;
	double maxptsb;
	double maxptsc;
	double absEps;
	double relEps;

	int maxStackDepth;

	std::vector< omxConfidenceInterval* > intervalList;
	void unpackConfidenceIntervals();
	void omxProcessConfidenceIntervals(SEXP intervalList, omxState *currentState);

	int computeCount; // protected by openmp atomic

	FreeVarGroup *findOrCreateVarGroup(int id);
	FreeVarGroup *findVarGroup(int id);

	// These lists exist only to free memory
	std::vector< omxCompute* > computeList;
	void omxProcessMxComputeEntities(SEXP rObj, omxState *currentState);

	std::vector< omxAlgebra* > algebraList;

	std::vector< std::string > bads;

	// Will need revision if multiple optimizers are running in parallel
	std::vector< omxCheckpoint* > checkpointList;
	FitContext *fc;

	omxGlobal();
	void deduplicateVarGroups();
	const char *getBads();
	void checkpointMessage(FitContext *fc, double *est, const char *fmt, ...) __attribute__((format (printf, 4, 5)));
	void checkpointPrefit(FitContext *fc, double *est, bool force);
	void checkpointPostfit(FitContext *fc);

	void cacheDependencies(omxState *os) {
		for (size_t vg=0; vg < freeGroup.size(); ++vg) {
			freeGroup[vg]->cacheDependencies(os);
		}
	};

	~omxGlobal();
};

// Use a pointer to ensure correct initialization and destruction
extern struct omxGlobal *Global;

// omxState is for stuff that must be duplicated for thread safety.
class omxState {
 private:
	void init();
	int wantStage;
 public:
	int getWantStage() const { return wantStage; }
	void setWantStage(int stage);

	std::vector< omxMatrix* > matrixList;
	std::vector< omxMatrix* > algebraList;
	std::vector< omxExpectation* > expectationList;
	std::vector< omxData* > dataList;

	int numConstraints;
	int ncnln;                                               // Number of linear and nonlinear constraints
	omxConstraint* conList;											// List of constraints

	long int currentRow; // only used for debugging

	omxState() { init(); };
	omxState(omxState *src);
	void omxProcessMxMatrixEntities(SEXP matList);
	void omxProcessMxAlgebraEntities(SEXP algList);
	void omxCompleteMxFitFunction(SEXP algList);
	void omxProcessConfidenceIntervals(SEXP intervalList);
	void omxProcessMxExpectationEntities(SEXP expList);
	void omxCompleteMxExpectationEntities();
	void omxProcessConstraints(SEXP constraints, FitContext *fc);
	void omxProcessMxDataEntities(SEXP data);
	omxData* omxNewDataFromMxData(SEXP dataObject, const char *name);
	void omxExportResults(MxRList *out);
	~omxState();
};

/* Initialize and Destroy */
omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element);

inline bool isErrorRaised() { return Global->bads.size() != 0; }
void omxRaiseError(const char* Rf_errorMsg); // DEPRECATED
void omxRaiseErrorf(const char* Rf_errorMsg, ...) __attribute__((format (printf, 1, 2)));

void omxStateNextRow(omxState *state);

void mxLog(const char* msg, ...) __attribute__((format (printf, 1, 2)));   // thread-safe
void mxLogBig(const std::string str);
std::string string_snprintf(const char *fmt, ...) __attribute__((format (printf, 1, 2)));
std::string string_vsnprintf(const char *fmt, va_list ap);

#endif /* _OMXSTATE_H_ */


