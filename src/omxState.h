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

#include "omxDefines.h"

#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <sys/types.h>

#include <time.h>
#include <unistd.h>
#include <string>
#include <stdarg.h>

/* Forward declarations for later includes */
typedef struct omxState omxState;
typedef struct omxFreeVar omxFreeVar;
struct ConfidenceInterval;

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxData.h"

struct omxFreeVarLocation {
	int matrix;
	int row, col;
};

class omxFreeVar {
	int numDeps;            // number of algebra/matrix dependencies
	int *depsPtr;           // indices of algebra/matrix dependencies
 public:
	int id;
	double lbound, ubound;
	std::vector<omxFreeVarLocation> locations;
	const char* name;
	
	// Be aware that a free variable might be assigned to more
	// than 1 location in the same matrix. This API just returns
	// the first matching location.
	const omxFreeVarLocation *getLocation(int matrix) const;
	const omxFreeVarLocation *getLocation(omxMatrix *mat) const;

	const omxFreeVarLocation *getOnlyOneLocation(int matrix, bool &moreThanOne) const;
	const omxFreeVarLocation *getOnlyOneLocation(omxMatrix *mat, bool &moreThanOne) const;

	void setDeps(int _numDeps, int *_deps) {
		numDeps = _numDeps;
		depsPtr = _deps;
	}
	const Eigen::Map< Eigen::VectorXi > getDeps() {
		Eigen::Map< Eigen::VectorXi > map(depsPtr, numDeps);
		return map;
	}
	// Warning: copyToState does not mark matrices dirty
	void copyToState(struct omxState *os, double val);
	void markDirty(omxState *os);
};

#define FREEVARGROUP_ALL      0
#define FREEVARGROUP_NONE    -1
#define FREEVARGROUP_INVALID -2

struct FreeVarGroup {
	std::vector<int> id;              // see omxGlobal::deduplicateVarGroups
	std::vector< omxFreeVar* > vars;

	// see cacheDependencies
	std::vector<bool> dependencies;
	std::vector<bool> locations;

	int lookupVar(const char *name);  // index or -1 if not found
	int lookupVar(int matrix, int row, int col);
	int lookupVar(omxMatrix *matrix, int row, int col);
	//int lookupVar(int id);
	void cacheDependencies(omxState *os);
	void markDirty(omxState *os);
	void log(omxState *os);
	bool hasSameVars(FreeVarGroup *g2);
	bool isDisjoint(FreeVarGroup *other);
};

// These were set in view of NPSOL option Infinite Bound Size. It
// probably makes sense to use std::numeric_limits<double>::max() and
// min().
#define NEG_INF -2e20
#define INF 2e20

class omxConstraint {
 public:
	enum Type {
		LESS_THAN=0,
		EQUALITY,
		GREATER_THAN
	};

	const char *name;
	int size, nrows, ncols;
	enum Type opCode;
	int linear;
	omxMatrix* jacobian;
	std::vector<int> jacMap;

	//Constraints created by backend for CIs use this, the base-class constructor:
        omxConstraint(const char *name) : name(name), linear(0), jacobian(NULL) {};
	virtual ~omxConstraint() {};
	void refreshAndGrab(FitContext *fc, double *out)
	{ refreshAndGrab(fc, opCode, out); };
	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out) = 0;
	virtual omxConstraint *duplicate(omxState *dest)=0;
	virtual void prep(FitContext *fc) {};
};

class UserConstraint : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxMatrix *pad;
	void refresh(FitContext *fc);
	UserConstraint(const char *name) : super(name) {};

 public:
 	//Constraints created from frontend MxConstraints use this, the derived-class constructor:
	UserConstraint(FitContext *fc, const char *name, omxMatrix *arg1, omxMatrix *arg2, omxMatrix *jac, int lin);
	virtual ~UserConstraint();
	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out);
	virtual omxConstraint *duplicate(omxState *dest);
	virtual void prep(FitContext *fc);
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

	void omxWriteCheckpointHeader();

 public:
	omxCheckpointType type;
	time_t timePerCheckpoint;
	int iterPerCheckpoint;
	int evalsPerCheckpoint;
	FILE* file;

	omxCheckpoint();
	void message(FitContext *fc, double *est, const char *msg);
	void postfit(const char *callerName, FitContext *fc, double *est, bool force);
	~omxCheckpoint();
};

struct ConfidenceInterval {
	enum { Lower=0, Upper=1 };
	std::string name;
	int matrixNumber;
	int row, col;		// Location of element to calculate
	bool boundAdj;          // Hu & Neale (2012)
	int varIndex;
	Eigen::Array<double,2,1> bound;		// distance from reference fit
	Eigen::Array<double,2,1> val;		// parameter value at bound
	Eigen::Array<int,2,1>    code;		// optimizer code at bound
	ConfidenceInterval();
	bool isWholeAlgebra() const { return row == -1 && col == -1; }
	omxMatrix *getMatrix(omxState *st) const;
	bool cmpBoundAndType(const ConfidenceInterval &other) {
		return (bound != other.bound).any() || boundAdj != other.boundAdj;
	}
};

// omxGlobal is for state that is read-only during parallel sections.
class omxGlobal {
	bool unpackedConfidenceIntervals;
	std::vector< FreeVarGroup* > freeGroup;
	time_t lastProgressReport;
	int previousReportLength;
	int previousComputeCount;
	double previousReportFit;
	void reportProgressStr(const char *msg);

 public:
	bool silent;
	int numThreads;
	int parallelDiag;
	int analyticGradients;
	double llScale;
	int debugProtectStack;
	int anonAlgebra;
	bool rowLikelihoodsWarning;
	double feasibilityTolerance;
	double optimalityTolerance;
	int majorIterations;
	bool intervals;
	double gradientTolerance;
	int dataTypeWarningCount;

	int maxOrdinalPerBlock;
	double maxptsa;
	double maxptsb;
	double maxptsc;
	int calcNumIntegrationPoints(int numVars) { return maxptsa + maxptsb * numVars + maxptsc * numVars * numVars; };
	double absEps;
	double relEps;

	int RAMInverseOpt;
	int RAMMaxDepth;

	int maxStackDepth;

	std::vector< ConfidenceInterval* > intervalList;
	void unpackConfidenceIntervals(omxState *currentState);
	void omxProcessConfidenceIntervals(SEXP intervalList, omxState *currentState);

	FreeVarGroup *findOrCreateVarGroup(int id);
	FreeVarGroup *findVarGroup(int id);
	bool boundsUpdated;

	// These lists exist only to free memory
	std::vector< omxCompute* > computeList;
	void omxProcessMxComputeEntities(SEXP rObj, omxState *currentState);

	std::vector< std::string > bads;

	// Will need revision if multiple optimizers are running in parallel
	std::vector< omxCheckpoint* > checkpointList;
	FitContext *topFc;

	omxGlobal();
	void deduplicateVarGroups();
	const char *getBads();
	void checkpointMessage(FitContext *fc, double *est, const char *fmt, ...) __attribute__((format (printf, 4, 5)));
	void checkpointPostfit(const char *callerName, FitContext *fc, double *est, bool force);
	double getGradientThreshold(double fit) { return std::max(fabs(fit) * gradientTolerance, .01); }

	void cacheDependencies(omxState *os) {
		for (size_t vg=0; vg < freeGroup.size(); ++vg) {
			freeGroup[vg]->cacheDependencies(os);
		}
	};

	~omxGlobal();
	void reportProgress(const char *context, FitContext *fc);
};

// Use a pointer to ensure correct initialization and destruction
extern struct omxGlobal *Global;

// omxState is for stuff that must be duplicated for thread safety.
class omxState {
 private:
	void init();
	static int nextId;
	int stateId;
	int wantStage; // hack because omxRecompute doesn't take 'want' as a parameter TODO
	bool clone;
 public:
	int getWantStage() const { return wantStage; }
	void setWantStage(int stage);
	int getId() const { return stateId; }
	bool isClone() const { return clone; }

	std::vector< omxMatrix* > matrixList;
	std::vector< omxMatrix* > algebraList;
	std::vector< omxExpectation* > expectationList;
	std::vector< omxData* > dataList;
	std::vector< omxConstraint* > conListX;

 	omxState() : clone(false) { init(); };
	omxState(omxState *src);
	void initialRecalc(FitContext *fc);
	void omxProcessMxMatrixEntities(SEXP matList);
	void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues);
	void omxProcessMxAlgebraEntities(SEXP algList);
	void omxCompleteMxFitFunction(SEXP algList, FitContext *fc);
	void omxProcessConfidenceIntervals(SEXP intervalList);
	void omxProcessMxExpectationEntities(SEXP expList);
	void omxCompleteMxExpectationEntities();
	void omxInitialMatrixAlgebraCompute(FitContext *fc);
	void omxProcessConstraints(SEXP constraints, FitContext *fc);
	void omxProcessMxDataEntities(SEXP data, SEXP defvars);
	omxData* omxNewDataFromMxData(SEXP dataObject, const char *name);
	void loadDefinitionVariables(bool start);
	void omxExportResults(MxRList *out, FitContext *fc);
	void invalidateCache();
	~omxState();

	omxMatrix *lookupDuplicate(omxMatrix *element) const;
	omxMatrix *getMatrixFromIndex(int matnum) const; // matrix (2s complement) or algebra
	omxMatrix *getMatrixFromIndex(omxMatrix *mat) const { return lookupDuplicate(mat); };
	const char *matrixToName(int matnum) const { return getMatrixFromIndex(matnum)->name(); };

	void countNonlinearConstraints(int &equality, int &inequality, bool distinguishLinear)
	{
		equality = 0;
		inequality = 0;
		for(int j = 0; j < int(conListX.size()); j++) {
			omxConstraint *cs = conListX[j];
			if(distinguishLinear && cs->linear){continue;}
			if (cs->opCode == omxConstraint::EQUALITY) {
				equality += cs->size;
			} else {
				inequality += cs->size;
			}
		}
	};
	void countLinearConstraints(int &l_equality, int &l_inequality)
	{
		l_equality = 0;
		l_inequality = 0;
		for(int j = 0; j < int(conListX.size()); j++) {
			omxConstraint *cs = conListX[j];
			if(!cs->linear){continue;}
			if (cs->opCode == omxConstraint::EQUALITY) {
				l_equality += cs->size;
			} else {
				l_inequality += cs->size;
			}
		}
	};
};

inline bool isErrorRaised() { return Global->bads.size() != 0; }
void omxRaiseError(const char* Rf_errorMsg); // DEPRECATED
void omxRaiseErrorf(const char* Rf_errorMsg, ...) __attribute__((format (printf, 1, 2)));

std::string string_vsnprintf(const char *fmt, va_list ap);

void diagParallel(int verbose, const char* msg, ...) __attribute__((format (printf, 2, 3)));
SEXP enableMxLog();

#endif /* _OMXSTATE_H_ */


