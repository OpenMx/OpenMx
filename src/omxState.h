/*
 *  Copyright 2007-2020 by the individuals mentioned in the source code history
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

#ifndef u_OMXSTATE_H_
#define u_OMXSTATE_H_

#include "omxDefines.h"

#include <R_ext/Rdynload.h>
#include <sys/types.h>
#include <functional>
#include <time.h>
#include <unistd.h>
#include <string>
#include <stdarg.h>

/* Forward declarations for later includes */
typedef class omxState omxState;
class omxFreeVar;
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

	void setDeps(int u_numDeps, int *u_deps) {
		numDeps = u_numDeps;
		depsPtr = u_deps;
	}
	const Eigen::Map< Eigen::VectorXi > getDeps() {
		Eigen::Map< Eigen::VectorXi > map(depsPtr, numDeps);
		return map;
	}
	// Warning: copyToState does not mark matrices dirty
	void copyToState(omxState *os, double val);
	double getCurValue(omxState *os);
	void markDirty(omxState *os);
};

#define FREEVARGROUP_ALL      0
#define FREEVARGROUP_NONE    -1
#define FREEVARGROUP_INVALID -2

struct FreeVarGroup {
	std::vector<int> id;              // see omxGlobal::deduplicateVarGroups
	std::vector< omxFreeVar* > vars;
  std::map< const char *, int, cstrCmp > byName;

	// see cacheDependencies
	std::vector<bool> dependencies;
	std::vector<bool> locations;

	int lookupVar(const char *name);  // index or -1 if not found
	int lookupVar(int matrix, int row, int col);
	int lookupVar(omxMatrix *matrix, int row, int col);
	//int lookupVar(int id);
	void reIndex();
	void cacheDependencies(omxState *os);
	void markDirty(omxState *os);
	void log(omxState *os);
	bool hasSameVars(FreeVarGroup *g2);
	bool isDisjoint(FreeVarGroup *other);
};

typedef std::function<void(int r, int c, double val)> MatrixStoreFn;

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
  int origSize;  // number of constraints
	int size;      // number of non-redundant constraints
  void setInitialSize(int sz);
  void recalcSize();
	enum Type opCode;
  // NPSOL is documented to have some special support for linear
  // constraints. Rob tried it out with NPSOL version 5 and it
  // didn't seem to work.
  const int linear;
  std::vector<bool> redundant;
  std::vector<bool> seenActive;
  Eigen::ArrayXXd initialJac;
  bool strict;

	//Constraints created by backend for CIs use this, the base-class constructor:
  omxConstraint(const char *name) : name(name), linear(0) {};
	virtual ~omxConstraint() {};
  virtual void getDim(int *rowsOut, int *colsOut) const = 0;
	virtual void refreshAndGrab(FitContext *fc, double *out) = 0;
	virtual void analyticJac(FitContext *fc, MatrixStoreFn out) {}
	virtual omxConstraint *duplicate(omxState *dest) const = 0;
	virtual void prep(FitContext *fc) {}; // only called once, not per-thread
	virtual void preeval(FitContext *fc) {};
  virtual int getVerbose() const { return 0; }
  virtual bool hasAnalyticJac(FitContext *fc) const = 0;
};

class UserConstraint : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxMatrix *pad;
	omxMatrix *jacobian;
	std::vector<int> jacMap;
  int verbose;
	void refresh(FitContext *fc);
	UserConstraint(const char *name) : super(name) {};

 public:
 	//Constraints created from frontend MxConstraints use this, the derived-class constructor:
	UserConstraint(FitContext *fc, const char *name, omxMatrix *arg1, omxMatrix *arg2, omxMatrix *jac, int verbose);
	virtual ~UserConstraint();
  virtual void getDim(int *rowsOut, int *colsOut) const override;
	virtual void refreshAndGrab(FitContext *fc, double *out) override;
	virtual void analyticJac(FitContext *fc, MatrixStoreFn out) override;
	virtual omxConstraint *duplicate(omxState *dest) const override;
	virtual void prep(FitContext *fc) override;
	virtual void preeval(FitContext *fc) override;
  virtual int getVerbose() const override;
  virtual bool hasAnalyticJac(FitContext *fc) const override { return jacobian; }
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
	void message(FitContext *fc, const char *msg);
	void postfit(const char *callerName, FitContext *fc, bool force);
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

enum GradientAlgorithm {
	GradientAlgorithm_Auto,
	GradientAlgorithm_Forward,
	GradientAlgorithm_Central
};

enum OptEngine {
	OptEngine_NPSOL,
	OptEngine_CSOLNP,
    OptEngine_NLOPT,
    OptEngine_SD
};

OptEngine nameToGradOptEngine(const char *engineName);

// omxGlobal is for state that is read-only during parallel sections.
class omxGlobal {
	bool unpackedConfidenceIntervals;
	std::vector< FreeVarGroup* > freeGroup;
	time_t lastProgressReport;
	std::string previousReport;
	int previousComputeCount;
	void reportProgressStr(std::string &str);

 public:
	bool RNGCheckedOut;
	ProtectAutoBalanceDoodad *mpi;
	bool silent;
	bool ComputePersist;
	int numThreads;
	int parallelDiag;
  OptEngine engine;
  GradientAlgorithm gradientAlgo;
  void setDefaultGradientAlgo();
  int gradientIter;
  double gradientStepSize;
	int analyticGradients;
	double llScale;
	int debugProtectStack;
	bool rowLikelihoodsWarning;
	double feasibilityTolerance;
	double optimalityTolerance;
	int majorIterations;
	time_t startTime;
	int maxSeconds;
	bool timedOut;
	bool intervals;
	double gradientTolerance;
	int dataTypeWarningCount;

  int NPSOL_HACK;
	int maxOrdinalPerBlock;
	double maxptsa;
	double maxptsb;
	double maxptsc;
	double maxptsd;
	double maxptse;
	int calcNumIntegrationPoints(int numVars) {
		int pts = (maxptsa + maxptsb * numVars + maxptsc * numVars * numVars +
			   exp(maxptsd + maxptse * numVars * log(relEps)));
		if (pts < 0) mxThrow("calcNumIntegrationPoints "
				     "%f + %f * %d + %f * %d * %d + "
				     "exp(%f + %f * %d * log(%f)) is too large (or non-positive)",
				      maxptsa, maxptsb, numVars, maxptsc, numVars, numVars,
				     maxptsd, maxptse, numVars, log(relEps));
		return pts;
	};
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

	// bundle computeLoop* into a structure?
	std::vector<const char *> computeLoopContext;
	std::vector<int> computeLoopIndex;
	std::vector<int> computeLoopIter;
	std::vector<int> computeLoopMax;
	int lastIndexDone;
	time_t lastIndexDoneTime;

	std::vector< std::string > checkpointColnames;
	std::vector< std::string > checkpointValues;
	std::vector< std::string > bads;
	bool userInterrupted;

	// Will need revision if multiple optimizers are running in parallel
	std::vector< omxCheckpoint* > checkpointList;
  Eigen::VectorXd startingValues;
	FitContext *topFc;

	omxGlobal();
	void deduplicateVarGroups();
	const char *getBads();
	void checkpointMessage(FitContext *fc, const char *fmt, ...) __attribute__((format (printf, 3, 4)));
	void checkpointPostfit(const char *callerName, FitContext *fc, bool force);
	double getGradientThreshold(double fit) {
		return( pow(optimalityTolerance, 1.0/3.0) * (1.0 + fabs(fit)) );
	}

	void cacheDependencies(omxState *os) {
		for (size_t vg=0; vg < freeGroup.size(); ++vg) {
			freeGroup[vg]->cacheDependencies(os);
		}
	};

	~omxGlobal();
	void reportProgress(const char *context, FitContext *fc);
	bool interrupted();
	void reportProgress1(const char *context, std::string detail);
	void throwOnUserInterrupted() { if (interrupted()) mxThrow("User interrupt"); };
};

// Use a pointer to ensure correct initialization and destruction
extern class omxGlobal *Global;

void diagParallel(int verbose, const char* msg, ...) __attribute__((format (printf, 2, 3)));

// omxState is for stuff that must be duplicated for thread safety.
class omxState {
 private:
	void init();
	static int nextId;
	int stateId;
	int wantStage; // hack because omxRecompute doesn't take 'want' as a parameter TODO
	omxState *parent;     // for read-only access to shared state
	omxState *workBoss; // for OpenMP when multiple omxState want to cooperate
	bool hasFakeParam;
 public:
	int getWantStage() const { return wantStage; }
	void setWantStage(int stage);
	int getId() const { return stateId; }
	bool isClone() const { return workBoss != 0; } // rename to isWorkBoss TODO
  bool isTopState() const { return parent == 0; }

	std::vector< omxMatrix* > matrixList;
	std::vector< omxMatrix* > algebraList;
	std::vector< omxExpectation* > expectationList;
	std::vector< omxData* > dataList;
	std::vector< omxConstraint* > conListX;
  void constraintPreeval(FitContext *fc) { for (auto c1 : conListX) c1->preeval(fc); }

	omxState() : wantStage(0), parent(0), workBoss(0), hasFakeParam(false) { init(); };
	omxState(omxState *src, bool isTeam);
	void initialRecalc(FitContext *fc);
	void omxProcessMxMatrixEntities(SEXP matList);
	void omxProcessFreeVarList(SEXP varList);
	void omxProcessMxAlgebraEntities(SEXP algList);
	void omxCompleteMxFitFunction(SEXP algList, FitContext *fc);
	void omxProcessConfidenceIntervals(SEXP intervalList);
	void omxProcessMxExpectationEntities(SEXP expList);
	void omxCompleteMxExpectationEntities();
	void omxInitialMatrixAlgebraCompute(FitContext *fc);
	void omxProcessConstraints(SEXP constraints, FitContext *fc);
  void hideBadConstraints(FitContext *fc);
	void omxProcessMxDataEntities(SEXP data, SEXP defvars);
	omxData* omxNewDataFromMxData(SEXP dataObject, const char *name);
	void loadDefinitionVariables(bool start);
	void omxExportResults(MxRList *out, FitContext *fc);
  void reportConstraints(MxRList &out);
	void invalidateCache();
	void connectToData();
	~omxState();

	omxExpectation *getParent(omxExpectation *element) const;
	omxExpectation *lookupDuplicate(omxExpectation *element) const;
	omxMatrix *lookupDuplicate(omxMatrix *element) const;
	omxMatrix *getMatrixFromIndex(int matnum) const; // matrix (2s complement) or algebra
	omxMatrix *getMatrixFromIndex(omxMatrix *mat) const { return lookupDuplicate(mat); };
	const char *matrixToName(int matnum) const { return getMatrixFromIndex(matnum)->name(); };

	bool isFakeParamSet() const { return hasFakeParam; }

	template <typename T1>
	void setFakeParam(Eigen::MatrixBase<T1> &point)
	{
		if (hasFakeParam) mxThrow("already has fake parameters loaded");
		hasFakeParam = true;

		auto varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
		size_t numParam = varGroup->vars.size();
		point.derived().resize(numParam);

		for(size_t k = 0; k < numParam; k++) {
			omxFreeVar* freeVar = varGroup->vars[k];
			point[k] = freeVar->getCurValue(this);
			freeVar->copyToState(this, 1.0);
		}
	}

	void restoreParam(const Eigen::Ref<const Eigen::VectorXd> point);

};

#include "autoTune.h"
#include "finiteDifferences.h"

// NOTE: All non-linear constraints are applied regardless of free
// variable group.
class ConstraintVec {
public:
  typedef std::function<bool(const omxConstraint&)> ClassifyFn;
  int verbose;
  bool verifyJac;

private:
  const char *name;
  ClassifyFn cf;
  int count;
  bool ineqAlwaysActive;
  bool anyAnalyticJacCache;
  std::unique_ptr<class AutoTune<class JacobianGadget>> jacTool;

public:
  ConstraintVec(FitContext *fc, const char *u_name, ClassifyFn u_cf);
  int getCount() const { return count; }
  void setIneqAlwaysActive() { ineqAlwaysActive = true; }
  void markUselessConstraints(FitContext *fc);
  void allocJacTool(FitContext *fc);
  void setAlgoOptions(GradientAlgorithm u_algo,	int u_numIter, double u_eps)
  { jacTool->work().setAlgoOptions(u_algo, u_numIter, u_eps); }
  // constrOut is a vector of size getCount
  // jacOut is a matrix of dimension getCount by numFree
  void eval(FitContext *fc, double *constrOut, double *jacOut);
  void eval(FitContext *fc, double *constrOut) { eval(fc, constrOut, 0); }
  bool anyAnalyticJac() const { return anyAnalyticJacCache; }
};

inline bool isErrorRaised() { return Global->bads.size() != 0 || Global->userInterrupted || Global->timedOut; }
inline bool isErrorRaisedIgnTime() { return Global->bads.size() != 0 || Global->userInterrupted; }

// rename from "Raise" to "Record" since no exception is thrown
void omxRaiseError(const char* msg); // DEPRECATED
void omxRaiseErrorf(const char* fmt, ...) __attribute__((format (printf, 1, 2)));

void string_vsnprintf(const char *fmt, va_list orig_ap, std::string &dest);

SEXP enableMxLog();

struct StateInvalidator {
	omxState &st;
	StateInvalidator(omxState &u_st) : st(u_st) {};
	virtual void doData();
	virtual void doMatrix();
	virtual void doExpectation();
	virtual void doAlgebra();
	void operator()()
	{
		doData();
		doMatrix();
		doExpectation();
		doAlgebra();
	}
};

#endif /* u_OMXSTATE_H_ */
