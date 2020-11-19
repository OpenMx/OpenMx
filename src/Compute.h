/*
 *  Copyright 2013-2019 by the individuals mentioned in the source code history
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

#ifndef _OMX_COMPUTE_H_
#define _OMX_COMPUTE_H_

#include "omxDefines.h"
#include <map>
#include <Eigen/SparseCore>
#include "glue.h"
#include "omxState.h"

// See R/MxRunHelperFunctions.R optimizerMessages
// Also see the NPSOL manual, section 6 "Error Indicators and Warnings"
// These are ordered from good to bad so we can use max() on a set
// of inform results to obtain a bound on convergence status.
typedef int ComputeInform;
#define INFORM_UNINITIALIZED NA_INTEGER
#define INFORM_CONVERGED_OPTIMUM 0
#define INFORM_UNCONVERGED_OPTIMUM 1
	// The final iterate satisfies the optimality conditions to the accuracy requested,
	// but the sequence of iterates has not yet converged.
	// Optimizer was terminated because no further improvement
	// could be made in the merit function (Mx status GREEN).
#define INFORM_LINEAR_CONSTRAINTS_INFEASIBLE 2
	// The linear constraints and bounds could not be satisfied.
	// The problem has no feasible solution.
#define INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE 3
	// The nonlinear constraints and bounds could not be satisfied.
	// The problem may have no feasible solution.
#define INFORM_ITERATION_LIMIT 4
	// The major iteration limit was reached (Mx status BLUE).
#define INFORM_NOT_CONVEX 5
        // Hessian is not positive definite (Mx status RED)
#define INFORM_NOT_AT_OPTIMUM 6
	// The model does not satisfy the first-order optimality conditions (i.e. gradient close to zero)
	// to the required accuracy, and no improved point for the
	// merit function could be found during the final linesearch (Mx status RED)
#define INFORM_BAD_DERIVATIVES 7
	// The function derivates returned by funcon or funobj
	// appear to be incorrect.
#define INFORM_INVALID_PARAM 9
	// An input parameter was invalid'
#define INFORM_STARTING_VALUES_INFEASIBLE 10

SEXP allocInformVector(int size);

enum ComputeInfoMethod {
	INFO_METHOD_DEFAULT,
	INFO_METHOD_HESSIAN,
	INFO_METHOD_SANDWICH,
	INFO_METHOD_BREAD,
	INFO_METHOD_MEAT
};

struct HessianBlock {
	//private:
	Eigen::MatrixXd mmat;   // including subblocks
	std::vector< HessianBlock* > subBlocks;
	bool merge;
	int useId;
	void addSubBlocks();

	//public:
	std::vector<int> vars;  // global freeVar ID in order
	Eigen::MatrixXd mat;    // vars * vars, only upper triangle referenced
	Eigen::MatrixXd imat;

	HessianBlock() : merge(false), useId(0) {}
	HessianBlock *clone();
	int estNonZero() const;
};

struct CIobjective {
	const ConfidenceInterval *CI;

	enum Diagnostic {
		DIAG_SUCCESS=1,
		DIAG_ALPHA_LEVEL,
		DIAG_BA_D1, DIAG_BA_D2,
		DIAG_BN_D1, DIAG_BN_D2,
		DIAG_BOUND_INFEASIBLE,
		DIAG_BOXED
	};

  CIobjective(const ConfidenceInterval *_CI) : CI(_CI) {}
  virtual std::unique_ptr<CIobjective> clone() const = 0;
	virtual bool gradientKnown() { return false; };
	virtual void gradient(FitContext *fc) {};
	virtual void evalIneq(FitContext *fc, omxMatrix *fitMat, double *out) {};
	virtual void evalEq(FitContext *fc, omxMatrix *fitMat, double *out) {};
	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc);
	virtual void checkSolution(FitContext *fc);
	virtual Diagnostic getDiag() = 0;
};

// The idea of FitContext is to eventually enable fitting from
// multiple starting values in parallel. There is one FitContext
// per thread.

class FitContext {
	static omxFitFunction *RFitFunction;

	FitContext *parent;

	std::vector<HessianBlock*> allBlocks;
	std::vector<HessianBlock*> mergeBlocks;
	std::vector<HessianBlock*> blockByVar;

	bool haveSparseHess;
	Eigen::SparseMatrix<double> sparseHess;
	bool haveSparseIHess;
	Eigen::SparseMatrix<double> sparseIHess;
	int estNonZero;
	int minBlockSize;
	int maxBlockSize;

	bool haveDenseHess;
	bool haveDenseIHess;

	void init();
	void analyzeHessian();
	void analyzeHessianBlock(HessianBlock *hb);
	void testMerge();
  void createChildren1();

	std::string IterationError;
	double ordinalRelativeError;
	int computeCount;
	ComputeInform inform;
	double previousReportFit;
  int _numFree;

 public:
	FreeVarGroup *varGroup;
	omxState *state;
	omxState *getParentState() const { return parent->state; };
	bool isClone() const;
	int numParam;
  int getNumFree() const { return _numFree; }
	std::vector<int> mapToParent;
	double mac;
	double fit;
	FitStatisticUnits fitUnits;
	int skippedRows;
  std::map< const char *, int, cstrCmp > freeToIndexMap;
  std::vector<int> freeToParamMap;  // length=numFree
	std::vector<bool> profiledOutZ;
	void calcNumFree();
  Eigen::ArrayXd est;  // length=numParam
  std::unique_ptr<class JacobianGadget> numericalGradTool;
	Eigen::VectorXd gradZ;  // length=numFree
	void initGrad() {
    int pars = getNumFree();
		gradZ = Eigen::VectorXd::Zero(pars);
	};
	int infoDefinite;
	double infoCondNum;
	Eigen::MatrixXd hess;
	Eigen::MatrixXd ihess;
	Eigen::VectorXd stderrs;   // plural to distinguish from stdio's stderr
	enum ComputeInfoMethod infoMethod;
	double *infoA; // sandwich, the bread
	double *infoB; // sandwich, the meat
	int iterations;
	int wanted;
	std::vector< class FitContext* > childList;

	//Outputs from gradient-based optimizers:
	Eigen::VectorXd constraintFunVals;
	Eigen::MatrixXd constraintJacobian;
	Eigen::VectorXd LagrMultipliers;
	Eigen::VectorXi constraintStates;
	Eigen::MatrixXd LagrHessian;

	//Constraint-related:
	void solEqBFun(bool wantAJ, int verbose);
	void myineqFun(bool wantAJ, bool alwaysActive);
	bool isUsingAnalyticJacobian(){ return state->usingAnalyticJacobian; }
	Eigen::MatrixXd analyticEqJacTmp; //<--temporarily holds analytic Jacobian (if present) for an equality constraint
	Eigen::MatrixXd analyticIneqJacTmp; //<--temporarily holds analytic Jacobian (if present) for an inequality constraint
	Eigen::VectorXd equality;
	Eigen::VectorXd inequality;
  bool isUnconstrained();
	void allConstraintsF(bool wantAJ, int verbose, int ineqType, bool ineqAlwaysActive, bool maskInactive);
	Eigen::MatrixXd vcov; //<--Repeated-sampling covariance matrix of the MLEs.
	int redundantEqualities;

	// for confidence intervals
  std::unique_ptr<CIobjective> ciobj;

	FitContext(omxState *_state);
	FitContext(FitContext *parent, FreeVarGroup *group);
	bool openmpUser;  // whether some fitfunction/expectation uses OpenMP
  bool permitParallel; // whether openmpUser is permitted
  void createChildren(omxMatrix *alg=0, bool permitParallel=false);
  int numOptimizerThreads() { return (childList.size() && !openmpUser)? childList.size() : 1; }
	void destroyChildren();
	void calcStderrs();
	void ensureParamWithinBox(bool nudge);
	void copyParamToModel();
	void copyParamToModelClean();
	omxMatrix *lookupDuplicate(omxMatrix* element);
	void maybeCopyParamToModel(omxState* os);
	void updateParent();
	void updateParentAndFree();
	template <typename T> void moveInsideBounds(std::vector<T> &prevEst);
	void log(int what);
	void resetToOriginalStarts();
	void setInform(int _in) { inform = _in; };
	int getInform() { return inform; };
	int wrapInform() {
		if (inform == INFORM_UNINITIALIZED) return NA_INTEGER;
		return 1 + inform;
	};
	bool haveReferenceFit(omxMatrix *fitMat) {
		if (std::isfinite(fit)) return true;
		if (inform == INFORM_UNINITIALIZED) {
			omxRecompute(fitMat, this);
			fit = omxMatrixElement(fitMat, 0, 0);
			if (std::isfinite(fit)) return true;
			setInform(INFORM_STARTING_VALUES_INFEASIBLE);
		}
		if (inform != INFORM_CONVERGED_OPTIMUM &&
		    inform != INFORM_UNCONVERGED_OPTIMUM) {
			return false;
		}
		mxThrow("%s: reference fit is not finite", fitMat->name());
	};
	~FitContext();

	// deriv related
	void clearHessian();
	void negateHessian();
	void queue(HessianBlock *hb);
	void refreshDenseHess();
	void refreshDenseIHess();
	void refreshSparseHess();
	bool refreshSparseIHess(); // NOTE: produces an ihess with filtered eigenvalues
	Eigen::VectorXd ihessGradProd();
	double *getDenseHessUninitialized();
	double *getDenseIHessUninitialized();
	double *getDenseHessianish();  // either a Hessian or inverse Hessian, remove TODO
	int getDenseHessianishSize();
	void copyDenseHess(double *dest);
	void copyDenseIHess(double *dest);
	Eigen::VectorXd ihessDiag();
	void preInfo();
	void postInfo();
	void resetIterationError();
	void recordIterationError(const char* msg, ...) __attribute__((format (printf, 2, 3)));
	void recordOrdinalRelativeError(double re) {
		// Could obtain NaN if density is exactly zero
		if (!std::isfinite(re) || re < ordinalRelativeError) return;
		ordinalRelativeError = re;
	};
	void resetOrdinalRelativeError();
	double getOrdinalRelativeError() const { return ordinalRelativeError; }

	int getGlobalComputeCount(); //approximate
	int getLocalComputeCount(); //approximate
	void incrComputeCount() { ++computeCount; };

	// If !std::isfinite(fit) then IterationError.size() should be nonzero but not all of
	// the code is audited to ensure that this condition is true.
	bool outsideFeasibleSet() const { return !std::isfinite(fit) || IterationError.size() > 0; }
	// Only check at the end of optimization
	bool insideFeasibleSet() const { return !outsideFeasibleSet() && skippedRows == 0; }

	std::string getIterationError();

	static void cacheFreeVarDependencies();
	static void setRFitFunction(omxFitFunction *rff);

	// profiledOut parameters are not subject to optimization
	template<typename T> void copyEstToOptimizer(T &out) {
		for (int vx=0; vx < getNumFree(); ++vx) {
			out[vx] = est[ freeToParamMap[vx] ];
		}
	};
	template<typename T> void copyGradToOptimizer(T &out) {
		for (int vx=0; vx < getNumFree(); ++vx) {
			out[vx] = gradZ[vx];
		}
	};
	template<typename T> void copyGradFromOptimizer(const T &in) {
		for (int vx=0; vx < getNumFree(); ++vx) {
			gradZ[vx] = in[vx];
		}
	};
	void setParamFromOptimizer(int px, double v1) {
    est[ freeToParamMap[px] ] = v1;
		copyParamToModel(); // could be smarter TODO
	}
	template<typename T> void setEstFromOptimizer(const T &in) {
		for (int vx=0; vx < getNumFree(); ++vx) {
			est[ freeToParamMap[vx] ] = in[vx];
		}
		copyParamToModel();
	};
	template<typename T> void setEstFromOptimizerClean(const T &in) {
		for (int vx=0; vx < getNumFree(); ++vx) {
			est[ freeToParamMap[vx] ] = in[vx];
		}
		copyParamToModelClean();
	};
	template<typename T1, typename T2>
	void setEstGradFromOptimizer(const T1 &ein, const T2 &gin) {
    setEstFromOptimizer(ein);
    if (gradZ.size()) copyGradFromOptimizer(gin);
		copyParamToModel();
	};
	template<typename T1, typename T2>
	void copyBoxConstraintToOptimizer(Eigen::MatrixBase<T1> &lb, Eigen::MatrixBase<T2> &ub)
	{
    // lb.size can be larger than numFree (see NPSOL)
    if (lb.size() < getNumFree()) OOPS;
		for (int vx=0; vx < getNumFree(); ++vx) {
      int px = freeToParamMap[vx];
			lb[vx] = varGroup->vars[px]->lbound;
			if (!std::isfinite(lb[vx])) lb[vx] = NEG_INF;
			ub[vx] = varGroup->vars[px]->ubound;
			if (!std::isfinite(ub[vx])) ub[vx] = INF;
		}
	};
  bool hasActiveBoxConstraint(int skip) const;
	std::string asProgressReport();
};

void copyParamToModelInternal(FreeVarGroup *varGroup, omxState *os, double *at);

typedef std::vector< std::pair<int, MxRList*> > LocalComputeResult;

struct allconstraints_functional {
	FitContext &fc;
	int verbose;

	allconstraints_functional(FitContext &_fc, int _verbose) : fc(_fc), verbose(_verbose) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		fc.setEstFromOptimizer(x.derived().data());
		fc.allConstraintsF(false, verbose, omxConstraint::LESS_THAN, false, true);
		result = fc.constraintFunVals;
	}

	template <typename T1, typename T2, typename T3>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result, Eigen::MatrixBase<T3> &jacobian) const {
		fc.setEstFromOptimizer(x.derived().data());
		fc.allConstraintsF(true, verbose, omxConstraint::LESS_THAN, false, true);
		result = fc.constraintFunVals;
		jacobian = fc.constraintJacobian;
	}
};

class omxCompute {
	int computeId;
	bool dotPersist;
 protected:
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *glob) {};
	void collectResultsHelper(FitContext *fc, std::vector< omxCompute* > &clist,
				  LocalComputeResult *lcr, MxRList *out);
	static enum ComputeInfoMethod stringToInfoMethod(const char *iMethod);
	void complainNoFreeParam();
 public:

	const char *name;
	FreeVarGroup *varGroup;
	omxCompute();
	virtual bool accumulateInform() { return true; };
        virtual void initFromFrontend(omxState *, SEXP rObj);
        void compute(FitContext *fc);
	void computeWithVarGroup(FitContext *fc);
        virtual void computeImpl(FitContext *fc) {}
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
        virtual ~omxCompute();
	void reportProgress(FitContext *fc) { Global->reportProgress(name, fc); }
	bool isPersist() { return dotPersist; };
};

class PushLoopIndex {
	void init(const char *name, int ix, int ii, int last)
	{
		Global->computeLoopContext.push_back(name);
		Global->computeLoopIndex.push_back(ix);
		Global->computeLoopIter.push_back(ii);
		Global->computeLoopMax.push_back(last);
	}
public:
	PushLoopIndex(const char *name, int ix, int last)
	{ init(name, ix, ix, last); }
	PushLoopIndex(const char *name, int ix, int ii, int last)
	{ init(name, ix, ii, last); }
	PushLoopIndex(const char *name)
	{ init(name, NA_INTEGER, 0, 0); }
	~PushLoopIndex() {
		Global->computeLoopContext.pop_back();
		Global->computeLoopIndex.pop_back();
		Global->computeLoopIter.pop_back();
		Global->computeLoopMax.pop_back();
	}
};

class RNGStateManager {
 protected:
	void checkOut()
	{
		if (Global->RNGCheckedOut) {
			mxThrow("Attempt to check out RNG but already checked out");
		}
		GetRNGstate();
		Global->RNGCheckedOut = true;
	}
	void checkIn()
	{
		if (!Global->RNGCheckedOut) {
			mxThrow("Attempt to return RNG but already returned");
		}
		PutRNGstate();
		Global->RNGCheckedOut = false;
	}
};

struct BorrowRNGState : public RNGStateManager {
	BorrowRNGState() { checkOut(); }
	~BorrowRNGState() { checkIn(); }
};

struct ReturnRNGState : public RNGStateManager {
	ReturnRNGState() { checkIn(); }
	~ReturnRNGState() { checkOut(); }
};

omxCompute *omxNewCompute(omxState* os, const char *type);

omxCompute *newComputeGradientDescent();
omxCompute *newComputeNumericDeriv();
omxCompute *newComputeNewtonRaphson();
omxCompute *newComputeConfidenceInterval();
omxCompute *newComputeTryHard();
omxCompute *newComputeNelderMead();
omxCompute *newComputeGenSA();

void omxApproxInvertPosDefTriangular(int dim, double *hess, double *ihess, double *stress);
void omxApproxInvertPackedPosDefTriangular(int dim, int *mask, double *packedHess, double *stress);
SEXP sparseInvert_wrapper(SEXP mat);

inline double addSkippedRowPenalty(double orig, int skipped) // orig does not have * Global->llScale
{
	orig -= 745 * skipped;  // a bit more than log(4.94066e-324) per row
	//
	// This is a good approximation of infinity because the
	// log of any number represented as a double will be
	// of smaller magnitude.
	//
	// http://www.cplusplus.com/forum/general/53760/
	//
	orig *= (1+skipped);    // add some extra badness
	return orig;
}

template <typename T>
void printSparse(Eigen::SparseMatrixBase<T> &sm) {
	typedef typename T::Index Index;
	typedef typename T::Storage Storage;
	// assume column major
	std::string buf;
	const Index *nzp = sm.derived().innerNonZeroPtr();
	//const Scalar *vp = sm.derived().valuePtr();
	//const Index *iip = sm.derived().innerIndexPtr();
	const Index *oip = sm.derived().outerIndexPtr();
	const Storage &m_data = sm.derived().data();
	if (!nzp) buf += "compressed ";
	buf += string_snprintf("%dx%d\n", sm.innerSize(), sm.outerSize());
	for (int rx=0; rx < sm.innerSize(); ++rx) {
		for (int cx=0; cx < sm.outerSize(); ++cx) {
			Index start = oip[cx];
			Index end = nzp ? oip[cx] + nzp[cx] : oip[cx+1];
			if (end <= start) {
				buf += " ***";
			} else {
				const Index p = m_data.searchLowerIndex(start,end-1,rx);
				if ((p<end) && (m_data.index(p)==rx)) {
					double v = m_data.value(p);
					if (v < 0) {
						buf += string_snprintf("%2.1f", v);
					} else {
						buf += string_snprintf(" %2.1f", v);
					}
				}
				else
					buf += " ***";
			}
			if (cx < sm.outerSize() - 1) buf += " ";
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

void AddLoadDataProvider(double version, unsigned int, class LoadDataProviderBase2 *ldp);

#endif
