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

/*File created 2013*/

#ifndef _OMX_COMPUTE_H_
#define _OMX_COMPUTE_H_

#include "omxDefines.h"
#include <Eigen/SparseCore>
#include "glue.h"
#include "omxState.h"

enum GradientAlgorithm {
	GradientAlgorithm_Forward,
	GradientAlgorithm_Central
};

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
	ConfidenceInterval *CI;

	enum Diagnostic {
		DIAG_SUCCESS=1,
		DIAG_ALPHA_LEVEL,
		DIAG_BA_D1, DIAG_BA_D2,
		DIAG_BN_D1, DIAG_BN_D2,
		DIAG_BOUND_INFEASIBLE,
		DIAG_BOXED
	};

	virtual bool gradientKnown() { return false; };
	virtual void gradient(FitContext *fc, double *gradOut) {};
	virtual void evalIneq(FitContext *fc, omxMatrix *fitMat, double *out) {};
	virtual void evalEq(FitContext *fc, omxMatrix *fitMat, double *out) {};
	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc);
	virtual void checkSolution(FitContext *fc);
	virtual Diagnostic getDiag() = 0;
};

// The idea of FitContext is to eventually enable fitting from
// multiple starting values in parallel.

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
	Eigen::MatrixXd hess;
	bool haveDenseIHess;
	Eigen::MatrixXd ihess;

	void init();
	void analyzeHessian();
	void analyzeHessianBlock(HessianBlock *hb);
	void testMerge();

	std::string IterationError;
	int computeCount;
	ComputeInform inform;

 public:
	FreeVarGroup *varGroup;
	omxState *state;
	omxState *getParentState() const { return parent->state; };
	bool isClone() const;
	size_t numParam;               // change to int type TODO
	std::vector<int> mapToParent;
	double mac;
	double fit;
	FitStatisticUnits fitUnits;
	int skippedRows;
	double *est;
	Eigen::Map< Eigen::VectorXd > getEst() {
		Eigen::Map< Eigen::VectorXd > vec(est, numParam);
		return vec;
	}
	std::vector<bool> profiledOut;
	Eigen::VectorXd grad;
	int infoDefinite;
	double infoCondNum;
	double *stderrs;   // plural to distinguish from stdio's stderr
	enum ComputeInfoMethod infoMethod;
	double *infoA; // sandwich, the bread
	double *infoB; // sandwich, the meat
	int iterations;
	int wanted;
	std::vector< class FitContext* > childList;
	Eigen::VectorXd constraintFunVals;
	Eigen::MatrixXd constraintJacobian;
	Eigen::VectorXd LagrMultipliers;
	Eigen::VectorXi constraintStates;
	Eigen::MatrixXd LagrHessian;

	// for confidence intervals
	CIobjective *ciobj;

	FitContext(omxState *_state, std::vector<double> &startingValues);
	FitContext(FitContext *parent, FreeVarGroup *group);
	bool openmpUser;  // whether some fitfunction/expectation uses OpenMP
	void createChildren(omxMatrix *alg);
	void destroyChildren();
	void allocStderrs();
	void ensureParamWithinBox(bool nudge);
	void copyParamToModel();
	void copyParamToModelClean();
	double *take(int want);
	omxMatrix *lookupDuplicate(omxMatrix* element);
	void maybeCopyParamToModel(omxState* os);
	void updateParent();
	void updateParentAndFree();
	template <typename T> void moveInsideBounds(std::vector<T> &prevEst);
	void log(int what);
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
		Rf_error("%s: reference fit is not finite", fitMat->name());
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
	void copyDenseHess(double *dest);
	void copyDenseIHess(double *dest);
	Eigen::VectorXd ihessDiag();
	void preInfo();
	void postInfo();
	void resetIterationError();
	void recordIterationError(const char* msg, ...) __attribute__((format (printf, 2, 3)));
	int getGlobalComputeCount(); //approximate
	int getLocalComputeCount(); //approximate
	void incrComputeCount() { ++computeCount; };

	// If !std::isfinite(fit) then IterationError.size() should be nonzero but not all of
	// the code is audited to ensure that this condition is true.
	bool outsideFeasibleSet() const { return !std::isfinite(fit) || IterationError.size() > 0; }

	std::string getIterationError();

	static void cacheFreeVarDependencies();
	static void setRFitFunction(omxFitFunction *rff);
};

void copyParamToModelInternal(FreeVarGroup *varGroup, omxState *os, double *at);

typedef std::vector< std::pair<int, MxRList*> > LocalComputeResult;

class omxCompute {
	int computeId;
 protected:
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *glob) {};
	void collectResultsHelper(FitContext *fc, std::vector< omxCompute* > &clist,
				  LocalComputeResult *lcr, MxRList *out);
	static enum ComputeInfoMethod stringToInfoMethod(const char *iMethod);
 public:
	const char *name;
	FreeVarGroup *varGroup;
	omxCompute();
        virtual void initFromFrontend(omxState *, SEXP rObj);
        void compute(FitContext *fc);
	void computeWithVarGroup(FitContext *fc);
        virtual void computeImpl(FitContext *fc) {}
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
        virtual ~omxCompute();
	void reportProgress(FitContext *fc) { Global->reportProgress(name, fc); }
};

omxCompute *omxNewCompute(omxState* os, const char *type);

omxCompute *newComputeGradientDescent();
omxCompute *newComputeNumericDeriv();
omxCompute *newComputeNewtonRaphson();
omxCompute *newComputeConfidenceInterval();
omxCompute *newComputeTryHard();
omxCompute *newComputeNelderMead();

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

#endif
