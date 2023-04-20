/*
 *  Copyright 2017-2021 by the individuals mentioned in the source code history
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

#ifndef u_OMX_COMPUTEGD_H_
#define u_OMX_COMPUTEGD_H_

#include "Compute.h"
#include "finiteDifferences.h"
#include "autoTune.h"

// The GradientOptimizerContext can manage multiple threads
// in parallel. Per-thread specific data should be located
// in the FitContext.
class GradientOptimizerContext {
 private:
	void copyBounds();
	int countNumFree() const { return fc->getNumFree(); }

	// We need to hide this from the optimizer because
	// some parameters might be profiled out and should
	// not be subject to optimization.
	FitContext *fc;
	std::string optName;
	const char *computeName;

 public:
	const int verbose;
	const int numFree;    // how many parameters are not profiled out
	const char *getOptName() const { return optName.c_str(); };
	void setEngineName(const char *engine) {
		optName = computeName;
		optName += "(";
		optName += engine;
		optName += ")";
	}

	bool feasible;
	void *extraData;
	omxMatrix *fitMatrix;
	int numOptimizerThreads;
	int maxMajorIterations;

	int ControlMinorLimit;
	double ControlRho;
	double ControlTolerance;
	bool warmStart;

	Eigen::VectorXd solLB;
	Eigen::VectorXd solUB;

	// NPSOL has bugs and can return the wrong fit & estimates
	// even when optimization proceeds correctly.
	double bestFit;
	Eigen::VectorXd est;    //like fc->est but omitting profiled out params
	Eigen::VectorXd bestEst;
	Eigen::VectorXd grad;
  ConstraintVec AllC;
  ConstraintVec IneqC;
  ConstraintVec EqC;
  void evalAllC(double *constrOut, double *jacOut=0) { AllC.eval(fc, constrOut, jacOut); }
  void evalIneq(double *constrOut, double *jacOut=0) { IneqC.eval(fc, constrOut, jacOut); }
  void evalEq(double *constrOut, double *jacOut=0) { EqC.eval(fc, constrOut, jacOut); }
	double eqNorm, ineqNorm;

	// output
	int informOut;
	Eigen::VectorXd gradOut;
	Eigen::MatrixXd hessOut;  // in-out for warmstart

	GradientOptimizerContext(FitContext *fc, int verbose,
                           omxCompute *owner);
	void reset();

	void setupSimpleBounds();          // NLOPT style
  void setupAllBounds(); //used with NPSOL.
  bool isUnconstrained();

	Eigen::VectorXd constraintFunValsOut;
	Eigen::MatrixXd constraintJacobianOut;
	Eigen::VectorXd LagrMultipliersOut;
	Eigen::VectorXi constraintStatesOut;
	Eigen::MatrixXd LagrHessianOut;

	double solFun(double *myPars, int* mode);
	double recordFit(double *myPars, int* mode);
	template <typename T1> void checkActiveBoxConstraints(Eigen::MatrixBase<T1> &nextEst);
	template <typename T1> void linearConstraintCoefficients(Eigen::MatrixBase<T1> &lcc);
	void useBestFit();
	void copyToOptimizer(double *myPars);
	void copyFromOptimizer(const double *myPars, FitContext *fc2);
	void copyFromOptimizer(const double *myPars) { copyFromOptimizer(myPars, fc); };
	void finish();
	double getFit() const { return fc->getUnscaledFit(); };
	// fc->iterations is a global counter that includes multiple optimizer runs
	//int getIteration() const { return fc->iterations; };
	int iterations;
	int getWanted() const { return fc->wanted; };
	void setWanted(int nw) { fc->wanted = nw; };
	omxState *getState() const { return fc->state; };
	bool doingCI(){
		if(fc->ciobj){return(true);}
		else{return(false);}
	};
};

template <typename T1>
double median(Eigen::MatrixBase<T1> &vec)
{
	if (vec.size() < 2) {
		return vec.array().sum() / vec.size();
	}

	std::vector<int> ind;
	ind.resize(vec.size());
	for (int xx=0; xx < int(vec.size()); ++xx) ind[xx] = xx;
	std::sort(ind.begin(), ind.end(), [&](int ii, int jj) {
			return vec[ii] < vec[jj];
		});
	if (vec.size() % 2 == 0) {
		int mid = vec.size() / 2 - 1;
		return (vec[ind[mid]] + vec[ind[mid+1]]) / 2.0;
	} else {
		int mid = vec.size() / 2;
		return vec[ind[mid]];
	}
}

template <typename T1>
void GradientOptimizerContext::checkActiveBoxConstraints(Eigen::MatrixBase<T1> &nextEst)
{
	if(verbose < 4) return;

	for (int index = 0; index < int(fc->numParam); index++) {
		if(nextEst[index] == solLB[index])
			mxLog("paramter %i hit lower bound %f", index, solLB[index]);
		if(nextEst[index] == solUB[index])
			mxLog("paramter %i hit upper bound %f", index, solUB[index]);
	}
}

typedef void (*GradientOptimizerType)(double *, GradientOptimizerContext &);

#endif
