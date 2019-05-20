/*
 *  Copyright 2017-2019 by the individuals mentioned in the source code history
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

#ifndef _OMX_COMPUTEGD_H_
#define _OMX_COMPUTEGD_H_

#include "Compute.h"
#include "finiteDifferences.h"

class GradientOptimizerContext {
 private:
	void copyBounds();
	int countNumFree();

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

	// Maybe the optimizer should not have information about
	// how the gradient is approximated?
	enum GradientAlgorithm gradientAlgo;
	int gradientIterations;
	double gradientStepSize;

	bool feasible;
	void *extraData;
	omxMatrix *fitMatrix;
	int numOptimizerThreads;
	int maxMajorIterations;

	int ControlMinorLimit;
	double ControlRho;
	double ControlTolerance;
	bool warmStart;
	bool useGradient;
	int ineqType;

	Eigen::VectorXd solLB;
	Eigen::VectorXd solUB;

	// TODO remove, better to pass as a parameter so we can avoid copies
	Eigen::VectorXd& equality;
	Eigen::VectorXd& inequality;
	bool CSOLNP_HACK;

	// NPSOL has bugs and can return the wrong fit & estimates
	// even when optimization proceeds correctly.
	double bestFit;
	Eigen::VectorXd est;    //like fc->est but omitting profiled out params
	Eigen::VectorXd bestEst;
	Eigen::VectorXd grad;
	double eqNorm, ineqNorm;

	// output
	int informOut;
	Eigen::VectorXd gradOut;
	Eigen::MatrixXd hessOut;  // in-out for warmstart

	GradientOptimizerContext(FitContext *fc, int verbose,
				 enum GradientAlgorithm _gradientAlgo,
				 int _gradientIterations,
				 double _gradientStepSize,
				 omxCompute *owner);
	void reset();

	void setupSimpleBounds();          // NLOPT style
	void setupIneqConstraintBounds();  // CSOLNP style
	void setupAllBounds();             // NPSOL style
	
	Eigen::VectorXd constraintFunValsOut;
	Eigen::MatrixXd constraintJacobianOut;
	Eigen::VectorXd LagrMultipliersOut;
	Eigen::VectorXi constraintStatesOut;
	Eigen::MatrixXd LagrHessianOut;

	double solFun(double *myPars, int* mode);
	double recordFit(double *myPars, int* mode);
	void solEqBFun(bool wantAJ);
	void myineqFun(bool wantAJ);
	template <typename T1, typename T2, typename T3> void allConstraintsFun(
			Eigen::MatrixBase<T1> &constraintOut, Eigen::MatrixBase<T2> &jacobianOut, Eigen::MatrixBase<T3> &needcIn, int mode);
	template <typename T1> void checkActiveBoxConstraints(Eigen::MatrixBase<T1> &nextEst);
	template <typename T1> void linearConstraintCoefficients(Eigen::MatrixBase<T1> &lcc);
	bool isUsingAnalyticJacobian(){ return fc->state->usingAnalyticJacobian; }
	Eigen::MatrixXd& analyticEqJacTmp; //<--temporarily holds analytic Jacobian (if present) for an equality constraint
	Eigen::MatrixXd& analyticIneqJacTmp; //<--temporarily holds analytic Jacobian (if present) for an inequality constraint
	void useBestFit();
	void copyToOptimizer(double *myPars);
	void copyFromOptimizer(double *myPars, FitContext *fc2);
	void copyFromOptimizer(double *myPars) { copyFromOptimizer(myPars, fc); };
	void finish();
	double getFit() const { return fc->fit; };
	// fc->iterations is a global counter that includes multiple optimizer runs
	//int getIteration() const { return fc->iterations; };
	int iterations;
	int getWanted() const { return fc->wanted; };
	void setWanted(int nw) { fc->wanted = nw; };
	bool hasKnownGradient() const;
	template <typename T1>
	void setKnownGradient(Eigen::MatrixBase<T1> &gradOut) {
		fc->ciobj->gradient(fc, gradOut.derived().data());
	};
	omxState *getState() const { return fc->state; };
	bool doingCI(){ 
		if(fc->ciobj){return(true);}
		else{return(false);}
	};

	GradientWithRef gwrContext;

	template <typename T1>
	void numericalGradientWithRef(Eigen::MatrixBase<T1> &Epoint);
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
void GradientOptimizerContext::numericalGradientWithRef(Eigen::MatrixBase<T1> &Epoint)
{
	if (getWanted() & FF_COMPUTE_GRADIENT) {
		return;
	} else if (hasKnownGradient()) {
		setKnownGradient(grad);
		return;
	}

	// fc assumed to hold the reference fit
	double refFit = fc->fit;

	gwrContext([&](double *myPars, int thrId)->double{
			FitContext *fc2 = thrId >= 0? fc->childList[thrId] : fc;
			Eigen::Map< Eigen::VectorXd > Est(myPars, fc2->numParam);
			// Only 1 parameter is different so we could
			// update only that parameter instead of all
			// of them.
			copyFromOptimizer(myPars, fc2);
			int want = FF_COMPUTE_FIT;
			ComputeFit(getOptName(), fc2->lookupDuplicate(fitMatrix), want, fc2);
			double fit = fc2->fit;
			if (fc2->outsideFeasibleSet()) {
				fit = nan("infeasible");
			}
			return fit;
		}, refFit, Epoint, grad);

	if (true) {
		Eigen::VectorXd absGrad = grad.array().abs();
		double m1 = std::max(median(absGrad), 1.0);
		double big = 1e4 * m1;
		int adj=0;
		for (int gx=0; gx < grad.size(); ++gx) {
			if (absGrad[gx] < big) continue;
			bool neg = grad[gx] < 0;
			double gg = m1;
			if (neg) gg = -gg;
			grad[gx] = gg;
			++adj;
		}
		if (false && adj) {
			mxLog("%d grad outlier", adj);
			mxPrintMat("absGrad", absGrad);
			mxPrintMat("robust grad", grad);
		}
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
