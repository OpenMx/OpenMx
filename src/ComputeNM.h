/*
 *  Copyright 2017 The OpenMx Project
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

#include "Compute.h"
#include <Eigen/Core>
#include <Eigen/Dense>



class omxComputeNM : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMatrix;
	int verbose;
	bool nudge;
	int maxIter;
	bool defaultMaxIter;
	std::vector<int> excludeVars;
	double alpha, betao, betai, gamma, sigma;
	double bignum;
	int iniSimplexType; //regular=1, right=2, smartRight=3, random=4
	double iniSimplexEdge;
	Eigen::MatrixXd iniSimplexMtx;
	bool greedyMinimize, altContraction;
	double degenLimit;
	Eigen::Vector2i stagnationCtrl;
	bool validationRestart;
	double xTolProx, fTolProx, xTolRelChange, fTolRelChange;
	bool doPseudoHessian;
	
public:
	omxComputeNM();
	virtual void initFromFrontend(omxState *, SEXP rObj);
	//virtual void computeImpl(FitContext *fc);
	//virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};
