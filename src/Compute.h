/*
 *  Copyright 2013 The OpenMx Project
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

#include <R.h>
#include <Rinternals.h>

#include "types.h"

// The idea of FitContext is to eventually enable fitting from
// multiple starting values in parallel.

class FitContext {
	static omxFitFunction *RFitFunction;

	FitContext *parent;

 public:
	FreeVarGroup *varGroup;
	double fit;
	double *est;
	int *flavor;
	//	double *denom;
	double *grad;
	double *hess;
	double *ihess;

	void init();
	FitContext(std::vector<double> &startingValues);
	FitContext(FitContext *parent, FreeVarGroup *group);
	void fixHessianSymmetry(int want);
	void copyParamToModel(omxState* os, double *at);
	void copyParamToModel(omxState *os);
	void copyParamToModel(omxMatrix *mat, double *at);
	void copyParamToModel(omxMatrix *mat);
	void updateParentAndFree();
	void log(const char *where, int what);
	~FitContext();
	
	static void cacheFreeVarDependencies();
	static void setRFitFunction(omxFitFunction *rff);
};

class omxCompute {
 public:
	FreeVarGroup *varGroup;
	omxCompute();
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(FitContext *fc) = 0;
        virtual void reportResults(FitContext *fc, MxRList *out) = 0;
	virtual double getOptimizerStatus() { return NA_REAL; }  // backward compatibility
        virtual ~omxCompute();
};

class omxCompute *omxNewCompute(omxState* os, const char *type);

class omxCompute *newComputeGradientDescent();
class omxCompute *newComputeEstimatedHessian();
class omxCompute *newComputeNewtonRaphson();

void omxApproxInvertPosDefTriangular(int dim, double *hess, double *ihess, double *stress);
void omxApproxInvertPackedPosDefTriangular(int dim, int *mask, double *packedHess, double *stress);

#endif
