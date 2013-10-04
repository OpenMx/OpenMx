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

struct matrixVectorProdTerm {
	int hentry;
	int gentry;
	int dest;
	matrixVectorProdTerm() {}
	matrixVectorProdTerm(int he, int ge, int de) {
		hentry=he;
		gentry=ge;
		dest=de;
	}
	bool operator< (const matrixVectorProdTerm &j) const {
	  if (hentry < j.hentry) return true;
	  if (hentry == j.hentry) {
	    if (gentry < j.gentry) return true;
	    if (gentry == j.gentry && dest < j.dest) return true;
	  }
	  return false;
	};
	bool operator==(const matrixVectorProdTerm &i) const {
		return i.hentry == hentry && i.gentry == gentry && i.dest == dest;
	}
};

// The idea of FitContext is to eventually enable fitting from
// multiple starting values in parallel.

class FitContext {
	static omxFitFunction *RFitFunction;

	FitContext *parent;

 public:
	FreeVarGroup *varGroup;
	double mac;
	double fit;
	double *est;
	int *flavor;
	//	double *denom;
	bool forwardDeriv;
	double *grad;
	double *hess;
	double *ihess;
	std::vector< matrixVectorProdTerm > hgProd;
	double caution;

	int wanted;
	bool changedEstimates; // only used for FF_COMPUTE_PREOPTIMIZE

	void init();
	FitContext(std::vector<double> &startingValues);
	FitContext(FitContext *parent, FreeVarGroup *group);
	void fixHessianSymmetry(int want);
	void copyParamToModel(omxState* os, double *at);
	void copyParamToModel(omxState *os);
	void copyParamToModel(omxMatrix *mat, double *at);
	void copyParamToModel(omxMatrix *mat);
	void maybeCopyParamToModel(omxState* os);
	void updateParentAndFree();
	void log(const char *where);
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
