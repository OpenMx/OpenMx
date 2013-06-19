/*
 *  Copyright 2007-2013 The OpenMx Project
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

/**
 * Based on:
 *
 * Paul Gilbert and Ravi Varadhan (2012). numDeriv: Accurate Numerical Derivatives. R package
 * version 2012.9-1. http://CRAN.R-project.org/package=numDeriv
 *
 **/

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxDefines.h"
#include "npsolWrap.h"
#include "omxState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxNPSOLSpecific.h"
#include "omxOptimizer.h"
#include "omxOpenmpWrap.h"
#include "omxExportBackendState.h"
#include "Compute.h"

class omxComputeEstimatedHessian : public omxCompute {
	double stepSize;
	int numIter;
	bool wantSE;

	omxMatrix *fitMat;
	double minimum;
	int numParams;
	double *optima;
	double *gradient;
	double *hessian;

	SEXP calculatedHessian;
	SEXP stdErrors;

	void init();
	void omxPopulateHessianWork(struct hess_struct *hess_work, omxState* state);
	void omxEstimateHessianOnDiagonal(int i, struct hess_struct* hess_work);
	void omxEstimateHessianOffDiagonal(int i, int l, struct hess_struct* hess_work);
	void doHessianCalculation(int numChildren, struct hess_struct *hess_work);

 public:
	omxComputeEstimatedHessian();
        virtual void initFromFrontend(SEXP rObj);
        virtual void compute(double *at);
        virtual void reportResults(MxRList *out);
	virtual double getFit() { return 0; }
	virtual double *getEstimate() { return optima; }
};

struct hess_struct {
	double* freeParams;
	double* Haprox;
	double* Gaprox;
	omxMatrix* fitMatrix;
};

void omxComputeEstimatedHessian::omxPopulateHessianWork(struct hess_struct *hess_work, omxState* state)
{
	double *freeParams = (double*) Calloc(numParams, double);

	hess_work->Haprox = (double*) Calloc(numIter, double);		// Hessian Workspace
	hess_work->Gaprox = (double*) Calloc(numIter, double);		// Gradient Workspace
	hess_work->freeParams = freeParams;
	for(int i = 0; i < numParams; i++) {
		freeParams[i] = optima[i];
	}

	hess_work->fitMatrix = omxLookupDuplicateElement(state, fitMat);

	handleFreeVarListHelper(state, freeParams, numParams);
}

/**
  @params i              parameter number
  @params hess_work      local copy
  @params optima         shared read-only variable
  @params gradient       shared write-only variable
  @params hessian        shared write-only variable
 */
void omxComputeEstimatedHessian::omxEstimateHessianOnDiagonal(int i, struct hess_struct* hess_work)
{
	static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
	static const double eps = 1E-4;	// Kept here for access purposes.

	double *Haprox             = hess_work->Haprox;
	double *Gaprox             = hess_work->Gaprox;
	double *freeParams         = hess_work->freeParams;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 

	/* Part the first: Gradient and diagonal */
	double iOffset = fabs(stepSize * optima[i]);
	if(fabs(iOffset) < eps) iOffset += eps;
	if(OMX_DEBUG) {mxLog("Hessian estimation: iOffset: %f.", iOffset);}
	for(int k = 0; k < numIter; k++) {			// Decreasing step size, starting at k == 0
		if(OMX_DEBUG) {mxLog("Hessian estimation: Parameter %d at refinement level %d (%f). One Step Forward.", i, k, iOffset);}
		freeParams[i] = optima[i] + iOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		if(OMX_DEBUG) {mxLog("Hessian estimation: One Step Back.");}

		freeParams[i] = optima[i] - iOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Gaprox[k] = (f1 - f2) / (2.0*iOffset); 						// This is for the gradient
		Haprox[k] = (f1 - 2.0 * minimum + f2) / (iOffset * iOffset);		// This is second derivative
		freeParams[i] = optima[i];									// Reset parameter value
		iOffset /= v;
		if(OMX_DEBUG) {mxLog("Hessian estimation: (%d, %d)--Calculating F1: %f F2: %f, Haprox: %f...", i, i, f1, f2, Haprox[k]);}
	}

	for(int m = 1; m < numIter; m++) {						// Richardson Step
		for(int k = 0; k < (numIter - m); k++) {
			Gaprox[k] = (Gaprox[k+1] * pow(4.0, m) - Gaprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
		}
	}

	if(OMX_DEBUG) { mxLog("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d) with value %f...", hessian, i, i, i*numParams+i, Haprox[0]); }
	gradient[i] = Gaprox[0];						// NPSOL reports a gradient that's fine.  Why report two?
	hessian[i*numParams + i] = Haprox[0];

	if(OMX_DEBUG) {mxLog("Done with parameter %d.", i);}

}

void omxComputeEstimatedHessian::omxEstimateHessianOffDiagonal(int i, int l, struct hess_struct* hess_work)
{
    static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
    static const double eps = 1E-4; // Kept here for access purposes.

	double *Haprox             = hess_work->Haprox;
	double *freeParams         = hess_work->freeParams;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 

	double iOffset = fabs(stepSize*optima[i]);
	if(fabs(iOffset) < eps) iOffset += eps;
	double lOffset = fabs(stepSize*optima[l]);
	if(fabs(lOffset) < eps) lOffset += eps;

	for(int k = 0; k < numIter; k++) {
		freeParams[i] = optima[i] + iOffset;
		freeParams[l] = optima[l] + lOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		if(OMX_DEBUG) {mxLog("Hessian estimation: One Step Back.");}

		freeParams[i] = optima[i] - iOffset;
		freeParams[l] = optima[l] - lOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Haprox[k] = (f1 - 2.0 * minimum + f2 - hessian[i*numParams+i]*iOffset*iOffset -
						hessian[l*numParams+l]*lOffset*lOffset)/(2.0*iOffset*lOffset);
		if(OMX_DEBUG) {mxLog("Hessian first off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f).", Haprox[k], iOffset, lOffset, f1, hessian[i*numParams+i], hessian[l*numParams+l], v, k, pow(v, k), stepSize*optima[i], stepSize*optima[l]);}

		freeParams[i] = optima[i];				// Reset parameter values
		freeParams[l] = optima[l];

		iOffset = iOffset / v;					//  And shrink step
		lOffset = lOffset / v;
	}

	for(int m = 1; m < numIter; m++) {						// Richardson Step
		for(int k = 0; k < (numIter - m); k++) {
			if(OMX_DEBUG) {mxLog("Hessian off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f, and %f).", Haprox[k], iOffset, lOffset, stepSize, optima[i], optima[l], v, m, pow(4.0, m), stepSize*optima[i], stepSize*optima[l], k);}
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k]) / (pow(4.0, m)-1);
		}
	}

	if(OMX_DEBUG) {mxLog("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d and %d) with value %f...", hessian, i, l, i*numParams+l, l*numParams+i, Haprox[0]);}
	hessian[i*numParams+l] = Haprox[0];
	hessian[l*numParams+i] = Haprox[0];

}

void omxComputeEstimatedHessian::doHessianCalculation(int numChildren, struct hess_struct *hess_work)
{
	int i,j;

	int numOffDiagonal = (numParams * (numParams - 1)) / 2;
	int *diags = Calloc(numOffDiagonal, int);
	int *offDiags = Calloc(numOffDiagonal, int);
	int offset = 0;
	// gcc does not detect the usage of the following variable
	// in the omp parallel pragma, and marks the variable as
	// unused, so the attribute is placed to silence the warning.
    int __attribute__((unused)) parallelism = (numChildren == 0) ? 1 : numChildren;

	// There must be a way to avoid constructing the
	// diags and offDiags arrays and replace them with functions
	// that produce these values given the input
	/// {0, 1, ..., numOffDiagonal - 1} -- M. Spiegel
	for(i = 0; i < numParams; i++) {
		for(j = i - 1; j >= 0; j--) {
			diags[offset] = i;
			offDiags[offset] = j;
			offset++;
		}
	}

	#pragma omp parallel for num_threads(parallelism) 
	for(i = 0; i < numParams; i++) {
		int threadId = (numChildren < 2) ? 0 : omx_absolute_thread_num();
		omxEstimateHessianOnDiagonal(i, hess_work + threadId);
	}

	#pragma omp parallel for num_threads(parallelism) 
	for(offset = 0; offset < numOffDiagonal; offset++) {
		int threadId = (numChildren < 2) ? 0 : omx_absolute_thread_num();
		omxEstimateHessianOffDiagonal(diags[offset], offDiags[offset],
			hess_work + threadId);
	}

	Free(diags);
	Free(offDiags);
}

void omxComputeEstimatedHessian::init()
{
	stepSize = .0001;
	numIter = 4;
	stdErrors = NULL;
	optima = NULL;
}

omxComputeEstimatedHessian::omxComputeEstimatedHessian()
{
	init();
}

void omxComputeEstimatedHessian::initFromFrontend(SEXP rObj)
{
	numParams = globalState->numFreeParams;
	if (numParams <= 0) error("Model has no free parameters");

	fitMat = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("se")));
	wantSE = asLogical(slotValue);
	UNPROTECT(1);
}

void omxComputeEstimatedHessian::compute(double *at)
{
	numParams = globalState->numFreeParams;

	optima = at;

	PROTECT(calculatedHessian = allocMatrix(REALSXP, numParams, numParams));

	// TODO: Check for nonlinear constraints and adjust algorithm accordingly.
	// TODO: Allow more than one hessian value for calculation

	int numChildren = globalState->numChildren;

	omxRecompute(fitMat);
	minimum = omxMatrixElement(fitMat, 0, 0);

	struct hess_struct* hess_work;
	if (numChildren < 2) {
		hess_work = Calloc(1, struct hess_struct);
		omxPopulateHessianWork(hess_work, globalState);
	} else {
		hess_work = Calloc(numChildren, struct hess_struct);
		for(int i = 0; i < numChildren; i++) {
			omxPopulateHessianWork(hess_work + i, globalState->childList[i]);
		}
	}
	if(OMX_DEBUG) mxLog("Hessian Calculation using %d children", numChildren);

	hessian = REAL(calculatedHessian);

	gradient = (double*) R_alloc(numParams, sizeof(double));
  
	doHessianCalculation(numChildren, hess_work);

	if(OMX_DEBUG) {mxLog("Hessian Computation complete.");}

	if (numChildren < 2) {
		Free(hess_work->Haprox);
		Free(hess_work->Gaprox);
		Free(hess_work->freeParams);
	    Free(hess_work);
	} else {
		for(int i = 0; i < numChildren; i++) {
			Free((hess_work + i)->Haprox);
			Free((hess_work + i)->Gaprox);
			Free((hess_work + i)->freeParams);
		}
		Free(hess_work);
	}

	if (wantSE) {
		// This function calculates the standard errors from the hessian matrix
		// sqrt(diag(solve(hessian)))

		const double scale = 2;
		double* workspace = (double *) Calloc(numParams * numParams, double);
	
		for(int i = 0; i < numParams; i++)
			for(int j = 0; j <= i; j++)
				workspace[i*numParams+j] = hessian[i*numParams+j];		// Populate upper triangle
	
		char u = 'U';
		std::vector<int> ipiv(numParams);
		int lwork = -1;
		double temp;
		int info = 0;
	
		F77_CALL(dsytrf)(&u, &numParams, workspace, &numParams, ipiv.data(), &temp, &lwork, &info);
	
		lwork = (temp > numParams?temp:numParams);
	
		double* work = (double*) Calloc(lwork, double);
	
		F77_CALL(dsytrf)(&u, &numParams, workspace, &numParams, ipiv.data(), work, &lwork, &info);
	
		if(info != 0) {
			// report error TODO
		} else {
			F77_CALL(dsytri)(&u, &numParams, workspace, &numParams, ipiv.data(), work, &info);
	
			if(info != 0) {
				// report error TODO
			} else {
				PROTECT(stdErrors = allocMatrix(REALSXP, numParams, 1));
				double* stdErr = REAL(stdErrors);
				for(int i = 0; i < numParams; i++) {
					stdErr[i] = sqrt(scale) * sqrt(workspace[i * numParams + i]);
				}
			}
		}
	
		Free(workspace);
		Free(work);
	}
}

void omxComputeEstimatedHessian::reportResults(MxRList *result)
{
	result->push_back(std::make_pair(mkChar("calculatedHessian"), calculatedHessian));

	if (stdErrors) {
		result->push_back(std::make_pair(mkChar("standardErrors"), stdErrors));
	}
}

omxCompute *newComputeEstimatedHessian()
{
	if (globalState->numConstraints != 0) {
		error("Cannot compute estimated Hessian with constraints (%d constraints found)",
		      globalState->numConstraints);
	}
	return new omxComputeEstimatedHessian;
}

