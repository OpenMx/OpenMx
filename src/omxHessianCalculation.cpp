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
#include "omxHessianCalculation.h"
#include "omxGlobalState.h"
#include "omxExportBackendState.h"

struct hess_struct {
	int     numParams;
	double* freeParams;
	double* Haprox;
	double* Gaprox;
	omxMatrix* fitMatrix;
	double  f0;
	double  functionPrecision;
	int r;
};

void omxPopulateHessianWork(struct hess_struct *hess_work, 
	double functionPrecision, int r, omxState* state,
	omxState *parentState) {

	omxFitFunction* oo = state->fitMatrix->fitFunction;
	int numParams = state->numFreeParams;

	double* optima = parentState->optimalValues;
	double *freeParams = (double*) Calloc(numParams, double);

	hess_work->numParams = numParams;
	hess_work->functionPrecision = functionPrecision;
	hess_work->r = r;
	hess_work->Haprox = (double*) Calloc(r, double);		// Hessian Workspace
	hess_work->Gaprox = (double*) Calloc(r, double);		// Gradient Workspace
	hess_work->freeParams = freeParams;
	for(int i = 0; i < numParams; i++) {
		freeParams[i] = optima[i];
	}

	omxMatrix *fitMatrix = oo->matrix;
	hess_work->fitMatrix = fitMatrix;

	handleFreeVarListHelper(state, freeParams, numParams);

	omxRecompute(fitMatrix);		// Initial recompute in case it matters.	
	hess_work->f0 = omxMatrixElement(fitMatrix, 0, 0);
}

/**
  @params i              parameter number
  @params hess_work      local copy
  @params optima         shared read-only variable
  @params gradient       shared write-only variable
  @params hessian        shared write-only variable
 */
void omxEstimateHessianOnDiagonal(int i, struct hess_struct* hess_work, 
	double *optima, double *gradient, double *hessian) {

	static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
	static const double eps = 1E-4;	// Kept here for access purposes.

	int     numParams          = hess_work->numParams;         // read-only
	double *Haprox             = hess_work->Haprox;
	double *Gaprox             = hess_work->Gaprox;
	double *freeParams         = hess_work->freeParams;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 
	double functionPrecision   = hess_work->functionPrecision; // read-only
	double f0                  = hess_work->f0;                // read-only
	int    r                   = hess_work->r;                 // read-only


	/* Part the first: Gradient and diagonal */
	double iOffset = fabs(functionPrecision * optima[i]);
	if(fabs(iOffset) < eps) iOffset += eps;
	if(OMX_DEBUG) {Rprintf("Hessian estimation: iOffset: %f.\n", iOffset);}
	for(int k = 0; k < r; k++) {			// Decreasing step size, starting at k == 0
		if(OMX_DEBUG) {Rprintf("Hessian estimation: Parameter %d at refinement level %d (%f). One Step Forward.\n", i, k, iOffset);}
		freeParams[i] = optima[i] + iOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		if(OMX_DEBUG) {Rprintf("Hessian estimation: One Step Back.\n");}

		freeParams[i] = optima[i] - iOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Gaprox[k] = (f1 - f2) / (2.0*iOffset); 						// This is for the gradient
		Haprox[k] = (f1 - 2.0 * f0 + f2) / (iOffset * iOffset);		// This is second derivative
		freeParams[i] = optima[i];									// Reset parameter value
		iOffset /= v;
		if(OMX_DEBUG) {Rprintf("Hessian estimation: (%d, %d)--Calculating F1: %f F2: %f, Haprox: %f...\n", i, i, f1, f2, Haprox[k]);}
	}

	for(int m = 1; m < r; m++) {						// Richardson Step
		for(int k = 0; k < (r - m); k++) {
			Gaprox[k] = (Gaprox[k+1] * pow(4.0, m) - Gaprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
		}
	}

	if(OMX_DEBUG) { Rprintf("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d) with value %f...\n", hessian, i, i, i*numParams+i, Haprox[0]); }
	gradient[i] = Gaprox[0];						// NPSOL reports a gradient that's fine.  Why report two?
	hessian[i*numParams + i] = Haprox[0];

	if(OMX_DEBUG) {Rprintf("Done with parameter %d.\n", i);}

}

void omxEstimateHessianOffDiagonal(int i, int l, struct hess_struct* hess_work, 
	double *optima, double *gradient, double *hessian) {

    static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
    static const double eps = 1E-4; // Kept here for access purposes.

	int     numParams          = hess_work->numParams;         // read-only
	double *Haprox             = hess_work->Haprox;
	double *freeParams         = hess_work->freeParams;
	omxMatrix* fitMatrix = hess_work->fitMatrix; 
	double functionPrecision   = hess_work->functionPrecision; // read-only
	double f0                  = hess_work->f0;                // read-only
	int    r                   = hess_work->r;                 // read-only

	double iOffset = fabs(functionPrecision*optima[i]);
	if(fabs(iOffset) < eps) iOffset += eps;
	double lOffset = fabs(functionPrecision*optima[l]);
	if(fabs(lOffset) < eps) lOffset += eps;

	for(int k = 0; k < r; k++) {
		freeParams[i] = optima[i] + iOffset;
		freeParams[l] = optima[l] + lOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f1 = omxMatrixElement(fitMatrix, 0, 0);

		if(OMX_DEBUG) {Rprintf("Hessian estimation: One Step Back.\n");}

		freeParams[i] = optima[i] - iOffset;
		freeParams[l] = optima[l] - lOffset;

		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, freeParams, numParams);

		omxRecompute(fitMatrix);
		double f2 = omxMatrixElement(fitMatrix, 0, 0);

		Haprox[k] = (f1 - 2.0 * f0 + f2 - hessian[i*numParams+i]*iOffset*iOffset -
						hessian[l*numParams+l]*lOffset*lOffset)/(2.0*iOffset*lOffset);
		if(OMX_DEBUG) {Rprintf("Hessian first off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f).\n", Haprox[k], iOffset, lOffset, f1, hessian[i*numParams+i], hessian[l*numParams+l], v, k, pow(v, k), functionPrecision*optima[i], functionPrecision*optima[l]);}

		freeParams[i] = optima[i];				// Reset parameter values
		freeParams[l] = optima[l];

		iOffset = iOffset / v;					//  And shrink step
		lOffset = lOffset / v;
	}

	for(int m = 1; m < r; m++) {						// Richardson Step
		for(int k = 0; k < (r - m); k++) {
			if(OMX_DEBUG) {Rprintf("Hessian off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f, and %f).\n", Haprox[k], iOffset, lOffset, functionPrecision, optima[i], optima[l], v, m, pow(4.0, m), functionPrecision*optima[i], functionPrecision*optima[l], k);}
			Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k]) / (pow(4.0, m)-1);
		}
	}

	if(OMX_DEBUG) {Rprintf("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d and %d) with value %f...", hessian, i, l, i*numParams+l, l*numParams+i, Haprox[0]);}
	hessian[i*numParams+l] = Haprox[0];
	hessian[l*numParams+i] = Haprox[0];

}

void doHessianCalculation(int numParams, int numChildren, 
	struct hess_struct *hess_work, omxState* parentState) {

	int i,j;

	omxFitFunction* parent_oo = parentState->fitMatrix->fitFunction;
	double* parent_gradient = parent_oo->gradient;
	double* parent_hessian = parent_oo->hessian;
	double* parent_optima = parentState->optimalValues;

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
		omxEstimateHessianOnDiagonal(i, hess_work + threadId, 
			parent_optima, parent_gradient, parent_hessian);
	}

	#pragma omp parallel for num_threads(parallelism) 
	for(offset = 0; offset < numOffDiagonal; offset++) {
		int threadId = (numChildren < 2) ? 0 : omx_absolute_thread_num();
		omxEstimateHessianOffDiagonal(diags[offset], offDiags[offset],
			hess_work + threadId, parent_optima, parent_gradient,
			parent_hessian);
	}

	Free(diags);
	Free(offDiags);

}

/*************************************************************************************
 *
 *   omxEstimateHessian
 *
 *  author: tbrick, 2010-02-04
 *
 *  Based on code in the numDeriv library for R <citation needed> // TODO: Find NumDeriv Citation
 *
 *  @params functionPrecision	functional precision for the calculation
 *  @params r					number of repetitions for Richardson approximation
 *  @params currentState		the current omxState
 *
 ************************************************************************************/
void omxEstimateHessian(double functionPrecision, int r)
{
	omxState* parentState = globalState;

	// TODO: Check for nonlinear constraints and adjust algorithm accordingly.
	// TODO: Allow more than one hessian value for calculation

	int numChildren = parentState->numChildren;
	int numParams = parentState->numFreeParams;
	int i;

    struct hess_struct* hess_work;
	if (numChildren < 2) {
		hess_work = Calloc(1, struct hess_struct);
		omxPopulateHessianWork(hess_work, functionPrecision, r, parentState, parentState);
	} else {
		hess_work = Calloc(numChildren, struct hess_struct);
		for(int i = 0; i < numChildren; i++) {
			omxPopulateHessianWork(hess_work + i, functionPrecision, r, parentState->childList[i], parentState);
		}
	}
	if(OMX_DEBUG) Rprintf("Hessian Calculation using %d children\n", numChildren);

	omxFitFunction* parent_oo = parentState->fitMatrix->fitFunction;

	if(parent_oo->hessian == NULL) {
		parent_oo->hessian = (double*) R_alloc(numParams * numParams, sizeof(double));
		if(OMX_DEBUG) {Rprintf("Generated hessian memory, (%d x %d), at 0x%x.\n", numParams, numParams, parent_oo->hessian);}
	}

	if(parent_oo->gradient == NULL) {
		parent_oo->gradient = (double*) R_alloc(numParams, sizeof(double));
		if(OMX_DEBUG) {Rprintf("Generated gradient memory, (%d), at 0x%x.\n", numParams, parent_oo->gradient);}
	}
  
	doHessianCalculation(numParams, numChildren, hess_work, parentState);

	if(OMX_DEBUG) {Rprintf("Hessian Computation complete.\n");}

	if (numChildren < 2) {
		Free(hess_work->Haprox);
		Free(hess_work->Gaprox);
		Free(hess_work->freeParams);
	    Free(hess_work);
	} else {
		for(i = 0; i < numChildren; i++) {
			Free((hess_work + i)->Haprox);
			Free((hess_work + i)->Gaprox);
			Free((hess_work + i)->freeParams);
		}
		Free(hess_work);
	}
}

class omxCompute *newComputeEstimateHessian()
{
	return new omxComputeEstimateHessian;
}

void omxComputeEstimateHessian::compute()
{
	int n = globalState->numFreeParams;

	PROTECT(calculatedHessian = allocMatrix(REALSXP, n, n));
	PROTECT(stdErrors = allocMatrix(REALSXP, n, 1));

	omxEstimateHessian(.0001, 4);
	if(globalState->calculateStdErrors) {
		if(OMX_DEBUG) { Rprintf("Calculating Standard Errors for Fit Function.\n");}
		omxFitFunction* oo = globalState->fitMatrix->fitFunction;
		omxCalculateStdErrorFromHessian(2.0, oo);
	}

	omxPopulateHessians(globalState->numHessians, globalState->fitMatrix, 
			    calculatedHessian, stdErrors, globalState->calculateStdErrors, n);
}

void omxComputeEstimateHessian::reportResults(MxRList *result)
{
	result->push_back(std::make_pair(mkChar("calculatedHessian"), calculatedHessian));

	if (globalState->calculateStdErrors) {
		result->push_back(std::make_pair(mkChar("standardErrors"), stdErrors));
	}
}
