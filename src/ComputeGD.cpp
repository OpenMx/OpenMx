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

#include "omxState.h"
#include "omxFitFunction.h"
#include "omxNPSOLSpecific.h"
#include "omxExportBackendState.h"
#include "Compute.h"

class omxComputeGD : public omxComputeOperation {
	typedef omxComputeOperation super;
	omxMatrix *fitMatrix;

	SEXP intervals, intervalCodes; // move to FitContext? TODO
	int inform, iter;

public:
	omxComputeGD();
	virtual void initFromFrontend(SEXP rObj);
	virtual void compute(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus() { return inform; }  // backward compatibility
};

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

omxComputeGD::omxComputeGD()
{
	intervals = 0;
	intervalCodes = 0;
	inform = 0;
	iter = 0;
}

void omxComputeGD::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, Global->freeGroup[paramGroup]);
}

void omxComputeGD::compute(FitContext *fc)
{
	if (fc->numParam <= 0) {
		error("Model has no free parameters");
		return;
	}

	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState);

	omxInvokeNPSOL(fitMatrix, fc, &inform, &iter);

	omxFreeChildStates(globalState);

	if (Global->numIntervals) {
		if (!(inform == 0 || inform == 1 || inform == 6)) {
			// TODO: Throw a warning, allow force()
			warning("Not calculating confidence intervals because of NPSOL status %d", inform);
		} else {
			PROTECT(intervals = allocMatrix(REALSXP, Global->numIntervals, 2));
			PROTECT(intervalCodes = allocMatrix(INTSXP, Global->numIntervals, 2));

			omxNPSOLConfidenceIntervals(fitMatrix, fc);
			omxPopulateConfidenceIntervals(intervals, intervalCodes); // TODO move code here
		}
	}  
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	int numFree = fc->numParam;

	SEXP minimum, estimate, gradient, hessian;
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, numFree));
	PROTECT(gradient = allocVector(REALSXP, numFree));
	PROTECT(hessian = allocMatrix(REALSXP, numFree, numFree));

	REAL(minimum)[0] = fc->fit;
	memcpy(REAL(estimate), fc->est, sizeof(double) * numFree);
	memcpy(REAL(gradient), fc->grad, sizeof(double) * numFree);
	memcpy(REAL(hessian), fc->hess, sizeof(double) * numFree * numFree);

	out->push_back(std::make_pair(mkChar("minimum"), minimum));
	out->push_back(std::make_pair(mkChar("estimate"), estimate));
	out->push_back(std::make_pair(mkChar("gradient"), gradient));
	out->push_back(std::make_pair(mkChar("hessianCholesky"), hessian));

	if (intervals && intervalCodes) {
		out->push_back(std::make_pair(mkChar("confidenceIntervals"), intervals));
		out->push_back(std::make_pair(mkChar("confidenceIntervalCodes"), intervalCodes));
	}

	SEXP code, iterations;

	PROTECT(code = NEW_NUMERIC(1));
	REAL(code)[0] = inform;
	out->push_back(std::make_pair(mkChar("npsol.code"), code));

	PROTECT(iterations = NEW_NUMERIC(1));
	REAL(iterations)[0] = iter;
	out->push_back(std::make_pair(mkChar("npsol.iterations"), iterations));
	out->push_back(std::make_pair(mkChar("iterations"), iterations)); // backward compatibility
}
