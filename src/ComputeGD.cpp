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
#include "omxOptimizer.h"
#include "omxNPSOLSpecific.h"
#include "omxExportBackendState.h"
#include "Compute.h"

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

void omxComputeGD::init()
{
	intervals = 0;
	intervalCodes = 0;
	inform = 0;
	iter = 0;
}

void omxComputeGD::initFromFrontend(SEXP rObj)
{
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");

	numFree = Global.numFreeParams;
	if (numFree <= 0) {
		error("Model has no free parameters");
		return;
	}

	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, numFree));
	PROTECT(gradient = allocVector(REALSXP, numFree));
	PROTECT(hessian = allocMatrix(REALSXP, numFree, numFree));
}

void omxComputeGD::compute(double *startVals)
{
	memcpy(REAL(estimate), startVals, sizeof(double)*numFree);

	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState);

	omxInvokeNPSOL(fitMatrix, REAL(minimum), REAL(estimate),
		       REAL(gradient), REAL(hessian), &inform, &iter);

	omxFreeChildStates(globalState);

	if (Global.numIntervals) {
		if (!(inform == 0 || inform == 1 || inform == 6)) {
			// TODO: Throw a warning, allow force()
			warning("Not calculating confidence intervals because of NPSOL status %d", inform);
		} else {
			PROTECT(intervals = allocMatrix(REALSXP, Global.numIntervals, 2));
			PROTECT(intervalCodes = allocMatrix(INTSXP, Global.numIntervals, 2));

			omxNPSOLConfidenceIntervals(fitMatrix, getFit(),
						    getEstimate(), Global.ciMaxIterations);
			omxPopulateConfidenceIntervals(intervals, intervalCodes);
		}
	}  
}

void omxComputeGD::reportResults(MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

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
