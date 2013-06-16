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

#include "omxGlobalState.h"
#include "omxFitFunction.h"
#include "omxOptimizer.h"
#include "omxNPSOLSpecific.h"
#include "omxHessianCalculation.h"
#include "omxExportBackendState.h"

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

void omxComputeGD::initFromFrontend(SEXP rObj)
{
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
}

void omxComputeGD::setStartValues(SEXP startVals)
{
	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState, globalState->numThreads);

	int numFree = globalState->numFreeParams;

	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, numFree));
	PROTECT(gradient = allocVector(REALSXP, numFree));
	PROTECT(hessian = allocMatrix(REALSXP, numFree, numFree));

	if (numFree>0) { memcpy(REAL(estimate), REAL(startVals), sizeof(double)*numFree); }
}

void omxComputeGD::compute(bool disableOpt)  // remove flag TODO
{
	omxInvokeNPSOL(fitMatrix, REAL(minimum), REAL(estimate),
		       REAL(gradient), REAL(hessian), disableOpt);
}

void omxComputeGD::reportResults(MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	out->push_back(std::make_pair(mkChar("minimum"), minimum));
	out->push_back(std::make_pair(mkChar("estimate"), estimate));
	out->push_back(std::make_pair(mkChar("gradient"), gradient));
	out->push_back(std::make_pair(mkChar("hessianCholesky"), hessian));
}

