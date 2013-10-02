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
#include "omxCsolnp.h"
#include "Compute.h"
#include "npsolswitch.h"

enum OptEngine {
	OptEngine_NPSOL,
	OptEngine_CSOLNP
};

class omxComputeGD : public omxCompute {
	typedef omxCompute super;
	enum OptEngine engine;
	omxMatrix *fitMatrix;
	bool useGradient;
	int verbose;
    
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
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);
    
	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("useGradient")));
	if (length(slotValue)) {
		useGradient = asLogical(slotValue);
	} else {
		useGradient = Global->analyticGradients;
	}
    
	PROTECT(slotValue = GET_SLOT(rObj, install("verbose")));
	verbose = asInteger(slotValue);
    
	PROTECT(slotValue = GET_SLOT(rObj, install("engine")));
	const char *engine_name = CHAR(asChar(slotValue));
	if (strcmp(engine_name, "CSOLNP")==0) {
		engine = OptEngine_CSOLNP;
	} else if (strcmp(engine_name, "NPSOL")==0) {
#if HAS_NPSOL
		engine = OptEngine_NPSOL;
#else
		error("NPSOL is not available in this build");
#endif
	} else {
		error("MxComputeGradientDescent engine %s unknown", engine_name);
	}
}

void omxComputeGD::compute(FitContext *fc)
{
	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		error("Model has no free parameters");
		return;
	}
    
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
	fc->maybeCopyParamToModel(globalState);
    
	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState);
    
	switch (engine) {
        case OptEngine_NPSOL:
#if HAS_NPSOL
            omxInvokeNPSOL(fitMatrix, fc, &inform, &iter, useGradient, varGroup, verbose);
#endif
            break;
        case OptEngine_CSOLNP:
            omxInvokeCSOLNP(fitMatrix, fc, &inform, &iter, useGradient, varGroup, verbose);
            break;
        default: error("huh?");
	}
    
	omxFreeChildStates(globalState);
    
	if (Global->numIntervals && engine == OptEngine_NPSOL) {
		if (!(inform == 0 || inform == 1 || inform == 6)) {
			// TODO: allow forcing
			warning("Not calculating confidence intervals because of NPSOL status %d", inform);
		} else {
			PROTECT(intervals = allocMatrix(REALSXP, Global->numIntervals, 2));
			PROTECT(intervalCodes = allocMatrix(INTSXP, Global->numIntervals, 2));
            
#if HAS_NPSOL
			omxNPSOLConfidenceIntervals(fitMatrix, fc);
#endif
			omxPopulateConfidenceIntervals(intervals, intervalCodes); // TODO move code here
		}
	}
    
	omxMarkDirty(fitMatrix); // not sure why it needs to be dirty
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);
    
	size_t numFree = varGroup->vars.size();
    
	SEXP estimate, gradient, hessian;
	PROTECT(estimate = allocVector(REALSXP, numFree));
	PROTECT(gradient = allocVector(REALSXP, numFree));
	PROTECT(hessian = allocMatrix(REALSXP, numFree, numFree));
    
	memcpy(REAL(estimate), fc->est, sizeof(double) * numFree);
	memcpy(REAL(gradient), fc->grad, sizeof(double) * numFree);
	memcpy(REAL(hessian), fc->hess, sizeof(double) * numFree * numFree);
    
	out->push_back(std::make_pair(mkChar("minimum"), ScalarReal(fc->fit)));
	out->push_back(std::make_pair(mkChar("Minus2LogLikelihood"), ScalarReal(fc->fit)));
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
