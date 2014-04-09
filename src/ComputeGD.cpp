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

#include "omxDefines.h"
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
    
	SEXP hessChol;
	SEXP intervals, intervalCodes; // move to FitContext? TODO
    
public:
	omxComputeGD();
	virtual void initFromFrontend(SEXP rObj);
	virtual omxFitFunction *getFitFunction();
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

omxComputeGD::omxComputeGD()
{
	intervals = 0;
	intervalCodes = 0;
	hessChol = NULL;
}

void omxComputeGD::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);
    
	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("useGradient")));
	if (Rf_length(slotValue)) {
		useGradient = Rf_asLogical(slotValue);
	} else {
		useGradient = Global->analyticGradients;
	}
    
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
    
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("engine")));
	const char *engine_name = CHAR(Rf_asChar(slotValue));
	if (strcmp(engine_name, "CSOLNP")==0) {
		engine = OptEngine_CSOLNP;
	} else if (strcmp(engine_name, "NPSOL")==0) {
#if HAS_NPSOL
		engine = OptEngine_NPSOL;
#else
		Rf_error("NPSOL is not available in this build");
#endif
	} else {
		Rf_error("MxComputeGradientDescent engine %s unknown", engine_name);
	}
}

omxFitFunction *omxComputeGD::getFitFunction()
{ return fitMatrix->fitFunction; }

void omxComputeGD::computeImpl(FitContext *fc)
{
    size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		Rf_error("Model has no free parameters");
		return;
	}
    
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState);
    
	switch (engine) {
        case OptEngine_NPSOL:
#if HAS_NPSOL
		Rf_protect(hessChol = Rf_allocMatrix(REALSXP, numParam, numParam));
		omxInvokeNPSOL(fitMatrix, fc, &fc->inform, &fc->iterations, useGradient, varGroup, verbose,
			       REAL(hessChol));
#endif
            break;
        case OptEngine_CSOLNP:
            omxInvokeCSOLNP(fitMatrix, fc, &fc->inform, &fc->iterations, varGroup, verbose, fc->getDenseHessUninitialized());
	    fc->wanted |= FF_COMPUTE_HESSIAN;
            break;
        default: Rf_error("huh?");
	}
	fc->wanted |= FF_COMPUTE_GRADIENT;
    
	if (!std::isfinite(fc->fit) || fc->fit == 1e24) {  // remove magic number 1e24 TODO
		std::string diag = fc->getIterationError();
		omxRaiseErrorf(0, "MxComputeGradientDescent: fitfunction %s evaluated to %f (%s)",
			       fitMatrix->name, fc->fit, diag.c_str());
		return;
	}

	double fitCopy = fc->fit;
	omxFreeChildStates(globalState);
    
	if (Global->numIntervals) {
		// I'm not sure why INFORM_NOT_AT_OPTIMUM is okay, but that's how it was.
		if (fc->inform >= INFORM_LINEAR_CONSTRAINTS_INFEASIBLE && fc->inform != INFORM_NOT_AT_OPTIMUM) {
			// TODO: allow forcing
			Rf_warning("Not calculating confidence intervals because of optimizer status %d", fc->inform);
		} else {
			Rf_protect(intervals = Rf_allocMatrix(REALSXP, Global->numIntervals, 2));
			Rf_protect(intervalCodes = Rf_allocMatrix(INTSXP, Global->numIntervals, 2));
			if (engine == OptEngine_NPSOL) {
#if HAS_NPSOL
				omxNPSOLConfidenceIntervals(fitMatrix, fc);
#endif
			}
			else if (engine == OptEngine_CSOLNP) {
				omxCSOLNPConfidenceIntervals(fitMatrix, fc, verbose);
			}
			omxPopulateConfidenceIntervals(intervals, intervalCodes);
			fc->copyParamToModel(globalState);
			fc->fit = fitCopy;
		}
	}

	fc->wanted |= FF_COMPUTE_BESTFIT;
    /*printf("fc->hess in computeGD\n");
    printf("%2f", fc->hess[0]); putchar('\n');
    printf("%2f", fc->hess[1]); putchar('\n');
    printf("%2f", fc->hess[2]); putchar('\n');
    */
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);
    
    /*printf("fc->hess in computeGD:report results\n");
    printf("%2f", fc->hess[0]); putchar('\n');
    printf("%2f", fc->hess[1]); putchar('\n');
    printf("%2f", fc->hess[2]); putchar('\n');
*/
	size_t numFree = varGroup->vars.size();
    
	if (engine == OptEngine_NPSOL) {
		out->add("hessianCholesky", hessChol);
	}
    
	if (intervals && intervalCodes) {
		out->add("confidenceIntervals", intervals);
		out->add("confidenceIntervalCodes", intervalCodes);
	}
}
