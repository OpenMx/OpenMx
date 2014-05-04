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

class ComputeGDBase : public omxCompute {
protected:
	typedef omxCompute super;
	enum OptEngine engine;
	omxMatrix *fitMatrix;
	int verbose;
	double optimalityTolerance;

	virtual void initFromFrontend(SEXP rObj);
};

class omxComputeGD : public ComputeGDBase {
	typedef ComputeGDBase super;
	bool useGradient;
	SEXP hessChol;
    
public:
	omxComputeGD();
	virtual void initFromFrontend(SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

class ComputeCI : public ComputeGDBase {
	typedef ComputeGDBase super;
	SEXP intervals, intervalCodes;

public:
	ComputeCI();
	virtual void initFromFrontend(SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

omxCompute *newComputeConfidenceInterval()
{
	return new ComputeCI();
}

omxComputeGD::omxComputeGD()
{
	hessChol = NULL;
}

void ComputeGDBase::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	SEXP slotValue;
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
    
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("tolerance")));
	optimalityTolerance = Rf_asReal(slotValue);
    
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
		Rf_error("%s: engine %s unknown", name, engine_name);
	}
}

void omxComputeGD::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);
    
	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("useGradient")));
	if (Rf_length(slotValue)) {
		useGradient = Rf_asLogical(slotValue);
	} else {
		useGradient = Global->analyticGradients;
	}
}

void omxComputeGD::computeImpl(FitContext *fc)
{
    size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		omxRaiseErrorf("%s: model has no free parameters", name);
		return;
	}
    
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

	if (fitMatrix->fitFunction && fitMatrix->fitFunction->usesChildModels)
		omxFitFunctionCreateChildren(globalState);
    
	switch (engine) {
        case OptEngine_NPSOL:{
#if HAS_NPSOL
		Rf_protect(hessChol = Rf_allocMatrix(REALSXP, numParam, numParam));
		omxInvokeNPSOL(fitMatrix, fc, &fc->inform, useGradient, varGroup, verbose,
			       REAL(hessChol), optimalityTolerance);
		Eigen::Map<Eigen::MatrixXd> hc(REAL(hessChol), numParam, numParam);
		Eigen::MatrixXd hcT = hc.transpose();
		Eigen::Map<Eigen::MatrixXd> dest(fc->getDenseHessUninitialized(), numParam, numParam);
		dest.noalias() = hcT * hc;
#endif
		break;}
        case OptEngine_CSOLNP:
            omxInvokeCSOLNP(fitMatrix, fc, &fc->inform, varGroup, verbose,
			    fc->getDenseHessUninitialized(), optimalityTolerance);
            break;
        default: Rf_error("huh?");
	}
	fc->wanted |= FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN;
    
	if (!std::isfinite(fc->fit) || fc->fit == 1e24) {  // remove magic number 1e24 TODO
		std::string diag = fc->getIterationError();
		omxRaiseErrorf("MxComputeGradientDescent: fitfunction %s evaluated to %f (%s)",
			       fitMatrix->name, fc->fit, diag.c_str());
		return;
	}

	omxFreeChildStates(globalState);
    
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
	if (engine == OptEngine_NPSOL) {
		out->add("hessianCholesky", hessChol);
	}
}

ComputeCI::ComputeCI()
{
	intervals = 0;
	intervalCodes = 0;
}

void ComputeCI::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);
}

void ComputeCI::computeImpl(FitContext *fc)
{
	int numInts = Global->numIntervals;
	if (verbose >= 1) mxLog("%s: starting work on %d intervals", name, numInts);
	if (!numInts) return;

	// I'm not sure why INFORM_NOT_AT_OPTIMUM is okay, but that's how it was.
	if (fc->inform >= INFORM_LINEAR_CONSTRAINTS_INFEASIBLE && fc->inform != INFORM_NOT_AT_OPTIMUM) {
		// TODO: allow forcing
		Rf_warning("Not calculating confidence intervals because of optimizer status %d", fc->inform);
		return;
	}

	Rf_protect(intervals = Rf_allocMatrix(REALSXP, numInts, 2));
	Rf_protect(intervalCodes = Rf_allocMatrix(INTSXP, numInts, 2));

	switch (engine) {
	case OptEngine_NPSOL:
#if HAS_NPSOL
		omxNPSOLConfidenceIntervals(fitMatrix, fc, optimalityTolerance);
#endif
		break;
	case OptEngine_CSOLNP:
		omxCSOLNPConfidenceIntervals(fitMatrix, fc, verbose, optimalityTolerance);
		break;
	default:
		Rf_error("huh?");
	}

	if(OMX_DEBUG) { mxLog("Populating CIs for %d fit functions.", numInts); }

	double* interval = REAL(intervals);
	int* intervalCode = INTEGER(intervalCodes);
	for(int j = 0; j < numInts; j++) {
		omxConfidenceInterval *oCI = Global->intervalList + j;
		interval[j] = oCI->min;
		interval[j + numInts] = oCI->max;
		intervalCode[j] = oCI->lCode;
		intervalCode[j + numInts] = oCI->uCode;
	}

	fc->copyParamToModel(globalState);
}

void ComputeCI::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	if (intervals && intervalCodes) {
		out->add("confidenceIntervals", intervals);
		out->add("confidenceIntervalCodes", intervalCodes);
	}
}
