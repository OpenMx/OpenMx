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
#include "nloptcpp.h"
#include "Compute.h"
#include "npsolswitch.h"
#include "glue.h"

enum OptEngine {
	OptEngine_NPSOL,
	OptEngine_CSOLNP,
    OptEngine_NLOPT
};

class ComputeGDBase : public omxCompute {
protected:
	typedef omxCompute super;
	enum OptEngine engine;
	omxMatrix *fitMatrix;
	int verbose;
	double optimalityTolerance;

	virtual void initFromFrontend(omxState *, SEXP rObj);
};

class omxComputeGD : public ComputeGDBase {
	typedef ComputeGDBase super;
	bool useGradient;
	SEXP hessChol;
	bool nudge;

	int warmStartSize;
	double *warmStart;
    
public:
	omxComputeGD();
	virtual void initFromFrontend(omxState *, SEXP rObj);
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
	virtual void initFromFrontend(omxState *, SEXP rObj);
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
	warmStart = NULL;
}

void ComputeGDBase::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
    
	ScopedProtect p2(slotValue, R_do_slot(rObj, Rf_install("tolerance")));
	optimalityTolerance = Rf_asReal(slotValue);
    
	ScopedProtect p3(slotValue, R_do_slot(rObj, Rf_install("engine")));
	const char *engine_name = CHAR(Rf_asChar(slotValue));
	if (strEQ(engine_name, "CSOLNP")) {
		engine = OptEngine_CSOLNP;
	} else if (strEQ(engine_name, "NLOPT")) {
#ifdef HAS_NLOPT
		engine = OptEngine_NLOPT;
#else
		Rf_error("NLOPT is not available in this build");
#endif
	} else if (strEQ(engine_name, "NPSOL")) {
#if HAS_NPSOL
		engine = OptEngine_NPSOL;
#else
		Rf_error("NPSOL is not available in this build");
#endif
	} else {
		Rf_error("%s: engine %s unknown", name, engine_name);
	}
}

void omxComputeGD::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);
    
	SEXP slotValue;
	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("useGradient")));
	if (Rf_length(slotValue)) {
		useGradient = Rf_asLogical(slotValue);
	} else {
		useGradient = Global->analyticGradients;
	}

	ScopedProtect p4(slotValue, R_do_slot(rObj, Rf_install("nudgeZeroStarts")));
	nudge = Rf_asLogical(slotValue);
    
	ScopedProtect p2(slotValue, R_do_slot(rObj, Rf_install("warmStart")));
	if (!Rf_isNull(slotValue)) {
		SEXP matrixDims;
		Rf_protect(matrixDims = Rf_getAttrib(slotValue, R_DimSymbol));
		int *dimList = INTEGER(matrixDims);
		int rows = dimList[0];
		int cols = dimList[1];
		if (rows != cols) Rf_error("%s: warmStart matrix must be square", name);

		warmStartSize = rows;
		warmStart = REAL(slotValue);
	}
}

void omxComputeGD::computeImpl(FitContext *fc)
{
    size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		omxRaiseErrorf("%s: model has no free parameters", name);
		return;
	}
    
	for (int px = 0; px < int(numParam); ++px) {
		omxFreeVar *fv = varGroup->vars[px];
		if (nudge && fc->est[px] == 0.0) {
			fc->est[px] += 0.1;
		}
		if (fv->lbound > fc->est[px]) {
			fc->est[px] = fv->lbound + 1.0e-6;
		}
		if (fv->ubound < fc->est[px]) {
			fc->est[px] = fv->ubound - 1.0e-6;
		}
        }
    
	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

	fc->createChildren();
    
	int beforeEval = Global->computeCount;

	switch (engine) {
        case OptEngine_NPSOL:{
#if HAS_NPSOL
		if (!hessChol) {
			Rf_protect(hessChol = Rf_allocMatrix(REALSXP, numParam, numParam));
		}
		bool doWarm = false;
		if (warmStart) {
			if (warmStartSize != int(numParam)) {
				Rf_warning("%s: warmStart size %d does not match number of free parameters %d (ignored)",
					   warmStartSize, numParam);
			} else {
				memcpy(REAL(hessChol), warmStart, sizeof(double) * numParam * numParam);
				doWarm = true;
			}
		}
		omxInvokeNPSOL(fitMatrix, fc, &fc->inform, useGradient, varGroup, verbose,
			       REAL(hessChol), optimalityTolerance, doWarm);
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
#ifdef HAS_NLOPT
        case OptEngine_NLOPT:
            omxInvokeNLOPTorSANN(fitMatrix, fc, &fc->inform, varGroup, verbose,
                fc->getDenseHessUninitialized(), optimalityTolerance);
        break;
#endif
        default: Rf_error("Optimizer %d is not available", engine);
	}
	fc->wanted |= FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN;
    
	if (fc->inform <= 0 && Global->computeCount - beforeEval == 1) {
		fc->inform = INFORM_STARTING_VALUES_INFEASIBLE;
	}

	// Optimizers can terminate with inconsistent fit and parameters
	fc->copyParamToModel();
	ComputeFit("ComputeGD", fitMatrix, FF_COMPUTE_FIT, fc);

	if (verbose >= 1) {
		mxLog("%s: final fit is %2f", name, fc->fit);
		fc->log(FF_COMPUTE_ESTIMATE);
	}

	if (fitMatrix->rows == 1) {
		if (!std::isfinite(fc->fit) || fc->fit == 1e24) {  // remove magic number 1e24 TODO
			std::string diag = fc->getIterationError();
			omxRaiseErrorf("MxComputeGradientDescent: fitfunction %s is not finite (%s)",
				       fitMatrix->name, diag.c_str());
			return;
		}
	}

	fc->wanted |= FF_COMPUTE_BESTFIT;
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);
    
	if (engine == OptEngine_NPSOL) {
		out->add("hessianCholesky", hessChol);
	}
}

ComputeCI::ComputeCI()
{
	intervals = 0;
	intervalCodes = 0;
}

void ComputeCI::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);
}

extern "C" { void F77_SUB(npoptn)(char* string, int Rf_length); };

void ComputeCI::computeImpl(FitContext *mle)
{
	Global->unpackConfidenceIntervals();

	int numInts = (int) Global->intervalList.size();
	if (verbose >= 1) mxLog("%s: starting work on %d intervals", name, numInts);
	if (!numInts) return;

	// I'm not sure why INFORM_NOT_AT_OPTIMUM is okay, but that's how it was.
	if (mle->inform >= INFORM_LINEAR_CONSTRAINTS_INFEASIBLE && mle->inform != INFORM_NOT_AT_OPTIMUM) {
		// TODO: allow forcing
		Rf_warning("Not calculating confidence intervals because of optimizer status %d", mle->inform);
		return;
	}

	Rf_protect(intervals = Rf_allocMatrix(REALSXP, numInts, 3));
	Rf_protect(intervalCodes = Rf_allocMatrix(INTSXP, numInts, 2));

	// Could be smarter about setting upper & lower once instead of every attempt TODO
	ConfidenceIntervalFit cif(0);
	cif.fitMatrix = fitMatrix;
	cif.ControlTolerance = std::isfinite(optimalityTolerance)? optimalityTolerance : 1.0e-16;

	GradientOptimizerType go = NULL;
	switch (engine) {
#if HAS_NPSOL
	case OptEngine_NPSOL:{
		std::string option = string_snprintf("Cold start");
		F77_CALL(npoptn)((char*) option.c_str(), option.size());
		//omxNPSOLConfidenceIntervals(fitMatrix, mle, optimalityTolerance);
		go = omxNPSOL;
		break;}
#endif
	case OptEngine_CSOLNP:
		//omxCSOLNPConfidenceIntervals(fitMatrix, mle, verbose, optimalityTolerance);
		go = omxCSOLNP;
		break;
#ifdef HAS_NLOPT
    case OptEngine_NLOPT:
	    //omxNLOPTorSANNConfidenceIntervals(fitMatrix, mle, optimalityTolerance);
        break;
#endif
	default:
		Rf_error("huh?");
	}

	{
		const int ciMaxIterations = Global->ciMaxIterations;
		FitContext fc(mle, mle->varGroup);
		fc.createChildren();
		cif.fc = &fc;
		FreeVarGroup *freeVarGroup = fc.varGroup;
    
		const int n = int(freeVarGroup->vars.size());
    
		if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }
    
		const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
    
		for(int i = 0; i < (int) Global->intervalList.size(); i++) {
			omxConfidenceInterval *currentCI = Global->intervalList[i];
        
			const char *matName = "anonymous matrix";
			if (currentCI->matrix->name) {
				matName = currentCI->matrix->name;
			}
        
			currentCI->lbound += mle->fit;          // Convert from offsets to targets
			currentCI->ubound += mle->fit;          // Convert from offsets to targets
        
			for (int lower=0; lower <= 1; ++lower) {
				if (lower  && !std::isfinite(currentCI->lbound)) continue;
				if (!lower && !std::isfinite(currentCI->ubound)) continue;

				memcpy(fc.est, mle->est, n * sizeof(double)); // Reset to previous optimum
        
				int tries = 0;
				int inform = -1;
				double bestFit = std::numeric_limits<double>::max();
            
				while (inform!= 0 && ++tries <= ciMaxIterations) {
					Global->checkpointMessage(mle, mle->est, "%s[%d, %d] %s CI (try %d)",
								  matName, currentCI->row + 1, currentCI->col + 1,
								  lower? "lower" : "upper", tries);

					cif.bestFit = std::numeric_limits<double>::max();
					cif.currentInterval = i;
					cif.calcLower = lower;
					go(fc.est, cif);

					fc.copyParamToModel();
					const double fitOut = fc.fit;

					if (fitOut < bestFit) {
						omxRecompute(currentCI->matrix, &fc);
						double val = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
						if (lower) currentCI->min = val;
						else       currentCI->max = val;
						bestFit = fitOut;
					}

					inform = cif.informOut;
					if (lower) currentCI->lCode = inform;
					else       currentCI->uCode = inform;
					if(verbose>=1) { mxLog("inform(%d,%d) is: %d", i, lower, inform);}
					if(inform == 0) break;

					bool jitter = TRUE;
					for(int j = 0; j < n; j++) {
						if(fabs(fc.est[j] - mle->est[j]) > objDiff) {
							jitter = FALSE;
							break;
						}
					}
					if(jitter) {
						for(int j = 0; j < n; j++) {
							double sign = 2 * (tries % 2) - 1;
							fc.est[j] = mle->est[j] + sign * objDiff * tries;
						}
					}
				}
			}
		}
	}

	mle->copyParamToModel();

	Eigen::Map< Eigen::ArrayXXd > interval(REAL(intervals), numInts, 3);
	interval.fill(NA_REAL);
	int* intervalCode = INTEGER(intervalCodes);
	for(int j = 0; j < numInts; j++) {
		omxConfidenceInterval *oCI = Global->intervalList[j];
		omxRecompute(oCI->matrix, mle);
		interval(j, 1) = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);
		interval(j, 0) = std::min(oCI->min, interval(j, 1));
		interval(j, 2) = std::max(oCI->max, interval(j, 1));
		intervalCode[j] = oCI->lCode;
		intervalCode[j + numInts] = oCI->uCode;
	}
}

void ComputeCI::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	if (!intervals) return;

	int numInt = (int) Global->intervalList.size();

	SEXP dimnames;
	SEXP names;
	Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
	Rf_protect(names = Rf_allocVector(STRSXP, 3));
	SET_STRING_ELT(names, 0, Rf_mkChar("lbound"));
	SET_STRING_ELT(names, 1, Rf_mkChar("estimate"));
	SET_STRING_ELT(names, 2, Rf_mkChar("ubound"));
	SET_VECTOR_ELT(dimnames, 1, names);

	Rf_protect(names = Rf_allocVector(STRSXP, numInt)); //shared between the two matrices
	for (int nx=0; nx < numInt; ++nx) {
		omxConfidenceInterval *ci = Global->intervalList[nx];
		SET_STRING_ELT(names, nx, Rf_mkChar(ci->name));
	}
	SET_VECTOR_ELT(dimnames, 0, names);

	Rf_setAttrib(intervals, R_DimNamesSymbol, dimnames);

	out->add("confidenceIntervals", intervals);

	Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, names);

	Rf_protect(names = Rf_allocVector(STRSXP, 2));
	SET_STRING_ELT(names, 0, Rf_mkChar("lbound"));
	SET_STRING_ELT(names, 1, Rf_mkChar("ubound"));
	SET_VECTOR_ELT(dimnames, 1, names);

	Rf_setAttrib(intervalCodes, R_DimNamesSymbol, dimnames);

	out->add("confidenceIntervalCodes", intervalCodes);
}
