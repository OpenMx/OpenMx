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
#include "ComputeSD.h"

#pragma GCC diagnostic warning "-Wshadow"

enum OptEngine {
	OptEngine_NPSOL,
	OptEngine_CSOLNP,
    OptEngine_NLOPT,
    OptEngine_SD
};

class omxComputeGD : public omxCompute {
	typedef omxCompute super;
	const char *engineName;
	enum OptEngine engine;
	const char *gradientAlgoName;
	enum GradientAlgorithm gradientAlgo;
	int gradientIterations;
	double gradientStepSize;
	omxMatrix *fitMatrix;
	int verbose;
	double optimalityTolerance;
	int maxIter;

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

omxComputeGD::omxComputeGD()
{
	hessChol = NULL;
	warmStart = NULL;
}

void omxComputeGD::initFromFrontend(omxState *globalState, SEXP rObj)
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
	if (!std::isfinite(optimalityTolerance)) {
		optimalityTolerance = Global->optimalityTolerance;
	}

	ScopedProtect p3(slotValue, R_do_slot(rObj, Rf_install("engine")));
	engineName = CHAR(Rf_asChar(slotValue));
	if (strEQ(engineName, "CSOLNP")) {
		engine = OptEngine_CSOLNP;
	} else if (strEQ(engineName, "SLSQP")) {
		engine = OptEngine_NLOPT;
	} else if (strEQ(engineName, "NPSOL")) {
#if HAS_NPSOL
		engine = OptEngine_NPSOL;
#else
		Rf_error("NPSOL is not available in this build. See ?omxGetNPSOL() to download this optimizer");
#endif
	} else if(strEQ(engineName, "SD")){
		engine = OptEngine_SD;
	} else {
		Rf_error("%s: engine %s unknown", name, engineName);
	}

	ScopedProtect p5(slotValue, R_do_slot(rObj, Rf_install("useGradient")));
	if (Rf_length(slotValue)) {
		useGradient = Rf_asLogical(slotValue);
	} else {
		useGradient = Global->analyticGradients;
	}

	ScopedProtect p4(slotValue, R_do_slot(rObj, Rf_install("nudgeZeroStarts")));
	nudge = Rf_asLogical(slotValue);

	ScopedProtect p6(slotValue, R_do_slot(rObj, Rf_install("warmStart")));
	if (!Rf_isNull(slotValue)) {
		SEXP matrixDims;
		ScopedProtect pws(matrixDims, Rf_getAttrib(slotValue, R_DimSymbol));
		int *dimList = INTEGER(matrixDims);
		int rows = dimList[0];
		int cols = dimList[1];
		if (rows != cols) Rf_error("%s: warmStart matrix must be square", name);

		warmStartSize = rows;
		warmStart = REAL(slotValue);
	}

	ScopedProtect p7(slotValue, R_do_slot(rObj, Rf_install("maxMajorIter")));
	if (Rf_length(slotValue)) {
		maxIter = Rf_asInteger(slotValue);
	} else {
		maxIter = -1; // different engines have different defaults
	}

	ScopedProtect p8(slotValue, R_do_slot(rObj, Rf_install("gradientAlgo")));
	gradientAlgoName = CHAR(Rf_asChar(slotValue));
	if (strEQ(gradientAlgoName, "forward")) {
		gradientAlgo = GradientAlgorithm_Forward;
	} else if (strEQ(gradientAlgoName, "central")) {
		gradientAlgo = GradientAlgorithm_Central;
	} else {
		Rf_error("%s: gradient algorithm '%s' unknown", name, gradientAlgoName);
	}

	ScopedProtect p9(slotValue, R_do_slot(rObj, Rf_install("gradientIterations")));
	gradientIterations = std::max(Rf_asInteger(slotValue), 1);

	ScopedProtect p10(slotValue, R_do_slot(rObj, Rf_install("gradientStepSize")));
	gradientStepSize = Rf_asReal(slotValue);
}

void omxComputeGD::computeImpl(FitContext *fc)
{
	omxFitFunctionPreoptimize(fitMatrix->fitFunction, fc);
	if (isErrorRaised()) return;

	size_t numParam = fc->varGroup->vars.size();
	if (fc->profiledOut.size()) {
		if (fc->profiledOut.size() != fc->numParam) Rf_error("Fail");
		for (size_t vx=0; vx < fc->varGroup->vars.size(); ++vx) {
			if (fc->profiledOut[vx]) --numParam;
		}
	}

	if (numParam <= 0) {
		omxRaiseErrorf("%s: model has no free parameters", name);
		return;
	}

	fc->ensureParamWithinBox(nudge);
	fc->createChildren();

	int beforeEval = fc->getComputeCount();

	if (verbose >= 1) mxLog("%s: engine %s (ID %d) gradient=%s tol=%g",
				name, engineName, engine, gradientAlgoName, optimalityTolerance);

	//if (fc->CI) verbose=3;
	GradientOptimizerContext rf(fc, verbose);
	rf.fitMatrix = fitMatrix;
	rf.ControlTolerance = optimalityTolerance;
	rf.useGradient = useGradient;
	rf.gradientAlgo = gradientAlgo;
	rf.gradientIterations = gradientIterations;
	rf.gradientStepSize = gradientStepSize;
	if (maxIter == -1) {
		rf.maxMajorIterations = -1;
	} else {
		rf.maxMajorIterations = fc->iterations + maxIter;
	}
	if (warmStart) {
		if (warmStartSize != int(numParam)) {
			Rf_warning("%s: warmStart size %d does not match number of free parameters %d (ignored)",
				   warmStartSize, numParam);
		} else {
			// Not sure if this code path works, need test TODO
			Eigen::Map< Eigen::MatrixXd > hessWrap(warmStart, numParam, numParam);
			rf.hessOut = hessWrap;
			rf.warmStart = true;
		}
	}

	switch (engine) {
        case OptEngine_NPSOL:{
#if HAS_NPSOL
		omxNPSOL(rf);
		rf.finish();
		fc->wanted |= FF_COMPUTE_GRADIENT;
		if (rf.hessOut.size() && fitMatrix->currentState->conList.size() == 0) {
			if (!hessChol) {
				Rf_protect(hessChol = Rf_allocMatrix(REALSXP, numParam, numParam));
			}
			Eigen::Map<Eigen::MatrixXd> hc(REAL(hessChol), numParam, numParam);
			hc = rf.hessOut;
			Eigen::Map<Eigen::MatrixXd> dest(fc->getDenseHessUninitialized(), numParam, numParam);
			dest.noalias() = rf.hessOut.transpose() * rf.hessOut;
			fc->wanted |= FF_COMPUTE_HESSIAN;
		}
#endif
		break;}
        case OptEngine_CSOLNP:
		if (rf.maxMajorIterations == -1) rf.maxMajorIterations = Global->majorIterations;
		rf.avoidRedundentEvals = true;
		rf.CSOLNP_HACK = true;
		omxCSOLNP(rf);
		rf.finish();
		if (rf.gradOut.size()) {
			fc->grad = rf.gradOut.tail(numParam);
			Eigen::Map< Eigen::MatrixXd > hess(fc->getDenseHessUninitialized(), numParam, numParam);
			hess = rf.hessOut.bottomRightCorner(numParam, numParam);
			fc->wanted |= FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN;
		}
		break;
        case OptEngine_NLOPT:
		rf.gradientStepSize = rf.gradientStepSize * GRADIENT_FUDGE_FACTOR(2.0);
		if (rf.maxMajorIterations == -1) rf.maxMajorIterations = Global->majorIterations;
		omxInvokeNLOPT(rf);
		rf.finish();
		fc->wanted |= FF_COMPUTE_GRADIENT;
		break;
        case OptEngine_SD:{
		fc->copyParamToModel();
		rf.setupSimpleBounds();
		rf.setupIneqConstraintBounds();
		rf.solEqBFun();
		rf.myineqFun();
		if(rf.inequality.size() == 0 && rf.equality.size() == 0) {
			omxSD(rf);   // unconstrained problems
			rf.finish();
		} else {
			Rf_error("Constrained problems are not implemented");
		}
		fc->wanted |= FF_COMPUTE_GRADIENT;
		break;}
        default: Rf_error("Optimizer %d is not available", engine);
	}

	if (!std::isfinite(fc->fit)) {
		fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);
	} else {
		fc->setInform(rf.informOut);
	}

	if (verbose >= 1) {
		mxLog("%s: engine %s done, iter=%d inform=%d",
		      name, engineName, fc->getComputeCount() - beforeEval, fc->getInform());
	}

	if (isErrorRaised()) return;

	// Optimizers can terminate with inconsistent fit and parameters
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, fc);

	if (verbose >= 2) {
		mxLog("%s: final fit is %2f", name, fc->fit);
		fc->log(FF_COMPUTE_ESTIMATE);
	}

	fc->wanted |= FF_COMPUTE_BESTFIT;
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	if (engine == OptEngine_NPSOL && hessChol) {
		out->add("hessianCholesky", hessChol);
	}
}

// -----------------------------------------------------------------------

class ComputeCI : public omxCompute {
	typedef omxCompute super;
	omxCompute *plan;
	omxMatrix *fitMatrix;
	int verbose;
	SEXP intervals, intervalCodes, detail;
	const char *ctypeName;
	bool useInequality;
	bool useEquality;

	enum Method {
		NEALE_MILLER_1997=1,
		WU_NEALE_2012=2
	};

	void regularCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int lower,
		       double &val, bool &better);
	void regularCI2(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow);
	void boundAdjCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow);
	void recordCI(Method meth, ConfidenceInterval *currentCI, int lower, FitContext &fc,
		      int &detailRow, double val, bool better);
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

ComputeCI::ComputeCI()
{
	intervals = 0;
	intervalCodes = 0;
	detail = 0;
	useInequality = false;
	useEquality = false;
}

void ComputeCI::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	{
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(slotValue);
	}
	{
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("constraintType")));
		ctypeName = CHAR(Rf_asChar(slotValue));
		if (strEQ(ctypeName, "ineq")) {
			useInequality = true;
		} else if (strEQ(ctypeName, "eq")) {
			useEquality = true;
		} else if (strEQ(ctypeName, "both")) {
			useEquality = true;
			useInequality = true;
		} else if (strEQ(ctypeName, "none")) {
			// OK
		} else {
			Rf_error("%s: unknown constraintType='%s'", name, ctypeName);
		}
	}

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = omxNewCompute(globalState, CHAR(s4class));
	plan->initFromFrontend(globalState, slotValue);
}

extern "C" { void F77_SUB(npoptn)(char* string, int Rf_length); };

class ciConstraint : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxState *state;
public:
	omxMatrix *fitMat;
	ciConstraint() : super("CI"), state(0) {};
	virtual ~ciConstraint() {
		pop();
	};
	void push(omxState *_state) {
		state = _state;
		state->conList.push_back(this);
	};
	void pop() {
		if (!state) return;
		size_t sz = state->conList.size();
		if (sz && state->conList[sz-1] == this) {
			state->conList.pop_back();
		}
		state = 0;
	};
};

class ciConstraintIneq : public ciConstraint {
 private:
	typedef ciConstraint super;
 public:
	ciConstraintIneq()
	{ size=1; opCode = LESS_THAN; };

	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out) {
		*out = fc->ciobj->evalIneq(fc, fitMat);
		if (ineqType != opCode) *out = - *out;
		//mxLog("fit %f diff %f", fit, diff);
	};
};

class ciConstraintEq : public ciConstraint {
 private:
	typedef ciConstraint super;
 public:
	ciConstraintEq()
	{ size=1; opCode = EQUALITY; };

	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out) {
		*out = fc->ciobj->evalEq(fc, fitMat);
		//mxLog("fit %f diff %f", fit, diff);
	};
};

static SEXP makeFactor(SEXP vec, int levels, const char **labels)
{
	SEXP classes;
	Rf_protect(classes = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("factor"));
	Rf_setAttrib(vec, R_ClassSymbol, classes);

	SEXP Rlev;
	Rf_protect(Rlev = Rf_allocVector(STRSXP, levels));
	for (int lx=0; lx < levels; ++lx) {
		SET_STRING_ELT(Rlev, lx, Rf_mkChar(labels[lx]));
	}

	Rf_setAttrib(vec, Rf_install("levels"), Rlev);
	return vec;
}

void ComputeCI::recordCI(Method meth, ConfidenceInterval *currentCI, int lower, FitContext &fc,
			 int &detailRow, double val, bool better)
{
	omxMatrix *ciMatrix = currentCI->getMatrix(fitMatrix->currentState);
	std::string &matName = ciMatrix->nameStr;

	if (better) {
		currentCI->val[!lower] = val;
		currentCI->code[!lower] = fc.getInform();
	}

	if(verbose >= 1) {
		mxLog("CI[%s,%s] %s[%d,%d] val=%f fit=%f accepted=%d",
		      currentCI->name.c_str(), (lower?"lower":"upper"), matName.c_str(),
		      1+currentCI->row, 1+currentCI->col, val, fc.fit, better);
	}

	SET_STRING_ELT(VECTOR_ELT(detail, 0), detailRow, Rf_mkChar(currentCI->name.c_str()));
	REAL(VECTOR_ELT(detail, 1))[detailRow] = val;
	INTEGER(VECTOR_ELT(detail, 2))[detailRow] = 1+lower;
	REAL(VECTOR_ELT(detail, 3))[detailRow] = fc.fit;
	Eigen::Map< Eigen::VectorXd > Est(fc.est, fc.numParam);
	for (int px=0; px < int(fc.numParam); ++px) {
		REAL(VECTOR_ELT(detail, 4+px))[detailRow] = Est[px];
	}
	INTEGER(VECTOR_ELT(detail, 4+fc.numParam))[detailRow] = meth;
	++detailRow;
}

struct regularCIobj : CIobjective {
	double targetFit;
	bool lowerBound;
	bool compositeCIFunction;

	virtual bool gradientKnown()
	{
		return CI->varIndex >= 0;
	}

	virtual void gradient(FitContext *fc, double *gradOut) {
		Eigen::Map< Eigen::VectorXd > Egrad(gradOut, fc->numParam);
		Egrad.setZero();
		Egrad[ CI->varIndex ] = lowerBound? 1 : -1;
	}

	virtual double evalIneq(FitContext *fc, omxMatrix *fitMat)
	{
		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		double diff = std::max(fit - targetFit, 0.0);
		if (diff > 100) diff = nan("infeasible");
		return diff;
	}

	virtual double evalEq(FitContext *fc, omxMatrix *fitMat)
	{
		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		double diff = fit - targetFit;
		diff *= diff;
		if (fabs(diff) > 100000) diff = nan("infeasible");
		return diff;
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc)
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			Rf_error("Not implemented yet");
		}

		// We need to compute the fit here because that's the only way to
		// check our soft feasibility constraints. If parameters don't
		// change between here and the constraint evaluation then we
		// should avoid recomputing the fit again in the constraint. TODO

		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);
		omxResizeMatrix(fitMat, 1, 1);

		if (!std::isfinite(fit) || !std::isfinite(CIElement)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fitMat->data[0] = nan("infeasible");
			return;
		}

		if (want & FF_COMPUTE_FIT) {
			double param = (lowerBound? CIElement : -CIElement);
			if (compositeCIFunction) {
				double diff = targetFit - fit;
				diff *= diff;
				if (diff > 1e2) {
					// Ensure there aren't any creative solutions
					diff = nan("infeasible");
					return;
				}
				fitMat->data[0] = diff + param;
			} else {
				fitMat->data[0] = param;
			}
			//mxLog("param at %f", fitMat->data[0]);
		}
		if (want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
			// add deriv adjustments here TODO
		}
	}
};

struct bound1CIobj : CIobjective {
	double bound;

	virtual bool gradientKnown() { return false; }
	
	virtual double evalEq(FitContext *fc, omxMatrix *fitMat)
	{
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);
		double diff = CIElement - bound;
		diff *= diff;
		if (fabs(diff) > 100000) diff = nan("infeasible");
		return diff;
	}
};

struct boundAwayCIobj : CIobjective {
	double logAlpha, sqrtCrit;
	double boundLL, bestLL;
	int lower;
	bool constrained;

	virtual bool gradientKnown() { return false; }
	
	template <typename T1>
	void computeConstraint(double fit, Eigen::ArrayBase<T1> &v1)
	{
		double d1 = sqrt(std::max(fit - boundLL,0.0));
		double d2 = sqrt(std::max(fit - bestLL, 0.0));
		//Rf_pnorm(d1,0,1,0,0) + Rf_pnorm(d2,0,1,0,0) > alpha
		//d1 < sqrtCrit
		//d2 < sqrtCrit
		double pA = Rf_pnorm5(d1,0,1,0,0) + Rf_pnorm5(d2,0,1,0,0);
		v1 << std::max(d1 - sqrtCrit, 0.0), std::max(d2 - sqrtCrit, 0.0),
			std::max(logAlpha - log(pA), 0.0);
		// mxPrintMat("v1", v1);
		// mxLog("fit %g sqrtCrit %g d1 %g d2 %g alpha %g pA %g",
		//       fit, sqrtCrit, d1, d2, exp(logAlpha), pA);
	}

	virtual double evalIneq(FitContext *fc, omxMatrix *fitMat)
	{
		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		Eigen::Array<double,3,1> v1;
		computeConstraint(fit, v1);
		double diff = v1.sum();
		if (diff > 100) diff = nan("infeasible");
		return diff;
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc)
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			Rf_error("Not implemented yet");
		}

		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);
		omxResizeMatrix(fitMat, 1, 1);

		if (!std::isfinite(fit) || !std::isfinite(CIElement)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fitMat->data[0] = nan("infeasible");
			return;
		}

		double cval = 0.0;
		if (constrained) {
			Eigen::Array<double,3,1> v1;
			computeConstraint(fit, v1);
			cval = v1.sum()*v1.sum();
			if (cval > 10000) cval = nan("infeasible");
		}
			     
		double param = (lower? CIElement : -CIElement);
		fitMat->data[0] = param + cval;
		//mxLog("param at %f", fitMat->data[0]);
	}
};

struct boundNearCIobj : CIobjective {
	double d0, sqrtCrit90, sqrtCrit95, logAlpha;
	double boundLL, bestLL;
	int lower;

	virtual bool gradientKnown() { return false; }
	
	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc)
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			Rf_error("Not implemented yet");
		}

		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);
		omxResizeMatrix(fitMat, 1, 1);

		if (!std::isfinite(fit) || !std::isfinite(CIElement)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fitMat->data[0] = nan("infeasible");
			return;
		}

		double lbd = std::max(d0/2, sqrtCrit90); //precompute TODO
		double ubd = std::min(d0, sqrtCrit95); //precompute TODO
		double dd = sqrt(std::max(fit - bestLL, 0.0));
		double pN = Rf_pnorm5(dd,0,1,0,0) + Rf_pnorm5((d0-dd)/2.0+dd*dd/(2*std::max(d0-dd,.001*dd*dd)),0,1,0,0);
		// pN > alpha
		// lbd < dd < ubd
		Eigen::Array<double,3,1> v1;
		v1 << std::max(lbd - dd, 0.0), std::max(dd - ubd, 0.0), std::max(logAlpha - log(pN), 0.0);
		//mxPrintMat("v1", v1);
		//mxLog("fit %g CIElement %g dd %g/%g/%g alpha %g pN %g",
		//fit, CIElement, lbd, dd, ubd, exp(logAlpha), pN);
			     
		double param = (lower? CIElement : -CIElement);
		fitMat->data[0] = param + v1.sum()*v1.sum();
		//mxLog("param at %f", fitMat->data[0]);
	}
};

// Must be free parameter (not algebra) and only 1 bound set
// Enable by default? modest performance cost
// Check against Hau Wu's code
// Symmetry test with different critical values
// For paper, growth curve factor variance, ACE
// Pek, J. & Wu, H. (in press). Profile likelihood-based confidence intervals and regions for structural equation models.
// \emph{Psychometrica.}
// Also look at 13 May 2015 email to pek@yorku.ca
void ComputeCI::boundAdjCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow)
{
	// need to skip some stuff if one bound is not wanted TODO
	omxState *state = fitMatrix->currentState;
	omxMatrix *ciMatrix = currentCI->getMatrix(state);
	std::string &matName = ciMatrix->nameStr;
	omxFreeVar *fv = mle->varGroup->vars[currentCI->varIndex];
	double bound;
	int side;
	if (fv->lbound == NEG_INF) {
		bound = fv->ubound;
		side = ConfidenceInterval::Upper;
	} else {
		bound = fv->lbound;
		side = ConfidenceInterval::Lower;
	}

	// Reset to previous optimum
	Eigen::Map< Eigen::VectorXd > Mle(mle->est, mle->numParam);
	Eigen::Map< Eigen::VectorXd > Est(fc.est, fc.numParam);
	Est = Mle;

	Global->checkpointMessage(mle, mle->est, "%s[%d, %d] at-bound CI",
				  matName.c_str(), currentCI->row + 1, currentCI->col + 1);

	ciConstraintEq constrEq;
	constrEq.fitMat = fitMatrix;
	constrEq.push(state);
	bound1CIobj ciobj;
	ciobj.CI = currentCI;
	ciobj.bound = bound;
	fc.ciobj = &ciobj;
	plan->compute(&fc);
	constrEq.pop();

	// check fc.inform and bail if something went wrong TODO
	double boundLL = fc.fit;
	if (false && boundLL - mle->fit > currentCI->bound[side]) {		// TODO remove?
		regularCI2(mle, fc, currentCI, detailRow);
		return;
	}

	if (currentCI->bound[!side]) {	// ------------------------------ away from bound side --
		//mxLog("boundLL %g bestLL %g", boundLL, mle->fit);
		Global->checkpointMessage(mle, mle->est, "%s[%d, %d] away side CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
		ciConstraintIneq constr;
		constr.fitMat = fitMatrix;
		constr.push(state);
		boundAwayCIobj baobj;
		baobj.CI = currentCI;
		baobj.boundLL = boundLL;
		baobj.bestLL = mle->fit;
		baobj.logAlpha = log(1.0 - Rf_pchisq(currentCI->bound[!side], 1, 1, 0));
		baobj.sqrtCrit = sqrt(currentCI->bound[!side]);
		baobj.lower = side;
		baobj.constrained = true;
		Est = Mle;
		fc.ciobj = &baobj;
		plan->compute(&fc);
		constr.pop();

		// check fc.inform and bail if something went wrong TODO
		omxRecompute(ciMatrix, &fc);
		double val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

		fc.ciobj = 0;
		ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);
		recordCI(WU_NEALE_2012, currentCI, side, fc, detailRow, val, true);
	}

	if (currentCI->bound[side]) {     // ------------------------------ near to bound side --
		Global->checkpointMessage(mle, mle->est, "%s[%d, %d] near side CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
		double sqrtCrit95 = sqrt(currentCI->bound[side]);
		double sqrtCrit90 = sqrt(Rf_qchisq(1-(2.0 * (1-Rf_pchisq(currentCI->bound[side], 1, 1, 0))),1,1,0));
		double d0 = sqrt(std::max(boundLL - mle->fit, 0.0));
		if (d0 < sqrtCrit90) {
			mxLog("near side is too close to bound");
			fc.fit = boundLL;
			recordCI(WU_NEALE_2012, currentCI, !side, fc, detailRow, bound, true); //too close to bound
			return;
		}
	
		if (sqrtCrit95 <= d0/2.0) {
			mxLog("near side is far enough away that no correction is needed");
			bool better;
			double val;
			regularCI(mle, fc, currentCI, !side, val, better);
			recordCI(NEALE_MILLER_1997, currentCI, !side, fc, detailRow, val, true);
			return;
		}
	
		mxLog("near side getting fancy adjustment");
		boundNearCIobj bnobj;
		bnobj.CI = currentCI;
		bnobj.d0 = d0;
		bnobj.sqrtCrit90 = sqrtCrit90;
		bnobj.sqrtCrit95 = sqrtCrit95;
		bnobj.boundLL = boundLL;
		bnobj.bestLL = mle->fit;
		bnobj.logAlpha = log(1.0 - Rf_pchisq(currentCI->bound[side], 1, 1, 0));
		bnobj.lower = !side;
		Est = Mle;
		fc.ciobj = &bnobj;
		plan->compute(&fc);

		// check fc.inform and bail if something went wrong TODO
		omxRecompute(ciMatrix, &fc);
		double val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

		fc.ciobj = 0;
		ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);
		recordCI(WU_NEALE_2012, currentCI, !side, fc, detailRow, val, true);
	}
}

void ComputeCI::regularCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int lower,
			  double &val, bool &better)
{
	omxState *state = fitMatrix->currentState;

	ciConstraintEq constrEq;
	ciConstraintIneq constrIneq;
	ciConstraint *constr = 0;
	bool constrained = useInequality || useEquality;
	if (constrained) {
		constr = useInequality? (ciConstraint*)&constrIneq : (ciConstraint*)&constrEq;
		constr->fitMat = fitMatrix;
		constr->push(state);
	}
	
	// Reset to previous optimum
	Eigen::Map< Eigen::VectorXd > Mle(mle->est, mle->numParam);
	Eigen::Map< Eigen::VectorXd > Est(fc.est, fc.numParam);
	Est = Mle;

	regularCIobj ciobj;
	ciobj.CI = currentCI;
	ciobj.compositeCIFunction = !constrained;
	ciobj.lowerBound = lower;
	ciobj.targetFit = currentCI->bound[!lower] + mle->fit;
	fc.ciobj = &ciobj;
	//mxLog("Set target fit to %f (MLE %f)", fc->targetFit, fc->fit);

	plan->compute(&fc);

	omxMatrix *ciMatrix = currentCI->getMatrix(fitMatrix->currentState);
	omxRecompute(ciMatrix, &fc);
	val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

	// We check the fit again so we can report it
	// in the detail data.frame.
	fc.ciobj = 0;
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);

	double dist = currentCI->bound[!lower];
	better = (!constrained || fabs(fc.fit - ciobj.targetFit) < (dist * .05));
}

void ComputeCI::regularCI2(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow)
{
	omxState *state = fitMatrix->currentState;
	omxMatrix *ciMatrix = currentCI->getMatrix(state);
	std::string &matName = ciMatrix->nameStr;

	for (int upper=0; upper <= 1; ++upper) {
		int lower = 1-upper;
		if (!(currentCI->bound[upper])) continue;

		Global->checkpointMessage(mle, mle->est, "%s[%d, %d] %s CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1,
					  upper? "upper" : "lower");
		double val;
		bool better;
		regularCI(mle, fc, currentCI, lower, val, better);
		recordCI(NEALE_MILLER_1997, currentCI, lower, fc, detailRow, val, better);
	}
}

void ComputeCI::computeImpl(FitContext *mle)
{
	if (intervals) Rf_error("Can only compute CIs once");
	if (!Global->intervals) {
		if (verbose >= 1) mxLog(name, "%s: mxRun(..., intervals=FALSE), skipping");
		return;
	}

	Global->unpackConfidenceIntervals(mle->state);

	// Not strictly necessary, but makes it easier to run
	// mxComputeConfidenceInterval alone without other compute
	// steps.
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, mle);

	int numInts = (int) Global->intervalList.size();
	if (verbose >= 1) mxLog("%s: %d intervals of '%s' (ref fit %f %s)",
				name, numInts, fitMatrix->name(), mle->fit, ctypeName);
	if (!numInts) return;

	if (!mle->haveReferenceFit(fitMatrix)) return;

	// I'm not sure why INFORM_NOT_AT_OPTIMUM is okay, but that's how it was.
	if (mle->getInform() >= INFORM_LINEAR_CONSTRAINTS_INFEASIBLE && mle->getInform() != INFORM_NOT_AT_OPTIMUM) {
		// TODO: allow forcing
		Rf_warning("Not calculating confidence intervals because of optimizer status %d", mle->getInform());
		return;
	}

	Rf_protect(intervals = Rf_allocMatrix(REALSXP, numInts, 3));
	Rf_protect(intervalCodes = Rf_allocMatrix(INTSXP, numInts, 2));

	int totalIntervals = 0;
	for(int j = 0; j < numInts; j++) {
		ConfidenceInterval *oCI = Global->intervalList[j];
		totalIntervals += (oCI->bound != 0.0).count();
	}

	Rf_protect(detail = Rf_allocVector(VECSXP, 5 + mle->numParam));
	SET_VECTOR_ELT(detail, 0, Rf_allocVector(STRSXP, totalIntervals));
	SET_VECTOR_ELT(detail, 1, Rf_allocVector(REALSXP, totalIntervals));
	const char *sideLabels[] = { "upper", "lower" };
	SET_VECTOR_ELT(detail, 2, makeFactor(Rf_allocVector(INTSXP, totalIntervals), 2, sideLabels));
	for (int cx=0; cx < 1+int(mle->numParam); ++cx) {
		SET_VECTOR_ELT(detail, 3+cx, Rf_allocVector(REALSXP, totalIntervals));
	}
	const char *methodLabels[] = { "neale-miller-1997", "wu-neale-2012" };
	SET_VECTOR_ELT(detail, 4+mle->numParam,
		       makeFactor(Rf_allocVector(INTSXP, totalIntervals), 2, methodLabels));

	SEXP detailCols;
	Rf_protect(detailCols = Rf_allocVector(STRSXP, 5 + mle->numParam));
	Rf_setAttrib(detail, R_NamesSymbol, detailCols);
	SET_STRING_ELT(detailCols, 0, Rf_mkChar("parameter"));
	SET_STRING_ELT(detailCols, 1, Rf_mkChar("value"));
	SET_STRING_ELT(detailCols, 2, Rf_mkChar("side"));
	SET_STRING_ELT(detailCols, 3, Rf_mkChar("fit"));
	for (int nx=0; nx < int(mle->numParam); ++nx) {
		SET_STRING_ELT(detailCols, 4+nx, Rf_mkChar(mle->varGroup->vars[nx]->name));
	}
	SET_STRING_ELT(detailCols, 4 + mle->numParam, Rf_mkChar("method"));

	markAsDataFrame(detail, totalIntervals);

	omxState *state = fitMatrix->currentState;
	FitContext fc(mle, mle->varGroup);
	FreeVarGroup *freeVarGroup = fc.varGroup;

	int detailRow = 0;
	for(int ii = 0; ii < (int) Global->intervalList.size(); ii++) {
		ConfidenceInterval *currentCI = Global->intervalList[ii];
		omxMatrix *ciMatrix = currentCI->getMatrix(state);

		currentCI->varIndex = freeVarGroup->lookupVar(ciMatrix, currentCI->row, currentCI->col);

		if (currentCI->boundAdj && currentCI->varIndex < 0) {
			currentCI->boundAdj = false;
		}
		if (currentCI->boundAdj) {
			omxFreeVar *fv = freeVarGroup->vars[currentCI->varIndex];
			if (!((fv->lbound == NEG_INF) ^ (fv->ubound == INF))) {
				currentCI->boundAdj = false;
			}
		}
		if (currentCI->boundAdj) {
			boundAdjCI(mle, fc, currentCI, detailRow);
		} else {
			regularCI2(mle, fc, currentCI, detailRow);
		}
	}

	mle->copyParamToModel();

	Eigen::Map< Eigen::ArrayXXd > interval(REAL(intervals), numInts, 3);
	interval.fill(NA_REAL);
	int* intervalCode = INTEGER(intervalCodes);
	for(int j = 0; j < numInts; j++) {
		ConfidenceInterval *oCI = Global->intervalList[j];
		omxMatrix *ciMat = oCI->getMatrix(state);
		omxRecompute(ciMat, mle);
		interval(j, 1) = omxMatrixElement(ciMat, oCI->row, oCI->col);
		interval(j, 0) = std::min(oCI->val[ConfidenceInterval::Lower], interval(j, 1));
		interval(j, 2) = std::max(oCI->val[ConfidenceInterval::Upper], interval(j, 1));
		intervalCode[j] = oCI->code[ConfidenceInterval::Lower];
		intervalCode[j + numInts] = oCI->code[ConfidenceInterval::Upper];
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
		ConfidenceInterval *ci = Global->intervalList[nx];
		SET_STRING_ELT(names, nx, Rf_mkChar(ci->name.c_str()));
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

	MxRList output;
	output.add("detail", detail);
	slots->add("output", output.asR());
}
