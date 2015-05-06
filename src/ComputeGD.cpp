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
#include "ComputeSD_AL.h"

enum OptEngine {
	OptEngine_NPSOL,
	OptEngine_CSOLNP,
    OptEngine_NLOPT,
    OptEngine_SD
};

class omxComputeGD : public omxCompute {
	typedef omxCompute super;
	enum OptEngine engine;
	omxMatrix *fitMatrix;
	int verbose;
	double optimalityTolerance;

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

	ScopedProtect p3(slotValue, R_do_slot(rObj, Rf_install("engine")));
	const char *engine_name = CHAR(Rf_asChar(slotValue));
	if (strEQ(engine_name, "CSOLNP")) {
		engine = OptEngine_CSOLNP;
	} else if (strEQ(engine_name, "SLSQP")) {
		engine = OptEngine_NLOPT;
	} else if (strEQ(engine_name, "NPSOL")) {
#if HAS_NPSOL
		engine = OptEngine_NPSOL;
#else
		Rf_error("NPSOL is not available in this build");
#endif
	} else if(strEQ(engine_name, "SD")){
		engine = OptEngine_SD;
	} else {
		Rf_error("%s: engine %s unknown", name, engine_name);
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

	omxFitFunctionPreoptimize(fitMatrix->fitFunction, fc);

	fc->createChildren();

	int beforeEval = Global->computeCount;

	//if (fc->CI) verbose=3;
	GradientOptimizerContext rf(verbose);
	rf.fc = fc;
	rf.fitMatrix = fitMatrix;
	rf.ControlTolerance = optimalityTolerance;
	rf.useGradient = useGradient;
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
		omxNPSOL(fc->est, rf);
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
		rf.avoidRedundentEvals = true;
		omxCSOLNP(fc->est, rf);
		if (rf.gradOut.size()) {
			fc->grad = rf.gradOut.tail(numParam);
			Eigen::Map< Eigen::MatrixXd > hess(fc->getDenseHessUninitialized(), numParam, numParam);
			hess = rf.hessOut.bottomRightCorner(numParam, numParam);
			fc->wanted |= FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN;
		}
		break;
        case OptEngine_NLOPT:
		omxInvokeNLOPT(fc->est, rf);
		fc->wanted |= FF_COMPUTE_GRADIENT;
		break;
        case OptEngine_SD:{
		rf.fc->copyParamToModel();
		rf.setupSimpleBounds();
		rf.setupIneqConstraintBounds();
		rf.solEqBFun();
		rf.myineqFun();
		int MAXIT = 50000;
		if(rf.inequality.size() == 0 && rf.equality.size() == 0)
			{
				omxSD(rf, MAXIT);   // unconstrained problems
			} else {
			omxSD_AL(rf);       // constrained problems
		}
		fc->wanted |= FF_COMPUTE_GRADIENT;
		break;}
        default: Rf_error("Optimizer %d is not available", engine);
	}

	fc->inform = rf.informOut;
	if (fc->inform <= 0 && Global->computeCount - beforeEval == 1) {
		fc->inform = INFORM_STARTING_VALUES_INFEASIBLE;
	}

	// Optimizers can terminate with inconsistent fit and parameters
	fc->copyParamToModel();
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
}

void ComputeCI::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	{
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(slotValue);
	}

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, Rf_install("class")), 0));
	plan = omxNewCompute(globalState, CHAR(s4class));
	plan->initFromFrontend(globalState, slotValue);
}

extern "C" { void F77_SUB(npoptn)(char* string, int Rf_length); };

class ciConstraintIneq : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxMatrix *fitMat;
 public:
	ciConstraintIneq(omxMatrix *fitMat) : super("CI"), fitMat(fitMat)
	{ size=1; opCode = LESS_THAN; };

	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out) {
		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		double diff = std::max(fit - fc->targetFit, 0.0);
		if (diff > 100) diff = nan("infeasible");
		if (ineqType != opCode) diff = -diff;
		//mxLog("fit %f diff %f", fit, diff);
		out[0] = diff;
	};
};

class ciConstraintEq : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxMatrix *fitMat;
 public:
	ciConstraintEq(omxMatrix *fitMat) : super("CI"), fitMat(fitMat)
	{ size=1; opCode = EQUALITY; };

	virtual void refreshAndGrab(FitContext *fc, Type ineqType, double *out) {
		omxFitFunctionCompute(fitMat->fitFunction, FF_COMPUTE_FIT, fc);
		const double fit = totalLogLikelihood(fitMat);
		double diff = fit - fc->targetFit;
		diff *= diff;
		//mxLog("fit %f diff %f", fit, diff);
		out[0] = diff;
	};
};

void markAsDataFrame(SEXP list)
{
	SEXP classes;
	Rf_protect(classes = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
	Rf_setAttrib(list, R_ClassSymbol, classes);
}

// Optimization: For profile CIs of free parameters, the gradient for
// the fit is trivial.  Many evaluations could be saved.

void ComputeCI::computeImpl(FitContext *mle)
{
	if (intervals) Rf_error("Can only compute CIs once");

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

	int totalIntervals = 0;
	for(int j = 0; j < numInts; j++) {
		omxConfidenceInterval *oCI = Global->intervalList[j];
		totalIntervals += std::isfinite(oCI->lbound) + std::isfinite(oCI->ubound);
	}

	Rf_protect(detail = Rf_allocVector(VECSXP, 3 + mle->numParam));
	SET_VECTOR_ELT(detail, 0, Rf_allocVector(STRSXP, totalIntervals));
	SET_VECTOR_ELT(detail, 1, Rf_allocVector(INTSXP, totalIntervals));
	for (int cx=0; cx < 1+int(mle->numParam); ++cx) {
		SET_VECTOR_ELT(detail, 2+cx, Rf_allocVector(REALSXP, totalIntervals));
	}

	SEXP detailCols;
	Rf_protect(detailCols = Rf_allocVector(STRSXP, 3 + mle->numParam));
	Rf_setAttrib(detail, R_NamesSymbol, detailCols);
	SET_STRING_ELT(detailCols, 0, Rf_mkChar("parameter"));
	SET_STRING_ELT(detailCols, 1, Rf_mkChar("lower"));
	SET_STRING_ELT(detailCols, 2, Rf_mkChar("fit"));
	for (int nx=0; nx < int(mle->numParam); ++nx) {
		SET_STRING_ELT(detailCols, 3+nx, Rf_mkChar(mle->varGroup->vars[nx]->name));
	}

	SEXP detailRowNames = Rf_allocVector(INTSXP, totalIntervals);
	Rf_setAttrib(detail, R_RowNamesSymbol, detailRowNames);
	markAsDataFrame(detail);

	const int ciMaxIterations = Global->ciMaxIterations;
	FitContext fc(mle, mle->varGroup);
	FreeVarGroup *freeVarGroup = fc.varGroup;

	const int n = int(freeVarGroup->vars.size());
	Eigen::Map< Eigen::VectorXd > Mle(mle->est, n);

	if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }

	ciConstraintIneq constr(fitMatrix);

	int detailRow = 0;
	for(int i = 0; i < (int) Global->intervalList.size(); i++) {
		omxConfidenceInterval *currentCI = Global->intervalList[i];

		const char *matName = "anonymous matrix";
		if (currentCI->matrix->name) {
			matName = currentCI->matrix->name;
		}

		for (int lower=0; lower <= 1; ++lower) {
			if (lower  && !std::isfinite(currentCI->lbound)) continue;
			if (!lower && !std::isfinite(currentCI->ubound)) continue;

			// Reset to previous optimum
			Eigen::Map< Eigen::VectorXd > Est(fc.est, n);
			Est = Mle;

			int tries = 0;
			int inform = INFORM_UNINITIALIZED;
			double *store = lower? &currentCI->min : &currentCI->max;
			Eigen::VectorXd bestEst = Mle;
			double bestFit = nan("uninit");

			while (++tries <= ciMaxIterations) {
				Global->checkpointMessage(mle, mle->est, "%s[%d, %d] %s CI (try %d)",
							  matName, currentCI->row + 1, currentCI->col + 1,
							  lower? "lower" : "upper", tries);

				mle->state->conList.push_back(&constr);
				fc.CI = currentCI;
				fc.lowerBound = lower;
				fc.fit = mle->fit;
				plan->compute(&fc);
				mle->state->conList.pop_back();

				omxRecompute(currentCI->matrix, &fc);
				double val = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);

				// We need to check the fit again because a random perturbation
				// could have moved us into a bad part of the space where we
				// got sucked into a local minimum. This could be addressed
				// by choosing a better random starting value.
				fc.CI = NULL;
				ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);

				bool better = (fc.inform != INFORM_STARTING_VALUES_INFEASIBLE &&
					       fc.fit - 1e-2 < fc.targetFit &&
					       ((!std::isfinite(*store) ||
						 (lower && val < *store) || (!lower && val > *store))));

				if(verbose>=2) {
					mxLog("CI[%d,%d] inform=%d best=%f val=%f fit-target=%f better=%d",
					      i, lower, inform, *store, val, fc.fit - fc.targetFit, better);
				}

				if (better) {
					*store = val;
					bestEst = Est;
					bestFit = fc.fit;

					inform = fc.inform;
					if (lower) currentCI->lCode = inform;
					else       currentCI->uCode = inform;
				}

				// Try different starting values, but try to stay close enough to
				// the MLE that we don't leave the feasible set.
				Eigen::VectorXd dir(fc.numParam);
				dir.setRandom();
				Est = Mle + ((Est - Mle).array() * dir.array()).matrix();
			}

			if (verbose >= 1) {
				mxLog("%s[%d, %d] %s CI %f fit %f", matName, currentCI->row + 1, currentCI->col + 1,
				      lower? "lower" : "upper", *store, bestFit);
			}

			INTEGER(detailRowNames)[detailRow] = 1 + detailRow;
			SET_STRING_ELT(VECTOR_ELT(detail, 0), detailRow, Rf_mkChar(currentCI->name));
			INTEGER(VECTOR_ELT(detail, 1))[detailRow] = lower;
			REAL(VECTOR_ELT(detail, 2))[detailRow] = bestFit;
			for (int px=0; px < int(fc.numParam); ++px) {
				REAL(VECTOR_ELT(detail, 3+px))[detailRow] = bestEst[px];
			}

			++detailRow;
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
		if (1) {
			interval(j, 0) = std::min(oCI->min, interval(j, 1));
			interval(j, 2) = std::max(oCI->max, interval(j, 1));
		} else {
			interval(j, 0) = oCI->min;
			interval(j, 2) = oCI->max;
		}
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

	MxRList output;
	output.add("detail", detail);
	slots->add("output", output.asR());
}
