/*
 *  Copyright 2013-2021 by the individuals mentioned in the source code history
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
#include "omxNLopt.h"
#include "Compute.h"
#include "npsolswitch.h"
#include "glue.h"
#include "ComputeSD.h"
#include "ComputeGD.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <RcppEigenCholmod.h>
#include <RcppEigenWrap.h>
#include "asa.h"

#include "EnableWarnings.h"

void GradientOptimizerContext::copyBounds()
{
	fc->copyBoxConstraintToOptimizer(solLB, solUB);
}

void GradientOptimizerContext::setupSimpleBounds()
{
	solLB.resize(numFree);
	solUB.resize(numFree);
	copyBounds();
}

bool GradientOptimizerContext::isUnconstrained()
{ return fc->isUnconstrained(); }

void GradientOptimizerContext::setupAllBounds() //used with NPSOL.
{
	omxState *st = fc->state;
	int n = (int) numFree;

	// Treat all constraints as non-linear.
	// Does the special handling of linear constraints work in NPSOL version 6? TODO
	int ncnln = AllC.getCount();
	solLB.resize(n + ncnln);
	solUB.resize(n + ncnln);

	copyBounds();

	int index = n;
	for(int constraintIndex = 0; constraintIndex < int(st->conListX.size()); constraintIndex++) {
		omxConstraint &cs = *st->conListX[constraintIndex];
		omxConstraint::Type type = cs.opCode;
		switch(type) {
		case omxConstraint::LESS_THAN:
		case omxConstraint::GREATER_THAN:
			for(int offset = 0; offset < cs.size; offset++) {
        if (cs.redundant[offset]) continue;
				solLB[index] = NEG_INF;
				solUB[index] = -0.0;
				index++;
			}
			break;
		case omxConstraint::EQUALITY:
			for(int offset = 0; offset < cs.size; offset++) {
        if (cs.redundant[offset]) continue;
				solLB[index] = -0.0;
				solUB[index] = 0.0;
				index++;
			}
			break;
		default:
			mxThrow("Unknown constraint type %d", type);
		}
	}
}

void GradientOptimizerContext::reset()
{
	feasible = false;
	bestFit = std::numeric_limits<double>::max();
	eqNorm = 0;
	ineqNorm = 0;
	iterations = 0;
}

GradientOptimizerContext::GradientOptimizerContext(FitContext *u_fc, int u_verbose,
						   omxCompute *owner)
	: fc(u_fc), verbose(u_verbose), numFree(countNumFree()),
	  numOptimizerThreads(u_fc->numOptimizerThreads()),
    AllC(u_fc, "constraint",
         [](const omxConstraint &con){ return true; }),
    IneqC(u_fc, "ineq",
          [](const omxConstraint &con){ return con.opCode != omxConstraint::EQUALITY; }),
    EqC(u_fc, "eq",
        [](const omxConstraint &con){ return con.opCode == omxConstraint::EQUALITY; })
{
	computeName = owner->name;
	fitMatrix = NULL;
	ControlMinorLimit = 800;
	ControlRho = 1.0;
	ControlTolerance = nan("uninit");
	warmStart = false;
	est.resize(numFree);
	grad.resize(numFree);
	copyToOptimizer(est.data());
	reset();
}

double GradientOptimizerContext::recordFit(double *myPars, int* mode)
{
	double fit = solFun(myPars, mode);
	if (std::isfinite(fit) && fit < bestFit && !fc->skippedRows) {
		bestFit = fit;
		Eigen::Map< Eigen::VectorXd > pvec(myPars, fc->numParam);
		bestEst = pvec;
	}
	return fit;
}

void GradientOptimizerContext::copyToOptimizer(double *myPars)
{
	Eigen::Map<Eigen::VectorXd> vec(myPars, numFree);
	fc->copyEstToOptimizer(vec);
}

void GradientOptimizerContext::useBestFit()
{
	fc->setUnscaledFit(bestFit);
	est = bestEst;
	// restore gradient too? TODO
}

void GradientOptimizerContext::copyFromOptimizer(const double *myPars, FitContext *fc2)
{
	Eigen::Map<const Eigen::VectorXd> vec(myPars, numFree);
	fc2->setEstFromOptimizer(vec);
}

void GradientOptimizerContext::finish()
{
	fc->setEstGradFromOptimizer(est, grad);
}

double GradientOptimizerContext::solFun(double *myPars, int* mode)
{
	Eigen::Map< Eigen::VectorXd > Est(myPars, fc->numParam);
	if (*mode == 1) {
		fc->iterations += 1;
		fc->resetOrdinalRelativeError();
	}
	copyFromOptimizer(myPars, fc);

	int want = FF_COMPUTE_FIT;
	if (*mode > 0) want |= FF_COMPUTE_GRADIENT;

	ComputeFit(getOptName(), fitMatrix, want, fc);

	if (*mode == 1) Global->reportProgress(getOptName(), fc);

	if (fc->outsideFeasibleSet() || isErrorRaised()) {
		*mode = -1;
	} else {
    if (!feasible) {
      feasible = true;
      if (verbose >= 2) mxLog("%s: Congratulations! Starting values are feasible!", getOptName());
    }
		if (want & FF_COMPUTE_GRADIENT) {
			fc->copyGradToOptimizer(grad);
		}
	}

	if (verbose >= 3) {
		mxLog("fit %.12f (mode %d)", fc->getUnscaledFit(), *mode);
	}

	return fc->getUnscaledFit();
}

// ------------------------------------------------------------

class omxComputeGD : public omxCompute {
	typedef omxCompute super;
	enum OptEngine engine;
	enum GradientAlgorithm gradientAlgo;
	int gradientIterations;
	double gradientStepSize;
	omxMatrix *fitMatrix;
	int verbose;
	double optimalityTolerance;
	int maxIter;

	int nudge;

	Eigen::MatrixXd hessChol;  // in-out for warmstart
	int warmStartSize;
	double *warmStart;
	int threads;

public:
	omxComputeGD();
	virtual void initFromFrontend(omxState *, SEXP rObj) override;
	virtual void computeImpl(FitContext *fc) override;
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out) override;
};

class omxCompute *newComputeGradientDescent()
{
	return new omxComputeGD();
}

omxComputeGD::omxComputeGD()
{
	warmStart = NULL;
}

void omxComputeGD::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	omxCompleteFitFunction(fitMatrix);

	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);

	ScopedProtect p2(slotValue, R_do_slot(rObj, Rf_install("tolerance")));
	optimalityTolerance = Rf_asReal(slotValue);
	if (!std::isfinite(optimalityTolerance)) {
		optimalityTolerance = Global->optimalityTolerance;
	}

	ScopedProtect p3(slotValue, R_do_slot(rObj, Rf_install("engine")));
  engine = nameToGradOptEngine(CHAR(Rf_asChar(slotValue)));

	ScopedProtect p4(slotValue, R_do_slot(rObj, Rf_install("nudgeZeroStarts")));
	nudge = false;
	friendlyStringToLogical("nudgeZeroStarts", slotValue, &nudge);

	ScopedProtect p6(slotValue, R_do_slot(rObj, Rf_install("warmStart")));
	if (!Rf_isNull(slotValue)) {
		SEXP matrixDims;
		ScopedProtect pws(matrixDims, Rf_getAttrib(slotValue, R_DimSymbol));
		int *dimList = INTEGER(matrixDims);
		int rows = dimList[0];
		int cols = dimList[1];
		if (rows != cols) mxThrow("%s: warmStart matrix must be square", name);

		warmStartSize = rows;
		warmStart = REAL(slotValue);
	}

	ScopedProtect p7(slotValue, R_do_slot(rObj, Rf_install("maxMajorIter")));
	if (Rf_length(slotValue)) {
		maxIter = Rf_asInteger(slotValue);
	} else {
		maxIter = -1; // different engines have different defaults
	}
}

void omxComputeGD::computeImpl(FitContext *fc)
{
	omxAlgebraPreeval(fitMatrix, fc);

	int numParam = fc->getNumFree();

	if (numParam <= 0) { complainNoFreeParam(); return; }

	fc->ensureParamWithinBox(nudge);
	fc->createChildren(fitMatrix, true);

	int beforeEval = fc->getLocalComputeCount();

	if (verbose >= 1) {
		int numConstr = fitMatrix->currentState->conListX.size();
		mxLog("%s: engine %d #P=%d tol=%g constraints=%d",
		      name, engine, int(numParam), optimalityTolerance,
		      numConstr);
	}

	//if (fc->ciobj) verbose=2;
	GradientOptimizerContext rf(fc, verbose, this);
	threads = rf.numOptimizerThreads;
	rf.fitMatrix = fitMatrix;
	rf.ControlTolerance = optimalityTolerance;
	if (maxIter == -1) {
		rf.maxMajorIterations = -1;
	} else {
		rf.maxMajorIterations = fc->iterations + maxIter;
	}
	/*Arguably, this effort involving a warm start should be conditioned on use of NPSOL.
	But hopefully, we'll support warm starts for other optimizers someday.*/
	if (warmStart) {
		if (warmStartSize != int(numParam)) {
			Rf_warning("%s: warmStart size %d does not match number of free parameters %d (ignored)",
				   name, warmStartSize, numParam);
		} else {
			Eigen::Map< Eigen::MatrixXd > hessWrap(warmStart, numParam, numParam);
			rf.hessOut = hessWrap;
			rf.warmStart = true;
		}
	}
	else{
		if(fc->wanted & FF_COMPUTE_HESSIAN){
			rf.hessOut.setZero(numParam,numParam);
			Eigen::LLT< Eigen::MatrixXd > chol4WS(numParam);
			fc->refreshDenseHess();
			fc->copyDenseHess(rf.hessOut.data());
			chol4WS.compute(rf.hessOut);
			if(chol4WS.info() == Eigen::Success){
				rf.hessOut = (Eigen::MatrixXd)(chol4WS.matrixU());
				rf.warmStart = true;
			}
			else{
				if(rf.verbose >= 1){
					mxLog("Hessian not positive-definite at initial values");
				}
			}
		}
		//else if(fc->wanted & FF_COMPUTE_IHESSIAN)
	}

	switch (engine) {
        case OptEngine_NPSOL:{
#if HAS_NPSOL
		omxNPSOL(rf);
    fc->initGrad();
		rf.finish();
		fc->wanted |= FF_COMPUTE_GRADIENT;
		if (rf.hessOut.size() ){
			hessChol = rf.hessOut;
			Eigen::Map<Eigen::MatrixXd> dest(fc->getDenseHessUninitialized(), numParam, numParam);
			dest.noalias() = rf.hessOut.transpose() * rf.hessOut;
			fc->wanted |= FF_COMPUTE_HESSIAN;
		}
		fc->constraintFunVals = rf.constraintFunValsOut;
		fc->constraintJacobian = rf.constraintJacobianOut;
		fc->LagrMultipliers = rf.LagrMultipliersOut;
		fc->constraintStates = rf.constraintStatesOut;
		//LagrHessian?
#endif
		break;}
        case OptEngine_CSOLNP:
		if (rf.maxMajorIterations == -1) rf.maxMajorIterations = Global->majorIterations;
		rf.IneqC.setIneqAlwaysActive();
		omxCSOLNP(rf);
		rf.finish();
		if (rf.gradOut.size()) {
			fc->gradZ = rf.gradOut.tail(numParam);
			Eigen::Map< Eigen::MatrixXd > hess(fc->getDenseHessUninitialized(), numParam, numParam);
			hess = rf.hessOut.bottomRightCorner(numParam, numParam);
			fc->wanted |= FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN;
		}
		fc->constraintFunVals = rf.constraintFunValsOut;
		fc->LagrMultipliers = rf.LagrMultipliersOut;
		fc->constraintJacobian = rf.constraintJacobianOut;
		fc->LagrHessian = rf.LagrHessianOut;
		break;
        case OptEngine_NLOPT:
		if (rf.maxMajorIterations == -1) rf.maxMajorIterations = Global->majorIterations;
		omxInvokeNLOPT(rf);
		rf.finish();
		fc->wanted |= FF_COMPUTE_GRADIENT;
		fc->constraintFunVals = rf.constraintFunValsOut;
		fc->LagrMultipliers = rf.LagrMultipliersOut;
		fc->constraintJacobian = rf.constraintJacobianOut;
		fc->LagrHessian = rf.LagrHessianOut;
		break;
        case OptEngine_SD:{
		rf.setupSimpleBounds();
		if(rf.isUnconstrained()) {
			omxSD(rf);
			rf.finish();
		} else {
			mxThrow("Constrained problems are not implemented");
		}
		fc->wanted |= FF_COMPUTE_GRADIENT;
		break;}
        default: mxThrow("Optimizer %d is not available", engine);
	}

	if (isErrorRaised()) return;

	if (!fc->insideFeasibleSet()) {
		fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);
	} else {
		fc->setInform(rf.informOut);
		if (fc->ciobj) fc->ciobj->checkSolution(fc);
	}

	if (verbose >= 1) {
		mxLog("%s: done, iter=%d inform=%d",
		      name, fc->getLocalComputeCount() - beforeEval, fc->getInform());
	}

	// Optimizers can terminate with inconsistent fit and parameters
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT | FF_COMPUTE_BESTFIT, fc);

	if (verbose >= 2) {
		mxLog("%s: final fit is %2f", name, fc->getFit());
		fc->log(FF_COMPUTE_ESTIMATE);
	}

	fc->wanted |= FF_COMPUTE_BESTFIT;
}

void omxComputeGD::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	omxPopulateFitFunction(fitMatrix, out);

	MxRList output;
	SEXP pn, cv, cjac, lambdas, cstates, lagrhess;
	size_t i=0;

	output.add("maxThreads", Rf_ScalarInteger(threads));
	if( fc->varGroup->vars.size() ){
		Rf_protect(pn = Rf_allocVector( STRSXP, fc->varGroup->vars.size() ));
		for(i=0; i < fc->varGroup->vars.size(); i++){
			SET_STRING_ELT( pn, i, Rf_mkChar(fc->varGroup->vars[i]->name) );
		}
		output.add("paramNames", pn);
	}
  fc->state->reportConstraints(output);
	if( fc->constraintFunVals.size() ){
		Rf_protect(cv = Rf_allocVector( REALSXP, fc->constraintFunVals.size() ));
		memcpy( REAL(cv), fc->constraintFunVals.data(), sizeof(double) * fc->constraintFunVals.size() );
		output.add("constraintFunctionValues", cv);
	}
	if( fc->constraintJacobian.size() ){
		Rf_protect(cjac = Rf_allocMatrix( REALSXP, fc->constraintJacobian.rows(), fc->constraintJacobian.cols() ));
		memcpy( REAL(cjac), fc->constraintJacobian.data(), sizeof(double) * fc->constraintJacobian.rows() * fc->constraintJacobian.cols() );
		output.add("constraintJacobian", cjac);
	}
	if( fc->LagrMultipliers.size() ){
		Rf_protect(lambdas = Rf_allocVector( REALSXP, fc->LagrMultipliers.size() ));
		memcpy( REAL(lambdas), fc->LagrMultipliers.data(), sizeof(double) * fc->LagrMultipliers.size() );
		output.add("LagrangeMultipliers", lambdas);
	}
	if( fc->constraintStates.size() ){
		Rf_protect(cstates = Rf_allocVector( INTSXP, fc->constraintStates.size() ));
		memcpy( INTEGER(cstates), fc->constraintStates.data(), sizeof(int) * fc->constraintStates.size() );
		output.add("istate", cstates); //<--Not sure if CSOLNP and SLSQP have their own constraint state codes.
	}
	if( fc->LagrHessian.size() ){
		Rf_protect(lagrhess = Rf_allocMatrix( REALSXP, fc->LagrHessian.rows(), fc->LagrHessian.cols() ));
		memcpy( REAL(lagrhess), fc->LagrHessian.data(), sizeof(double) * fc->LagrHessian.rows() * fc->LagrHessian.cols() );
		output.add("LagrHessian", lagrhess);
	}
	slots->add("output", output.asR());

	if (engine == OptEngine_NPSOL && hessChol.size()) {
		out->add("hessianCholesky", Rcpp::wrap(hessChol));
	}
}

// -----------------------------------------------------------------------

class ComputeCI : public omxCompute {
	typedef omxCompute super;
	typedef CIobjective::Diagnostic Diagnostic;
  std::unique_ptr<omxCompute> plan;
	omxMatrix *fitMatrix;
	int verbose;
	int totalIntervals;
	SEXP intervals, intervalCodes, detail;
	const char *ctypeName;
	bool useInequality;

	enum Method {
		NEALE_MILLER_1997=1,
		WU_NEALE_2012
	};

	void regularCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int lower,
		       double &val, Diagnostic &better);
	void regularCI2(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow);
	void boundAdjCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow);
	void recordCI(Method meth, ConfidenceInterval *currentCI, int lower, FitContext &fc,
		      int &detailRow, double val, Diagnostic diag);
	void checkOtherBoxConstraints(FitContext &fc, ConfidenceInterval *currentCI, Diagnostic &diag);
	void checkBoxConstraints(FitContext &fc, int skip, Diagnostic &diag);
	void runPlan(FitContext *fc);
public:
	ComputeCI();
	virtual void initFromFrontend(omxState *, SEXP rObj) override;
	virtual void computeImpl(FitContext *fc) override;
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out) override;
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out) override;
};

omxCompute *newComputeConfidenceInterval()
{
	return new ComputeCI();
}

ComputeCI::ComputeCI()
{
	plan = 0;
	intervals = 0;
	intervalCodes = 0;
	detail = 0;
	useInequality = false;
}

void ComputeCI::runPlan(FitContext *fc)
{
	plan->compute(fc);

	// any derivs are only relevent to the current problem
	fc->wanted &= ~(FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN);
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
		// really only constrained and unconstrained TODO
		if (strEQ(ctypeName, "ineq")) {
			useInequality = true;
		} else if (strEQ(ctypeName, "eq")) {
			useInequality = true;
		} else if (strEQ(ctypeName, "both")) {
			useInequality = true;
		} else if (strEQ(ctypeName, "none")) {
			// OK
		} else {
			mxThrow("%s: unknown constraintType='%s'", name, ctypeName);
		}
	}

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	omxCompleteFitFunction(fitMatrix);

	PushLoopIndex pli(name);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = std::unique_ptr<omxCompute>(omxNewCompute(globalState, CHAR(s4class)));
	plan->initFromFrontend(globalState, slotValue);
}

extern "C" { void F77_SUB(npoptn)(char* string, int Rf_length); }

class ciConstraint : public omxConstraint {
 private:
	typedef omxConstraint super;
	omxState *state;
public:
	omxMatrix *fitMat;
	ciConstraint(omxState *u_state) : super("CI"), state(u_state) {}
  virtual void getDim(int *rowsOut, int *colsOut) const override
  { *rowsOut = size; *colsOut = 1; }
  void push()
  {
		state->conListX.push_back(this);
  }
  void pop()
  {
		size_t sz = state->conListX.size();
		if (!sz || state->conListX[sz-1] != this)
      mxThrow("Error destroying ciConstraint");
    state->conListX.pop_back();
    state = 0;
  }
	virtual omxConstraint *duplicate(omxState *dest) const override = 0;
  virtual bool hasAnalyticJac(FitContext *fc) const override
  { return fc->ciobj->hasAnalyticJac(); }
  virtual int getVerbose() const override { return 0; }
};

class ciConstraintIneq : public ciConstraint {
 private:
	typedef ciConstraint super;
 public:
	ciConstraintIneq(omxState *u_state, int u_size) : super(u_state)
	{ size=u_size; opCode = LESS_THAN; setInitialSize(size); };

	virtual void refreshAndGrab(FitContext *fc, double *out) override
  {
		fc->ciobj->evalIneq(fc, fitMat, out);
		Eigen::Map< Eigen::ArrayXd > Eout(out, size);
		//mxLog("fit %f diff %f", fit, diff);
	};
	virtual void analyticJac(FitContext *fc, MatrixStoreFn out) override
  { fc->ciobj->ineqAnalyticJac(fc, fitMat, out); }
	virtual omxConstraint *duplicate(omxState *dest) const override
  {
    auto *ptr = new ciConstraintIneq(dest, size);
    ptr->fitMat = dest->lookupDuplicate(fitMat);
    return ptr;
  }
};

class ciConstraintEq : public ciConstraint {
 private:
	typedef ciConstraint super;
 public:
	ciConstraintEq(omxState *u_state, int u_size) : super(u_state)
	{ size=u_size; opCode = EQUALITY; setInitialSize(size); };

	virtual void refreshAndGrab(FitContext *fc, double *out) override {
		fc->ciobj->evalEq(fc, fitMat, out);
		//mxLog("fit %f diff %f", fit, diff);
	};
	virtual void analyticJac(FitContext *fc, MatrixStoreFn out) override
  { fc->ciobj->eqAnalyticJac(fc, fitMat, out); }
	virtual omxConstraint *duplicate(omxState *dest) const override
  {
    auto *ptr = new ciConstraintEq(dest, size);
    ptr->fitMat = dest->lookupDuplicate(fitMat);
    return ptr;
  }
};

void ComputeCI::recordCI(Method meth, ConfidenceInterval *currentCI, int lower, FitContext &fc,
			 int &detailRow, double val, Diagnostic diag)
{
	omxMatrix *ciMatrix = currentCI->getMatrix(fitMatrix->currentState);
	std::string &matName = ciMatrix->nameStr;

	if (diag == CIobjective::DIAG_SUCCESS) {
		currentCI->val[!lower] = val;
		currentCI->code[!lower] = fc.getInform();
	}

	if(verbose >= 1) {
		mxLog("CI[%s,%s] %s[%d,%d] val=%f fit=%f status=%d accepted=%d",
		      currentCI->name.c_str(), (lower?"lower":"upper"), matName.c_str(),
		      1+currentCI->row, 1+currentCI->col, val, fc.getFit(), fc.getInform(), diag);
	}

	SET_STRING_ELT(VECTOR_ELT(detail, 0), detailRow, Rf_mkChar(currentCI->name.c_str()));
	INTEGER(VECTOR_ELT(detail, 1))[detailRow] = 1+lower;
	REAL(VECTOR_ELT(detail, 2))[detailRow] = val;
	REAL(VECTOR_ELT(detail, 3))[detailRow] = fc.getFit();
	INTEGER(VECTOR_ELT(detail, 4))[detailRow] = diag;
	INTEGER(VECTOR_ELT(detail, 5))[detailRow] = fc.wrapInform();
	INTEGER(VECTOR_ELT(detail, 6))[detailRow] = meth;
	for (int px=0; px < int(fc.numParam); ++px) {
		REAL(VECTOR_ELT(detail, 7+px))[detailRow] = fc.est[px];
	}
	++detailRow;
}

struct regularCIobj : CIobjective {
	double targetFit;
	double diff;

  regularCIobj(const ConfidenceInterval *u_CI,
               bool u_constrained, bool u_lowerBound, double u_targetFit) :
    CIobjective(u_CI, u_constrained, u_lowerBound), targetFit(u_targetFit) {}

  virtual std::unique_ptr<CIobjective> clone() const override
  { return std::make_unique<regularCIobj>(*this); }

	void computeConstraint(double fit)
	{
		diff = fit - targetFit;
	}

	virtual void evalIneq(FitContext *fc, omxMatrix *fitMat, double *out) override
	{
    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		const double fit = fc->getFit();
		computeConstraint(fit);
		*out = std::max(diff, 0.0);
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc) override
	{
    omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			mxThrow("Not implemented yet");
		}
		if (want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
      mxThrow("Not implemeneted");
		}

    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });

		const double fit = fc->getFit();

		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);

		if (!std::isfinite(fit)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fc->setFit(nan("infeasible"));
			return;
		}

		if (want & FF_COMPUTE_FIT) {
			double param = (lowerBound? CIElement : -CIElement);
			computeConstraint(fit);
			if (fabs(diff) > 100) param = nan("infeasible");
			if (constrained) {
				fc->setFit(diff*diff + param);
			} else {
				fc->setFit(param);
			}
		}
    if (want & FF_COMPUTE_GRADIENT) setGrad(fc);
	}

	virtual Diagnostic getDiag() override
	{
		Diagnostic diag = DIAG_SUCCESS;
		if (fabs(diff) > 1e-1) {
			diag = DIAG_ALPHA_LEVEL;
		}
		return diag;
	}

  virtual bool hasAnalyticJac() override { return !Global->NPSOL_HACK; }

	virtual void ineqAnalyticJac(FitContext *fc, omxMatrix *fitMat, MatrixStoreFn out) override
  {
    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_GRADIENT, fc); });

    double scale = fc->getFitScale();
    auto &gr = fc->gradZ;
    for (int gx=0; gx < gr.size(); ++gx) out(0, gx, gr[gx] * scale);
  }

	virtual void eqAnalyticJac(FitContext *fc, omxMatrix *fitMat, MatrixStoreFn out) override
  { mxThrow("Not implemented"); }
};

struct bound1CIobj : CIobjective {
  double bound;
	Eigen::Array<double,1,1> eq;

  bound1CIobj(const ConfidenceInterval *u_CI,
              double u_bound, bool u_constrained) :
    CIobjective(u_CI, u_constrained, true), bound(u_bound) {}

  virtual std::unique_ptr<CIobjective> clone() const override
  { return std::make_unique<bound1CIobj>(*this); }

	template <typename T1>
	void computeConstraint(FitContext *fc, omxMatrix *fitMat, double fit, Eigen::ArrayBase<T1> &v1)
	{
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);
		v1 << CIElement - bound;
		eq = v1;
		//mxPrintMat("v1", v1);
	}

	virtual void evalEq(FitContext *fc, omxMatrix *fitMat, double *out) override
	{
    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		const double fit = fc->getFit();
		Eigen::Array<double,1,1> v1;
		computeConstraint(fc, fitMat, fit, v1);
		*out = v1(0);
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc) override
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			mxThrow("Not implemented yet");
		}

    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		double fit = fc->getFit();

		if (!std::isfinite(fit)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fc->setFit(nan("infeasible"));
			return;
		}

		double cval = 0.0;
		Eigen::Array<double,1,1> v1;
		computeConstraint(fc, fitMat, fit, v1);
		if (fabs(v1(0)) > 100) fit = nan("infeasible");
		if (!constrained) {
			cval = v1(0) * v1(0);
		}

		fc->setFit(fit + cval);

    if (want & FF_COMPUTE_GRADIENT) setGrad(fc);
	}

	virtual Diagnostic getDiag() override
	{
		Diagnostic diag = DIAG_SUCCESS;
		if (fabs(eq(0)) > 1e-3) {
			diag = DIAG_BOUND_INFEASIBLE;
		}
		return diag;
	}
};

struct boundAwayCIobj : CIobjective {
	double logAlpha, sqrtCrit;
	double unboundedLL, bestLL;
	Eigen::Array<double, 3, 1> ineq;

  boundAwayCIobj(const ConfidenceInterval *u_CI,
                 double u_la, double u_sq,
                 double u_unboundedLL, double u_bestLL,
                 int u_lower, bool u_constrained) :
    CIobjective(u_CI, u_constrained, u_lower), logAlpha(u_la), sqrtCrit(u_sq),
    unboundedLL(u_unboundedLL), bestLL(u_bestLL) {}

  virtual std::unique_ptr<CIobjective> clone() const override
  { return std::make_unique<boundAwayCIobj>(*this); }

	template <typename T1>
	void computeConstraint(double fit, Eigen::ArrayBase<T1> &v1)
	{
		double d1 = sqrt(std::max(fit - bestLL, 0.0));
		double d2 = sqrt(std::max(fit - unboundedLL,0.0));
		//Rf_pnorm(d1,0,1,0,0) + Rf_pnorm(d2,0,1,0,0) > alpha
		//d1 < sqrtCrit
		//d2 > sqrtCrit
		double pA = Rf_pnorm5(d1,0,1,0,0) + Rf_pnorm5(d2,0,1,0,0);
		v1 << std::max(d1 - sqrtCrit, 0.0), std::max(sqrtCrit - d2, 0.0),
			std::max(logAlpha - log(pA), 0.0);
		ineq = v1;
		//mxPrintMat("v1", v1);
		//mxLog("fit %g sqrtCrit %g d1 %g d2 %g alpha %g pA %g",
		//fit, sqrtCrit, d1, d2, exp(logAlpha), pA);
	}

	virtual void evalIneq(FitContext *fc, omxMatrix *fitMat, double *out) override
	{
    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		double fit = fc->getFit();
		Eigen::Map< Eigen::Array<double,3,1> >v1(out);
		computeConstraint(fit, v1);
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc) override
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			mxThrow("Not implemented yet");
		}

    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		double fit = fc->getFit();
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);

		if (!std::isfinite(fit)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fc->setFit(nan("infeasible"));
			return;
		}

		double param = (lowerBound? CIElement : -CIElement);
		double cval = 0.0;
		Eigen::Array<double,3,1> v1;
		computeConstraint(fit, v1);
		if ((v1 > 10).any()) param = nan("infeasible");
		if (!constrained) {
			cval = v1.sum()*v1.sum();
		}

		fc->setFit(param + cval);
		//mxLog("param at %f", fitMat->data[0]);

    if (want & FF_COMPUTE_GRADIENT) setGrad(fc);
	}

	virtual Diagnostic getDiag() override
	{
		Diagnostic diag = DIAG_SUCCESS;
		if (ineq[0] > 1e-3) {
			diag = DIAG_BA_D1;
		} else if (ineq[1] > 1e-3) {
			diag = DIAG_BA_D2;
		} else if (ineq[2] > 1e-1) {
			diag = DIAG_ALPHA_LEVEL;
		}
		return diag;
	}
};

struct boundNearCIobj : CIobjective {
	double d0, logAlpha;
	double boundLL, bestLL;
	double pN;
	double lbd, ubd;
	Eigen::Array<double,3,1> ineq;

  boundNearCIobj(const ConfidenceInterval *u_CI,
                 double u_d0, double u_logAlpha,
                 double u_boundLL, double u_bestLL,
                 int u_lower, bool u_constrained,
                 double u_lbd, double u_ubd)
    : CIobjective(u_CI, u_constrained, u_lower),
      d0(u_d0), logAlpha(u_logAlpha), boundLL(u_boundLL), bestLL(u_bestLL),
      lbd(u_lbd), ubd(u_ubd) {}

  virtual std::unique_ptr<CIobjective> clone() const override
  { return std::make_unique<boundNearCIobj>(*this); }

	template <typename T1>
	void computeConstraint(double fit, Eigen::ArrayBase<T1> &v1)
	{
		double dd = sqrt(std::max(fit - bestLL, 0.0));
		double pN1 = Rf_pnorm5(dd,0,1,0,0);
		double pN2 = Rf_pnorm5((d0-dd)/2.0+dd*dd/(2*std::max(d0-dd,.001*dd*dd)),0,1,0,0);
		pN = pN1 + pN2;
		// pN > alpha
		// lbd < dd < ubd
		v1 << std::max(lbd - dd, 0.0), std::max(dd - ubd, 0.0), std::max(logAlpha - log(pN), 0.0);
		ineq = v1;
		// mxLog("pN1 %.6g pN2 %.6g d0-dd %.5g dd*dd %.5g expr %.5g",
		//       pN1, pN2, d0-dd, dd*dd, dd*dd/(2*std::max(d0-dd,.001*dd*dd)));
		// mxPrintMat("v1", v1);
		// mxLog("fit %g dd %g/%g/%g alpha %g pN %.6g",
		//       fit, lbd, dd, ubd, alpha, pN);
	}

	virtual void evalIneq(FitContext *fc, omxMatrix *fitMat, double *out) override
	{
    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		double fit = fc->getFit();
		Eigen::Map< Eigen::Array<double,3,1> > v1(out);
		computeConstraint(fit, v1);
	}

	virtual void evalFit(omxFitFunction *ff, int want, FitContext *fc) override
	{
		omxMatrix *fitMat = ff->matrix;

		if (!(want & FF_COMPUTE_FIT)) {
			if (want & (FF_COMPUTE_PREOPTIMIZE | FF_COMPUTE_INITIAL_FIT)) return;
			mxThrow("Not implemented yet");
		}

    fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
		double fit = fc->getFit();
		omxMatrix *ciMatrix = CI->getMatrix(fitMat->currentState);
		omxRecompute(ciMatrix, fc);
		double CIElement = omxMatrixElement(ciMatrix, CI->row, CI->col);

		if (!std::isfinite(fit) || !std::isfinite(CIElement)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			fc->setFit(nan("infeasible"));
			return;
		}

		double param = (lowerBound? CIElement : -CIElement);
		double cval = 0.0;
		Eigen::Array<double,3,1> v1;
		computeConstraint(fit, v1);
		if ((v1 > 10).any()) param = nan("infeasible");
		if (!constrained) {
			cval = v1.sum() * v1.sum();
		}

		fc->setFit(param  + cval);

    if (want & FF_COMPUTE_GRADIENT) setGrad(fc);
	}

	virtual Diagnostic getDiag() override
	{
		Diagnostic diag = DIAG_SUCCESS;
		if (ineq[0] > 1e-3) {
			diag = DIAG_BN_D1;
		} else if (ineq[1] > 1e-2) {
			diag = DIAG_BN_D2;
		} else if (fabs(pN - exp(logAlpha)) > 1e-3) {
			diag = DIAG_ALPHA_LEVEL;
		}
		return diag;
	}
};

void ComputeCI::boundAdjCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow)
{
	omxState *state = fitMatrix->currentState;
	omxMatrix *ciMatrix = currentCI->getMatrix(state);
	std::string &matName = ciMatrix->nameStr;
	omxFreeVar *fv = mle->varGroup->vars[currentCI->varIndex];
	double &nearBox = fv->lbound == NEG_INF? fv->ubound : fv->lbound;
	double &farBox = fv->lbound == NEG_INF? fv->lbound : fv->ubound;
	int side = fv->lbound == NEG_INF? ConfidenceInterval::Upper : ConfidenceInterval::Lower;

	Eigen::VectorXd Mle = mle->est;
	fc.est = Mle;

	bool boundActive = fabs(Mle[currentCI->varIndex] - nearBox) < sqrt(std::numeric_limits<double>::epsilon());
	if (currentCI->bound[!side] > 0.0) {	// ------------------------------ away from bound side --
		PushLoopIndex pli(name, detailRow, totalIntervals);
		if (!boundActive) {
			Diagnostic diag;
			double val;
			regularCI(mle, fc, currentCI, side, val, diag);
			recordCI(NEALE_MILLER_1997, currentCI, side, fc, detailRow, val, diag);
			goto part2;
		}

		Global->checkpointMessage(mle, "%s[%d, %d] unbounded fit",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
		{
			double boxSave = nearBox;
			nearBox = NA_REAL;
			fc.est = Mle;
			runPlan(&fc);
			nearBox = boxSave;
			if (verbose >= 2) {
				omxRecompute(ciMatrix, &fc);
				double val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);
				mxLog("without bound, fit %.8g val %.8g", fc.getFit(), val);
			}
		}

		double unboundedLL = fc.getFit();
		//mxLog("unboundedLL %g bestLL %g", unboundedLL, mle->fit);

		Global->checkpointMessage(mle, "%s[%d, %d] adjusted away CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
		bool useConstr = useInequality;
		ciConstraintIneq constr(state, 3);
		constr.fitMat = fitMatrix;
		if (useConstr) {
      constr.push();
      fc.calcNumFree();
    }
		// Could set farBox to the regular UB95, but we'd need to optimize to get it
		fc.est = Mle;
		fc.ciobj = std::make_unique<boundAwayCIobj>
      (currentCI,
       log(1.0 - Rf_pchisq(currentCI->bound[!side], 1, 1, 0)),
       sqrt(currentCI->bound[!side]), unboundedLL, mle->getFit(), side, useConstr);
		runPlan(&fc);
		if (useConstr) constr.pop();

		omxRecompute(ciMatrix, &fc);
		double val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

		ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);
		Diagnostic diag = fc.ciobj->getDiag();
		fc.ciobj.reset();
		checkOtherBoxConstraints(fc, currentCI, diag);
		recordCI(WU_NEALE_2012, currentCI, side, fc, detailRow, val, diag);
	}

 part2:
	if (currentCI->bound[side] > 0.0) {     // ------------------------------ near to bound side --
		PushLoopIndex pli(name, detailRow, totalIntervals);
		double boundLL = NA_REAL;
		double sqrtCrit95 = sqrt(currentCI->bound[side]);
		if (!boundActive) {
			Global->checkpointMessage(mle, "%s[%d, %d] at-bound CI",
						  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
			fc.est = Mle;
			fc.est[currentCI->varIndex] = nearBox; // might be infeasible
			fc.profiledOutZ[currentCI->varIndex] = true;
      fc.calcNumFree();
			runPlan(&fc);
			fc.profiledOutZ[currentCI->varIndex] = false;
      fc.calcNumFree();
			if (fc.getInform() == 0) {
				boundLL = fc.getFit();
			}

			if (verbose >= 2) {
				mxLog("%s[%d, %d]=%.2f (at bound) fit %.2f status %d",
				      matName.c_str(), currentCI->row + 1, currentCI->col + 1,
				      nearBox, fc.getFit(), fc.getInform());
			}
			if (!std::isfinite(boundLL)) {
				// Might work if simple approach failed
				ciConstraintEq constr(state, 1);
				constr.fitMat = fitMatrix;
				constr.push();
        fc.calcNumFree();
				fc.ciobj = std::make_unique<bound1CIobj>
          (currentCI, nearBox, useInequality);
				fc.est = Mle;
				runPlan(&fc);
				constr.pop();
				boundLL = fc.getFit();
				Diagnostic diag = fc.ciobj->getDiag();
				if (diag != CIobjective::DIAG_SUCCESS) {
					recordCI(WU_NEALE_2012, currentCI, !side, fc, detailRow,
						 NA_REAL, diag);
					return;
				}
				//mxLog("got2 %.4g", fc.fit);
			}
			//mxPrintMat("bound1", ciobj.ineq);
			//mxLog("mle %g boundLL %g", mle->fit, boundLL);
		} else {
			boundLL = mle->getFit();
		}

		Global->checkpointMessage(mle, "%s[%d, %d] near side CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1);
		double sqrtCrit90 = sqrt(Rf_qchisq(1-(2.0 * (1-Rf_pchisq(currentCI->bound[side], 1, 1, 0))),1,1,0));
		double d0 = sqrt(std::max(boundLL - mle->getFit(), 0.0));
		//mxLog("d0 %.4g sqrtCrit90 %.4g d0/2 %.4g", d0, sqrtCrit90, d0/2.0);
		if (d0 < sqrtCrit90) {
			fc.setFit(boundLL);
			Diagnostic diag = CIobjective::DIAG_SUCCESS;
			checkOtherBoxConstraints(fc, currentCI, diag);
			recordCI(WU_NEALE_2012, currentCI, !side, fc, detailRow, nearBox, diag);
			return;
		}

		if (sqrtCrit95 <= d0/2.0) {
			Diagnostic better;
			double val;
			regularCI(mle, fc, currentCI, !side, val, better);
			recordCI(NEALE_MILLER_1997, currentCI, !side, fc, detailRow, val, better);
			return;
		}

		double alphalevel = 1.0 - Rf_pchisq(currentCI->bound[side], 1, 1, 0);
		bool useConstr = useInequality;
		ciConstraintIneq constr(state, 3);
		constr.fitMat = fitMatrix;
		if (useConstr) {
      constr.push();
      fc.calcNumFree();
    }
		double boxSave = farBox;
		farBox = Mle[currentCI->varIndex];
		fc.est = Mle;
		// Perspective helps? Optimizer seems to like to start further away
		fc.est[currentCI->varIndex] = (9*Mle[currentCI->varIndex] + nearBox) / 10.0;
		fc.ciobj = std::make_unique<boundNearCIobj>
      (currentCI, d0, log(alphalevel), boundLL, mle->getFit(), !side,
       useConstr, std::max(d0/2, sqrtCrit90), std::min(d0, sqrtCrit95));
		runPlan(&fc);
		farBox = boxSave;
		if (useConstr) constr.pop();

		omxRecompute(ciMatrix, &fc);
		double val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

		//mxLog("val=%g", val);
		//mxPrintMat("bn", bnobj.ineq);
		ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);
		Diagnostic diag = fc.ciobj->getDiag();
		fc.ciobj.reset();
		checkOtherBoxConstraints(fc, currentCI, diag);
		recordCI(WU_NEALE_2012, currentCI, !side, fc, detailRow, val, diag);
	}
}

void ComputeCI::checkOtherBoxConstraints(FitContext &fc, ConfidenceInterval *currentCI,
					 Diagnostic &diag)
{
	checkBoxConstraints(fc, currentCI->varIndex, diag);
}

void ComputeCI::checkBoxConstraints(FitContext &fc, int skip, Diagnostic &diag)
{
	if (diag != CIobjective::DIAG_SUCCESS) return;
  if (fc.hasActiveBoxConstraint(skip)) diag = CIobjective::DIAG_BOXED;
}

void ComputeCI::regularCI(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int lower,
			  double &val, Diagnostic &diag)
{
	omxState *state = fitMatrix->currentState;

	ciConstraintIneq constr(state, 1);
	bool constrained = useInequality;
	if (constrained) {
		constr.fitMat = fitMatrix;
		constr.push();
    fc.calcNumFree();
	}

	fc.est = mle->est;
	fc.ciobj = std::make_unique<regularCIobj>
    (currentCI, !constrained, lower, currentCI->bound[!lower] + mle->getFit());
	//mxLog("Set target fit to %f (MLE %f)", fc->targetFit, fc->fit);

	runPlan(&fc);
	if (constrained) constr.pop();

	omxMatrix *ciMatrix = currentCI->getMatrix(fitMatrix->currentState);
	omxRecompute(ciMatrix, &fc);
	val = omxMatrixElement(ciMatrix, currentCI->row, currentCI->col);

	diag = fc.ciobj->getDiag();
	fc.ciobj.reset();

	// We check the fit again so we can report it
	// in the detail data.frame.
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, &fc);

	checkBoxConstraints(fc, -1, diag);
}

void ComputeCI::regularCI2(FitContext *mle, FitContext &fc, ConfidenceInterval *currentCI, int &detailRow)
{
	omxState *state = fitMatrix->currentState;
	omxMatrix *ciMatrix = currentCI->getMatrix(state);
	std::string &matName = ciMatrix->nameStr;

	for (int upper=0; upper <= 1; ++upper) {
		int lower = 1-upper;
		if (!(currentCI->bound[upper])) continue;

		PushLoopIndex pli(name, detailRow, totalIntervals);
		Global->checkpointMessage(mle, "%s[%d, %d] %s CI",
					  matName.c_str(), currentCI->row + 1, currentCI->col + 1,
					  upper? "upper" : "lower");
		double val;
		Diagnostic better;
		regularCI(mle, fc, currentCI, lower, val, better);
		recordCI(NEALE_MILLER_1997, currentCI, lower, fc, detailRow, val, better);
	}
}

void ComputeCI::computeImpl(FitContext *mle)
{
	if (intervals) mxThrow("Can only compute CIs once");
	if (!Global->intervals) {
		if (verbose >= 1) mxLog(name, "%s: mxRun(..., intervals=FALSE), skipping");
		return;
	}

	omxState *state = fitMatrix->currentState;
	Global->unpackConfidenceIntervals(state);

	// Not strictly necessary, but makes it easier to run
	// mxComputeConfidenceInterval alone without other compute
	// steps.
	omxAlgebraPreeval(fitMatrix, mle);
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, mle);
  auto *ff = fitMatrix->fitFunction;
  if (!fitUnitsIsChiSq(ff->units)) {
    mxThrow("%s: units '%s' are not supported", name, fitUnitsToName(ff->units));
  }

	int numInts = (int) Global->intervalList.size();
	if (verbose >= 1) mxLog("%s: %d intervals of '%s' (ref fit %f %s)",
                          name, numInts, fitMatrix->name(), mle->getFit(), ctypeName);
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

	totalIntervals = 0;
	{
		Eigen::Map< Eigen::ArrayXXd > interval(REAL(intervals), numInts, 3);
		interval.fill(NA_REAL);
		for(int j = 0; j < numInts; j++) {
			ConfidenceInterval *oCI = Global->intervalList[j];
			omxMatrix *ciMat = oCI->getMatrix(state);
			omxRecompute(ciMat, mle);
			if (!ciMat->isValidElem(oCI->row, oCI->col)) {
				mxThrow("%s: attempt to find confidence interval of "
					 "nonexistent element (%d,%d) in %dx%d matrix '%s'",
					 name, 1+oCI->row, 1+oCI->col, ciMat->rows, ciMat->cols, ciMat->name());
			}
			interval(j, 1) = omxMatrixElement(ciMat, oCI->row, oCI->col);
			totalIntervals += (oCI->bound != 0.0).count();
		}
	}

	int numDetailCols = 7 + mle->numParam;
	Rf_protect(detail = Rf_allocVector(VECSXP, numDetailCols));
	SET_VECTOR_ELT(detail, 0, Rf_allocVector(STRSXP, totalIntervals));
	const char *sideLabels[] = { "upper", "lower" };
	SET_VECTOR_ELT(detail, 1, makeFactor(Rf_allocVector(INTSXP, totalIntervals), 2, sideLabels));
	SET_VECTOR_ELT(detail, 2, Rf_allocVector(REALSXP, totalIntervals));
	SET_VECTOR_ELT(detail, 3, Rf_allocVector(REALSXP, totalIntervals));
	const char *diagLabels[] = {
		"success", "alpha level not reached",
		"bound-away mle distance", "bound-away unbounded distance",
		"bound-near lower distance", "bound-near upper distance",
		"bound infeasible", "active box constraint"
	};
	SET_VECTOR_ELT(detail, 4,
		       makeFactor(Rf_allocVector(INTSXP, totalIntervals),
				  OMX_STATIC_ARRAY_SIZE(diagLabels), diagLabels));
	SET_VECTOR_ELT(detail, 5, allocInformVector(totalIntervals));
	const char *methodLabels[] = { "neale-miller-1997", "wu-neale-2012" };
	SET_VECTOR_ELT(detail, 6,
		       makeFactor(Rf_allocVector(INTSXP, totalIntervals),
				  OMX_STATIC_ARRAY_SIZE(methodLabels), methodLabels));
	for (int cx=0; cx < int(mle->numParam); ++cx) {
		SET_VECTOR_ELT(detail, 7+cx, Rf_allocVector(REALSXP, totalIntervals));
	}

	SEXP detailCols;
	Rf_protect(detailCols = Rf_allocVector(STRSXP, numDetailCols));
	Rf_setAttrib(detail, R_NamesSymbol, detailCols);
	SET_STRING_ELT(detailCols, 0, Rf_mkChar("parameter"));
	SET_STRING_ELT(detailCols, 1, Rf_mkChar("side"));
	SET_STRING_ELT(detailCols, 2, Rf_mkChar("value"));
	SET_STRING_ELT(detailCols, 3, Rf_mkChar("fit"));
	SET_STRING_ELT(detailCols, 4, Rf_mkChar("diagnostic"));
	SET_STRING_ELT(detailCols, 5, Rf_mkChar("statusCode"));
	SET_STRING_ELT(detailCols, 6, Rf_mkChar("method"));
	for (int nx=0; nx < int(mle->numParam); ++nx) {
		SET_STRING_ELT(detailCols, 7+nx, Rf_mkChar(mle->varGroup->vars[nx]->name));
	}

	markAsDataFrame(detail, totalIntervals);

	FitContext fc(mle, mle->varGroup);
  fc.calcNumFree();
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
	int* intervalCode = INTEGER(intervalCodes);
	for(int j = 0; j < numInts; j++) {
		ConfidenceInterval *oCI = Global->intervalList[j];
		interval(j, 0) = std::min(oCI->val[ConfidenceInterval::Lower], interval(j, 1));
		interval(j, 2) = std::max(oCI->val[ConfidenceInterval::Upper], interval(j, 1));
		intervalCode[j] = oCI->code[ConfidenceInterval::Lower];
		intervalCode[j + numInts] = oCI->code[ConfidenceInterval::Upper];
	}
}

void ComputeCI::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);
	std::vector< omxCompute* > clist(1);
	clist[0] = plan.get();
	collectResultsHelper(fc, clist, lcr, out);
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

// ---------------------------------------------------------------

class ComputeTryH : public omxCompute {
	typedef omxCompute super;
	std::unique_ptr< omxCompute > plan;
	int numFree;
	int verbose;
	double loc;
	double scale;
	int maxRetries;
	int invocations;
	int numRetries;
	Eigen::ArrayXd bestEst;
	int bestStatus;
	double bestFit;
	Eigen::VectorXd solLB;
	Eigen::VectorXd solUB;

	static bool satisfied(FitContext *fc);
public:
	virtual void initFromFrontend(omxState *, SEXP rObj) override;
	virtual void computeImpl(FitContext *fc) override;
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out) override;
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out) override;
};

omxCompute *newComputeTryHard()
{ return new ComputeTryH(); }

void ComputeTryH::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	{
		ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(Rverbose);

		ProtectedSEXP Rloc(R_do_slot(rObj, Rf_install("location")));
		loc = Rf_asReal(Rloc);
		ProtectedSEXP Rscale(R_do_slot(rObj, Rf_install("scale")));
		scale = Rf_asReal(Rscale);
		ProtectedSEXP Rretries(R_do_slot(rObj, Rf_install("maxRetries")));
		maxRetries = Rf_asReal(Rretries);
	}

	invocations = 0;
	numRetries = 0;

	PushLoopIndex pli(name);

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = std::unique_ptr< omxCompute >(omxNewCompute(globalState, CHAR(s4class)));
	plan->initFromFrontend(globalState, slotValue);
}

bool ComputeTryH::satisfied(FitContext *fc)
{
	return (fc->getInform() == INFORM_CONVERGED_OPTIMUM ||
		fc->getInform() == INFORM_UNCONVERGED_OPTIMUM ||
		fc->getInform() == INFORM_ITERATION_LIMIT);
}

void ComputeTryH::computeImpl(FitContext *fc)
{
	using Eigen::Map;
	using Eigen::ArrayXd;
  fc->calcNumFree();
	numFree = fc->getNumFree();
	ArrayXd origStart(numFree);
  ArrayXd attempt(numFree);
  fc->copyEstToOptimizer(origStart);
	bestEst = origStart;

	solLB.resize(numFree);
	solUB.resize(numFree);
	fc->copyBoxConstraintToOptimizer(solLB, solUB);

	++invocations;

	int retriesRemain = maxRetries - 1;
	if (verbose >= 1) {
		mxLog("%s: at most %d attempts (Welcome)", name, retriesRemain);
	}

	{
		PushLoopIndex pli(name, 1, 0);
		bestStatus = INFORM_UNINITIALIZED;
		bestFit = NA_REAL;
		fc->setInform(INFORM_UNINITIALIZED);
		plan->compute(fc);
		if (fc->getInform() != INFORM_UNINITIALIZED &&
		    fc->getInform() != INFORM_STARTING_VALUES_INFEASIBLE) {
			bestStatus = fc->getInform();
			fc->copyEstToOptimizer(bestEst);
			bestFit = fc->getFit();
		}
	}

	while (!satisfied(fc) && retriesRemain > 0) {
		if (verbose >= 2) {
			mxLog("%s: fit %.2f inform %d, %d retries remain", name, fc->getFit(),
			      fc->getInform(), retriesRemain);
		}

    attempt = origStart;
		{
			BorrowRNGState grs;
			for (int vx=0; vx < numFree; ++vx) {
				double adj1 = loc + unif_rand() * 2.0 * scale - scale;
				double adj2 = 0.0 + unif_rand() * 2.0 * scale - scale;
				if (verbose >= 3) {
					mxLog("%d %g %g", vx, adj1, adj2);
				}
				attempt[vx] = attempt[vx] * adj1 + adj2;
				if(attempt[vx] < solLB[vx]){attempt[vx] = solLB[vx];}
				if(attempt[vx] > solUB[vx]){attempt[vx] = solUB[vx];}
			}
      fc->setEstFromOptimizer(attempt);
		}

		--retriesRemain;
		PushLoopIndex pli(name, maxRetries - retriesRemain, 0);

		fc->setInform(INFORM_UNINITIALIZED);
		fc->wanted &= ~(FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN);
		plan->compute(fc);
		if (fc->getInform() != INFORM_UNINITIALIZED && fc->getInform() != INFORM_STARTING_VALUES_INFEASIBLE &&
		    (bestStatus == INFORM_UNINITIALIZED || fc->getInform() < bestStatus)) {
			bestStatus = fc->getInform();
			fc->copyEstToOptimizer(bestEst);
			bestFit = fc->getFit();
		}
	}

	fc->setInform(bestStatus);
  fc->setEstFromOptimizer(bestEst);
	fc->setFit(bestFit);

	numRetries += maxRetries - retriesRemain;
	if (verbose >= 1) {
		mxLog("%s: fit %.2f inform %d after %d attempt(s)", name, fc->getFit(), fc->getInform(),
		      maxRetries - retriesRemain);
	}
}

void ComputeTryH::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);
	std::vector< omxCompute* > clist(1);
	clist[0] = plan.get();
	collectResultsHelper(fc, clist, lcr, out);
}

void ComputeTryH::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	MxRList info;
	info.add("invocations", Rf_ScalarInteger(invocations));
	info.add("retries", Rf_ScalarInteger(numRetries));
	slots->add("debug", info.asR());
}

// ---------------------------------------------------------------

class ComputeGenSA : public omxCompute {
	typedef omxCompute super;
  std::unique_ptr<omxCompute> plan;
	static const char *optName;
	const char *methodName;
	std::string contextStr;
	int numFree;
	Eigen::ArrayXd equality;
	Eigen::ArrayXd inequality;
	omxMatrix *fitMatrix;
	int verbose;
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;
	Eigen::VectorXd range;
	Eigen::VectorXd xMini;
	Eigen::VectorXd curBest;
	double curBestFit;
	double curBestPenalty;
	double qv;
	double qaInit;
	double lambda;
	double temSta;
	double temEnd;
	int stepsPerTemp;
	double visita(double temp);
  std::unique_ptr<ConstraintVec> cvec;
	double getConstraintPenalty(FitContext *fc);
	std::string asaOut;
	Eigen::VectorXd quenchParamScale;
	Eigen::VectorXd quenchCostScale;
	USER_DEFINES *OPTIONS;
	enum algo {
		ALGO_TSALLIS1996,
		ALGO_INGBER2012,
	} method;
	FitContext *u_fc;
	void tsallis1996(FitContext *fc);
	void ingber2012(FitContext *fc);

 public:
	virtual ~ComputeGenSA();
        virtual void initFromFrontend(omxState *, SEXP rObj) override;
        virtual void computeImpl(FitContext *fc) override;
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out) override;

	double asa_cost(double *x, int *cost_flag, int *exit_code, USER_DEFINES *opt);
};

const char *ComputeGenSA::optName = "GenSA";

class omxCompute *newComputeGenSA()
{ return new ComputeGenSA(); }

ComputeGenSA::~ComputeGenSA()
{
	delete OPTIONS;
}

void ComputeGenSA::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	// ingber2012 defaults
	OPTIONS = new USER_DEFINES;
	OMXZERO(OPTIONS, 1);
	OPTIONS->Limit_Acceptances = 10000;
	OPTIONS->Limit_Generated = 99999;
	OPTIONS->Limit_Invalid_Generated_States = 1000;
	OPTIONS->Accepted_To_Generated_Ratio = 1.0E-6;
	OPTIONS->Maximum_Cost_Repeat = 5;
	OPTIONS->Number_Cost_Samples = 5;
	OPTIONS->Temperature_Ratio_Scale = 1.0E-5;
	OPTIONS->Temperature_Anneal_Scale = 100.;
	OPTIONS->Cost_Parameter_Scale_Ratio = 1;
	OPTIONS->Sequential_Parameters = -1;
	OPTIONS->Initial_Parameter_Temperature = 1.0;
	OPTIONS->Acceptance_Frequency_Modulus = 100;
	OPTIONS->Generated_Frequency_Modulus = 10000;
	OPTIONS->Reanneal_Cost = 1;
	OPTIONS->Reanneal_Parameters = TRUE;
	asaOut = "/dev/null";
	OPTIONS->Asa_Out_File = asaOut.c_str();

	// tsallis1996 defaults
	qv=2.62;
	qaInit=-3.;
	lambda=0.85;
	temSta=5230.;
	temEnd=0.1;
	stepsPerTemp=1;
	{
		fitMatrix = omxNewMatrixFromSlot(rObj, state, "fitfunction");
		omxCompleteFitFunction(fitMatrix);

		ProtectedSEXP Rdgss(R_do_slot(rObj, Rf_install("defaultGradientStepSize")));
		OPTIONS->Delta_X = Rf_asReal(Rdgss);

		ProtectedSEXP Rdfp(R_do_slot(rObj, Rf_install("defaultFunctionPrecision")));
		OPTIONS->Cost_Precision = Rf_asReal(Rdfp);

		ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(Rverbose);

		ProtectedSEXP Rmethod(R_do_slot(rObj, Rf_install("method")));
		methodName = R_CHAR(STRING_ELT(Rmethod, 0));
		if (strEQ(methodName, "tsallis1996")) method = ALGO_TSALLIS1996;
		else if (strEQ(methodName, "ingber2012")) method = ALGO_INGBER2012;
		else mxThrow("%s: unknown method '%s'", name, methodName);
		contextStr = string_snprintf("%s(%s)", name, methodName);

		ProtectedSEXP Rcontrol(R_do_slot(rObj, Rf_install("control")));
		ProtectedSEXP RcontrolName(Rf_getAttrib(Rcontrol, R_NamesSymbol));
		for (int ax=0; ax < Rf_length(Rcontrol); ++ax) {
			const char *key = R_CHAR(STRING_ELT(RcontrolName, ax));
			ProtectedSEXP Rval(VECTOR_ELT(Rcontrol, ax));
			if (method == ALGO_TSALLIS1996) {
				if (strEQ(key, "qv")) {
					qv = Rf_asReal(Rval);
				} else if (strEQ(key, "qaInit")) {
					qaInit = Rf_asReal(Rval);
				} else if (strEQ(key, "lambda")) {
					lambda = Rf_asReal(Rval);
				} else if (strEQ(key, "tempStart")) {
					temSta = Rf_asReal(Rval);
				} else if (strEQ(key, "tempEnd")) {
					temEnd = Rf_asReal(Rval);
				} else if (strEQ(key, "stepsPerTemp")) {
					stepsPerTemp = Rf_asInteger(Rval);
				} else {
					Rf_warning("%s: unknown key '%s' for method '%s'",
						   name, key, methodName);
				}
			} else if (method == ALGO_INGBER2012) {
				if (strEQ(key, "Limit_Acceptances")) {
					OPTIONS->Limit_Acceptances = Rf_asInteger(Rval);
				} else if (strEQ(key, "Limit_Generated")) {
					OPTIONS->Limit_Generated = Rf_asInteger(Rval);
				} else if (strEQ(key, "Limit_Invalid_Generated_States")) {
					OPTIONS->Limit_Invalid_Generated_States = Rf_asInteger(Rval);
				} else if (strEQ(key, "Accepted_To_Generated_Ratio")) {
					OPTIONS->Accepted_To_Generated_Ratio = Rf_asReal(Rval);
				} else if (strEQ(key, "Cost_Precision")) {
					OPTIONS->Cost_Precision = Rf_asReal(Rval);
				} else if (strEQ(key, "Maximum_Cost_Repeat")) {
					OPTIONS->Maximum_Cost_Repeat = Rf_asInteger(Rval);
				} else if (strEQ(key, "Number_Cost_Samples")) {
					OPTIONS->Number_Cost_Samples = Rf_asInteger(Rval);
				} else if (strEQ(key, "Temperature_Ratio_Scale")) {
					OPTIONS->Temperature_Ratio_Scale = Rf_asReal(Rval);
				} else if (strEQ(key, "Temperature_Anneal_Scale")) {
					OPTIONS->Temperature_Anneal_Scale = Rf_asReal(Rval);
				} else if (strEQ(key, "Cost_Parameter_Scale_Ratio")) {
					OPTIONS->Cost_Parameter_Scale_Ratio = Rf_asReal(Rval);
				} else if (strEQ(key, "Initial_Parameter_Temperature")) {
					OPTIONS->Initial_Parameter_Temperature = Rf_asReal(Rval);
				} else if (strEQ(key, "Acceptance_Frequency_Modulus")) {
					OPTIONS->Acceptance_Frequency_Modulus = Rf_asInteger(Rval);
				} else if (strEQ(key, "Generated_Frequency_Modulus")) {
					OPTIONS->Generated_Frequency_Modulus = Rf_asInteger(Rval);
				} else if (strEQ(key, "Reanneal_Cost")) {
					OPTIONS->Reanneal_Cost = Rf_asInteger(Rval);
				} else if (strEQ(key, "Reanneal_Parameters")) {
					OPTIONS->Reanneal_Parameters = Rf_asLogical(Rval);
				} else if (strEQ(key, "Delta_X")) {
					OPTIONS->Delta_X = Rf_asReal(Rval);
				} else if (strEQ(key, "Asa_Out_File")) {
					OPTIONS->Asa_Out_File = R_CHAR(STRING_ELT(Rval, 0));
				} else if (strEQ(key, "User_Quench_Param_Scale")) {
					const auto &ev(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rval));
					quenchParamScale = ev;
				} else if (strEQ(key, "User_Quench_Cost_Scale")) {
					const auto &ev(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Rval));
					quenchCostScale = ev;
				} else {
					Rf_warning("%s: unknown key '%s' for method '%s'",
						   name, key, methodName);
				}
			} else mxThrow("%s: method %d unimplemented", name, method);
		}
	}

	PushLoopIndex pli(name);

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = std::unique_ptr<omxCompute>(omxNewCompute(state, CHAR(s4class)));
	plan->initFromFrontend(state, slotValue);
}

double ComputeGenSA::visita(double temp)
{
	// This code is copied from the Appendix of Tsallis and Stariolo (1996)
	const double pi = M_PI;
	double fator1 = exp(log(temp) / (qv - 1.));
	double fator2 = exp((4. - qv) * log(qv - 1.));
	double fator3 = exp((2. - qv) * M_LN2 / (qv - 1.));
	double fator4 = M_SQRT_PI * fator1 * fator2 / (fator3 * (3. - qv));
	double fator5 = 1. / (qv - 1.) - .5;
	// calculates the gamma function using the reflection formula for
	// 0<arg<1
	double fator6 = pi * (1. - fator5) / sin(pi * (1. - fator5))
		/ fabs(R::gammafn(2.-fator5));
	double sigmax = exp(-(qv - 1.) * log(fator6 / fator4) / (3. - qv));
	double x = sigmax * norm_rand();
	double y = norm_rand();
	double den = exp((qv - 1.) * log((fabs(y))) / (3. - qv));
	double ret_val = x / den;
	if (ret_val > INF) {
		return INF * unif_rand();
	} else if (ret_val < NEG_INF) {
		return NEG_INF * unif_rand();
	}
	return ret_val;
}

double ComputeGenSA::getConstraintPenalty(FitContext *fc)
{
  if (!cvec) {
    cvec =  std::make_unique<ConstraintVec>(fc, "constraint",
                                            [](const omxConstraint &con){ return true; });
  }

  Eigen::ArrayXd constr(cvec->getCount());
  cvec->eval(fc, constr.data());

	double penalty = constr.abs().sum();
	return penalty;
}

static double
asa_cost_function_stub(double *x, double *lower, double *upper,
		       double *tangents, double *curvature,
		       ALLOC_INT *par_dim, int *par_int_real,
		       int *cost_flag, int *exit_code, USER_DEFINES *opt)
{
	return opt->Asa_Data_Ptr->asa_cost(x, cost_flag, exit_code, opt);
}

double ComputeGenSA::asa_cost(double *x, int *cost_flag, int *exit_code, USER_DEFINES *opt)
{
	using Eigen::Map;
	using Eigen::VectorXd;

	FitContext *fc = u_fc;
	Map< VectorXd > proposal(x, numFree);

	{
		ReturnRNGState rrs;
		PushLoopIndex pli(name, opt->N_Generated, opt->Limit_Generated);
		fc->setInform(INFORM_UNINITIALIZED);
    fc->setEstFromOptimizer(proposal);
		fc->wanted = FF_COMPUTE_FIT;
		plan->compute(fc);
	}

	if (Global->interrupted()) return nan("abort");
	if (fc->outsideFeasibleSet()) {
		return std::numeric_limits<double>::max();
	}
	double penalty = getConstraintPenalty(fc);
	fc->setUnscaledFit(fc->getUnscaledFit() + penalty * round(opt->N_Generated / 100)); //OK? TODO
	Global->reportProgress1(contextStr.c_str(), fc->asProgressReport());
	return fc->getUnscaledFit();
}

static double
asa_random_generator(LONG_INT *)
{ return unif_rand(); }

void ComputeGenSA::ingber2012(FitContext *fc)
{
	u_fc = fc;
	LONG_INT seed = 0; //unused
	ALLOC_INT num_parameters = numFree;
	Eigen::VectorXd tangents(numFree);
	tangents.setZero();
	Eigen::VectorXi paramType(numFree);
	paramType.setConstant(REAL_TYPE); // OpenMx doesn't offer integral parameters
	int valid_state_generated_flag=0;
	int exit_status=0;

	if (quenchParamScale.size() == 0) {
		quenchParamScale.resize(numFree);
		quenchParamScale.setConstant(1.);
	}
	if (quenchParamScale.size() != numFree) {
		mxThrow("%s: quenchParamScale must have %d entries instead of %d",
			name, numFree, quenchParamScale.size());
	}
	OPTIONS->User_Quench_Param_Scale = quenchParamScale.data();

	if (quenchCostScale.size() == 0) {
		quenchCostScale.resize(numFree);
		quenchCostScale.setConstant(1.);
	}
	if (quenchCostScale.size() != numFree) {
		mxThrow("%s: quenchCostScale must have %d entries instead of %d",
			name, numFree, int(quenchCostScale.size()));
	}
	OPTIONS->User_Quench_Cost_Scale = quenchCostScale.data();

	OPTIONS->User_Initial_Parameters = 1;
	OPTIONS->Curvature_0 = 1; // can always use ComputeNumericDeriv later
	OPTIONS->Asa_Data_Dim_Ptr = 1;
	OPTIONS->Asa_Data_Ptr = this;

	{
    Eigen::VectorXd curEst(numFree);
    fc->copyEstToOptimizer(curEst);
		BorrowRNGState grs;
		(void) asa(asa_cost_function_stub, asa_random_generator, &seed,
               curEst.data(), lbound.data(), ubound.data(), tangents.data(),
			   0, &num_parameters, paramType.data(), &valid_state_generated_flag,
			   &exit_status, OPTIONS);
	}

	if (!valid_state_generated_flag && verbose) mxLog("invalid state generated");

	switch (exit_status) {
	case NORMAL_EXIT: break; //OK
	case P_TEMP_TOO_SMALL:
		if (verbose >= 1) mxLog("%s: P_TEMP_TOO_SMALL", name);
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case C_TEMP_TOO_SMALL:
		if (verbose >= 1) mxLog("%s: C_TEMP_TOO_SMALL", name);
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case COST_REPEATING:
		if (verbose >= 1) mxLog("%s: COST_REPEATING", name);
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case TOO_MANY_INVALID_STATES:
		if (verbose >= 1) mxLog("%s: TOO_MANY_INVALID_STATES", name);
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case IMMEDIATE_EXIT:
		if (verbose >= 1) mxLog("%s: IMMEDIATE_EXIT", name);
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case INVALID_USER_INPUT:
	case INVALID_COST_FUNCTION:
	case INVALID_COST_FUNCTION_DERIV:
		mxThrow("%s: ASA error %d", name, exit_status);
		break;
	case CALLOC_FAILED:
		mxThrow("%s: out of memory", name);
		break;
	default:
		Rf_warning("%s: unknown exit_status %d", name, exit_status);
		break;
	}
}

void ComputeGenSA::tsallis1996(FitContext *fc)
{
	using Eigen::Map;
	using Eigen::VectorXd;

	int markovLength = stepsPerTemp * numFree;
	curBestFit = fc->getUnscaledFit();
	curBestPenalty = getConstraintPenalty(fc);
  curBest.resize(numFree);
  fc->copyEstToOptimizer(curBest);
  xMini.resize(numFree);
  fc->copyEstToOptimizer(xMini);
	double curFit = curBestFit;
	double curPenalty = curBestPenalty;

	BorrowRNGState grs;
	const double t1 = exp((qv - 1.) * M_LN2) - 1.;

	// Tsallis & Stariolo (1996) start at t=1 (see around Eqn 4)
	for (int tt = 1; !isErrorRaised(); ++tt) {
		// Equation 14' from Tsallis & Stariolo (1996)
		double t2 = exp((qv - 1.) * log(tt + 1.));
		double tem = temSta * t1 / t2;
		if (tem < temEnd) break;

		for (int jj = 0; jj < markovLength && !isErrorRaised(); ++jj) {
			int vx = jj % numFree;
			double va = visita(tem);
			double a = xMini[vx] + va;
			// fmod will truncate some significant digits
			// so we want to avoid it if possible.
			if (a > ubound[vx]) {
				a = ubound[vx] - fmod(a - ubound[vx], range[vx]);
			} else if (lbound[vx] > a) {
				a = fmod(lbound[vx] - a, range[vx]) + lbound[vx];
			}
      fc->setParamFromOptimizer(vx, a);

			{
				ReturnRNGState rrs;
				PushLoopIndex pli(name, jj, markovLength);
				fc->setInform(INFORM_UNINITIALIZED);
				fc->wanted = FF_COMPUTE_FIT;
				plan->compute(fc);
			}

			if (fc->outsideFeasibleSet()) {
        fc->setParamFromOptimizer(vx, xMini[vx]);
				continue;
			}
			double candidateFit = fc->getUnscaledFit();
			double candidatePenalty = getConstraintPenalty(fc);
			if (verbose >= 3) {
				mxLog("%s: [%d] raw penalty %f",
				      name, tt, candidatePenalty);
			}
			double candidateMini = candidateFit + tt * candidatePenalty;
			double eMini = curFit + tt * curPenalty;
			double curBestMini = curBestFit + tt * curBestPenalty;

			// Equation 5 from Tsallis & Stariolo (1996)
			if (candidateMini < eMini) {
        fc->copyEstToOptimizer(xMini);
				curFit = candidateFit;
				curPenalty = candidatePenalty;
				if (verbose >= 2) mxLog("%s: temp %f downhill to %f",
							name, tem, candidateMini);
			} else {
				double worse1 = candidateMini - eMini;
				double qa = qaInit - lambda * tt;
				double worse2 = (1. + (qa-1.) * worse1);
				if (worse2 > 0) {
					double thresh = pow(worse2 / tem, -1./(qa-1.));
					if (thresh >= 1.0 || thresh > unif_rand()) {
            fc->copyEstToOptimizer(xMini);
						curFit = candidateFit;
						curPenalty = candidatePenalty;
						if (verbose >= 2) mxLog("%s: temp %f uphill to %f", name, tem, eMini);
					}
				}
			}
			if (eMini < curBestMini) {
				curBest = xMini;
				curBestFit = curFit;
				curBestPenalty = curPenalty;
			}
			if (curFit != candidateFit) {
				// candidate rejected
        fc->setParamFromOptimizer(vx, xMini[vx]);
			}
			Global->reportProgress(contextStr.c_str(), fc);
		}
	}

  fc->setEstFromOptimizer(curBest);
}

void ComputeGenSA::computeImpl(FitContext *fc)
{
	omxAlgebraPreeval(fitMatrix, fc); // should enable parallel TODO

	using Eigen::Map;
	using Eigen::VectorXd;

	numFree = fc->getNumFree();
	if (numFree <= 0) { complainNoFreeParam(); return; }

	VectorXd attempt(numFree);
  fc->copyEstToOptimizer(attempt);

	lbound.resize(numFree);
	ubound.resize(numFree);
	fc->copyBoxConstraintToOptimizer(lbound, ubound);
	range = ubound - lbound;

	if (verbose >= 1) {
		mxLog("Welcome to %s/%s (%d param)", name, methodName, numFree);
	}

	// ----------

	ComputeFit(optName, fitMatrix, FF_COMPUTE_FIT, fc);
	{
		BorrowRNGState grs;
		int retries = 5;
		while (fc->outsideFeasibleSet() && retries-- > 0) {
			for (int vx=0; vx < attempt.size(); ++vx) {
				attempt[vx] = lbound[vx] + unif_rand() * range[vx];
			}
      fc->setEstFromOptimizer(attempt);
			ComputeFit(optName, fitMatrix, FF_COMPUTE_FIT, fc);
		}
	}

	if (fc->outsideFeasibleSet()) {
		fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);
		return;
	}

	switch (method) {
	case ALGO_TSALLIS1996: tsallis1996(fc); break;
	case ALGO_INGBER2012: ingber2012(fc); break;
	default: mxThrow("%s: unknown method %d", name, method);
	}

	fc->copyParamToModel();  // omit? TODO
	ComputeFit(name, fitMatrix, FF_COMPUTE_FIT, fc);

	if (fc->getInform() != INFORM_UNINITIALIZED || isErrorRaised()) return;
	fc->wanted |= FF_COMPUTE_BESTFIT;
	fc->setInform(INFORM_CONVERGED_OPTIMUM);
}

void ComputeGenSA::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
}
