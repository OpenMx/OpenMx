/*
 *  Copyright 2007-2017 The OpenMx Project
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
#include "omxExportBackendState.h"
#include "nloptcpp.h"
#include "Compute.h"
#include "glue.h"
#include "ComputeGD.h"
#include "ComputeNM.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include <Rmath.h>
#include <R_ext/Utils.h>

#include "EnableWarnings.h"

static const char engineName[] = "NldrMd";


class omxCompute *newComputeNelderMead()
{
	return new omxComputeNM();
}


omxComputeNM::omxComputeNM()
{

}


void omxComputeNM::initFromFrontend(omxState *globalState, SEXP rObj){
	super::initFromFrontend(globalState, rObj);
	
	const double myPI	=	3.141592653589793238462643383280;
	
	SEXP slotValue;
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);
	
	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
	if(OMX_DEBUG){
		mxLog("omxComputeNM member 'verbose' is %d", verbose);
	}
	
	ScopedProtect p2(slotValue, R_do_slot(rObj, Rf_install("nudgeZeroStarts")));
	nudge = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'nudgeZeroStarts' is %d", nudge);
	}
	
	ScopedProtect p3(slotValue, R_do_slot(rObj, Rf_install("defaultMaxIter")));
	defaultMaxIter = Rf_asLogical(slotValue);
	
	ScopedProtect p4(slotValue, R_do_slot(rObj, Rf_install("maxIter")));
	if(defaultMaxIter){maxIter = Global->majorIterations * 5;}
	else{maxIter = Rf_asInteger(slotValue);}
	if(verbose){
		mxLog("omxComputeNM member 'maxIter' is %d", maxIter);
	}
	
	ScopedProtect p5(slotValue, R_do_slot(rObj, Rf_install("alpha")));
	alpha = Rf_asReal(slotValue);
	if(alpha<=0){Rf_error("reflection coefficient 'alpha' must be positive");}
	if(verbose){
		mxLog("omxComputeNM member 'alpha' is %f", alpha);
	}
	
	ScopedProtect p6(slotValue, R_do_slot(rObj, Rf_install("betao")));
	betao = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'betao' is %f", betao);
	}
	
	ScopedProtect p7(slotValue, R_do_slot(rObj, Rf_install("betai")));
	betai = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'betai' is %f", betai);
	}
	if(betao<=0 || betao>=1 || betai<=0 || betai>=1){
		Rf_error("contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)");
	}
	
	ScopedProtect p8(slotValue, R_do_slot(rObj, Rf_install("gamma")));
	gamma = Rf_asReal(slotValue);
	if(gamma>0 && gamma<=alpha){
		Rf_error("if positive, expansion coefficient 'gamma' must be greater than reflection coefficient 'alpha'");
	}
	if(verbose){
		mxLog("omxComputeNM member 'gamma' is %f", gamma);
	}
	
	ScopedProtect p9(slotValue, R_do_slot(rObj, Rf_install("sigma")));
	sigma = Rf_asReal(slotValue);
	if(sigma>=1){Rf_error("shrink coefficient 'sigma' must be less than 1.0");}
	if(verbose){
		mxLog("omxComputeNM member 'sigma' is %f", sigma);
	}
	
	ScopedProtect p10(slotValue, R_do_slot(rObj, Rf_install("bignum")));
	bignum = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'bignum' is %f", bignum);
	}
	
	ScopedProtect p11(slotValue, R_do_slot(rObj, Rf_install("iniSimplexType")));
	if(strEQ(CHAR(Rf_asChar(slotValue)),"regular")){iniSimplexType = 1;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"right")){iniSimplexType = 2;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"smartRight")){iniSimplexType = 3;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"random")){iniSimplexType = 4;}
	else{Rf_error("unrecognized character string provided for Nelder-Mead 'iniSimplexType'");}
	if(verbose){
		mxLog("omxComputeNM member 'iniSimplexType' is %d", iniSimplexType);
	}
	
	ScopedProtect p12(slotValue, R_do_slot(rObj, Rf_install("iniSimplexEdge")));
	iniSimplexEdge = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'iniSimplexEdge' is %f", iniSimplexEdge);
	}
	
	ScopedProtect p13(slotValue, R_do_slot(rObj, Rf_install("iniSimplexMat")));
	if (Rf_length(slotValue)) {
		SEXP matrixDims;
		ScopedProtect pipm(matrixDims, Rf_getAttrib(slotValue, R_DimSymbol));
		int *dimList = INTEGER(matrixDims);
		int rows = dimList[0];
		int cols = dimList[1];
		iniSimplexMat = Eigen::Map< Eigen::MatrixXd >(REAL(slotValue), rows, cols);
	}
	
	ScopedProtect p26(slotValue, R_do_slot(rObj, Rf_install(".iniSimplexColnames")));
	int cnameslen = Rf_length(slotValue);
	if(cnameslen){
		iniSimplexColnames.resize(cnameslen);
		int i;
		for(i=0; i<cnameslen; i++){
			SEXP elem;
			{
				ScopedProtect p27(elem, STRING_ELT(slotValue, i));
				iniSimplexColnames[i] = CHAR(elem);
			}
		}
	}
	
	ScopedProtect p14(slotValue, R_do_slot(rObj, Rf_install("greedyMinimize")));
	greedyMinimize = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'greedyMinimize' is %d", greedyMinimize);
	}
	
	ScopedProtect p15(slotValue, R_do_slot(rObj, Rf_install("altContraction")));
	altContraction = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'altContraction' is %d", altContraction);
	}
	
	ScopedProtect p16(slotValue, R_do_slot(rObj, Rf_install("degenLimit")));
	degenLimit = Rf_asReal(slotValue);
	if(degenLimit<0 || degenLimit>myPI){Rf_error("'degenLimit' must ge within interval [0,pi]");}
	if(verbose){
		mxLog("omxComputeNM member 'degenLimit' is %f", degenLimit);
	}
	
	ScopedProtect p17(slotValue, R_do_slot(rObj, Rf_install("stagnCtrl")));
	if(Rf_length(slotValue)!=2){Rf_error("'stagnCtrl' must be an integer vector of length 2");}
	stagnCtrl[0] = INTEGER(slotValue)[0];
	stagnCtrl[1] = INTEGER(slotValue)[1];
	if(verbose){
		mxPrintMat("omxComputeNM member 'stagnCtrl':", stagnCtrl);
	}
	
	/*ScopedProtect p18(slotValue, R_do_slot(rObj, Rf_install("validationRestart")));
	validationRestart = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'validationRestart' is %d", validationRestart);
	}*/
	
	ScopedProtect p19(slotValue, R_do_slot(rObj, Rf_install("xTolProx")));
	xTolProx = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'xTolProx' is %f", xTolProx);
	}
	
	ScopedProtect p20(slotValue, R_do_slot(rObj, Rf_install("fTolProx")));
	fTolProx = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'fTolProx' is %f", fTolProx);
	}
	
	//Prevent user blunder w/r/t convergence criteria:
	if(xTolProx<=0 && fTolProx<=0){
		fTolProx = 1e-14;
		Rf_warning("both 'xTolProx' and 'fTolProx' are non-positive; 'fTolProx' will be assigned a value of 1e-14");
	}
	
	/*ScopedProtect p23(slotValue, R_do_slot(rObj, Rf_install("doPseudoHessian")));
	doPseudoHessian = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'doPseudoHessian' is %d", doPseudoHessian);
	}*/
	
	ScopedProtect p24(slotValue, R_do_slot(rObj, Rf_install("ineqConstraintMthd")));
	if(strEQ(CHAR(Rf_asChar(slotValue)),"soft")){ineqConstraintMthd = 0;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"eqMthd")){ineqConstraintMthd = 1;}
	else{Rf_error("unrecognized character string provided for Nelder-Mead 'ineqConstraintMthd'");}
	if(verbose){
		mxLog("omxComputeNM member 'ineqConstraintMthd' is %d", ineqConstraintMthd);
	}
	
	ScopedProtect p25(slotValue, R_do_slot(rObj, Rf_install("eqConstraintMthd")));
	if(strEQ(CHAR(Rf_asChar(slotValue)),"soft")){eqConstraintMthd = 1;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"backtrack")){eqConstraintMthd = 2;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"GDsearch")){eqConstraintMthd = 3;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"augLag")){eqConstraintMthd = 4;}
	else{Rf_error("unrecognized character string provided for Nelder-Mead 'eqConstraintMthd'");}
	if(verbose){
		mxLog("omxComputeNM member 'eqConstraintMthd' is %d", eqConstraintMthd);
	}
	
	feasTol = Global->feasibilityTolerance;
	
	ProtectedSEXP Rexclude(R_do_slot(rObj, Rf_install(".excludeVars")));
	excludeVars.reserve(Rf_length(Rexclude));
	for (int ex=0; ex < Rf_length(Rexclude); ++ex) {
		int got = varGroup->lookupVar(CHAR(STRING_ELT(Rexclude, ex)));
		if (got < 0) continue;
		excludeVars.push_back(got);
	}
}


void omxComputeNM::computeImpl(FitContext *fc){
	
	omxAlgebraPreeval(fitMatrix, fc);
	if (isErrorRaised()) return;
	
	size_t numParam = fc->varGroup->vars.size();
	if (excludeVars.size()) {
		fc->profiledOut.assign(fc->numParam, false);
		for (auto vx : excludeVars) fc->profiledOut[vx] = true;
	}
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
	fc->createChildren(fitMatrix);
	
	
	//int beforeEval = fc->getLocalComputeCount();
	
	if (verbose >= 1){
		//mxLog something here
	}
	
	NelderMeadOptimizerContext nmoc(fc, this);
	if(eqConstraintMthd==4){Rf_error("'augLag' Not Yet Implemented");}
	nmoc.verbose = verbose;
	nmoc.countConstraintsAndSetupBounds();
	nmoc.invokeNelderMead();
	
	//TODO: validation restart, pseudoHessian, etc.
	return;
}


void omxComputeNM::reportResults(FitContext *fc, MxRList *slots, MxRList *out){
	omxPopulateFitFunction(fitMatrix, out);
	//TODO: export final simplex and fit-value vector, pseudoHessian, etc.
	return;
}

//-------------------------------------------------------

NelderMeadOptimizerContext::NelderMeadOptimizerContext(FitContext* _fc, omxComputeNM* _nmo)
	: fc(_fc), NMobj(_nmo), numFree(countNumFree())
{
	est.resize(numFree);
	copyParamsFromFitContext(est.data());
	backtrackSteps=10; //<--Eventually should be made user-settable
}

void NelderMeadOptimizerContext::copyBounds()
{
	FreeVarGroup *varGroup = fc->varGroup;
	int px=0;
	for (size_t vx=0; vx < fc->profiledOut.size(); ++vx) {
		if (fc->profiledOut[vx]) continue;
		solLB[px] = varGroup->vars[vx]->lbound;
		if (!std::isfinite(solLB[px])) solLB[px] = NEG_INF;
		solUB[px] = varGroup->vars[vx]->ubound;
		if (!std::isfinite(solUB[px])) solUB[px] = INF;
		++px;
	}
}

void NelderMeadOptimizerContext::countConstraintsAndSetupBounds()
{
	solLB.resize(numFree);
	solUB.resize(numFree);
	copyBounds();
	
	omxState *globalState = fc->state;
	globalState->countNonlinearConstraints(numEqC, numIneqC, false);
	if(verbose){
		mxLog("counted %d equality constraints",numEqC);
		mxLog("counted %d inequality constraints",numIneqC);
	}
	//If there aren't any of one of the two constraint types, then the
	//method for handling them shouldn't matter.  But, switching the
	//method to the simplest setting helps simplify programming logic:
	if(!numIneqC){NMobj->ineqConstraintMthd = 0;}
	if(!numEqC){NMobj->eqConstraintMthd = 1;}
	equality.resize(numEqC);
	inequality.resize(numIneqC);
}

int NelderMeadOptimizerContext::countNumFree()
{
	int nf = 0;
	for (size_t vx=0; vx < fc->profiledOut.size(); ++vx) {
		if (fc->profiledOut[vx]) continue;
		++nf;
	}
	return nf;
}

void NelderMeadOptimizerContext::copyParamsFromFitContext(double *ocpars)
{
	int px=0;
	for (size_t vx=0; vx < fc->profiledOut.size(); ++vx) {
		if (fc->profiledOut[vx]) continue;
		ocpars[px] = fc->est[vx];
		++px;
	}
}

void NelderMeadOptimizerContext::copyParamsFromOptimizer(Eigen::VectorXd &x, FitContext* fc2)
{
	int px=0;
	for (size_t vx=0; vx < fc2->profiledOut.size(); ++vx) {
		if (fc2->profiledOut[vx]) continue;
		fc2->est[vx] = x[px];
		++px;
	}
	fc2->copyParamToModel();
}

//----------------------------------------------------------------------

void NelderMeadOptimizerContext::enforceBounds(Eigen::VectorXd &x){
	int i=0;
	for(i=0; i < x.size(); i++){
		if(x[i] < solLB[i]){x[i] = solLB[i];}
		if(x[i] > solUB[i]){x[i] = solUB[i];}
	}
}

bool NelderMeadOptimizerContext::checkBounds(Eigen::VectorXd &x){
	bool retval=true;
	int i=0;
	for(i=0; i < x.size(); i++){
		if(x[i] < solLB[i] && x[i] > solUB[i]){
			retval=false;
			break;
		}
	}
	return(retval);
}

void NelderMeadOptimizerContext::evalIneqC()
{
	if(!numIneqC){return;}
	
	omxState *st = fc->state;
	int ineqType = omxConstraint::LESS_THAN;
	
	int cur=0, j=0;
	for (j=0; j < int(st->conListX.size()); j++) {
		omxConstraint &con = *st->conListX[j];
		if (con.opCode == omxConstraint::EQUALITY) continue;
		con.refreshAndGrab(fc, (omxConstraint::Type) ineqType, &inequality(cur));
		//Nelder-Mead, of course, does not use constraint Jacobians...
		cur += con.size;
	}
	//Nelder-Mead will not care about the function values of inactive inequality constraints:
	inequality = inequality.array().max(0.0);
	
	if (NMobj->verbose >= 3) {
		mxPrintMat("inequality", inequality);
	}
	
}

void NelderMeadOptimizerContext::evalEqC()
{
	if(!numEqC){return;}
	
	omxState *st = fc->state;
	
	int cur=0, j=0;
	for(j = 0; j < int(st->conListX.size()); j++) {
		omxConstraint &con = *st->conListX[j];
		if (con.opCode != omxConstraint::EQUALITY) continue;
		con.refreshAndGrab(fc, &equality(cur));
		//Nelder-Mead, of course, does not use constraint Jacobians...
		cur += con.size;
	}
	
	if (NMobj->verbose >= 3) {
		mxPrintMat("equality", equality);
	}
}

double NelderMeadOptimizerContext::evalFit(Eigen::VectorXd &x)
{
	copyParamsFromOptimizer(x,fc);
	ComputeFit(engineName, NMobj->fitMatrix, FF_COMPUTE_FIT, fc);
	if(!std::isfinite(fc->fit)){return(NMobj->bignum);}
	else{
		double fv = fc->fit;
		if(NMobj->eqConstraintMthd==4){
			//TODO: add terms from augmented Lagrangian to fv
		}
		return(fv);
	}
}

void NelderMeadOptimizerContext::checkNewPointInfeas(Eigen::VectorXd &x, Eigen::Vector2i &ifcr)
{
	int i=0;
	double feasTol = NMobj->feasTol;
	ifcr.setZero(2);
	if(!numIneqC && !numEqC){return;}
	copyParamsFromOptimizer(x,fc);
	evalIneqC();
	evalEqC();
	if(numIneqC){
		for(i=0; i < inequality.size(); i++){
			if(inequality[i] > feasTol){
				ifcr[0] = 1;
				break;
			}
		}
	}
	if(numEqC){
		for(i=0; i < equality.size(); i++){
			if(fabs(equality[i]) > feasTol){
				ifcr[1] = 1;
				break;
			}
		}
	}
}


void NelderMeadOptimizerContext::evalFirstPoint(Eigen::VectorXd &x, double &fv, int &infeas)
{
	Eigen::Vector2i ifcr;
	int ineqConstraintMthd = NMobj->ineqConstraintMthd;
	int eqConstraintMthd = NMobj->eqConstraintMthd;
	double bignum = NMobj->bignum;
	enforceBounds(x);
	checkNewPointInfeas(x, ifcr);
	if(!ifcr.sum()){
		infeas = 0L;
		fv = (evalFit(x));
		//if(fv==NMobj->bignum){infeas=1L;}
		return;
	}
	else if(ifcr[1] || (ifcr[0] && ineqConstraintMthd)){
		switch(eqConstraintMthd){
		case 1:
			infeas = 1L;
			fv = bignum;
			return;
		case 2:
			//Can't backtrack to someplace else if it's the very first point.
			//Rf_error("starting values not feasible; re-specify the initial simplex, or consider mxTryHard()");
			infeas = 1L;
			fv = bignum;
			break;
		case 3:
			Rf_error("'GDsearch' Not Yet Implemented");
		case 4:
			if(ifcr[0]){
				fv = bignum;
				infeas = 1L;
			}
			else{
				fv = evalFit(x);
				infeas = 0L;
			}
			return;
		}
	}
	else if(ifcr[0]){
		fv = bignum;
		infeas = 1L;
		return;
	}
}

//oldpt is used for backtracking:
void NelderMeadOptimizerContext::evalNewPoint(Eigen::VectorXd &newpt, Eigen::VectorXd &oldpt, double &fv, int &newInfeas, int oldInfeas)
{
	Eigen::Vector2i ifcr;
	int ineqConstraintMthd = NMobj->ineqConstraintMthd;
	int eqConstraintMthd = NMobj->eqConstraintMthd;
	double bignum = NMobj->bignum;
	enforceBounds(newpt);
	checkNewPointInfeas(newpt, ifcr);
	if(!ifcr.sum()){
		newInfeas = 0L;
		fv = (evalFit(newpt));
		return;
	}
	else if(ifcr[1] || (ifcr[0] && ineqConstraintMthd)){
		switch(eqConstraintMthd){
		case 1:
			newInfeas = 1L;
			fv = bignum;
			return;
		case 2:
			//If old point is not feasible, there's no sense in backtracking toward it:
			if(oldInfeas){
				newInfeas = 1L;
				fv = bignum;
				return;
			}
			else{Rf_error("'backtrack' Not Yet Implemented");}
		case 3:
			Rf_error("'GDsearch' Not Yet Implemented");
		case 4:
			if(ifcr[0]){
				fv = bignum;
				newInfeas = 1L;
			}
			else{
				fv = evalFit(newpt);
				newInfeas = 0L;
			}
			return;
		}
	}
	else if(ifcr[0]){
		fv = bignum;
		newInfeas = 1L;
		return;
	}
}

void NelderMeadOptimizerContext::jiggleCoord(Eigen::VectorXd &xin, Eigen::VectorXd &xout){
	double a,b;
	int i;
	GetRNGstate();
	for(i=0; i < xin.size(); i++){
		b = Rf_runif(0.75,1.25);
		a = Rf_runif(-0.25,0.25);
		xout[i] = b*xin[i] + a;
	}
	PutRNGstate();
}

void NelderMeadOptimizerContext::printProblemState()
{
	int i=0;
	Eigen::MatrixXd tmpvrt(n+1,numFree);
	for(i=0; i<n+1; i++){tmpvrt.row(i) = vertices[i];}
	mxPrintMat("working simplex:",tmpvrt);
	mxPrintMat("fitfunction values:",fvals);
	mxPrintMat("infeasibility states:",vertexInfeas);
}

void NelderMeadOptimizerContext::printNewPoint(Eigen::VectorXd &x, double fv, int isbad)
{
	mxPrintMat("coordinates:",x);
	mxLog("fitfunction value: %f",fv);
	mxLog("infeasible?: %d",isbad);
}

void NelderMeadOptimizerContext::initializeSimplex(Eigen::VectorXd &startpt, double edgeLength, bool isRestart)
{
	if(verbose){mxLog("(re-)initializing simplex");}
	int i=0;
	Eigen::VectorXd xin, xout, newpt, oldpt;
	if(NMobj->iniSimplexMat.rows() && NMobj->iniSimplexMat.cols() && !isRestart){
		Eigen::MatrixXd SiniSupp, iniSimplexMat2;
		Eigen::VectorXi paramMap(numFree);
		if(NMobj->iniSimplexMat.cols() != numFree){
			Rf_error("'iniSimplexMat' has %d columns, but %d columns expected",NMobj->iniSimplexMat.cols(), numFree);
		}
		if( int(NMobj->iniSimplexColnames.size()) != numFree){
			Rf_error("'iniSimplexMat' has %d column names, but %d column names expected", NMobj->iniSimplexColnames.size(), numFree);
		}
		if(NMobj->iniSimplexMat.rows()>n+1){
			Rf_warning("'iniSimplexMat' has %d rows, but %d rows expected; extraneous rows will be ignored",NMobj->iniSimplexMat.rows(), n+1);
			NMobj->iniSimplexMat.conservativeResize(n+1,numFree);
		}
		iniSimplexMat2.resize(NMobj->iniSimplexMat.rows(), numFree);
		int gx=0;
		/*If there are no problems, then every time vx gets incremented, it should become equal to the current 
		 value of gx*/
		for (int vx=0; vx < int(fc->varGroup->vars.size()); ++vx) {
			for (int nx=0; nx < int(NMobj->iniSimplexColnames.size()); ++nx) {
				if (strEQ(NMobj->iniSimplexColnames[nx], fc->varGroup->vars[vx]->name)) {
					paramMap[gx] = vx;
					++gx;
					break;
				}
			}
		}
		if ( gx != int(NMobj->iniSimplexColnames.size()) ){
			Rf_error("error in mapping column names of 'iniSimplexMat' to free-parameter labels");
		}
		for(i=0; i < NMobj->iniSimplexMat.cols(); i++){
			iniSimplexMat2.col(paramMap[i]) = NMobj->iniSimplexMat.col(i);
		}
		if(NMobj->iniSimplexMat.rows()<n+1){
			Rf_warning("'iniSimplexMat' has %d rows, but %d rows expected; omitted rows will be generated randomly",NMobj->iniSimplexMat.rows(),n+1);
			SiniSupp.resize(n + 1 - NMobj->iniSimplexMat.rows(), numFree);
			xin=iniSimplexMat2.row(0);
			for(i=0; i<SiniSupp.rows(); i++){
				xout=SiniSupp.row(i);
				jiggleCoord(xin, xout);
				SiniSupp.row(i) = xout;
			}
		}
		for(i=0; i < NMobj->iniSimplexMat.rows(); i++){
			vertices[i] = iniSimplexMat2.row(i);
		}
		if(SiniSupp.rows()){
			for(i=0; i<SiniSupp.rows(); i++){
				vertices[i+NMobj->iniSimplexMat.rows()] = SiniSupp.row(i);
			}
		}
	}
	else{
		//TODO: regular simplex for other than unit edge-length:
		double k = (double) n;
		double shhp = (edgeLength/sqrt(k))*(1/k/sqrt(2))*(-1.0 + k + sqrt(1.0+k));
		double shhq = (edgeLength/sqrt(k))*(1/k/sqrt(2))*(sqrt(1.0+k)-1);
		//double shhp = (1/k/sqrt(2))*(-1.0 + k + sqrt(1.0+k));
		//double shhq = (1/k/sqrt(2))*(sqrt(1.0+k)-1);
		Eigen::VectorXd xu, xd;
		double fu=0, fd=0;
		int badu=0, badd=0;
		//TODO: None of these except case 4 are OK to use with equality MxConstraints
		switch(NMobj->iniSimplexType){
		case 1:
			vertices[0].setZero(numFree);
			for(i=1; i<n+1; i++){
				vertices[i].setConstant(numFree,shhq);
				vertices[i][i-1] = shhp;
			}
			for(i=0; i<n+1; i++){
				vertices[i] += startpt;
			}
			break;
		case 2:
			vertices[0] = startpt;
			for(i=1; i<n+1; i++){
				vertices[i] = vertices[0];
				vertices[i][i-1] += edgeLength;
			}
			break;
		case 3:
			//TODO: this could be even smarter if it also figured out different edge lengths 
			//to account for different scaling of the free parameters:
			//oldpt = vertices.row(0); //<--oldpt
			vertices[0] = startpt;
			evalFirstPoint(vertices[0], fvals[0], vertexInfeas[0]);
			for(i=0; i<n; i++){
				xu = vertices[0];
				xu[i] += edgeLength;
				xd = vertices[0];
				xd[i] -= edgeLength;
				evalNewPoint(xu, vertices[0], fu, badu, vertexInfeas[0]);
				evalNewPoint(xd, vertices[0], fd, badd, vertexInfeas[0]);
				vertices[i+1] = fu<fd ? xu : xd;
				fvals[i+1] = fu<fd ? fu : fd;
				vertexInfeas[i+1] = fu<fd ? badu : badd;
			}
			if(verbose){printProblemState();}
			return;
		case 4:
			vertices[0] = startpt;
			for(i=1; i<n+1; i++){
				vertices[i].setZero(numFree);
				jiggleCoord(vertices[0],vertices[i]);
			}
			break;
		}
	}
	//Now evaluate each vertex:
	evalFirstPoint(vertices[0], fvals[0], vertexInfeas[0]);
	for(i=1; i<n+1; i++){
		evalNewPoint(vertices[i], vertices[0], fvals[i], vertexInfeas[i], vertexInfeas[0]);
	}
	//if(verbose){printProblemState();}
}


void NelderMeadOptimizerContext::fullSort()
{
	int i=0;
	Eigen::VectorXi ind(n+1);
	for(i=0; i<=n; i++){
		ind[i] = i;
	}
	Eigen::VectorXi tmpVertexInfeas = vertexInfeas;
	std::vector<Eigen::VectorXd> tmpVertices = vertices;
	//If we don't care about tie-breaking rules:
	if( (fvals.tail(n).array() < fvals[0]).any() ){
		unchangedx0count = 0;
		rsort_with_index(fvals.data(), ind.data(), n+1);
		for(i=0; i<n+1; i++){
			vertices[i] = tmpVertices[ind[i]];
			vertexInfeas[i] = tmpVertexInfeas[ind[i]];
		}
	}
	else{
		unchangedx0count++;
		Eigen::VectorXi ind_tail = ind.tail(n);
		Eigen::VectorXd fvals_tail = fvals.tail(n);
		rsort_with_index(fvals_tail.data(), ind_tail.data(), n);
		for(i=1; i<n+1; i++){
			fvals[i] = fvals_tail[i-1];
			vertices[i] = tmpVertices[ind_tail[i-1]];
			vertexInfeas[i] = tmpVertexInfeas[ind_tail[i-1]];
		}
	}
	//Calculate centroids:
	subcentroid.setZero(numFree);
	eucentroidCurr.setZero(numFree);
	for(i=0; i<n+1; i++){
		eucentroidCurr += vertices[i] / (n+1);
		if(i<n){subcentroid += vertices[i] / n;}
	}
	Eigen::Vector2i scfcr;
	scfcr.setZero();
	checkNewPointInfeas(subcentroid, scfcr);
	badsc = (scfcr.sum()) ? 1 : 0;
	
	needFullSort = false;
	if(verbose){
		mxLog("full sort complete...");
		printProblemState();
	}
	return;
}

void NelderMeadOptimizerContext::fastSort()
{
	int i=0, j;
	Eigen::VectorXi tmpVertexInfeas = vertexInfeas;
	std::vector<Eigen::VectorXd> tmpVertices = vertices;
	Eigen::VectorXd tmpFvals = fvals;
	if(tmpFvals[n]<tmpFvals[0]){
		unchangedx0count = 0;
		fvals[0] = tmpFvals[n];
		vertices[0] = tmpVertices[n];
		vertexInfeas[0] = tmpVertexInfeas[n];
		for(i=1; i<=n; i++){
			fvals[i] = tmpFvals[i-1];
			vertices[i] = tmpVertices[i-1];
			vertexInfeas[i] = tmpVertexInfeas[i-1];
		}
	}
	else{
		unchangedx0count++;
		for(i=n-1; i>=0; i--){
			if(tmpFvals[i] > tmpFvals[n]){
				fvals[i+1] = tmpFvals[i];
				vertices[i+1] = tmpVertices[i];
				vertexInfeas[i+1] = tmpVertexInfeas[i];
			}
			else{
				fvals[i+1] = tmpFvals[n];
				vertices[i+1] = tmpVertices[n];
				vertexInfeas[i+1] = tmpVertexInfeas[n];
				break;
			}
		}
		for(j=i; j>=0; j--){
			fvals[j] = tmpFvals[j];
			vertices[j] = tmpVertices[j];
			vertexInfeas[j] = tmpVertexInfeas[j];
		}
	}
	//TODO: this could be made faster, since we do fastSort() when only one vertex of the simplex has changed:
	subcentroid.setZero(numFree);
	eucentroidCurr.setZero(numFree);
	for(i=0; i<n+1; i++){
		eucentroidCurr += vertices[i] / (n+1);
		if(i<n){subcentroid += vertices[i] / n;}
	}
	Eigen::Vector2i scfcr;
	scfcr.setZero();
	checkNewPointInfeas(subcentroid, scfcr);
	badsc = (scfcr.sum()) ? 1 : 0;
	
	if(verbose){
		mxLog("fast sort complete...");
		printProblemState();
	}
	return;
}


void NelderMeadOptimizerContext::simplexTransformation()
{
	failedContraction = false;
	oldWorstVertex = vertices[n];
	//Reflection transformation:
	xr = subcentroid + NMobj->alpha*(subcentroid - vertices[n]);
	evalNewPoint(xr, subcentroid, fr, badr, badsc);
	if(verbose){
		mxLog("reflection point...");
		printNewPoint(xr, fr, badr);
	}
	if(fr<fvals[n-1]){ //<--If fit at reflection point is better than second worst fit
		//If fit at reflection point is worse than best fit, or expansions are turned off, accept reflection point:
		if(fr>=fvals[0] || NMobj->gamma<0){
			fvals[n] = fr;
			vertices[n] = xr;
			vertexInfeas[n] = badr;
			needFullSort=false;
			if(verbose){mxLog("reflection point accepted");}
			return;
		}
		else{ //<--If fit at reflection point is better than best fit and expansions are turned on
			//Expansion transformation:
			xe = subcentroid + NMobj->gamma*(xr - subcentroid);
			evalNewPoint(xe, subcentroid, fe, bade, badsc);
			if(verbose){
				mxLog("expansion point...");
				printNewPoint(xe, fe, bade);
			}
			if(NMobj->greedyMinimize){ //<--If using greedy minimization
				//Accept the better of the reflection and expansion points:
				fvals[n] = (fr<fe) ? fr : fe;
				vertices[n] = (fr<fe) ? xr : xe;
				vertexInfeas[n] = (fr<fe) ? badr : bade;
				needFullSort=false;
				if(verbose){
					if(fr<fe){mxLog("reflection point accepted");}
					else{mxLog("expansion point accepted");}
				}
				return;
			}
			else{ //<--If using greedy expansion
				//Accept expansion point unless reflection point is strictly better:
				fvals[n] = (fe<fvals[0]) ? fe : fr;
				vertices[n] = (fe<fvals[0]) ? xe : xr;
				vertexInfeas[n] = (fe<fvals[0]) ? bade : badr;
				needFullSort=false;
				if(verbose){
					if(fe<fvals[0]){mxLog("expansion point accepted");}
					else{mxLog("reflection point accepted");}
				}
				return;
			}
		}
	}
	else{
		if(fr<fvals[n]){ //<--If fit at reflection point is at least better than the worst fit
			//Outside-contraction transformation:
			if(!NMobj->altContraction){
				xoc = subcentroid + NMobj->betao*(xr - subcentroid);
				evalNewPoint(xoc, subcentroid, foc, badoc, badsc);
			}
			else{
				//Eigen::VectorXd oldpt = vertices.row(0);
				xoc = vertices[0] + NMobj->betao*(xr - vertices[0]);
				evalNewPoint(xoc, vertices[0], foc, badoc, vertexInfeas[0]);
			}
			if(verbose){
				mxLog("outside contraction point...");
				printNewPoint(xoc, foc, badoc);
			}
			if(foc<=fr){ //<--If fit at xoc is no worse than fit at reflection point
				//Accept xoc:
				fvals[n] = foc;
				vertices[n] = xoc;
				vertexInfeas[n] = badoc;
				needFullSort=false;
				if(verbose){mxLog("outside contraction point accepted");}
				return;
			}
			else if(NMobj->sigma<=0){ //<--If fit at xoc is worse than fit at reflection point, and shrinks are turned off
				failedContraction = true;
				//Accept reflection point:
				fvals[n] = fr;
				vertices[n] = xr;
				vertexInfeas[n] = badr;
				needFullSort=false;
				if(verbose){mxLog("reflection point accepted");}
				return;
			}
		}
		else{ //<--If fit at reflection point is no better than worst fit
			//Inside-contraction transformation:
			if(!NMobj->altContraction){
				xic = subcentroid + NMobj->betai*(vertices[n] - subcentroid);
				evalNewPoint(xic, subcentroid, fic, badic, badsc);
			}
			else{
				xic = vertices[0] + NMobj->betai*(vertices[n] - vertices[0]);
				evalNewPoint(xic, vertices[0], fic, badic, vertexInfeas[0]);
			}
			if(verbose){
				mxLog("inside contraction point...");
				printNewPoint(xic, fic, badic);
			}
			if(fic<fvals[n]){ //<--If fit at xic is better than worst fit
				//Accept xic:
				fvals[n] = fic;
				vertices[n] = xic;
				vertexInfeas[n] = badic;
				needFullSort=false;
				if(verbose){mxLog("inside contraction point accepted");}
				return;
			}
			else if(NMobj->sigma<=0){
				failedContraction = true;
				if(verbose){mxLog("contraction failed and shrinks are switched off...");}
				return;
			}
		}
		//Shrink transformation:
		if(NMobj->sigma>0){
			int i=0;
			std::vector<Eigen::VectorXd> tmpVertices = vertices;
			Eigen::VectorXi tmpVertexInfeas = vertexInfeas;
			for(i=1; i<n+1; i++){
				vertices[i] = vertices[0] + NMobj->sigma*(vertices[i] - vertices[0]);
				evalNewPoint(vertices[i], tmpVertices[i], fvals[i], vertexInfeas[i], tmpVertexInfeas[i]);
			}
			needFullSort=true;
			if(verbose){mxLog("shrink transformation complete");}
			return;
		}
	}
}


//TODO:status codes
bool NelderMeadOptimizerContext::checkConvergence(){
	int i=0;
	Eigen::VectorXd xdiffs(n);
	Eigen::VectorXd fdiffs(n);
	double fprox, xprox;
	//Range-convergence test:
	if(NMobj->fTolProx > 0){
		for(i=0; i<n; i++){
			fdiffs[i] = fabs(fvals[i+1] - fvals[0]);
		}
		fprox = fdiffs.array().maxCoeff();
		if(verbose){mxLog("range proximity measure: %f",fprox);}
		if(fprox < NMobj->fTolProx){
			return(true);
		}
	}
	//Domain-convergence test:
	if(NMobj->fTolProx > 0){
		for(i=0; i<n; i++){
			xdiffs[i] = (vertices[i+1] - vertices[0]).array().abs().maxCoeff();
		}
		xprox = xdiffs.array().maxCoeff();
		if(verbose){mxLog("domain proximity measure: %f",xprox);}
		if(xprox < NMobj->xTolProx){
			return(true);
		}
	}
	if(itersElapsed >= NMobj->maxIter){
		//TODO: set code 4 (status blue)
		return(true);
	}
	return(false);
}


bool NelderMeadOptimizerContext::checkProgress(){
	const double myPI	=	3.141592653589793238462643383280;
	Eigen::VectorXd d1, d2;
	double t;
	int i, j, k;
	if(failedContraction && NMobj->sigma<=0){
		return(true);
	}
	if(NMobj->stagnCtrl[0]>0 && NMobj->stagnCtrl[1]>0 && 
    unchangedx0count>=NMobj->stagnCtrl[0] && NMobj->stagnCtrl[1]<restartsUsed){
		return(true);
	}
	if(NMobj->degenLimit>0){
		for(i=0; i<n+1; i++){
			for(j=0; j<n; j++){
				if(j==i){continue;}
				for(k=j+1; k<n+1; k++){
					d1 = vertices[i] - vertices[j];
					d2 = vertices[i] - vertices[k];
					t = acos( d1.dot(d2) / sqrt(d1.dot(d1)) / sqrt(d2.dot(d2)) );
					if(t < NMobj->degenLimit || myPI - t < NMobj->degenLimit){
						return(true);
					}
				}
			}
		}
	}
	return(false);
}


//TODO: pass optimization status codes to FitContext.
void NelderMeadOptimizerContext::invokeNelderMead(){
	n = numFree - numEqC;
	vertices.resize(n+1);
	fvals.resize(n+1);
	vertexInfeas.resize(n+1);
	subcentroid.resize(numFree);
	eucentroidCurr.resize(numFree);
	initializeSimplex(est, NMobj->iniSimplexEdge, false);
	if(vertexInfeas.sum()==n+1 || (fvals.array()==NMobj->bignum).all()){
		Rf_error("initial simplex is not feasible; specify it differently, try different start values, or use mxTryHard()");
	}
	fullSort();
	needFullSort=false;
	bool needRestart = false;
	bool stopflag=false;
	itersElapsed = 0;
	restartsUsed = 0;
	
	//Loop is: sort, check convergence, check progress, transform;
	do{
		if(itersElapsed){
			//Order the vertices by fit value:
			if(needFullSort){fullSort();}
			else{fastSort();}
			
			stopflag = checkConvergence();
			if(stopflag){break;}
			
			needRestart = checkProgress();
			if(needRestart){
				initializeSimplex(vertices[0], sqrt((vertices[0]-vertices[1]).dot(vertices[0]-vertices[1])), true);
				needRestart = false;
				needFullSort = true;
				itersElapsed++;
				continue;
			}
		}
		
		simplexTransformation();
		
		eucentroidPrev = eucentroidCurr;
		itersElapsed++;
	} while (!stopflag);
	
	if(verbose){mxPrintMat("solution?",vertices[0]);}
	
	//TODO: check to see if fit at centroids is any better than at best vertex
	fvals[0] = evalFit(vertices[0]); //<--This is to assign the parameter and fit values at the solution to the FitContext as part of evalFit() .
	fc->iterations = itersElapsed;
}


