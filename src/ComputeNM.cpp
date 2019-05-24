/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
#include "omxNLopt.h"
#include "Compute.h"
#include "glue.h"
#include "ComputeGD.h"
#include "ComputeNM.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include <Rmath.h>
#include <R_ext/Utils.h>

#include "nlopt.h"
#include "nlopt-internal.h"
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
	
	//TODO: use const defined in Rmath.h:
	const double myPI	=	3.141592653589793238462643383280;
	
	SEXP slotValue;
	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
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
	if(defaultMaxIter){maxIter = Global->majorIterations * 10;}
	else{maxIter = Rf_asInteger(slotValue);}
	if(verbose){
		mxLog("omxComputeNM member 'maxIter' is %d", maxIter);
	}
	
	ScopedProtect p5(slotValue, R_do_slot(rObj, Rf_install("alpha")));
	alpha = Rf_asReal(slotValue);
	if(alpha<=0){mxThrow("reflection coefficient 'alpha' must be positive");}
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
		mxThrow("contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)");
	}
	
	ScopedProtect p8(slotValue, R_do_slot(rObj, Rf_install("gamma")));
	gamma = Rf_asReal(slotValue);
	if(gamma>0 && gamma<=alpha){
		mxThrow("if positive, expansion coefficient 'gamma' must be greater than reflection coefficient 'alpha'");
	}
	if(verbose){
		mxLog("omxComputeNM member 'gamma' is %f", gamma);
	}
	
	ScopedProtect p9(slotValue, R_do_slot(rObj, Rf_install("sigma")));
	sigma = Rf_asReal(slotValue);
	if(sigma>=1){mxThrow("shrink coefficient 'sigma' must be less than 1.0");}
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
	else{mxThrow("unrecognized character string provided for Nelder-Mead 'iniSimplexType'");}
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
	if(degenLimit<0 || degenLimit>myPI){mxThrow("'degenLimit' must ge within interval [0,pi]");}
	if(verbose){
		mxLog("omxComputeNM member 'degenLimit' is %f", degenLimit);
	}
	
	ScopedProtect p17(slotValue, R_do_slot(rObj, Rf_install("stagnCtrl")));
	if(Rf_length(slotValue)!=2){mxThrow("'stagnCtrl' must be an integer vector of length 2");}
	stagnCtrl[0] = INTEGER(slotValue)[0];
	stagnCtrl[1] = INTEGER(slotValue)[1];
	if(verbose){
		mxPrintMat("omxComputeNM member 'stagnCtrl':", stagnCtrl);
	}
	
	ScopedProtect p18(slotValue, R_do_slot(rObj, Rf_install("validationRestart")));
	validationRestart = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'validationRestart' is %d", validationRestart);
	}
	
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
	
	ScopedProtect p30(slotValue, R_do_slot(rObj, Rf_install("doPseudoHessian")));
	doPseudoHessian = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'doPseudoHessian' is %d", doPseudoHessian);
	}
	
	ScopedProtect p24(slotValue, R_do_slot(rObj, Rf_install("ineqConstraintMthd")));
	if(strEQ(CHAR(Rf_asChar(slotValue)),"soft")){ineqConstraintMthd = 0;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"eqMthd")){ineqConstraintMthd = 1;}
	else{mxThrow("unrecognized character string provided for Nelder-Mead 'ineqConstraintMthd'");}
	if(verbose){
		mxLog("omxComputeNM member 'ineqConstraintMthd' is %d", ineqConstraintMthd);
	}
	
	ScopedProtect p25(slotValue, R_do_slot(rObj, Rf_install("eqConstraintMthd")));
	if(strEQ(CHAR(Rf_asChar(slotValue)),"soft")){eqConstraintMthd = 1;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"backtrack")){eqConstraintMthd = 2;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"GDsearch")){eqConstraintMthd = 3;}
	else if(strEQ(CHAR(Rf_asChar(slotValue)),"l1p")){eqConstraintMthd = 4;}
	else{mxThrow("unrecognized character string provided for Nelder-Mead 'eqConstraintMthd'");}
	if(verbose){
		mxLog("omxComputeNM member 'eqConstraintMthd' is %d", eqConstraintMthd);
	}
	
	ScopedProtect p28(slotValue, R_do_slot(rObj, Rf_install("backtrackCtrl1")));
	backtrackCtrl1 = Rf_asReal(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'backtrackCtrl1' is %f", backtrackCtrl1);
	}
	ScopedProtect p29(slotValue, R_do_slot(rObj, Rf_install("backtrackCtrl2")));
	backtrackCtrl2 = Rf_asInteger(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'backtrackCtrl2' is %d", backtrackCtrl2);
	}
	
	ScopedProtect p31(slotValue, R_do_slot(rObj, Rf_install("centerIniSimplex")));
	centerIniSimplex = Rf_asLogical(slotValue);
	if(verbose){
		mxLog("omxComputeNM member 'centerIniSimplex' is %d", centerIniSimplex);
	}
	
	feasTol = Global->feasibilityTolerance;
}


void omxComputeNM::computeImpl(FitContext *fc){
	
	omxAlgebraPreeval(fitMatrix, fc);
	if (isErrorRaised()) return;
	
	fc->ensureParamWithinBox(nudge);
	fc->createChildren(fitMatrix);
	
	NelderMeadOptimizerContext nmoc(fc, this);

	if (nmoc.numFree <= 0) { complainNoFreeParam(); return; }
	
	nmoc.verbose = verbose;
	nmoc.maxIter = maxIter;
	nmoc.iniSimplexType = iniSimplexType;
	nmoc.iniSimplexEdge = iniSimplexEdge;
	nmoc.centerIniSimplex = centerIniSimplex;
	nmoc.fit2beat = R_PosInf;
	nmoc.bignum = bignum;
	nmoc.iniSimplexMat = iniSimplexMat;
	nmoc.ineqConstraintMthd = ineqConstraintMthd;
	nmoc.eqConstraintMthd = eqConstraintMthd;
	nmoc.countConstraintsAndSetupBounds();
	if(nmoc.eqConstraintMthd==4 && (nmoc.numEqC || (nmoc.ineqConstraintMthd && nmoc.numIneqC))){
		if(verbose){mxLog("starting l1-penalty algorithm");}
		fc->iterations = 0; //<--Not sure about this
		nmoc.maxIter = maxIter/10;
		nmoc.addPenalty = true;
		int k;
		for(k=0; k<=10; k++){
			if(verbose){mxLog("l1p iteration %d",k);}
			if(k>0){
				if(nmoc.iniSimplexMat.rows() || nmoc.iniSimplexMat.cols()){nmoc.iniSimplexMat.resize(0,0);}
				if(nmoc.statuscode==10){break;}
				if( !nmoc.estInfeas && nmoc.statuscode==0 ){
					if(verbose){mxLog("l1p solution found");}
					break;
				}
				if(nmoc.estInfeas){
					nmoc.rho *= 10.0;
					if(verbose){mxLog("penalty factor rho = %f",nmoc.rho);}
					nmoc.iniSimplexEdge = iniSimplexEdge;
				}
				else{ //It's making progress w/r/t the constraints, so re-initialize the simplex with a small edge:
					nmoc.iniSimplexEdge = 
						sqrt((nmoc.vertices[nmoc.n] - nmoc.vertices[0]).dot(nmoc.vertices[nmoc.n] - nmoc.vertices[0]));
					//It's a good idea to reduce the penalty coefficient if the algorithm is making progress.
					//That helps prevent it from stopping at a non-optimal point:
					nmoc.rho /= 5.0;
					if(verbose){mxLog("penalty factor rho = %f",nmoc.rho);}
				}
				if(fc->iterations >= maxIter){
					nmoc.statuscode = 4;
					if(verbose){mxLog("l1p algorithm ended with status 4");}
					break;
				}
			}
			nmoc.invokeNelderMead();
			fc->iterations += nmoc.itersElapsed;
			if(verbose){mxLog("total Nelder-Mead iterations elapsed: %d",fc->iterations);}
		}
	}
	else{
		nmoc.invokeNelderMead();
		fc->iterations = nmoc.itersElapsed;
	}
	
	if(validationRestart && nmoc.statuscode==0){
		NelderMeadOptimizerContext nmoc2(fc, this);
		nmoc2.verbose = verbose;
		nmoc2.maxIter = 2 * nmoc.n;
		nmoc2.iniSimplexType = 1;
		nmoc2.iniSimplexEdge = 
			sqrt((nmoc.vertices[nmoc.n] - nmoc.vertices[0]).dot(nmoc.vertices[nmoc.n] - nmoc.vertices[0]));
		nmoc2.centerIniSimplex = true;
		nmoc2.fit2beat = nmoc.bestfit;
		nmoc2.bignum = nmoc.bignum;
		nmoc2.est = nmoc.est;
		nmoc2.rho = nmoc.rho;
		nmoc2.addPenalty = nmoc.addPenalty;
		nmoc2.eqConstraintMthd = nmoc.eqConstraintMthd;
		nmoc2.ineqConstraintMthd = nmoc.ineqConstraintMthd;
		nmoc2.countConstraintsAndSetupBounds();
		nmoc2.invokeNelderMead();
		if(nmoc2.statuscode==10){
			fc->resetIterationError();
		}
		
		if(nmoc2.bestfit < nmoc.bestfit && (nmoc2.statuscode==0 || nmoc2.statuscode==4)){
			nmoc.bestfit = nmoc2.bestfit;
			nmoc.est = nmoc2.est;
			nmoc.estInfeas = nmoc2.estInfeas;
			if(nmoc2.statuscode==0){
				nmoc.fvals = nmoc2.fvals;
				nmoc.vertices = nmoc2.vertices;
				nmoc.vertexInfeas = nmoc2.vertexInfeas;
				nmoc.subcentroid = nmoc2.subcentroid;
				nmoc.eucentroidPrev = nmoc2.eucentroidPrev;
				nmoc.equality = nmoc2.equality;
				nmoc.inequality = nmoc2.inequality;
			}
			else if(Global->timedOut){ //i.e., if time ran out during the validation restart
				nmoc.statuscode = 4;
			}
		}
		
		//Not sure about this:
		fc->iterations += nmoc2.itersElapsed;
	}
	
	if(doPseudoHessian && (nmoc.statuscode==0 || nmoc.statuscode==4) && !nmoc.vertexInfeas.sum() && !nmoc.numEqC && !nmoc.addPenalty){
		nmoc.calculatePseudoHessian();
	}
	
	if(nmoc.estInfeas && nmoc.statuscode!=10){nmoc.statuscode = 3;}
	
	switch(nmoc.statuscode){
	case -1:
		mxThrow("unknown Nelder-Mead optimizer error");
		break;
	case 0:
		fc->setInform(INFORM_CONVERGED_OPTIMUM);
		break;
	case 3:
		fc->setInform(INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE);
		break;
	case 4:
		fc->setInform(INFORM_ITERATION_LIMIT);
		break;
	case 10:
		fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);
		break;
	}
	
	size_t i=0;
	Eigen::VectorXd xdiffs(nmoc.n);
	Eigen::MatrixXd fdiffs(nmoc.n,1);
	Eigen::MatrixXd Q(nmoc.n, nmoc.n);
	verticesOut.resize(nmoc.vertices.size(), nmoc.vertices[0].size());
	for(i=0; i < nmoc.vertices.size(); i++){
		verticesOut.row(i) = nmoc.vertices[i];
	}
	fvalsOut = nmoc.fvals;
	vertexInfeasOut = nmoc.vertexInfeas;
	for(i=0; i < size_t(nmoc.n); i++){
		fdiffs(i,0) = fabs(fvalsOut[i+1] - fvalsOut[0]);
	}
	fproxOut = fdiffs.array().maxCoeff();
	for(i=0; i < size_t(nmoc.n); i++){
		if(!nmoc.numEqC){
			Q.col(i) = verticesOut.row(i+1) - verticesOut.row(0);
			xdiffs[i] = (Q.col(i)).array().abs().maxCoeff();
		}
		else{
			xdiffs[i] = (verticesOut.row(i+1) - verticesOut.row(0)).array().abs().maxCoeff();
		}
	}
	xproxOut = xdiffs.array().maxCoeff();
	if(!nmoc.vertexInfeas.sum() && !nmoc.numEqC && !nmoc.addPenalty){
		Eigen::FullPivLU< Eigen::MatrixXd > luq(Q);
		if(luq.isInvertible()){
			Eigen::MatrixXd Qinv(nmoc.n, nmoc.n);
			Qinv = luq.inverse();
			//This is the "simplex gradient" of Kelley (1999):
			simplexGradient = Qinv.transpose() * fdiffs;
			if(verbose){mxPrintMat("simplex gradient: ",simplexGradient);}
		}
	}
	
	nmoc.finalize();
	
	fc->wanted |= FF_COMPUTE_BESTFIT;
	
	return;
}


void omxComputeNM::reportResults(FitContext *fc, MxRList *slots, MxRList *out){
	omxPopulateFitFunction(fitMatrix, out);
	
	MxRList output;
	SEXP pn, cn, cr, cc, cv, vrt, fv, vinf, fpm, xpm, phess, sg, bf;
	size_t i=0;
	
	if( fc->varGroup->vars.size() ){
		Rf_protect(pn = Rf_allocVector( STRSXP, fc->varGroup->vars.size() ));
		for(i=0; i < fc->varGroup->vars.size(); i++){
			SET_STRING_ELT( pn, i, Rf_mkChar(fc->varGroup->vars[i]->name) );
		}
		output.add("paramNames", pn);
	}
	if( fc->state->conListX.size() ){
		Rf_protect(cn = Rf_allocVector( STRSXP, fc->state->conListX.size() ));
		Rf_protect(cr = Rf_allocVector( INTSXP, fc->state->conListX.size() ));
		Rf_protect(cc = Rf_allocVector( INTSXP, fc->state->conListX.size() ));
		for(i=0; i < fc->state->conListX.size(); i++){
			SET_STRING_ELT( cn, i, Rf_mkChar(fc->state->conListX[i]->name) );
			INTEGER(cr)[i] = fc->state->conListX[i]->nrows;
			INTEGER(cc)[i] = fc->state->conListX[i]->ncols;
		}
		output.add("constraintNames", cn);
		output.add("constraintRows", cr);
		output.add("constraintCols", cc);
	}
	if( fc->constraintFunVals.size() ){
		Rf_protect(cv = Rf_allocVector( REALSXP, fc->constraintFunVals.size() ));
		memcpy( REAL(cv), fc->constraintFunVals.data(), sizeof(double) * fc->constraintFunVals.size() );
		output.add("constraintFunctionValues", cv);
	}
	if( verticesOut.rows() && verticesOut.cols() ){
		Rf_protect(vrt = Rf_allocMatrix( REALSXP, verticesOut.rows(), verticesOut.cols() ));
		memcpy( REAL(vrt), verticesOut.data(), sizeof(double) * verticesOut.rows() * verticesOut.cols() );
		output.add("finalSimplexMat", vrt);
	}
	if( fvalsOut.size() ){
		Rf_protect(fv = Rf_allocVector( REALSXP, fvalsOut.size() ));
		memcpy( REAL(fv), fvalsOut.data(), sizeof(double) * fvalsOut.size() );
		output.add("finalFitValues", fv);
	}
	if( vertexInfeasOut.size() ){
		Rf_protect(vinf = Rf_allocVector( INTSXP, vertexInfeasOut.size() ));
		memcpy( INTEGER(vinf), vertexInfeasOut.data(), sizeof(int) * vertexInfeasOut.size() );
		output.add("finalVertexInfeas", vinf); 
	}
	if( pseudohess.rows() && pseudohess.cols() ){
		Rf_protect(phess = Rf_allocMatrix( REALSXP, pseudohess.rows(), pseudohess.cols() ));
		memcpy( REAL(phess), pseudohess.data(), sizeof(double) * pseudohess.rows() * pseudohess.cols() );
		output.add("pseudoHessian", phess); 
	}
	if( simplexGradient.rows() && simplexGradient.cols() ){
		Rf_protect(sg = Rf_allocVector( REALSXP, simplexGradient.rows() ));
		memcpy( REAL(sg), simplexGradient.data(), sizeof(double) * simplexGradient.rows() );
		output.add("simplexGradient", sg); 
	}
	
	Rf_protect(fpm = Rf_allocVector(REALSXP, 1));
	//it would also work to do 'REAL(fpm)[0] = fproxOut;':
	memcpy( REAL(fpm), &fproxOut, sizeof(double) );
	output.add("rangeProximityMeasure", fpm);
	
	Rf_protect(xpm = Rf_allocVector(REALSXP, 1));
	memcpy( REAL(xpm), &xproxOut, sizeof(double) );
	output.add("domainProximityMeasure", xpm);
	
	Rf_protect(bf = Rf_allocVector(REALSXP, 1));
	memcpy( REAL(bf), &bestfitOut, sizeof(double) );
	output.add("penalizedFit", bf);
	
	slots->add("output", output.asR());
	return;
}

//-------------------------------------------------------

NelderMeadOptimizerContext::NelderMeadOptimizerContext(FitContext* _fc, omxComputeNM* _nmo)
	: fc(_fc), NMobj(_nmo), numFree(_fc->calcNumFree()),
	  subsidiarygoc(GradientOptimizerContext(_fc, 0L, GradientAlgorithm_Forward, 1L, 1e-5, _nmo))
{
	est.resize(numFree);
	copyParamsFromFitContext(est.data());
	statuscode = -1;
	addPenalty = false;
	rho = 1;
	checkRedundantEqualities = true;
}

void NelderMeadOptimizerContext::copyBounds()
{
	fc->copyBoxConstraintToOptimizer(solLB, solUB);
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
	if(!numEqC && !ineqConstraintMthd){eqConstraintMthd = 1;}
	if(!numIneqC){ineqConstraintMthd = 0;}
	equality.resize(numEqC);
	inequality.resize(numIneqC);
	
	fc->equality.resize(numEqC);
	
	//Check for redundant equality constraints, and warn if found:
	if(numEqC > 1 && checkRedundantEqualities){
		NldrMd_equality_functional eqf(this, fc);
		Eigen::MatrixXd ej(numEqC, numFree);
		ej.setConstant(NA_REAL);
		eqf(est, equality);
		fd_jacobian<true>(
			GradientAlgorithm_Central, 4, 1.0e-7,
			eqf, equality, est, ej);
		Eigen::MatrixXd ejt = ej.transpose();
		//mxPrintMat("ej: ",ej);
		Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qrj;
		qrj.compute(ejt);
		if(qrj.rank() < numEqC){
			Rf_warning(
				"counted %d equality constraints, but equality-constraint Jacobian is apparently rank %d at the start values; " 
				"Nelder-Mead will not work correctly unless equality constraints are linearly independent "
				"(this warning may be spurious if there are non-smooth equality constraints)", numEqC, qrj.rank()
			);
			checkRedundantEqualities = false;
		}
	}
	
	if(numEqC + numIneqC || eqConstraintMthd==3){
		subsidiarygoc.setEngineName("SLSQP");
		subsidiarygoc.ControlTolerance = 2 * Global->optimalityTolerance;
		subsidiarygoc.useGradient = true;
		subsidiarygoc.maxMajorIterations = Global->majorIterations;
		subsidiarygoc.setupSimpleBounds();
		//mxThrow("so far, so good");
	}
}

void NelderMeadOptimizerContext::copyParamsFromFitContext(double *ocpars)
{
	Eigen::Map<Eigen::VectorXd> vec(ocpars, numFree);
	fc->copyEstToOptimizer(vec);
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
	if( fc->outsideFeasibleSet() ){
		return(bignum);
	}
	else{
		if(fc->fit > bignum){bignum = 10 * fc->fit;}
		double fv = fc->fit;
		if(eqConstraintMthd==4 && addPenalty){
			int i;
			for(i=0; i < equality.size(); i++){
				fv += rho * fabs(equality[i]);
			}
			if(ineqConstraintMthd){
				for(i=0; i < inequality.size(); i++){
					fv += rho * fabs(inequality[i]);
				}
			}
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
	enforceBounds(x);
	checkNewPointInfeas(x, ifcr);
	if(!ifcr.sum()){
		infeas = 0L;
		fv = evalFit(x);
		if(fv==bignum){infeas=1L;}
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
			infeas = 1L;
			fv = bignum;
			break;
		case 3:
			gdfsIter = 0;
			tentativpt = x;
			if (NMobj->verbose >= 3) {
				mxPrintMat("tentative point", tentativpt);
			}
			omxInvokeSLSQPfromNelderMead(this, x);
			if (NMobj->verbose >= 3) {
				mxPrintMat("replacement point", x);
			}
			checkNewPointInfeas(x, ifcr);
			if(!ifcr.sum()){
				infeas = 0L;
				fv = evalFit(x);
				if(fv==bignum){infeas=1L;}
				return;
			}
			else{
				fv = bignum;
				infeas = 1L;
				return;
			}
			//mxThrow("'GDsearch' Not Yet Implemented");
		case 4:
			fv = evalFit(x);
			infeas = 1L;
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
void NelderMeadOptimizerContext::evalNewPoint(Eigen::VectorXd &newpt, Eigen::VectorXd oldpt, double &fv, int &newInfeas, int oldInfeas)
{
	Eigen::Vector2i ifcr;
	enforceBounds(newpt);
	checkNewPointInfeas(newpt, ifcr);
	if(!ifcr.sum()){
		newInfeas = 0L;
		fv = (evalFit(newpt));
		if(fv==bignum){newInfeas=1L;}
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
			else{
				int i;
				for(i=1; i <= NMobj->backtrackCtrl2; i++){
					ifcr.setZero();
					newpt = oldpt + NMobj->backtrackCtrl1*(newpt - oldpt);
					enforceBounds(newpt);
					checkNewPointInfeas(newpt, ifcr);
					if(!ifcr.sum()){
						newInfeas = 0L;
						fv = evalFit(newpt);
						if(fv==bignum){continue;}
						return;
					}
				}
				fv = bignum;
				newInfeas = 1L;
				return;
			}
		case 3:
			gdfsIter = 0;
			tentativpt = newpt;
			if (NMobj->verbose >= 3) {
				mxPrintMat("tentative point", tentativpt);
			}
			omxInvokeSLSQPfromNelderMead(this, newpt);
			if (NMobj->verbose >= 3) {
				mxPrintMat("replacement point", newpt);
			}
			checkNewPointInfeas(newpt, ifcr);
			if(!ifcr.sum()){
				newInfeas = 0L;
				fv = evalFit(newpt);
				if(fv==bignum){newInfeas=1L;}
				return;
			}
			else{
				fv = bignum;
				newInfeas = 1L;
				return;
			}
			//mxThrow("'GDsearch' Not Yet Implemented");
		case 4:
			fv = evalFit(newpt);
			newInfeas = 1L;
			return;
		}
	}
	else if(ifcr[0]){
		fv = bignum;
		newInfeas = 1L;
		return;
	}
}

void NelderMeadOptimizerContext::jiggleCoord(Eigen::VectorXd &xin, Eigen::VectorXd &xout, double scal){
	double a,b;
	int i;
	BorrowRNGState grs;
	for(i=0; i < xin.size(); i++){
		b = Rf_runif(1.0-scal,1.0+scal);
		a = Rf_runif(0.0-scal,0.0+scal);
		xout[i] = b*xin[i] + a;
	}
}

//TODO: make the different parts of the printing subject to different verbose levels
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

//Want to pass startpt as value, not reference:
void NelderMeadOptimizerContext::initializeSimplex(Eigen::VectorXd startpt, double edgeLength, bool isRestart)
{
	if(verbose){mxLog("(re-)initializing simplex");}
	int i=0;
	Eigen::VectorXd xin, xout, newpt, oldpt;
	if(iniSimplexMat.rows() && iniSimplexMat.cols() && !isRestart){
		Eigen::MatrixXd SiniSupp, iniSimplexMat2;
		Eigen::VectorXi paramMap(numFree);
		if(iniSimplexMat.cols() != numFree){
			mxThrow("'iniSimplexMat' has %d columns, but %d columns expected",iniSimplexMat.cols(), numFree);
		}
		if( int(NMobj->iniSimplexColnames.size()) != numFree){
			mxThrow("'iniSimplexMat' has %d column names, but %d column names expected", int(NMobj->iniSimplexColnames.size()), numFree);
		}
		if(iniSimplexMat.rows()>n+1){
			Rf_warning("'iniSimplexMat' has %d rows, but %d rows expected; extraneous rows will be ignored",iniSimplexMat.rows(), n+1);
			iniSimplexMat.conservativeResize(n+1,numFree);
		}
		iniSimplexMat2.resize(iniSimplexMat.rows(), numFree);
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
			mxThrow("error in mapping column names of 'iniSimplexMat' to free-parameter labels");
		}
		for(i=0; i < iniSimplexMat.cols(); i++){
			iniSimplexMat2.col(paramMap[i]) = iniSimplexMat.col(i);
		}
		if(iniSimplexMat.rows()<n+1){
			Rf_warning("'iniSimplexMat' has %d rows, but %d rows expected; omitted rows will be generated randomly",iniSimplexMat.rows(),n+1);
			SiniSupp.resize(n + 1 - iniSimplexMat.rows(), numFree);
			xin=iniSimplexMat2.row(0);
			for(i=0; i<SiniSupp.rows(); i++){
				xout=SiniSupp.row(i);
				jiggleCoord(xin, xout, edgeLength/4.0);
				SiniSupp.row(i) = xout;
			}
		}
		for(i=0; i < iniSimplexMat.rows(); i++){
			vertices[i] = iniSimplexMat2.row(i);
		}
		if(SiniSupp.rows()){
			for(i=0; i<SiniSupp.rows(); i++){
				vertices[i+iniSimplexMat.rows()] = SiniSupp.row(i);
			}
		}
	}
	else{
		double k = (double) n;
		double shhp = edgeLength*(1.0/k/sqrt(2.0))*(-1.0 + k + sqrt(1.0+k));
		double shhq = edgeLength*(1.0/k/sqrt(2.0))*(sqrt(1.0+k)-1.0);
		Eigen::VectorXd xu, xd;
		double fu=0, fd=0;
		int badu=0, badd=0;
		switch(iniSimplexType){
		case 1:
			vertices[0] = startpt;
			if(n==numFree){
				for(i=1; i<n+1; i++){
					vertices[i].setConstant(numFree,shhq);
					vertices[i][i-1] = shhp;
					vertices[i] += startpt;
				}
			}
			else{
				for(i=1; i<n+1; i++){
					vertices[i].setConstant(numFree,shhq);
					vertices[i] += startpt;
				}
				int j=1;
				for(i=0; i<numFree; i++){
					vertices[j%(n+1)][i] += (shhp - shhq);
					j++;
					if(j==n+1){j = 1;}
				}
			}
			break;
		case 2:
			vertices[0] = startpt;
			if(n==numFree){
				for(i=1; i<n+1; i++){
					vertices[i] = startpt;
					vertices[i][i-1] += edgeLength;
				}
			}
			else{
				for(i=1; i<n+1; i++){
					vertices[i] = startpt;
				}
				int j=1;
				for(i=0; i<numFree; i++){
					vertices[j%(n+1)][i] += edgeLength;
					j++;
					if(j==n+1){j = 1;}
				}
			}
			break;
		case 3:
			//TODO: this could be even smarter if it also figured out different edge lengths 
			//to account for different scaling of the free parameters:
			if(n==numFree){
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
			}
			else{
				vertices[0] = startpt;
				evalFirstPoint(vertices[0], fvals[0], vertexInfeas[0]);
				for(i=1; i<n+1; i++){
					vertices[i] = startpt;
				}
				int j=1;
				for(i=0; i<numFree; i++){
					xu = vertices[j%(n+1)];
					xu[i] += edgeLength;
					xd = vertices[j%(n+1)];
					xd[i] -= edgeLength;
					evalNewPoint(xu, vertices[0], fu, badu, vertexInfeas[0]);
					evalNewPoint(xd, vertices[0], fd, badd, vertexInfeas[0]);
					vertices[j%(n+1)] = fu<fd ? xu : xd;
					fvals[j%(n+1)] = fu<fd ? fu : fd;
					vertexInfeas[j%(n+1)] = fu<fd ? badu : badd;
					j++;
					if(j==n+1){j = 1;}
				}
				if(verbose){printProblemState();}
				return;
			}
		case 4:
			vertices[0] = startpt;
			for(i=1; i<n+1; i++){
				vertices[i].setZero(numFree);
				jiggleCoord(vertices[0],vertices[i],edgeLength/4.0);
			}
			break;
		}
		if(centerIniSimplex && !isRestart){
			eucentroidCurr.setZero(numFree);
			for(i=0; i<n+1; i++){
				eucentroidCurr += vertices[i] / (n+1.0);
			}
			for(i=0; i<n+1; i++){
				vertices[i] += startpt - eucentroidCurr;
			}
		}
	}
	//Now evaluate each vertex:
	evalFirstPoint(vertices[0], fvals[0], vertexInfeas[0]);
	for(i=1; i<n+1; i++){
		evalNewPoint(vertices[i], vertices[0], fvals[i], vertexInfeas[i], vertexInfeas[0]);
	}
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
			evalNewPoint(xe, xr, fe, bade, badr);
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
				//This case is considered a failed contraction:
				failedContraction = true;
				if(verbose){mxLog("outside contraction failed and shrinks are switched off...");}
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
				if(verbose){mxLog("inside contraction failed and shrinks are switched off...");}
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
		if(fprox < NMobj->fTolProx && fvals[0] < fit2beat){
			statuscode = 0;
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
		if(xprox < NMobj->xTolProx && fvals[0] < fit2beat){
			statuscode = 0;
			return(true);
		}
	}
	if(itersElapsed >= maxIter || isErrorRaised()){
		statuscode = 4;
		return(true);
	}
	return(false);
}


bool NelderMeadOptimizerContext::checkProgress(){
	//TODO: use const defined in Rmath.h:
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


void NelderMeadOptimizerContext::invokeNelderMead(){
	n = numFree - numEqC;
	vertices.resize(n+1);
	fvals.resize(n+1);
	vertexInfeas.resize(n+1);
	subcentroid.resize(numFree);
	eucentroidCurr.resize(numFree);
	initializeSimplex(est, iniSimplexEdge, false);
	if( (vertexInfeas.sum()==n+1 && eqConstraintMthd != 4) || (fvals.array()==bignum).all()){
		fc->recordIterationError("initial simplex is not feasible; specify it differently, try different start values, or use mxTryHard()");
		statuscode = 10;
		return;
	}
	fullSort();
	needFullSort=false;
	bool needRestart = false;
	bool stopflag=false;
	itersElapsed = 0;
	restartsUsed = 0;
	
	//Loop is: sort, check convergence, check progress, transform;
	do{
		if(verbose){
			mxLog("Nelder-Mead iteration %d / %d",itersElapsed,maxIter);
		}
		fc->resetOrdinalRelativeError();
		if(itersElapsed){
			//Order the vertices by fit value:
			if(needFullSort){fullSort();}
			else{fastSort();}
			
			stopflag = checkConvergence();
			if(stopflag){
				break;
			}
			
			needRestart = checkProgress();
			if(needRestart){
				initializeSimplex(vertices[0], sqrt((vertices[0]-vertices[1]).dot(vertices[0]-vertices[1])), true);
				needRestart = false;
				restartsUsed++;
				needFullSort = true;
				itersElapsed++;
				continue;
			}
		}
		
		simplexTransformation();
		
		eucentroidPrev = eucentroidCurr;
		itersElapsed++;
		Global->reportProgress("MxComputeNelderMead", fc);
	} while (!stopflag);
	
	est = vertices[0];
	bestfit = fvals[0];
	estInfeas = vertexInfeas[0];
	
	double centFit;
	int centInfeas;
	evalNewPoint(subcentroid, vertices[0], centFit, centInfeas, vertexInfeas[0]);
	if(centFit < bestfit && !centInfeas){
		est = subcentroid;
		bestfit = centFit;
		estInfeas = 0;
		
	}
	evalNewPoint(eucentroidCurr, vertices[0], centFit, centInfeas, vertexInfeas[0]);
	if(centFit < bestfit && !centInfeas){
		est = eucentroidCurr;
		bestfit = centFit;
		estInfeas = 0;
	}	
	
	//if(estInfeas){statuscode = 3;}
	
	if(verbose){mxPrintMat("solution?",est);}
	
}


void NelderMeadOptimizerContext::calculatePseudoHessian()
{
	int numpts = (n+1)*(n+2)/2;
	bool canDoAnalyt=true;
	int i, j, k, pminInfeas;
	double a0, pminfit;
	NMobj->pseudohess.resize(n, n);
	NMobj->phpts.resize(numpts, n);
	NMobj->phFvals.resize(numpts, 1);
	NMobj->phInfeas.resize(numpts);
	Eigen::VectorXd currpt(n);
	Eigen::VectorXd currpt2(n);
	Eigen::VectorXi jvec(numpts);
	Eigen::VectorXi kvec(numpts);
	Eigen::VectorXd a(n), pmin(n);
	Eigen::MatrixXd B(n,n), Q(n, n);
	
	NMobj->pseudohess.setZero(n, n);
	NMobj->phpts.setZero(numpts, n);
	NMobj->phFvals.setZero(numpts, 1);
	NMobj->phInfeas.setZero(numpts);
	
	for(i=0; i<n; i++){
		Q.col(i) = vertices[i+1] - vertices[0];
	}
	Eigen::FullPivLU< Eigen::MatrixXd > luq(Q);
	
	
	for(i=0; i<n+1; i++){
		NMobj->phpts.row(i) = vertices[i];
		NMobj->phFvals(i,0) = fvals[i];
		NMobj->phInfeas[i] = 0; //<--Assuming that this function is not called if any vertices are infeasible.
		kvec[i] = -1;
		jvec[i] = -1;
	}
	
	i=n+1;
	for(j=0; j<n; j++){
		for(k=j+1; k<n+1; k++){
			jvec[i] = j;
			kvec[i] = k;
			currpt = (vertices[j] + vertices[k])/2;
			currpt2 = currpt;
			evalNewPoint(currpt, vertices[j], NMobj->phFvals(i,0), NMobj->phInfeas[i], 0);
			if(NMobj->phInfeas[i]){
				//TODO: export a message about the pseudohessian for the user
				NMobj->pseudohess.resize(0,0);
				NMobj->phpts.resize(0,0);
				NMobj->phFvals(0,0);
				NMobj->phInfeas.resize(0);
				return;
			}
			else if(NMobj->phFvals(i,0) < bestfit){
				est = currpt;
				bestfit = NMobj->phFvals(i,0);
				estInfeas = 0;
			}
			//We can't use Nelder & Mead's analytic solution if the midpoints of the edges aren't actually such:
			if( (currpt.array() != currpt2.array()).any() ){
				canDoAnalyt = false;
			}
			NMobj->phpts.row(i) = currpt;
			i++;
		}
	}
	
	if(canDoAnalyt && luq.isInvertible()){
		if(verbose){mxLog("analytically calculating pseudoHessian");}
		a0 = fvals[0];
		for(i=0; i<n; i++){
			a[i] = 2*NMobj->phFvals(i+(n+1),0) - (fvals[i+1] + 3*a0)/2;
			B(i,i) = 2*( fvals[i+1] + a0 - 2*NMobj->phFvals(i+(n+1),0) );
		}
		for(i=n+n+1; i<numpts; i++){
			if(jvec[i] == kvec[i]){continue;}
			B(jvec[i]-1,kvec[i]-1) = 2*( NMobj->phFvals(i,0) + a0 - NMobj->phFvals(jvec[i]+(n+1)-1, 0) - 
				NMobj->phFvals(kvec[i]+(n+1)-1, 0) );
			B(kvec[i]-1,jvec[i]-1) = B(jvec[i]-1,kvec[i]-1);
		}
		Eigen::FullPivLU< Eigen::MatrixXd > lub(B);
		if(lub.isInvertible()){
			pmin = vertices[0] - (Q * lub.inverse() * a);
			evalNewPoint(pmin, vertices[0], pminfit, pminInfeas, vertexInfeas[0]);
			if(pminfit<bestfit && !pminInfeas){
				est = pmin;
				bestfit = pminfit;
				estInfeas = 0;
			}
		}
		Eigen::MatrixXd Qinv = luq.inverse();
		//NMobj->pseudohess = luq.inverse().transpose() * B * luq.inverse();
		NMobj->pseudohess = Qinv.transpose() * B * Qinv;
	}
	else{
		if(verbose){mxLog("numerically calculating pseudoHessian");}
		Eigen::MatrixXd X(numpts, numpts), polynomb(numpts,1), Binv;
		for(i=0; i<numpts; i++){
			X(i,0) = 1;
		}
		for(i=0; i<n; i++){
			X.col(i+1) = NMobj->phpts.col(i);
		}
		i=n+1;
		for(j=0; j<n; j++){
			for(k=j; k<n; k++){
				X.col(i) = (NMobj->phpts.col(j).array() * NMobj->phpts.col(k).array());
				i++;
			}
		}
		polynomb.setZero(numpts,1);
		
		Eigen::ColPivHouseholderQR< Eigen::MatrixXd > qrx(X);
		if(qrx.info() != Eigen::Success){
			NMobj->pseudohess.resize(0,0);
			NMobj->phpts.resize(0,0);
			NMobj->phFvals(0,0);
			NMobj->phInfeas.resize(0);
			return;
		}
		polynomb = qrx.solve(NMobj->phFvals);
		if(verbose){mxPrintMat("polynomial coefficients:",polynomb);}
		
		i=n+1;
		for(j=0; j<n; j++){
			for(k=j; k<n; k++){
				NMobj->pseudohess(j,k) = polynomb(i,0);
				if(j != k){NMobj->pseudohess(k,j) = polynomb(i,0);}
				i++;
			}
		}
		Eigen::FullPivLU< Eigen::MatrixXd > lub(NMobj->pseudohess);
		if(lub.isInvertible()){
			Binv = lub.inverse();
			for(i=0; i<n; i++){
				a[i] = polynomb(i+1,0);
			}
			pmin = vertices[0] - (Binv * a);
			evalNewPoint(pmin, vertices[0], pminfit, pminInfeas, vertexInfeas[0]);
			if(pminfit<bestfit && !pminInfeas){
				est = pmin;
				bestfit = pminfit;
				estInfeas = pminInfeas;
			}
		}
		NMobj->Xout = X;
	}
	if(verbose){
		mxPrintMat("pseudoHessian is ", NMobj->pseudohess);
	}
	return;
}


void NelderMeadOptimizerContext::finalize()
{
	//The omxComputeNM object stows the possibly penalized fit value; the FitContext here recomputes the unpenalized fit value, at the
	//best parameter values:
	NMobj->bestfitOut = bestfit;
	copyParamsFromOptimizer(est,fc);
	ComputeFit(engineName, NMobj->fitMatrix, FF_COMPUTE_FIT, fc);
	/*Doing this here ensures (1) that the fit has just been freshly evaluated at the solution, (2) that this check is done as part of the
	MxComputeNelderMead step (necessary for bootstrapping), and (3) that Nelder-Mead reports status code 3 for solutions that violate 
	MxConstraints, and status code 10 for	all other kinds of infeasible solutions:*/
	if(!fc->insideFeasibleSet() && (statuscode==0 || statuscode==4)){fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);}
	
	omxState *st = fc->state;
	int ineqType = omxConstraint::LESS_THAN;
	int cur=0, j=0;
	Eigen::VectorXd cfv(numEqC + numIneqC);
	
	for (j=0; j < int(st->conListX.size()); j++) {
		omxConstraint &con = *st->conListX[j];
		con.refreshAndGrab(fc, (omxConstraint::Type) ineqType, &cfv(cur));
		cur += con.size;
	}
	
	fc->constraintFunVals = cfv;
}


double nmgdfso(unsigned n, const double *x, double *grad, void *f_data)
{
	NelderMeadOptimizerContext *nmoc = (NelderMeadOptimizerContext *) f_data;
	nlopt_opt opt = (nlopt_opt) nmoc->extraData;
	unsigned i;
	double ssq=0, currdiff=0;
	if(grad){
		if(nmoc->gdfsIter >= nmoc->subsidiarygoc.maxMajorIterations){
			nlopt_force_stop(opt);
		}
		(nmoc->gdfsIter)++;
	}
	for(i=0; i < n; i++){
		currdiff = x[i] - nmoc->tentativpt[i];
		if(grad){grad[i] = 2*currdiff;}
		currdiff *= currdiff;
		ssq += currdiff;
	}
	return(ssq);
}
