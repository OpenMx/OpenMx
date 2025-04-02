 /*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
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

#include "omxFitFunction.h"
#include "omxGREMLExpectation.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Rmath.h>
#include "Compute.h"
#include "EnableWarnings.h"
#include "finiteDifferences.h"

struct omxGREMLFitState : omxFitFunction {
	//TODO(?): Some of these members might be redundant with what's stored in the FitContext,
	//and could therefore be cut
	omxMatrix *y, *X, *cov, *invcov, *means, *residual;
	std::vector< omxMatrix* > dV, dyhat;
	/*`dVnames` provides the mapping from elements of `dV` to free parameters when `dyhat` is empty and `dV` is not; 
	`dVnames` also provides the mapping to free parameters from rows and rows/columns of the augmentation's gradient and Hessian;
	`dNames` provides the mapping from elements of `dyhat` and (if applicable) `dV` whenever `dyhat` has nonzero length:*/
	std::vector< const char* > dVnames, dNames;
	std::vector<bool> indyAlg, indyAlg2; //<--Will keep track of which algebras don't get marked dirty after dropping cases.
	std::vector<int> origdVdim, origdyhatdim;
	std::vector<bool> didUserGivedV, didUserGivedyhat;
	void dVupdate(FitContext *fc); //<--Also updates dyhat.
	void dVupdate_final(); //<--Also updates dyhat.
	int origVdim, numExplicitFreePar, origyhatdim;
	int dVnameslength, dyhatnameslength; //<--The lengths of frontend arguments `dV` and `dyhat`, as originally provided by the user.
	int dNamesLength; //<--Has length equal to the number of unique parameter labels present in the names of `dV` and/or `dyhat`.
	/*Each has length either equal to zero, if no such user-provided derivatives; or to the number of explicit free parameters,
	if any such user-provided derivatives:*/
	int dyhatlength, dVlength; 
	int parallelDerivScheme, derivType, oldWantHess, infoMatType, parallelDerivSchemeFromFrontend;
	bool usingGREMLExpectation, doREML, didUserSpecifyParallelDerivScheme;
	bool didUserProvideYhat; //<--Value of `true` means that doREML is FALSE *and* the user provided a non-empty, valid name for 'yhat'.
	double nll, REMLcorrection;
	Eigen::VectorXd gradient;
	Eigen::MatrixXd infoMat; //<--The Average Information matrix or the Expected Information matrix, as the case may be.
	FreeVarGroup *varGroup;
	std::vector<int> gradMap;
	void buildParamMap(FreeVarGroup *newVarGroup);
	std::vector< Eigen::VectorXi > rowbins, AIMelembins;
	void planParallelDerivs(int nThreadz, int wantHess, int Vrows);
	omxMatrix *aug, *augGrad, *augHess;
	std::vector<int> dAugMap;
	double pullAugVal(int thing, int row, int col);
	void recomputeAug(int thing, FitContext *fc);
	JacobianGadget jg;

	//Foundation for separate functions that compute fitfunction derivatives from dV:
	template <typename T1, typename T2, typename T3>
	void gradientAndAIM1(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv);

	template <typename T1, typename T2, typename T3>
	void gradientAndAIM2(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv);

	template <typename T1, typename T2, typename T3>
	void gradientAndAIM3(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv);

	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM1(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv);

	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM2(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv);

	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM3(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv);
	
	//Only implemented for `didUserProvideYhat==true` case; impractically slow; do not use:
	template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
	void gradientAndOIM1(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
			double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv, Eigen::MatrixBase<T5> &u_VinvResid,
			Eigen::MatrixBase<T6> &u_VinvResidResidT);
	
	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM1_yhat(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
			double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
			Eigen::MatrixBase<T4> &u_VinvResidResidT);
	
	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM2_yhat(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
			double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
			Eigen::MatrixBase<T4> &u_VinvResidResidT);
	
	template <typename T1, typename T2, typename T3, typename T4>
	void gradientAndEIM3_yhat(
			int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
			double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
			Eigen::MatrixBase<T4> &u_VinvResidResidT);

	template <typename T1, typename T2>
	void crude_numeric_dV(
		FitContext *u_fc, Eigen::MatrixBase<T1> &u_curEst, Eigen::MatrixBase<T2> &dV_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId);
	
	template <typename T1, typename T2>
	void crude_numeric_dyhat(
			FitContext *u_fc, Eigen::MatrixBase<T1> &u_curEst, Eigen::MatrixBase<T2> &dyhat_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId);

	omxGREMLFitState() : jg(1) {};
	virtual void init() override;
	virtual void compute2(int want, FitContext *fc) override;
	virtual void populateAttr(SEXP algebra) override;
};

omxFitFunction *omxInitGREMLFitFunction()
{ return new omxGREMLFitState; }

void omxGREMLFitState::init()
{
	auto *oo = this;
	auto *newObj = this;

  if(OMX_DEBUG) { mxLog("Initializing GREML fitfunction."); }

  jg.setAlgoOptions(GradientAlgorithm_Forward, 2, 1e-4);
  oo->units = FIT_UNITS_MINUS2LL;
  oo->canDuplicate = true;

  omxState* currentState = expectation->currentState;
  newObj->usingGREMLExpectation = (strcmp(expectation->name, "MxExpectationGREML")==0 ? true : false);
  if(!newObj->usingGREMLExpectation){
    //Maybe someday GREML fitfunction could be made compatible with another expectation, but not at present:
    mxThrow("GREML fitfunction is currently only compatible with GREML expectation");
  }
  omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);

  newObj->y = omxGetExpectationComponent(expectation, "y");
  newObj->cov = omxGetExpectationComponent(expectation, "cov");
  newObj->invcov = omxGetExpectationComponent(expectation, "invcov");
  newObj->X = omxGetExpectationComponent(expectation, "X");
  newObj->means = omxGetExpectationComponent(expectation, "means");
  newObj->origVdim = oge->origVdim;
  newObj->nll = 0;
  newObj->REMLcorrection = 0;
  newObj->parallelDerivScheme = 0;
  newObj->parallelDerivSchemeFromFrontend = 0;
  newObj->didUserSpecifyParallelDerivScheme = false;
  newObj->varGroup = NULL;
  newObj->augGrad = NULL;
  newObj->augHess = NULL;
  newObj->dVlength = 0;
  newObj->dVnameslength = 0;
  newObj->dNamesLength = 0;
  newObj->dyhatlength = 0;
  newObj->oldWantHess = 0;
  newObj->doREML = oge->doREML;
  newObj->didUserProvideYhat = oge->didUserProvideYhat;
  newObj->hessianAvailable = false;

	//autoDerivType:
  {
  	ProtectedSEXP adt(R_do_slot(rObj, Rf_install("autoDerivType")));
  	if(strEQ(CHAR(Rf_asChar(adt)),"semiAnalyt")){derivType = 1;}
  	else if(strEQ(CHAR(Rf_asChar(adt)),"numeric")){derivType = 0;}
  	else{mxThrow("unrecognized character string provided for GREML fitfunction 'autoDerivType'");}
  }
  /*if(derivType==1 && didUserProvideYhat){
  	Rf_warning("use of semi-analytic derivatives with 'yhat' is Not Yet Implemented; numeric derivatives will be used instead");
  }*/

  //infoMatType:
  {
  	ProtectedSEXP imt(R_do_slot(rObj, Rf_install("infoMatType")));
  	if(strEQ(CHAR(Rf_asChar(imt)),"expected")){infoMatType = 1;}
  	else if(strEQ(CHAR(Rf_asChar(imt)),"average")){infoMatType = 0;}
  	else{mxThrow("unrecognized character string provided for GREML fitfunction 'infoMatType'");}
  }

  //Augmentation:
  newObj->aug = 0;
  if (R_has_slot(rObj, Rf_install("aug"))) {
	  ProtectedSEXP Raug(R_do_slot(rObj, Rf_install("aug")));
	  if(Rf_length(Raug)){
		  int* augint = INTEGER(Raug);
		  newObj->aug = omxMatrixLookupFromStateByNumber(augint[0], currentState);
	  }
  }
  
  //Parallel-derivative scheme:
  if (R_has_slot(rObj, Rf_install(".parallelDerivScheme"))) {
  	ProtectedSEXP R_pds(R_do_slot(rObj, Rf_install(".parallelDerivScheme")));
  	parallelDerivSchemeFromFrontend = Rf_asInteger(R_pds);
  	if(parallelDerivSchemeFromFrontend!=0){
  		didUserSpecifyParallelDerivScheme = true;
  		//`parallelDerivSchemeFromFrontend` will get assigned to `parallelDerivScheme` during the call to `planParallelDerivs()`
  		if(parallelDerivSchemeFromFrontend<0 || parallelDerivSchemeFromFrontend>3){
  			Rf_warning("`.parallelDerivScheme` not equal to 1, 2, or 3");
  		}
  	}
  }

	/*
	 * Possible TODO: If we want GREML fitfunction to work properly in scenarios where the freeVarGroup changes after the
	 * omxGREMLFitState has been initialized, then all code from here to the end of this function would need to be moved
	 * to buildParamMap().
	*/
	
	//Stuff common to all derivatives:
	if (R_has_slot(rObj, Rf_install("dNames"))) {
		ProtectedSEXP RdNames(R_do_slot(rObj, Rf_install("dNames")));
		newObj->dNamesLength = Rf_length(RdNames);
		newObj->dNames.resize(newObj->dNamesLength);
		if(newObj->dNamesLength){
			if(OMX_DEBUG) { mxLog("Processing common derivative labels."); }
			for(int i=0; i < newObj->dNamesLength; i++){
				SEXP elem;
				{ScopedProtect p3(elem, STRING_ELT(RdNames, i));
					newObj->dNames[i] = CHAR(elem);}
			}
		}
	}
	
  //Derivatives of V:
  if (R_has_slot(rObj, Rf_install("dV"))) {
  	ProtectedSEXP RdV(R_do_slot(rObj, Rf_install("dV")));
  	ProtectedSEXP RdVnames(R_do_slot(rObj, Rf_install("dVnames")));
  	newObj->dVlength = Rf_length(RdV);
  	newObj->dVnameslength = Rf_length(RdVnames);
  	newObj->dV.resize(newObj->dVlength);
  	newObj->dVnames.resize(newObj->dVnameslength);
  	newObj->origdVdim.resize(newObj->dVlength);
  	if(newObj->dVlength){
  		/*if(!newObj->usingGREMLExpectation){
  			//Probably best not to allow use of dV if we aren't sure means will be calculated GREML-GLS way:
  			mxThrow("derivatives of 'V' matrix in GREML fitfunction only compatible with GREML expectation");
  		}*/
  		/*if(didUserProvideYhat){
  			Rf_warning("derivatives of 'V' matrix with 'yhat' are Not Yet Implemented; numeric derivatives will be used instead");
  		}*/
  		if(OMX_DEBUG) { mxLog("Processing derivatives of V."); }
  		int* dVint = INTEGER(RdV);
  		for(int i=0; i < newObj->dVlength; i++){
  			if(!ISNA(dVint[i])){ newObj->dV[i] = omxMatrixLookupFromStateByNumber(dVint[i], currentState); }
  		}
  		for(int i=0; i < newObj->dVnameslength; i++){
  			SEXP elem;
  			{ScopedProtect p3(elem, STRING_ELT(RdVnames, i));
  				newObj->dVnames[i] = CHAR(elem);}
  		}
  	}
  }
  
  /*Derivatives of yhat:*/
  if (R_has_slot(rObj, Rf_install("dyhat"))) {
  	ProtectedSEXP Rdyhat(R_do_slot(rObj, Rf_install("dyhat")));
  	ProtectedSEXP Rdyhatnames(R_do_slot(rObj, Rf_install("dyhatnames")));
  	newObj->dyhatlength = Rf_length(Rdyhat);
  	newObj->dyhat.resize(newObj->dyhatlength);
  	newObj->dyhatnameslength = Rf_length(Rdyhatnames);
  	newObj->origdyhatdim.resize(newObj->dyhatlength);
  	if(newObj->dyhatlength){
  		/*if(!newObj->usingGREMLExpectation){
  		 //Probably best not to allow use of dV if we aren't sure means will be calculated GREML-GLS way:
  		 mxThrow("derivatives of 'V' matrix in GREML fitfunction only compatible with GREML expectation");
  		 }*/
  		if(didUserProvideYhat){
  			//Rf_warning("derivatives of 'yhat' are Not Yet Implemented; numeric derivatives will be used instead");
  		}
  		if(OMX_DEBUG) { mxLog("Processing derivatives of yhat."); }
  		int* dyhatint = INTEGER(Rdyhat);
  		for(int i=0; i < newObj->dyhatlength; i++){
  			if(!ISNA(dyhatint[i])){ newObj->dyhat[i] = omxMatrixLookupFromStateByNumber(dyhatint[i], currentState); }
  			/*SEXP elem;
  			{ScopedProtect p3(elem, STRING_ELT(Rdyhatnames, i));
  				newObj->dyhatnames[i] = CHAR(elem);}*/
  		}
  	}
  }

  if(derivType==1 && !newObj->usingGREMLExpectation){
  	mxThrow("semi-analytic derivatives only compatible with GREML expectation");
  }
  
  newObj->rowbins.resize(Global->numThreads);
  newObj->AIMelembins.resize(Global->numThreads);

  if(newObj->dVlength || newObj->dyhatlength || derivType==1){
  	oo->hessianAvailable = true;
  	//^^^Gets changed to false in `buildParamMap()` if it turns out that derivType=0 and 0 < dVlength < numExplicitFreePar.
  	for(int i=0; i < newObj->dVlength; i++){
  		if(newObj->dV[i]){
  			/*Each dV must either (1) match the dimensions of V, OR (2) match the length of y if that is less than the
  			 dimension of V (implying downsizing due to missing observations):*/
  			if( ((newObj->dV[i]->rows == newObj->cov->rows)&&(newObj->dV[i]->cols == newObj->cov->cols)) ||
         ((newObj->y->cols < newObj->cov->rows)&&(newObj->dV[i]->rows == newObj->y->cols)&&
         (newObj->dV[i]->cols == newObj->y->cols)) ){
  				newObj->origdVdim[i] = newObj->dV[i]->rows;
  			}
  			else{
  				mxThrow("all derivatives of V must have the same dimensions as V");
  			}
  		}
  		else{ newObj->origdVdim[i] = 0; }
  	}
  	for(int i=0; i < newObj->dyhatlength; i++){
  		if(newObj->dyhat[i]){
  			/*Each dyhat must either (1) match the dimensions of yhat, OR (2) match the length of y if that is less than the
  			 dimension of yhat (implying downsizing due to missing observations):*/
  			if((newObj->dyhat[i]->size() == oge->yhatFromUser->size()) || 
        ((newObj->y->cols < oge->yhatFromUser->size())&&(newObj->dyhat[i]->size() == newObj->y->cols))){
  				newObj->origdyhatdim[i] = newObj->dyhat[i]->size();
  			}
  			else{
  				mxThrow("all derivatives of yhat must have the same length as yhat");
  			}
  		}
  		else{ newObj->origdyhatdim[i] = 0; }
  	}
  }
  	
  //Augmentation derivatives:
  if( (newObj->dVlength || derivType==1) && newObj->aug){
  	//^^^Ignore derivatives of aug unless aug itself and objective derivatives are supplied, or if aug is provided and derivType==1.
  	ProtectedSEXP RaugGrad(R_do_slot(rObj, Rf_install("augGrad")));
  	ProtectedSEXP RaugHess(R_do_slot(rObj, Rf_install("augHess")));
  	if(!Rf_length(RaugGrad)){
  		if(Rf_length(RaugHess)){
  			mxThrow("if argument 'augHess' has nonzero length, then argument 'augGrad' must as well");
  		}
  		else{
  			if(dVlength){
  				mxThrow("if arguments 'dV' and 'aug' have nonzero length, then 'augGrad' must as well");
  			}
  			if(dyhatlength){
  				mxThrow("if arguments 'dyhat' and 'aug' have nonzero length, then 'augGrad' must as well");	
  			}
  			mxThrow("if using semi-analytic derivatives and 'aug' has nonzero length, then 'augGrad' must as well");
  		}
  	}
  	else{
  		int* augGradint = INTEGER(RaugGrad);
  		newObj->augGrad = omxMatrixLookupFromStateByNumber(augGradint[0], currentState);
  		if(Rf_length(RaugHess)){
  			//Conformability of augGrad and augHess are checked later, during buildParamMap().
  			int* augHessint = INTEGER(RaugHess);
  			newObj->augHess = omxMatrixLookupFromStateByNumber(augHessint[0], currentState);
  		}
  		else{oo->hessianAvailable = false;}
  	}
  }
}

void omxGREMLFitState::compute2(int want, FitContext *fc)
 {
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;

 	//Recompute Expectation:
 	omxExpectationCompute(fc, expectation, NULL);

 	omxGREMLFitState *gff = this;
	auto *oo = this;

 	//Ensure that the pointer in the GREML fitfunction is directed at the right FreeVarGroup
 	//(not necessary for most compute plans):
 	if(fc && gff->varGroup != fc->varGroup){
 		gff->buildParamMap(fc->varGroup);
 	}

 	gff->recomputeAug(0, fc);

 	//Declare local variables used in more than one scope in this function:
 	const double Scale = fabs(Global->llScale); //<--absolute value of loglikelihood scale
 	const double NATLOG_2PI = 1.837877066409345483560659472811;	//<--log(2*pi)
 	Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(gff->y), gff->y->cols, 1);
 	Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(gff->invcov), gff->invcov->rows, gff->invcov->cols);
 	EigenMatrixAdaptor EigX(gff->X);
 	Eigen::MatrixXd P, Py, VinvResid;
 	double logdetV=0, logdetquadX=0, ytPy=0;

 	if(want & (FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 		if(gff->usingGREMLExpectation){
 			omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
 			
 			//Check that factorizations of V and the quadratic form in X succeeded:
 			if(oge->cholV_fail){
 				oo->matrix->data[0] = NA_REAL;
 				if (fc) fc->recordIterationError("expected covariance matrix is non-positive-definite");
 				return;
 			}
 			if(!didUserProvideYhat && oge->cholquadX_fail){
 				oo->matrix->data[0] = NA_REAL;
 				if (fc) fc->recordIterationError("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
 				return;
 			}
 			
 			//Log determinant of V:
 			logdetV = oge->logdetV_om->data[0];
 			
 			//Log determinant of quadX:
 			gff->REMLcorrection = 0;
 			if(doREML){
 				for(int i=0; i < gff->X->cols; i++){
 					logdetquadX += log(oge->cholquadX_vectorD[i]);
 				}
 				logdetquadX *= 2;
 				gff->REMLcorrection = Scale*0.5*logdetquadX;
 			}
 			
 			/*Finish computing fit (negative loglikelihood) if wanted.  P and Py will be needed later if analytic or semi-analytic derivatives in use;
 			 otherwise, extraneous calculations can be avoided:*/
 			if(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 				if(!didUserProvideYhat){ //<--Either doREML is true, or doREML is false but we're not using a user-supplied yhat.
 					P.setZero(gff->invcov->rows, gff->invcov->cols);
 					Eigen::MatrixXd Hatmat = EigX * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv;
 					subtractFromIdentityMatrixInPlace(Hatmat);
 					P.triangularView<Eigen::Lower>() = (Vinv.selfadjointView<Eigen::Lower>() * Hatmat).triangularView<Eigen::Lower>(); //P = Vinv * (I-Hatmat)
 					Py = P.selfadjointView<Eigen::Lower>() * Eigy;
 					P.triangularView<Eigen::Upper>() = P.triangularView<Eigen::Lower>().transpose();
 					if(want & FF_COMPUTE_FIT){
 						ytPy = (Eigy.transpose() * Py)(0,0);
 						if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
 						oo->matrix->data[0] = Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + ytPy) + Scale*gff->pullAugVal(0L,0,0);
 						oo->matrix->data[0] += doREML ? gff->REMLcorrection : 0;
 						gff->nll = oo->matrix->data[0];
 						if(OMX_DEBUG){mxLog("augmentation is %3.3f",gff->pullAugVal(0L,0,0));}
 					}
 				}
 				else{
 					VinvResid = Vinv.selfadjointView<Eigen::Lower>() * oge->residual;
 					if(want & FF_COMPUTE_FIT){
 						oo->matrix->data[0] = 
 							Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + (oge->residual.transpose() * VinvResid)(0,0) ) + 
 							Scale*gff->pullAugVal(0L,0,0);
 						gff->nll = oo->matrix->data[0];
 					}
 				}
 			}
 			/*ytPy can be calculated so that rate-limiting step is O(2kN^2), where k is the number of covariates,
 			 and N is the dimension of Vinv (and typically N>>k):*/
 			else{
 				if(!didUserProvideYhat){//<--Either doREML is true, or doREML is false but we're not using a user-supplied yhat.
 					ytPy = (( Eigy.transpose() * Vinv.selfadjointView<Eigen::Lower>() * Eigy ) -
 						( Eigy.transpose() * oge->XtVinv.transpose() * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv * Eigy ))(0,0);
 					if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
 					oo->matrix->data[0] = Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + ytPy) + Scale*gff->pullAugVal(0L,0,0);
 					oo->matrix->data[0] += doREML ? gff->REMLcorrection : 0;
 					gff->nll = oo->matrix->data[0];
 					if(OMX_DEBUG){mxLog("augmentation is %3.3f",gff->pullAugVal(0L,0,0));}
 					return; //<--Since only fit value is wanted.
 				}
 				else{
 					//Eigen::Map< Eigen::MatrixXd > Eigyhat(omxMatrixDataColumnMajor(gff->yhatFromUser), gff->yhatFromUser->size(), 1);
 					oo->matrix->data[0] = 
 						Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + (oge->residual.transpose() * Vinv.selfadjointView<Eigen::Lower>() * oge->residual)(0,0) ) +
 						Scale*gff->pullAugVal(0L,0,0);
 					gff->nll = oo->matrix->data[0];
 					return; //<--Since only fit value is wanted.
 				}
 			}
 		}
 		else{ //If not using GREML expectation, deal with means and cov in a general way to compute fit...
 			//(There is currently no codepath that can reach the code in this else conditional.)
 			//Declare locals:
 			EigenMatrixAdaptor yhat(gff->means);
 			EigenMatrixAdaptor EigV(gff->cov);
 			Eigen::MatrixXd Vinv2, quadX;
 			Eigen::LLT< Eigen::MatrixXd > cholV(gff->cov->rows);
 			Eigen::LLT< Eigen::MatrixXd > cholquadX(gff->X->cols);
 			Eigen::VectorXd cholV_vectorD, cholquadX_vectorD;

 			//Cholesky factorization of V:
 			cholV.compute(EigV);
 			if(cholV.info() != Eigen::Success){
 				omxRaiseErrorf("expected covariance matrix is non-positive-definite");
 				oo->matrix->data[0] = NA_REAL;
 				return;
 			}
 			//Log determinant of V:
 			cholV_vectorD = (( Eigen::MatrixXd )(cholV.matrixL())).diagonal();
 			for(int i=0; i < gff->X->rows; i++){
 				logdetV += log(cholV_vectorD[i]);
 			}
 			logdetV *= 2;

 			Vinv2 = cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse

 			quadX = EigX.transpose() * Vinv2 * EigX; //<--Quadratic form in X

 			cholquadX.compute(quadX); //<--Cholesky factorization of quadX
 			if(cholquadX.info() != Eigen::Success){
 				omxRaiseErrorf("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
 				oo->matrix->data[0] = NA_REAL;
 				return;
 			}
 			cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal();
 			for(int i=0; i < gff->X->cols; i++){
 				logdetquadX += log(cholquadX_vectorD[i]);
 			}
 			logdetquadX *= 2;
 			gff->REMLcorrection = Scale*0.5*logdetquadX;

 			//Finish computing fit:
 			oo->matrix->data[0] = gff->REMLcorrection + Scale*0.5*( ((double)gff->y->rows * NATLOG_2PI) + logdetV +
 				( Eigy.transpose() * Vinv2 * (Eigy - yhat) )(0,0));
 			gff->nll = oo->matrix->data[0];
 			return;
 		}
 	}

 	if(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 		//This part requires GREML expectation:
 		omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
 		
 		/*if(!doREML){
 			return;
 			//mxThrow("analytic and semi-analytic GREML derivatives are Not Yet Implemented for the `REML=FALSE` case");
 		}*/

 		//Recompute derivatives:
 		gff->dVupdate(fc);
 		gff->recomputeAug(1, fc);

 		//Declare local variables for this scope:
 		int nThreadz = Global->numThreads;
 		int wantHess = 0;
 		//GREMLSense sense(this, fc, oge->numcases2drop, oge->dropcase);

 		//Set up new HessianBlock:
 		HessianBlock *hb = new HessianBlock;
 		if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 			if(gff->dVlength < gff->numExplicitFreePar && derivType==0){
        if (gff->hessianAvailable) mxThrow("%s: sorry, i lied about whether analytic derivs were available", gff->name());
 				omxRaiseErrorf("when 'autoDerivType' is 'numeric', GREML fitfunction cannot compute information matrix without analytic derivatives of V with respect to EVERY free parameter");
 			}
 			hb->vars.resize(gff->numExplicitFreePar);
 			hb->mat.resize(gff->numExplicitFreePar, gff->numExplicitFreePar);
 			gff->recomputeAug(2, fc);
 			wantHess = 1;
 		}

 		if(oldWantHess!=wantHess){
 			parallelDerivScheme = 0;
 			oldWantHess = wantHess;
 		}

 		if(gff->parallelDerivScheme==0){gff->planParallelDerivs(nThreadz,wantHess,gff->cov->rows);}

		/*If datapoints need to be dropped due to missingness, and we want the analytic Hessian,
 		then we need to resize any derivatives of V that come from front-end MxAlgebras
		 ahead of time,	in order to assure thread-safety of the parallelized code for
 		evaluating the gradient and Hessian (AIM):*/
		if(oge->numcases2drop &&  wantHess){
			for(int i=0; i < numExplicitFreePar; i++){
				if(didUserGivedV[i] && gff->dV[i]->rows > Eigy.rows()){
					dropCasesFromAlgdV(gff->dV[i], oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[i]);
				}
			}
		}


 		//Begin parallelized evaluation of fitfunction derivatives:
 		if(didUserProvideYhat){
 			Eigen::MatrixXd VinvResidResidT = VinvResid * oge->residual.transpose();
 			switch(gff->parallelDerivScheme){
 			case 2: //bin by row
 				gradientAndEIM2_yhat(nThreadz, Eigy.rows(), fc, want, hb, oge, Scale, Eigy, Vinv, VinvResid, VinvResidResidT);
 				break;
 			case 3: //bin by cell
 				gradientAndEIM3_yhat(nThreadz, Eigy.rows(), fc, want, hb, oge, Scale, Eigy, Vinv, VinvResid, VinvResidResidT);
 				break;
 			default: //bin naively (which is perfectly adequate for gradient-only, or for a single thread)
 				gradientAndEIM1_yhat(nThreadz, Eigy.rows(), fc, want, hb, oge, Scale, Eigy, Vinv, VinvResid, VinvResidResidT);
 			break;
 			}
 		}
 		else{
 			switch(gff->parallelDerivScheme){
 			case 2: //bin by row
 				if(infoMatType==1){
 					gradientAndEIM2(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Eigy, Vinv);
 				}
 				else{
 					gradientAndAIM2(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Vinv);
 				}
 				break;
 			case 3: //bin by cell
 				if(infoMatType==1){
 					gradientAndEIM3(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Eigy, Vinv);
 				}
 				else{
 					gradientAndAIM3(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Vinv);
 				}
 				break;
 			default: //bin naively (which is perfectly adequate for gradient-only, or for a single thread)
 				if(infoMatType==1){
 					gradientAndEIM1(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Eigy, Vinv);
 				}
 				else{
 					gradientAndAIM1(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Vinv);
 				}
 				break;
 			}
 		}
 			//Assign upper triangle elements of infoMat to the HessianBlock:
 			if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 				for (size_t d1=0, h1=0; h1 < size_t(numExplicitFreePar); ++h1) {
 					for (size_t d2=0, h2=0; h2 <= h1; ++h2) {
 						hb->mat(d2,d1) = gff->infoMat(h2,h1);
 						++d2;
 					}
 					++d1;
 				}
 				fc->queue(hb);
 			}
			else{
				delete hb;
			}
 	}
 	return;
 }


//Bin "naively"; Hessian is average information matrix:
template <typename T1, typename T2, typename T3>
void omxGREMLFitState::gradientAndAIM1(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
    double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
		mxLog("Now beginning `gradientAndAIM1()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, j=0, t1=0, t2=0, a1=0, a2=0;
		double tr=0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = threadID * numExplicitFreePar / u_nThreadz;
		int iend = (threadID+1) * numExplicitFreePar / u_nThreadz;
		if(threadID == u_nThreadz-1){iend = numExplicitFreePar;}
		for(i=istart; i < iend; i++){
			tr=0;
			t1 = gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
				double *ptrToMatrix1=0;
				Eigen::MatrixXd filteredCopy1;
				a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[i] = t1;}
				if(didUserGivedV[t1]){
					if( u_oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
						dropCasesAndEigenizeSquareMatrix(dV[i], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[i], false);
					}
					else{
						ptrToMatrix1 = omxMatrixDataColumnMajor(dV[i]);
					}
				}
				else{
					filteredCopy1.setZero(Eigyrows, Eigyrows);
					crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1 = filteredCopy1.data();
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1( ptrToMatrix1, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter i.
				Eigen::MatrixXd ytPdV_dtheta1 = u_Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
				for(j=i; j < numExplicitFreePar; j++){
					if(j==i){
						tr = doREML ? trace_prod_symm(u_P,dV_dtheta1) : trace_prod_symm(u_Vinv,dV_dtheta1);
						gradient(t1) = u_Scale*0.5*(tr - (ytPdV_dtheta1 * u_Py)(0,0)) +
							u_Scale*pullAugVal(1,a1,0);
						if(u_want & FF_COMPUTE_GRADIENT){
							u_fc->gradZ(t1) += gradient(t1);
						}
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							infoMat(t1,t1) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
								u_Scale*pullAugVal(2,a1,a1);
						}
					}
					else{
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							Eigen::MatrixXd filteredCopy2;
							double *ptrToMatrix2=0;
							t2 = gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							a2 = dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
							if(didUserGivedV[t2]){
								if( u_oge->numcases2drop && (dV[j]->rows > Eigyrows) ){
									dropCasesAndEigenizeSquareMatrix(dV[j], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[j], false);
								}
								else{
									ptrToMatrix2 = omxMatrixDataColumnMajor(dV[j]);
								}
							}
							else{
								filteredCopy2.setZero(Eigyrows, Eigyrows);
								crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2 = filteredCopy2.data();
							}
							Eigen::Map< Eigen::MatrixXd > dV_dtheta2( ptrToMatrix2, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter j.
							infoMat(t1,t2) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() *
								dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Py)(0,0) + u_Scale*pullAugVal(2,a1,a2);
							infoMat(t2,t1) = infoMat(t1,t2);
						}
					}
				}
			}
			else{
				gradient(t1) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(t1) = NA_REAL;
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


//Bin by row; Hessian is average information matrix:
template <typename T1, typename T2, typename T3>
void omxGREMLFitState::gradientAndAIM2(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndAIM2()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, t1, t2;
		double tr=0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = 0;
		int iend = rowbins[threadID].size();
		for(i=istart; i < iend; i++){
			tr=0;
			hrn = rowbins[threadID](i); //Current row number of the AIM.
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
				Eigen::MatrixXd filteredCopy1;
				double *ptrToMatrix1=0;
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				if(didUserGivedV[t1]){
					if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows) ){
						dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
					}
					else{
						ptrToMatrix1 = omxMatrixDataColumnMajor(dV[hrn]);
					}
				}
				else{
					filteredCopy1.setZero(Eigyrows, Eigyrows);
					crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1 = filteredCopy1.data();
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1( ptrToMatrix1, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter hrn.
				Eigen::MatrixXd ytPdV_dtheta1 = u_Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
				for(hcn=hrn; hcn < numExplicitFreePar; hcn++){
					if(hcn==hrn){
						tr = doREML ? trace_prod_symm(u_P,dV_dtheta1) : trace_prod_symm(u_Vinv,dV_dtheta1);
						gradient(hrn) = u_Scale*0.5*(tr - (ytPdV_dtheta1 * u_Py)(0,0)) +
							u_Scale*pullAugVal(1,a1,0);
						if(u_want & FF_COMPUTE_GRADIENT){
							u_fc->gradZ(hrn) += gradient(hrn);
						}
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							infoMat(hrn,hrn) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
								u_Scale*pullAugVal(2,a1,a1);
						}
					}
					//I think it can be assumed at this point that the Hessian is wanted?:
					else{
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							Eigen::MatrixXd filteredCopy2;
							double *ptrToMatrix2=0;
							t2 = gradMap[hcn];
							if(t2 < 0){continue;}
							a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
							if(didUserGivedV[t2]){
								if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
									dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
								}
								else{
									ptrToMatrix2 = omxMatrixDataColumnMajor(dV[hcn]);
								}
							}
							else{
								filteredCopy2.setZero(Eigyrows, Eigyrows);
								crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2 = filteredCopy2.data();
							}
							Eigen::Map< Eigen::MatrixXd > dV_dtheta2( ptrToMatrix2, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter hcn.
							infoMat(hrn,hcn) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() *
								dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Py)(0,0) + u_Scale*pullAugVal(2,a1,a2);
							infoMat(hcn,hrn) = infoMat(hrn,hcn);
						}
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


//Bin by cells; Hessian is average information matrix
template <typename T1, typename T2, typename T3>
void omxGREMLFitState::gradientAndAIM3(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndAIM3()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, inielem=0, t1, t2;
		double tr=0;
		double *ptrToMatrix1=0;
		Eigen::MatrixXd filteredCopy1;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int iend = AIMelembins[threadID].size();
		if(iend){inielem = AIMelembins[threadID](0);}
		while(inielem > 0){
			hcn++;
			inielem--;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
			}
		}
		while(i < iend){
			tr=0;
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				if(hrn==hcn || i==0){
					if(didUserGivedV[t1]){
						if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
						}
						else{
							ptrToMatrix1 = omxMatrixDataColumnMajor(dV[hrn]);
						}
					}
					else{
						filteredCopy1.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1 = filteredCopy1.data();
					}
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1( ptrToMatrix1, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter hrn.
				Eigen::MatrixXd ytPdV_dtheta1 = u_Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
				if(hrn==hcn){
					tr = doREML ? trace_prod_symm(u_P,dV_dtheta1) : trace_prod_symm(u_Vinv,dV_dtheta1);
					gradient(hrn) = u_Scale*0.5*(tr - (ytPdV_dtheta1 * u_Py)(0,0)) +
						u_Scale*pullAugVal(1,a1,0);
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(hrn) += gradient(hrn);
					}
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						infoMat(hrn,hrn) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
							u_Scale*pullAugVal(2,a1,a1);
					}
				}
				else{
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						Eigen::MatrixXd filteredCopy2;
						double *ptrToMatrix2=0;
						t2 = gradMap[hcn];
						if(t2 < 0){continue;}
						a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
						if(didUserGivedV[t2]){
							if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
								dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
							}
							else{
								ptrToMatrix2 = omxMatrixDataColumnMajor(dV[hcn]);
							}
						}
						else{
							filteredCopy2.setZero(Eigyrows, Eigyrows);
							crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
							ptrToMatrix2 = filteredCopy2.data();
						}
						Eigen::Map< Eigen::MatrixXd > dV_dtheta2( ptrToMatrix2, Eigyrows, Eigyrows ); //<--Derivative of V w/r/t parameter hcn.
						infoMat(hrn,hcn) = u_Scale*0.5*(ytPdV_dtheta1 * u_P.template selfadjointView<Eigen::Lower>() *
							dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Py)(0,0) + u_Scale*pullAugVal(2,a1,a2);
						infoMat(hcn,hrn) = infoMat(hrn,hcn);
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
			hcn++;
			i++;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


//Bin "naively"; Hessian is expected information matrix:
template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM1(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM1()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, t1=0, a1=0;
		double tr1=0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = threadID * numExplicitFreePar / u_nThreadz;
		int iend = (threadID+1) * numExplicitFreePar / u_nThreadz;
		if(threadID == u_nThreadz-1){iend = numExplicitFreePar;}
		for(i=istart; i < iend; i++){
			tr1=0;
			t1 = gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			//If the Hessian isn't wanted, then we compute the gradient elements the fast way.
			//If the Hessian is wanted, then we compute the gradient in a way that is slower, but
			//which saves time when computing the expected information matrix:
			if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
				int j=0, t2=0, a2=0;
				double tr2=0;
				u_hb->vars[i] = t1;
				if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
					double *ptrToMatrix1=0;
					Eigen::MatrixXd filteredCopy1;
					a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
					if(didUserGivedV[t1]){
						if( u_oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[i], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[i], false);
						}
						else{
							ptrToMatrix1 = omxMatrixDataColumnMajor(dV[i]);
						}
					}
					else{
						filteredCopy1.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1 = filteredCopy1.data();
					}
					Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1, Eigyrows, Eigyrows);
					//Eigen does not define a multiplication operator for selfadjointView times selfadjointView:
					Eigen::MatrixXd dV_dtheta1P = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_P;
					Eigen::MatrixXd dV_dtheta1Vinv;
					//dV_dtheta1Vinv.setZero(1,1);
					if(doREML){
						tr1 = dV_dtheta1P.trace();
					}
					else{
						dV_dtheta1Vinv = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_Vinv;
						tr1 = dV_dtheta1Vinv.trace();
					}
					for(j=i; j < numExplicitFreePar; j++){
						tr2=0;
						if(j==i){
							gradient(t1) = u_Scale*0.5*(tr1 - (u_Py.transpose() * dV_dtheta1P * u_Eigy)(0,0)) +
								u_Scale*pullAugVal(1,a1,0);
							if(u_want & FF_COMPUTE_GRADIENT){
								u_fc->gradZ(t1) += gradient(t1);
							}
							tr2 = doREML ? trace_prod(dV_dtheta1P,dV_dtheta1P) : trace_prod(dV_dtheta1Vinv,dV_dtheta1Vinv);
							infoMat(t1,t1) = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a1);
						}
						else{
							Eigen::MatrixXd filteredCopy2;
							double *ptrToMatrix2=0;
							t2 = gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							a2 = dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
							if(didUserGivedV[t2]){
								if( u_oge->numcases2drop && (dV[j]->rows > Eigyrows) ){
									dropCasesAndEigenizeSquareMatrix(dV[j], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[j], false);
								}
								else{
									ptrToMatrix2 = omxMatrixDataColumnMajor(dV[j]);
								}
							}
							else{
								filteredCopy2.setZero(Eigyrows, Eigyrows);
								crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2 = filteredCopy2.data();
							}
							Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2, Eigyrows, Eigyrows);
							if(doREML){
								Eigen::MatrixXd dV_dtheta2P = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_P;
								tr2 = trace_prod(dV_dtheta1P,dV_dtheta2P);
							}
							else{
								Eigen::MatrixXd dV_dtheta2Vinv = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Vinv;
								tr2 = trace_prod(dV_dtheta1Vinv,dV_dtheta2Vinv);
							}
							infoMat(t1,t2) = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a2);
							infoMat(t2,t1) = infoMat(t1,t2);
						}
					}
				}
				else{
					gradient(t1) = NA_REAL;
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(t1) = NA_REAL;
					}
				}
			}

			else{
				if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
					double *ptrToMatrix1=0;
					Eigen::MatrixXd filteredCopy1;
					a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
					if(didUserGivedV[t1]){
						if( u_oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[i], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[i], false);
						}
						else{
							ptrToMatrix1 = omxMatrixDataColumnMajor(dV[i]);
						}
					}
					else{
						filteredCopy1.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1 = filteredCopy1.data();
					}
					Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1, Eigyrows, Eigyrows);
					Eigen::MatrixXd ytPdV_dtheta1 = u_Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
					tr1 = doREML ? trace_prod_symm(u_P,dV_dtheta1) : trace_prod_symm(u_Vinv,dV_dtheta1);
					gradient(t1) = u_Scale*0.5*(tr1 - (ytPdV_dtheta1 * u_Py)(0,0)) +
						u_Scale*pullAugVal(1,a1,0);
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(t1) += gradient(t1);
					}
				}
				else{
					gradient(t1) = NA_REAL;
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(t1) = NA_REAL;
					}
				}
			}
		}

	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


//Bin by row; Hessian is expected information matrix:
template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM2(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM2()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, t1, t2;
		double tr1=0, tr2=0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = 0;
		int iend = rowbins[threadID].size();
		for(i=istart; i < iend; i++){
			tr1=0;
			hrn = rowbins[threadID](i); //Current row number of the EIM.
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
				double *ptrToMatrix1=0;
				Eigen::MatrixXd filteredCopy1;
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				if(didUserGivedV[t1]){
					if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows)){
						dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
					}
					else{
						ptrToMatrix1 = omxMatrixDataColumnMajor(dV[hrn]);
					}
				}
				else{
					filteredCopy1.setZero(Eigyrows, Eigyrows);
					crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1 = filteredCopy1.data();
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1, Eigyrows, Eigyrows); //<--Derivative of V w/r/t parameter hrn.
				//Eigen does not define a multiplication operator for selfadjointView times selfadjointView:
				Eigen::MatrixXd dV_dtheta1P = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_P;
				Eigen::MatrixXd dV_dtheta1Vinv;
				if(doREML){
					tr1 = dV_dtheta1P.trace();
				}
				else{
					dV_dtheta1Vinv = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_Vinv;
					tr1 = dV_dtheta1Vinv.trace();
				}
				for(hcn=hrn; hcn < numExplicitFreePar; hcn++){
					tr2=0;
					if(hcn==hrn){
						gradient(hrn) = u_Scale*0.5*(tr1 - (u_Py.transpose() * dV_dtheta1P * u_Eigy)(0,0)) +
							u_Scale*pullAugVal(1,a1,0);
						if(u_want & FF_COMPUTE_GRADIENT){
							u_fc->gradZ(hrn) += gradient(hrn);
						}
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							tr2 = doREML ? trace_prod(dV_dtheta1P,dV_dtheta1P) : trace_prod(dV_dtheta1Vinv,dV_dtheta1Vinv);
							infoMat(hrn,hrn)  = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a1);
						}
					}
					//I think it can be assumed at this point that the Hessian is wanted?:
					else{
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							Eigen::MatrixXd filteredCopy2;
							double *ptrToMatrix2=0;
							t2 = gradMap[hcn];
							if(t2 < 0){continue;}
							a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
							if(didUserGivedV[t2]){
								if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
									dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
								}
								else{
									ptrToMatrix2 = omxMatrixDataColumnMajor(dV[hcn]);
								}
							}
							else{
								filteredCopy2.setZero(Eigyrows, Eigyrows);
								crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2 = filteredCopy2.data();
							}
							Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2, Eigyrows, Eigyrows); //<--Derivative of V w/r/t parameter hcn.
							if(doREML){
								Eigen::MatrixXd dV_dtheta2P = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_P;
								tr2 = trace_prod(dV_dtheta1P,dV_dtheta2P);
							}
							else{
								Eigen::MatrixXd dV_dtheta2Vinv = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Vinv;
								tr2 = trace_prod(dV_dtheta1Vinv,dV_dtheta2Vinv);
							}
							infoMat(hrn,hcn) = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a2);
							infoMat(hcn,hrn) = infoMat(hrn,hcn);
						}
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


//Bin by cell; Hessian is expected information matrix:
template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM3(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM3()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, inielem=0, t1, t2;
		double tr1=0, tr2=0;
		double *ptrToMatrix1=0;
		Eigen::MatrixXd filteredCopy1;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int iend = AIMelembins[threadID].size();
		if(iend){inielem = AIMelembins[threadID](0);}
		while(inielem > 0){
			hcn++;
			inielem--;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
			}
		}
		while(i < iend){
			tr2=0;
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( (didUserGivedV[t1] || derivType==1) && !didUserProvideYhat){
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				if(hrn==hcn || i==0){
					if(didUserGivedV[t1]){
						if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1, ptrToMatrix1, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
						}
						else{
							ptrToMatrix1 = omxMatrixDataColumnMajor(dV[hrn]);
						}
					}
					else{
						filteredCopy1.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1 = filteredCopy1.data();
					}
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1, Eigyrows, Eigyrows); //<--Derivative of V w/r/t parameter hrn.
				Eigen::MatrixXd dV_dtheta1P = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_P;
				Eigen::MatrixXd dV_dtheta1Vinv;
				if(doREML){
					tr1 = dV_dtheta1P.trace();
				}
				else{
					dV_dtheta1Vinv = dV_dtheta1.selfadjointView<Eigen::Lower>() * u_Vinv;
					tr1 = dV_dtheta1Vinv.trace();
				}
				if(hrn==hcn){
					gradient(hrn) = u_Scale*0.5*(tr1 - (u_Py.transpose() * dV_dtheta1P * u_Eigy)(0,0)) +
						u_Scale*pullAugVal(1,a1,0);
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(hrn) += gradient(hrn);
					}
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						tr2 = doREML ? trace_prod(dV_dtheta1P,dV_dtheta1P) : trace_prod(dV_dtheta1Vinv,dV_dtheta1Vinv);
						infoMat(hrn,hrn)  = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a1);
					}
				}
				else{
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						Eigen::MatrixXd filteredCopy2;
						double *ptrToMatrix2=0;
						t2 = gradMap[hcn];
						if(t2 < 0){continue;}
						a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
						if(didUserGivedV[t2]){
							if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
								dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2, ptrToMatrix2, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
							}
							else{
								ptrToMatrix2 = omxMatrixDataColumnMajor(dV[hcn]);
							}
						}
						else{
							filteredCopy2.setZero(Eigyrows, Eigyrows);
							crude_numeric_dV(u_fc, curEst, filteredCopy2, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
							ptrToMatrix2 = filteredCopy2.data();
						}
						Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2, Eigyrows, Eigyrows); //<--Derivative of V w/r/t parameter hcn.
						if(doREML){
							Eigen::MatrixXd dV_dtheta2P = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_P;
							tr2 = trace_prod(dV_dtheta1P,dV_dtheta2P);
						}
						else{
							Eigen::MatrixXd dV_dtheta2Vinv = dV_dtheta2.selfadjointView<Eigen::Lower>() * u_Vinv;
							tr2 = trace_prod(dV_dtheta1Vinv,dV_dtheta2Vinv);
						}
						infoMat(hrn,hcn) = u_Scale*0.5*tr2 + u_Scale*pullAugVal(2,a1,a2);
						infoMat(hcn,hrn) = infoMat(hrn,hcn);
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
			hcn++;
			i++;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
				tr1=0;
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void omxGREMLFitState::gradientAndOIM1(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, Eigen::MatrixBase<T1> &u_P,
		double u_Scale, Eigen::MatrixBase<T2> &u_Py, Eigen::MatrixBase<T3> &u_Eigy, Eigen::MatrixBase<T4> &u_Vinv, Eigen::MatrixBase<T5> &u_VinvResid,
		Eigen::MatrixBase<T6> &u_VinvResidResidT){
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		if(didUserProvideYhat){
			int i=0, t1=0, a1=0, j=0;
			double term1=0.0, term2=0.0;
			Eigen::VectorXd curEst(numExplicitFreePar);
			u_fc->copyEstToOptimizer(curEst);
			Eigen::VectorXd curEst1p(1);
			int threadID = omx_absolute_thread_num();
			int istart = threadID * numExplicitFreePar / u_nThreadz;
			int iend = (threadID+1) * numExplicitFreePar / u_nThreadz;
			if(threadID == u_nThreadz-1){iend = numExplicitFreePar;}
			for(i=istart; i < iend; i++){
				term1=0.0, term2=0.0;
				t1 = gradMap[i]; //<--Parameter number for parameter i.
				if(t1 < 0){continue;}
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[i] = t1;}
				if( derivType==1 || (didUserGivedV[t1] && didUserGivedyhat[t1]) ){
					double *ptrToMatrix1m=0, *ptrToMatrix1v=0;
					Eigen::MatrixXd filteredCopy1m, filteredCopy1v;
					a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
					
					if(didUserGivedV[t1] && dV[i]){
						if( u_oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[i], filteredCopy1v, ptrToMatrix1v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[i], false);
						}
						else{
							ptrToMatrix1v = omxMatrixDataColumnMajor(dV[i]);
						}
					}
					else{
						filteredCopy1v.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1v, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1v = filteredCopy1v.data();
					}
					Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1v, Eigyrows, Eigyrows);
					
					if(didUserGivedyhat[t1] && dyhat[i]){
						if( u_oge->numcases2drop && (dyhat[i]->size() > Eigyrows) ){
							dropCasesAndEigenizeColumnVector(dyhat[i], filteredCopy1m, ptrToMatrix1m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[i], false);
						}
						else{
							ptrToMatrix1m = omxMatrixDataColumnMajor(dyhat[i]);
						}
					}
					else{
						filteredCopy1m.setZero(Eigyrows, 1);
						crude_numeric_dyhat(u_fc, curEst, filteredCopy1m, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1m = filteredCopy1m.data();
					}
					Eigen::Map< Eigen::MatrixXd > dyhat_dtheta1(ptrToMatrix1m, Eigyrows, 1);
					
					Eigen::MatrixXd VinvdV_dtheta1 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta1;
					
					if(u_want & FF_COMPUTE_GRADIENT){
						term1 = VinvdV_dtheta1.trace() - trace_prod(VinvdV_dtheta1,u_VinvResidResidT);
						term2 = -2*(dyhat_dtheta1.transpose()*u_VinvResid)(0,0);
						gradient(t1) = u_Scale*0.5*(term1+term2) + u_Scale*pullAugVal(1,a1,0);
						u_fc->gradZ(t1) += gradient(t1);
						//mxLog("gradient element %d is %f", t1, gradient(t1));
					}
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						for(j=i; j < numExplicitFreePar; j++){
							double tt0=0.0, tt0_0=0.0, tt0_1=0.0, tt0_2=0.0, tt1=0.0, tt2=0.0, tt3=0.0, tt4=0.0, tt5=0.0;
							int t2=0, a2;
							t2 = gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							if(didUserGivedV[t1] || didUserGivedyhat[t1] || derivType==1){
								double *ptrToMatrix2m=0, *ptrToMatrix2v=0;
								Eigen::MatrixXd filteredCopy2m, filteredCopy2v;
								a2 = dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
								
								if(didUserGivedV[t2] && dV[j]){
									if( u_oge->numcases2drop && (dV[j]->rows > Eigyrows) ){
										dropCasesAndEigenizeSquareMatrix(dV[j], filteredCopy2v, ptrToMatrix2v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[j], false);
									}
									else{
										ptrToMatrix2v = omxMatrixDataColumnMajor(dV[j]);
									}
								}
								else{
									filteredCopy2v.setZero(Eigyrows, Eigyrows);
									crude_numeric_dV(u_fc, curEst, filteredCopy2v, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2v = filteredCopy2v.data();
								}
								Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2v, Eigyrows, Eigyrows);
								
								if(didUserGivedyhat[t2] && dyhat[j]){
									if( u_oge->numcases2drop && (dyhat[j]->size() > Eigyrows) ){
										dropCasesAndEigenizeColumnVector(dyhat[j], filteredCopy2m, ptrToMatrix2m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[j], false);
									}
									else{
										ptrToMatrix2m = omxMatrixDataColumnMajor(dyhat[j]);
									}
								}
								else{
									filteredCopy2m.setZero(Eigyrows, 1);
									crude_numeric_dyhat(u_fc, curEst, filteredCopy2m, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2m = filteredCopy2m.data();
								}
								Eigen::Map< Eigen::MatrixXd > dyhat_dtheta2(ptrToMatrix2m, Eigyrows, 1);
								
								Eigen::MatrixXd VinvdV_dtheta2 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta2;
								Eigen::MatrixXd VinvdV_dtheta2VinvdV_dtheta1 = VinvdV_dtheta2 * VinvdV_dtheta1;
								
								tt0_0 = 0.0 - VinvdV_dtheta2VinvdV_dtheta1.trace(); //Vinv2ndDer.trace() - trace_prod(VinvDer2,VinvDer)
								tt0_1 = -1.0*0.0; //trace_prod(Vinv2ndDer,u_VinvResidResidT)
								tt0_2 = trace_prod(VinvdV_dtheta2VinvdV_dtheta1,u_VinvResidResidT);
								tt0 = -0.5*(tt0_0+tt0_1+tt0_2);
								tt1 = -0.5*tt0_2;
								tt2 = -0.5*trace_prod((VinvdV_dtheta1*u_Vinv),2*(dyhat_dtheta2*u_oge->residual.transpose()));
								tt3 = 0.0; //dyhat_dtheta1dtheta2.transpose() * u_VinvResid;
								tt4 = (dyhat_dtheta1.transpose() * dV_dtheta2 * u_oge->residual)(0,0);
								tt5 = -1.0*(dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta2)(0,0);
								infoMat(t1,t2) = -1.0*u_Scale*(tt0 + tt1 + tt2 + tt3 + tt4 + tt5) + u_Scale*pullAugVal(2,a1,a2);
								infoMat(t2,t1) = infoMat(t1,t2);
							}
						}
					}
				}
				else{
					gradient(t1) = NA_REAL;
					if(u_want & FF_COMPUTE_GRADIENT){
						u_fc->gradZ(t1) = NA_REAL;
					}
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM1_yhat(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
		double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
		Eigen::MatrixBase<T4> &u_VinvResidResidT){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM1_yhat()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, t1=0, a1=0, j=0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = threadID * numExplicitFreePar / u_nThreadz;
		int iend = (threadID+1) * numExplicitFreePar / u_nThreadz;
		if(threadID == u_nThreadz-1){iend = numExplicitFreePar;}
		for(i=istart; i < iend; i++){
			t1 = gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			if( derivType==1 || (didUserGivedV[t1] && didUserGivedyhat[t1]) ){
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[i] = t1;}
				double *ptrToMatrix1m=0, *ptrToMatrix1v=0;
				Eigen::MatrixXd filteredCopy1m, filteredCopy1v;
				a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
				
				if(didUserGivedV[t1] && dV[i]){
					if( u_oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
						dropCasesAndEigenizeSquareMatrix(dV[i], filteredCopy1v, ptrToMatrix1v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[i], false);
					}
					else{
						ptrToMatrix1v = omxMatrixDataColumnMajor(dV[i]);
					}
				}
				else{
					filteredCopy1v.setZero(Eigyrows, Eigyrows);
					crude_numeric_dV(u_fc, curEst, filteredCopy1v, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1v = filteredCopy1v.data();
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1v, Eigyrows, Eigyrows);
				
				if(didUserGivedyhat[t1] && dyhat[i]){
					if( u_oge->numcases2drop && (dyhat[i]->size() > Eigyrows) ){
						dropCasesAndEigenizeColumnVector(dyhat[i], filteredCopy1m, ptrToMatrix1m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[i], false);
					}
					else{
						ptrToMatrix1m = omxMatrixDataColumnMajor(dyhat[i]);
					}
				}
				else{
					filteredCopy1m.setZero(Eigyrows, 1);
					crude_numeric_dyhat(u_fc, curEst, filteredCopy1m, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1m = filteredCopy1m.data();
				}
				Eigen::Map< Eigen::MatrixXd > dyhat_dtheta1(ptrToMatrix1m, Eigyrows, 1);
				
				Eigen::MatrixXd VinvdV_dtheta1 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta1;
				
				if(u_want & FF_COMPUTE_GRADIENT){
					double term1 = VinvdV_dtheta1.trace() - trace_prod(VinvdV_dtheta1,u_VinvResidResidT);
					double term2 = -2*(dyhat_dtheta1.transpose()*u_VinvResid)(0,0);
					gradient(t1) = u_Scale*0.5*(term1+term2) + u_Scale*pullAugVal(1,a1,0);
					u_fc->gradZ(t1) += gradient(t1);
					//mxLog("gradient element %d is %f", t1, gradient(t1));
				}
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
					for(j=i; j < numExplicitFreePar; j++){
						double tt1=0.0, tt2=0.0;
						int t2=0, a2=0;
						if(j==i){
							tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta1);
							tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta1)(0,0);
							infoMat(t1,t1) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a2);
						}
						else{
							t2 = gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							if( (didUserGivedV[t2] && didUserGivedyhat[t2]) || derivType==1 ){
								double *ptrToMatrix2m=0, *ptrToMatrix2v=0;
								Eigen::MatrixXd filteredCopy2m, filteredCopy2v;
								a2 = dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
								
								if(didUserGivedV[t2] && dV[j]){
									if( u_oge->numcases2drop && (dV[j]->rows > Eigyrows) ){
										dropCasesAndEigenizeSquareMatrix(dV[j], filteredCopy2v, ptrToMatrix2v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[j], false);
									}
									else{
										ptrToMatrix2v = omxMatrixDataColumnMajor(dV[j]);
									}
								}
								else{
									filteredCopy2v.setZero(Eigyrows, Eigyrows);
									crude_numeric_dV(u_fc, curEst, filteredCopy2v, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2v = filteredCopy2v.data();
								}
								Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2v, Eigyrows, Eigyrows);
								
								if(didUserGivedyhat[t2] && dyhat[j]){
									if( u_oge->numcases2drop && (dyhat[j]->size() > Eigyrows) ){
										dropCasesAndEigenizeColumnVector(dyhat[j], filteredCopy2m, ptrToMatrix2m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[j], false);
									}
									else{
										ptrToMatrix2m = omxMatrixDataColumnMajor(dyhat[j]);
									}
								}
								else{
									filteredCopy2m.setZero(Eigyrows, 1);
									crude_numeric_dyhat(u_fc, curEst, filteredCopy2m, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2m = filteredCopy2m.data();
								}
								Eigen::Map< Eigen::MatrixXd > dyhat_dtheta2(ptrToMatrix2m, Eigyrows, 1);
								
								Eigen::MatrixXd VinvdV_dtheta2 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta2;
								
								tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta2);
								tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta2)(0,0);
								infoMat(t1,t2) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a2);
								infoMat(t2,t1) = infoMat(t1,t2);
							}
						}
					}
				}
			}
			else{
				gradient(t1) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(t1) = NA_REAL;
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM2_yhat(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
		double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
		Eigen::MatrixBase<T4> &u_VinvResidResidT){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM2_yhat()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, t1;
		double term1=0.0, term2=0.0;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = 0;
		int iend = rowbins[threadID].size();
		for(i=istart; i < iend; i++){
			term1=0.0, term2=0.0;
			hrn = rowbins[threadID](i); //Current row number of the EIM.
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( derivType==1 || (didUserGivedV[t1] && didUserGivedyhat[t1]) ){
				double *ptrToMatrix1m=0, *ptrToMatrix1v=0;
				Eigen::MatrixXd filteredCopy1m, filteredCopy1v;
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter i.
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				
				if(didUserGivedV[t1] && dV[hrn]){
					if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows) ){
						dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1v, ptrToMatrix1v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
					}
					else{
						ptrToMatrix1v = omxMatrixDataColumnMajor(dV[hrn]);
					}
				}
				else{
					filteredCopy1v.setZero(Eigyrows, Eigyrows);
					crude_numeric_dV(u_fc, curEst, filteredCopy1v, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1v = filteredCopy1v.data();
				}
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1v, Eigyrows, Eigyrows);
				
				if(didUserGivedyhat[t1] && dyhat[hrn]){
					if( u_oge->numcases2drop && (dyhat[hrn]->size() > Eigyrows) ){
						dropCasesAndEigenizeColumnVector(dyhat[hrn], filteredCopy1m, ptrToMatrix1m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[hrn], false);
					}
					else{
						ptrToMatrix1m = omxMatrixDataColumnMajor(dyhat[hrn]);
					}
				}
				else{
					filteredCopy1m.setZero(Eigyrows, 1);
					crude_numeric_dyhat(u_fc, curEst, filteredCopy1m, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
					ptrToMatrix1m = filteredCopy1m.data();
				}
				Eigen::Map< Eigen::MatrixXd > dyhat_dtheta1(ptrToMatrix1m, Eigyrows, 1);
				
				Eigen::MatrixXd VinvdV_dtheta1 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta1;
				
				for(hcn=hrn; hcn < numExplicitFreePar; hcn++){
					double tt1=0.0, tt2=0.0;
					if(hcn==hrn){
						if(u_want & FF_COMPUTE_GRADIENT){
							term1 = VinvdV_dtheta1.trace() - trace_prod(VinvdV_dtheta1,u_VinvResidResidT);
							term2 = -2*(dyhat_dtheta1.transpose()*u_VinvResid)(0,0);
							gradient(hrn) = u_Scale*0.5*(term1+term2) + u_Scale*pullAugVal(1,a1,0);
							u_fc->gradZ(hrn) += gradient(hrn);
							//mxLog("gradient element %d is %f", t1, gradient(t1));
						}
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta1);
							tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta1)(0,0);
							infoMat(hrn,hrn) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a1);
						}
					}
					else{
						if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							int t2=0, a2;
							t2 = gradMap[hcn]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							if( (didUserGivedV[t2] && didUserGivedyhat[t2]) || derivType==1 ){
								double *ptrToMatrix2m=0, *ptrToMatrix2v=0;
								Eigen::MatrixXd filteredCopy2m, filteredCopy2v;
								a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter j.
								
								if(didUserGivedV[t2] && dV[hcn]){
									if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
										dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2v, ptrToMatrix2v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
									}
									else{
										ptrToMatrix2v = omxMatrixDataColumnMajor(dV[hcn]);
									}
								}
								else{
									filteredCopy2v.setZero(Eigyrows, Eigyrows);
									crude_numeric_dV(u_fc, curEst, filteredCopy2v, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2v = filteredCopy2v.data();
								}
								Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2v, Eigyrows, Eigyrows);
								
								if(didUserGivedyhat[t2] && dyhat[hcn]){
									if( u_oge->numcases2drop && (dyhat[hcn]->size() > Eigyrows) ){
										dropCasesAndEigenizeColumnVector(dyhat[hcn], filteredCopy2m, ptrToMatrix2m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[hcn], false);
									}
									else{
										ptrToMatrix2m = omxMatrixDataColumnMajor(dyhat[hcn]);
									}
								}
								else{
									filteredCopy2m.setZero(Eigyrows, 1);
									crude_numeric_dyhat(u_fc, curEst, filteredCopy2m, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
									ptrToMatrix2m = filteredCopy2m.data();
								}
								Eigen::Map< Eigen::MatrixXd > dyhat_dtheta2(ptrToMatrix2m, Eigyrows, 1);
								
								Eigen::MatrixXd VinvdV_dtheta2 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta2;
								
								tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta2);
								tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta2)(0,0);
								infoMat(hrn,hcn) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a2);
								infoMat(hcn,hrn) = infoMat(hrn,hcn);
							}
						}
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


template <typename T1, typename T2, typename T3, typename T4>
void omxGREMLFitState::gradientAndEIM3_yhat(
		int u_nThreadz, int Eigyrows, FitContext *u_fc, int u_want, HessianBlock *u_hb, omxGREMLExpectation *u_oge, 
		double u_Scale, Eigen::MatrixBase<T1> &u_Eigy, Eigen::MatrixBase<T2> &u_Vinv, Eigen::MatrixBase<T3> &u_VinvResid, 
		Eigen::MatrixBase<T4> &u_VinvResidResidT){
	//if(OMX_DEBUG_ALGEBRA){
	mxLog("Now beginning `gradientAndEIM3_yhat()`");
	//}
#pragma omp parallel num_threads(u_nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, inielem=0, t1, t2;
		double *ptrToMatrix1m=0, *ptrToMatrix1v=0;
		Eigen::MatrixXd filteredCopy1m, filteredCopy1v;
		Eigen::VectorXd curEst(numExplicitFreePar);
		u_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int iend = AIMelembins[threadID].size();
		if(iend){inielem = AIMelembins[threadID](0);}
		while(inielem > 0){
			hcn++;
			inielem--;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
			}
		}
		while(i < iend){
			t1 = gradMap[hrn];
			if(t1 < 0){continue;} //Check for negative parameter number.
			if( derivType==1 || (didUserGivedV[t1] && didUserGivedyhat[t1]) ){
				if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){u_hb->vars[hrn] = t1;}
				a1 = dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter i.
				if(hrn==hcn || i==0){
					if(didUserGivedV[t1] && dV[hrn]){
						if( u_oge->numcases2drop && (dV[hrn]->rows > Eigyrows) ){
							dropCasesAndEigenizeSquareMatrix(dV[hrn], filteredCopy1v, ptrToMatrix1v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hrn], false);
						}
						else{
							ptrToMatrix1v = omxMatrixDataColumnMajor(dV[hrn]);
						}
					}
					else{
						filteredCopy1v.setZero(Eigyrows, Eigyrows);
						crude_numeric_dV(u_fc, curEst, filteredCopy1v, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1v = filteredCopy1v.data();
					}
					
					if(didUserGivedyhat[t1] && dyhat[hrn]){
						if( u_oge->numcases2drop && (dyhat[hrn]->size() > Eigyrows) ){
							dropCasesAndEigenizeColumnVector(dyhat[hrn], filteredCopy1m, ptrToMatrix1m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[hrn], false);
						}
						else{
							ptrToMatrix1m = omxMatrixDataColumnMajor(dyhat[hrn]);
						}
					}
					else{
						filteredCopy1m.setZero(Eigyrows, 1);
						crude_numeric_dyhat(u_fc, curEst, filteredCopy1m, t1, u_oge, (u_nThreadz>1 ? threadID : -1));
						ptrToMatrix1m = filteredCopy1m.data();
					}
				}
				
				Eigen::Map< Eigen::MatrixXd > dV_dtheta1(ptrToMatrix1v, Eigyrows, Eigyrows);
				Eigen::Map< Eigen::MatrixXd > dyhat_dtheta1(ptrToMatrix1m, Eigyrows, 1);
				Eigen::MatrixXd VinvdV_dtheta1 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta1;
				
				if(hrn==hcn){
					if(u_want & FF_COMPUTE_GRADIENT){
						double term1 = VinvdV_dtheta1.trace() - trace_prod(VinvdV_dtheta1,u_VinvResidResidT);
						double term2 = -2*(dyhat_dtheta1.transpose()*u_VinvResid)(0,0);
						gradient(hrn) = u_Scale*0.5*(term1+term2) + u_Scale*pullAugVal(1,a1,0);
						u_fc->gradZ(hrn) += gradient(hrn);
						//mxLog("gradient element %d is %f", t1, gradient(t1));
					}
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						double tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta1);
						double tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta1)(0,0);
						infoMat(hrn,hrn) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a1);
					}
				}
				else{
					if(u_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						t2 = gradMap[hcn]; //<--Parameter number for parameter j.
						if(t2 < 0){continue;}
						if( (didUserGivedV[t2] && didUserGivedyhat[t2]) || derivType==1 ){
							double *ptrToMatrix2m=0, *ptrToMatrix2v=0;
							Eigen::MatrixXd filteredCopy2m, filteredCopy2v;
							a2 = dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter j.
							
							if(didUserGivedV[t2] && dV[hcn]){
								if( u_oge->numcases2drop && (dV[hcn]->rows > Eigyrows) ){
									dropCasesAndEigenizeSquareMatrix(dV[hcn], filteredCopy2v, ptrToMatrix2v, u_oge->numcases2drop, u_oge->dropcase, true, origdVdim[hcn], false);
								}
								else{
									ptrToMatrix2v = omxMatrixDataColumnMajor(dV[hcn]);
								}
							}
							else{
								filteredCopy2v.setZero(Eigyrows, Eigyrows);
								crude_numeric_dV(u_fc, curEst, filteredCopy2v, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2v = filteredCopy2v.data();
							}
							Eigen::Map< Eigen::MatrixXd > dV_dtheta2(ptrToMatrix2v, Eigyrows, Eigyrows);
							
							if(didUserGivedyhat[t2] && dyhat[hcn]){
								if( u_oge->numcases2drop && (dyhat[hcn]->size() > Eigyrows) ){
									dropCasesAndEigenizeColumnVector(dyhat[hcn], filteredCopy2m, ptrToMatrix2m, u_oge->numcases2drop, u_oge->dropcase, false, origdyhatdim[hcn], false);
								}
								else{
									ptrToMatrix2m = omxMatrixDataColumnMajor(dyhat[hcn]);
								}
							}
							else{
								filteredCopy2m.setZero(Eigyrows, 1);
								crude_numeric_dyhat(u_fc, curEst, filteredCopy2m, t2, u_oge, (u_nThreadz>1 ? threadID : -1));
								ptrToMatrix2m = filteredCopy2m.data();
							}
							Eigen::Map< Eigen::MatrixXd > dyhat_dtheta2(ptrToMatrix2m, Eigyrows, 1);
							
							Eigen::MatrixXd VinvdV_dtheta2 = u_Vinv.template selfadjointView<Eigen::Lower>() * dV_dtheta2;
							
							double tt1 = 0.5*trace_prod(VinvdV_dtheta1,VinvdV_dtheta2);
							double tt2 = (dyhat_dtheta1.transpose() * u_Vinv * dyhat_dtheta2)(0,0);
							infoMat(hrn,hcn) = u_Scale*(tt1 + tt2) + u_Scale*pullAugVal(2,a1,a2);
							infoMat(hcn,hrn) = infoMat(hrn,hcn);
						}
					}
				}
			}
			else{
				gradient(hrn) = NA_REAL;
				if(u_want & FF_COMPUTE_GRADIENT){
					u_fc->gradZ(hrn) = NA_REAL;
				}
			}
			hcn++;
			i++;
			if(hcn == numExplicitFreePar){
				hrn++;
				hcn=hrn;
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
	return;
}


template <typename T1, typename T2>
void omxGREMLFitState::crude_numeric_dV(
		FitContext *u_fc, Eigen::MatrixBase<T1> &u_curEst, Eigen::MatrixBase<T2> &dV_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId)
{
	int c, r;
	FitContext *fc2 = u_fc;
	double *ptrToMatrix1=0, *ptrToMatrix2=0;
	if(thrId >= 0){
		fc2 = u_fc->childList[thrId];
	}
	omxState *st = fc2->state;

	/*
	 * Store current elements of V:
	 */
	Eigen::MatrixXd Vcurr(cov->rows, cov->cols);
	//By this time, `cov` has most certainly been filtered for missing observations in y or X:
	Eigen::Map< Eigen::MatrixXd > EigVmap(omxMatrixDataColumnMajor(cov), cov->rows, cov->cols);
	for(c=0; c < EigVmap.cols(); c++){
		for(r=c; r < EigVmap.rows(); r++){
			Vcurr(r,c) = EigVmap(r,c);
		}
	}

	//Shift parameter of interest and compute V at new point in parameter space:
	u_curEst[Parnum] += 1e-4;
	fc2->setEstFromOptimizer(u_curEst);
	omxMatrix *mat = st->lookupDuplicate(cov);
	omxRecompute(mat, fc2);

	//dV_dtheta_tmp is only needed in the single-threaded case???:
	if(thrId >= 0){
		if( ge->numcases2drop && mat->rows > y->cols ){
			dropCasesAndEigenizeSquareMatrix(mat, dV_dtheta, ptrToMatrix1, ge->numcases2drop, ge->dropcase, true, mat->rows, true);
		}
		else{
			dV_dtheta = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols);
		}
	}
	else{
		{
			Eigen::MatrixXd dV_dtheta_tmp; //<--Not needed outside this scope.
			if( ge->numcases2drop && mat->rows > y->cols ){
				dropCasesAndEigenizeSquareMatrix(mat, dV_dtheta_tmp, ptrToMatrix2, ge->numcases2drop, ge->dropcase, true, mat->rows, true);
			}
			else{
				dV_dtheta_tmp = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols);
			}
			dV_dtheta.resize(dV_dtheta_tmp.rows(), dV_dtheta_tmp.rows());
			for(c=0; c < dV_dtheta.rows(); c++){
				for(r=c; r < dV_dtheta.rows(); r++){
					dV_dtheta(r,c) = dV_dtheta_tmp(r,c);
				}
			}
		}
	}
	//Does this actually save memory...?:
	dV_dtheta.template triangularView<Eigen::Lower>() = (dV_dtheta - Vcurr)/1e-4;

	//Put things back the way they were:
	u_curEst[Parnum] -= 1e-4;
	fc2->setEstFromOptimizer(u_curEst);
	Eigen::Map< Eigen::MatrixXd > EigV2(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols);
	for(c=0; c < mat->rows; c++){
		for(r=c; r < mat->rows; r++){
			EigV2(r,c) = Vcurr(r,c);
		}
	}
	return;
}


template <typename T1, typename T2>
void omxGREMLFitState::crude_numeric_dyhat(
		FitContext *u_fc, Eigen::MatrixBase<T1> &u_curEst, Eigen::MatrixBase<T2> &dyhat_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId)
{
	int r;
	FitContext *fc2 = u_fc;
	double *ptrToMatrix1=0, *ptrToMatrix2=0;
	if(thrId >= 0){
		fc2 = u_fc->childList[thrId];
	}
	omxState *st = fc2->state;
	
	/*
	 * Store current elements of yhat:
	 */
	Eigen::MatrixXd yhatcurr(ge->yhatFromUser->size(), 1);
	//By this time, `yhat` has most certainly been filtered for missing observations:
	Eigen::Map< Eigen::MatrixXd > Eigyhatmap(omxMatrixDataColumnMajor(ge->yhatFromUser), ge->yhatFromUser->size(), 1);
	for(r=0; r < Eigyhatmap.rows(); r++){
		yhatcurr(r,0) = Eigyhatmap(r,0);
	}

	//Shift parameter of interest and compute yhat at new point in parameter space:
	u_curEst[Parnum] += 1e-4;
	fc2->setEstFromOptimizer(u_curEst);
	omxMatrix *mat = st->lookupDuplicate(ge->yhatFromUser);
	omxRecompute(mat, fc2);
	
	//dyhat_dtheta_tmp is only needed in the single-threaded case???:
	if(thrId >= 0){
		if( ge->numcases2drop && mat->size() > y->cols ){
			dropCasesAndEigenizeColumnVector(mat, dyhat_dtheta, ptrToMatrix1, ge->numcases2drop, ge->dropcase, false, mat->size(), true);
		}
		else{
			dyhat_dtheta = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->size(), 1);
		}
	}
	else{
		{
			Eigen::MatrixXd dyhat_dtheta_tmp; //<--Not needed outside this scope.
			if( ge->numcases2drop && mat->size() > y->cols ){
				dropCasesAndEigenizeColumnVector(mat, dyhat_dtheta_tmp, ptrToMatrix2, ge->numcases2drop, ge->dropcase, false, mat->size(), true);
			}
			else{
				dyhat_dtheta_tmp = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->size(), 1);
			}
			dyhat_dtheta.resize(dyhat_dtheta_tmp.rows(), 1);
			for(r=0; r < dyhat_dtheta.rows(); r++){
				dyhat_dtheta(r,0) = dyhat_dtheta_tmp(r,0);
			}
		}
	}
	//Does this actually save memory...?:
	dyhat_dtheta = (dyhat_dtheta - yhatcurr)/1e-4;
	
	//Put things back the way they were:
	u_curEst[Parnum] -= 1e-4;
	fc2->setEstFromOptimizer(u_curEst);
	Eigen::Map< Eigen::MatrixXd > Eigyhat2(omxMatrixDataColumnMajor(mat), mat->size(), 1);
	for(r=0; r < mat->size(); r++){
		Eigyhat2(r,0) = yhatcurr(r,0);
	}
	return;
}



void omxGREMLFitState::populateAttr(SEXP algebra)
{
	auto *oo = this;
	omxGREMLFitState *gff = (omxGREMLFitState*)this;
	gff->dVupdate_final();
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }

  {
    //Tell the frontend fitfunction counterpart how many observations there are...:
    /*^^^^The end result is that number of observations will be reported as 1 in summary()...
      which is always correct with GREML.  This is a bit of a hack, since it is sneaking this
      negative numObs into the pre-backend fitfunction that summary() looks at...*/
    int userSuppliedDataNumObs = (int)(( (omxGREMLExpectation*)(oo->expectation) )->data2->numObs);
    ProtectedSEXP numObs(Rf_ScalarInteger(1L - userSuppliedDataNumObs));
    Rf_setAttrib(algebra, Rf_install("numObsAdjust"), numObs);
  }

	{
  SEXP mlfitval;
	ScopedProtect p1(mlfitval, Rf_allocVector(REALSXP, 1));
	REAL(mlfitval)[0] = doREML ? (gff->nll - gff->REMLcorrection) : gff->nll;
	Rf_setAttrib(algebra, Rf_install("MLfit"), mlfitval);
	//mxLog("gff->nll: %f",gff->nll);
	//mxLog("gff->REMLcorrection: %f",gff->REMLcorrection);
	}
}


void omxGREMLFitState::buildParamMap(FreeVarGroup *newVarGroup)
{
	if(OMX_DEBUG) { mxLog("Building parameter map for GREML fitfunction."); }
	varGroup = newVarGroup;
	numExplicitFreePar = int(varGroup->vars.size());
	//Now that numExplicitFreePar is known, check to see if derivType=0 and 0 < dVlength < numExplicitFreePar:
	if( (dVlength < numExplicitFreePar || (didUserProvideYhat && dyhatlength>0 && dyhatlength<numExplicitFreePar)) && derivType==0 ){
		//The fitfunction is not allowed to claim to be able to provide a Hessian in this case:
		hessianAvailable = false;
	}
	gradient.setZero(numExplicitFreePar);
	infoMat.setZero(numExplicitFreePar,numExplicitFreePar);
	didUserGivedV.resize(numExplicitFreePar);
	didUserGivedV.assign(size_t(numExplicitFreePar),false);
	didUserGivedyhat.resize(numExplicitFreePar);
	didUserGivedyhat.assign(size_t(numExplicitFreePar),false);
	gradMap.resize(numExplicitFreePar);
	//gradMap2.resize(numExplicitFreePar);
	dAugMap.resize(numExplicitFreePar);
	dAugMap.assign(size_t(numExplicitFreePar),-1);
	indyAlg.resize(numExplicitFreePar);
	indyAlg.assign(size_t(numExplicitFreePar),false);
	indyAlg2.resize(numExplicitFreePar);
	indyAlg2.assign(size_t(numExplicitFreePar),false);
	int gx=0;
	if(dVlength){
		if(!dyhatlength){
			/*^^^This is the case where `dV`, but not `dyhat`, has nonzero length; this block of code is how the parameter-mapping was done
			before `dyhat` was implemented.*/
			if(dVnameslength > numExplicitFreePar){
				mxThrow("length of argument 'dV' is greater than the number of explicit free parameters");
			}
			/*The pointers to the derivatives of V, their names, and their original dimensions get temporarily
			 copied here:*/
			std::vector< omxMatrix* > dV_temp = dV;
			std::vector< const char* > dVnames_temp = dVnames;
			std::vector<int> origdVdim_temp = origdVdim;
			dV.resize(numExplicitFreePar);
			dVnames.resize(numExplicitFreePar);
			origdVdim.resize(numExplicitFreePar);
			/*gx holds the write location for objects with length equal to numExplicitFreePar,
			 vx holds the read location for free parameters in varGroup:*/
			for (int vx=0; vx < numExplicitFreePar; ++vx) {
				//nx holds the read location for objects with length equal to dVlength:
				for (int nx=0; nx <= dVlength; ++nx) {
					if(nx==dVlength){
						gradMap[gx] = vx;
						dV[gx] = NULL;
						dVnames[gx] = NULL;
						origdVdim[gx] = 0;
						//Remember that didUserGivedV was set to all false just above, and dAugMap to all -1s.
						++gx;
						break;
					}
					if (strEQ(dVnames_temp[nx], varGroup->vars[vx]->name)) {
						gradMap[gx] = vx;
						dV[gx] = dV_temp[nx];
						dVnames[gx] = dVnames_temp[nx]; //<--Probably not strictly necessary...
						origdVdim[gx] = origdVdim_temp[nx];
						dAugMap[gx] = nx;
						indyAlg[gx] = ( dV_temp[nx]->algebra && !(dV_temp[nx]->dependsOnParameters()) ) ? true : false;
						didUserGivedV[gx] = true;
						++gx;
						break;
					}
				}
			}
			if (gx != numExplicitFreePar) mxThrow("Problem in dVnames mapping");
		}
		if(dyhatlength){
			if(dVnameslength > numExplicitFreePar){
				mxThrow("length of argument 'dV' is greater than the number of explicit free parameters");
			}
			if(dyhatnameslength > numExplicitFreePar){
				mxThrow("length of argument 'dyhat' is greater than the number of explicit free parameters");
			}
			/*The pointers to the derivatives of V, their names, and their original dimensions get temporarily
			 copied here:*/
			std::vector< omxMatrix* > dV_temp = dV;
			std::vector< omxMatrix* > dyhat_temp = dyhat;
			std::vector< const char* > dVnames_temp = dVnames;
			std::vector< const char* > dNames_temp = dNames;
			std::vector<int> origdVdim_temp = origdVdim;
			std::vector<int> origdyhatdim_temp = origdyhatdim;
			dV.resize(numExplicitFreePar);
			dyhat.resize(numExplicitFreePar);
			dVnames.resize(numExplicitFreePar);
			dNames.resize(numExplicitFreePar);
			origdVdim.resize(numExplicitFreePar);
			origdyhatdim.resize(numExplicitFreePar);
			bool breakFlag1, breakFlag2;
			/*gx holds the write location for objects with length equal to numExplicitFreePar,
			 vx holds the read location for free parameters in varGroup:*/
			for (int vx=0; vx < numExplicitFreePar; ++vx) {
				breakFlag1 = breakFlag2 = false;
				//nx holds the read location for objects with length equal to dVlength:
				for (int nx=0; nx <= dVlength; ++nx) {
					if(nx==dVlength){
						if(!breakFlag1){
							gradMap[gx] = vx;
							dV[gx] = NULL;
							dyhat[gx] = NULL;
							origdVdim[gx] = 0;
							origdyhatdim[gx] = 0;
						}
						if(!breakFlag2){ dVnames[gx] = NULL; }
						//Remember that didUserGivedV was set to all false just above, and dAugMap to all -1s.
						++gx;
						break;
					}
					if (strEQ(dNames_temp[nx], varGroup->vars[vx]->name)) {
						gradMap[gx] = vx;
						dV[gx] = dV_temp[nx];
						dyhat[gx] = dyhat_temp[nx];
						dNames[gx] = dNames_temp[nx]; //<--Probably not strictly necessary...
						origdVdim[gx] = origdVdim_temp[nx];
						origdyhatdim[gx] = origdyhatdim_temp[nx];
						if(dV_temp[nx]){ indyAlg[gx] = ( dV_temp[nx]->algebra && !(dV_temp[nx]->dependsOnParameters()) ) ? true : false; }
						if(dyhat_temp[nx]){ indyAlg2[gx] = ( dyhat_temp[nx]->algebra && !(dyhat_temp[nx]->dependsOnParameters()) ) ? true : false; }
						didUserGivedV[gx] = true;
						didUserGivedyhat[gx] = true;
						breakFlag1 = true;
					}
					if(nx < dVnameslength){
						if(strEQ(dVnames_temp[nx], varGroup->vars[vx]->name)){
							dVnames[gx] = dVnames_temp[nx]; //<--Probably not strictly necessary...
							dAugMap[gx] = nx;
							breakFlag2 = true;
						}
					}
					if(breakFlag1 && breakFlag2){
						++gx;
						break;
					}
				}
			}
			if (gx != numExplicitFreePar) mxThrow("Problem in dVnames mapping");
		}
		if(augGrad){
			int ngradelem = std::max(augGrad->rows, augGrad->cols);
			if(ngradelem != numExplicitFreePar){
				mxThrow("matrix referenced by 'augGrad' must have as many elements as there are explicit free parameters");
			}
			if(augHess){
				if (augHess->rows != augHess->cols) {
					mxThrow("matrix referenced by 'augHess' must be square (instead of %dx%d)",
             augHess->rows, augHess->cols);
				}
				if(augHess->rows != ngradelem){
					mxThrow("Augmentation derivatives non-conformable (gradient is size %d and Hessian is %dx%d)",
             ngradelem, augHess->rows, augHess->cols);
				}
			}
		}
	}
	/*else{
		for(gx=0; gx < numExplicitFreePar; ++gx){
			gradMap[gx] = gx;
		}
	}*/
	if(!dVlength){
		if(augGrad){
			mxThrow("if argument 'dV' has length zero, then so must argument 'augGrad'");
		}
		if(augHess){
			mxThrow("if argument 'dV' has length zero, then so must argument 'augHess'");
		}
		if(dyhatlength){
			if(dyhatnameslength > numExplicitFreePar){
				mxThrow("length of argument 'dyhat' is greater than the number of explicit free parameters");
			}
			gx=0;
			/*The pointers to the derivatives of yhat, their names, and their original dimensions get temporarily
			 copied here:*/
			std::vector< omxMatrix* > dyhat_temp = dyhat;
			std::vector< const char* > dNames_temp = dNames;
			std::vector<int> origdyhatdim_temp = origdyhatdim;
			dyhat.resize(numExplicitFreePar);
			dNames.resize(numExplicitFreePar);
			origdyhatdim.resize(numExplicitFreePar);
			/*gx holds the write location for objects with length equal to numExplicitFreePar,
			 vx holds the read location for free parameters in varGroup:*/
			for (int vx=0; vx < numExplicitFreePar; ++vx) {
				//nx holds the read location for objects with length equal to dyhatlength:
				for (int nx=0; nx <= dyhatlength; ++nx) {
					if(nx==dyhatlength){
						gradMap[gx] = vx;
						dyhat[gx] = NULL;
						dNames[gx] = NULL;
						origdyhatdim[gx] = 0;
						//Remember that didUserGivedyhat was set to all false above.
						++gx;
						break;
					}
					if (strEQ(dNames_temp[nx], varGroup->vars[vx]->name)) {
						gradMap[gx] = vx;
						dyhat[gx] = dyhat_temp[nx];
						dNames[gx] = dNames_temp[nx]; //<--Probably not strictly necessary...
						origdyhatdim[gx] = origdyhatdim_temp[nx];
						indyAlg2[gx] = ( dyhat_temp[nx]->algebra && !(dyhat_temp[nx]->dependsOnParameters()) ) ? true : false;
						didUserGivedyhat[gx] = true;
						++gx;
						break;
					}
				}
			}
			if (gx != numExplicitFreePar) mxThrow("Problem in dyhatnames mapping");
		}
		else{
			for(gx=0; gx < numExplicitFreePar; ++gx){
				gradMap[gx] = gx;
			}
		}
	}
}

void omxGREMLFitState::planParallelDerivs(int nThreadz, int wantHess, int Vrows){
	//Note: AIM = Average Information Matrix (Hessian)
	if( !didUserSpecifyParallelDerivScheme && (wantHess==0 || nThreadz<2 || numExplicitFreePar < 2) ){
		parallelDerivScheme = 1; //Divvy up parameters the old, naive way.
		return;
	}

	/*Under the AIM-row-binning scheme, each thread will calculate the gradient element and the row
	of the AIM (starting with its diagonal element) for each of its parameters.*/

	//Stuff for row binning:
	int i, j, minbin;
	std::vector<int> rownums(numExplicitFreePar,0);
	Eigen::VectorXi rowbinsums(nThreadz);
	rowbinsums.setZero(nThreadz);
	for(i=numExplicitFreePar; i>0; i--){rownums[numExplicitFreePar-i] = i;}

	//Greedy partitioning algorithm to bin rows of the AIM:
	for(i=0; i<numExplicitFreePar; i++){
		minbin=0;
		for(j=1; j<nThreadz; j++){
			if(rowbinsums(j) < rowbinsums(j-1)){minbin=j;}
		}
		rowbins[minbin].conservativeResize(rowbins[minbin].size() + 1);
		rowbins[minbin](rowbins[minbin].size()-1) = rownums[i]-1; //C array indexing starts at 0, not 1.
		rowbinsums(minbin) += rownums[i];
	}

	/*Alternately, we could partition the elements ("cells") of the upper triangle of the AIM among the
	threads.  The elements are 	numbered sequentially, starting with the upper left element,
	moving across each row, and starting again at the diagonal element of the next row.  Under
	this scheme, the elements of the AI matrix are divided among the threads as evenly as possible.*/

	//Stuff for AIM cell binning:
	int numcells = numExplicitFreePar*(numExplicitFreePar+1)/2;
	int targ = int(trunc(numcells/nThreadz));
	targ = (targ<1) ? 1 : targ;
	int remainder = numcells % nThreadz;
	int jlim, cellnum=0;

	//Bin AIM elements:
	if(numcells >= nThreadz){
		for(i=0; i<nThreadz && cellnum<numcells; i++){
			jlim = targ;
			if(remainder){
				jlim++;
				remainder--;
			}
			AIMelembins[i].resize(jlim);
			for(j=0; j<jlim && cellnum<numcells; j++){
				AIMelembins[i](j) = cellnum;
				cellnum++;
			}
		}
	}
	else{
		for(i=0; i<numcells; ++i){
			AIMelembins[i].resize(1);
			AIMelembins[i](0) = i;
		}
	}

	//Stuff for assessing slowest thread under row-binning scheme:
	double N = double(Vrows);
	double diagcost, offdiagcost, inicost;
	if(didUserProvideYhat){
		diagcost = 2*R_pow_di(N,3) + 3*R_pow_di(N,2) + 3*N;
		offdiagcost = R_pow_di(N,3) + 3*R_pow_di(N,2) + N;
	}
	else{
		/*The computational cost of computing a diagonal element includes the upfront cost of
		 computing ytPdV_dtheta, and the cost of computing the gradient element.
		 N^2 for ytPdV_dtheta
		 1.5*N^2 - 0.5*N to efficiently calculate trace of PdV_dtheta
		 N to finish gradient element
		 (N^2) + N for diagonal element:*/
		diagcost = (infoMatType==0) ? 3.5*R_pow_di(N,2) + 1.5*N : R_pow_di(N,3) + 2*R_pow_di(N,2) + N;
		offdiagcost = (infoMatType==0) ? 2*R_pow_di(N,2) + N : R_pow_di(N,3) + R_pow_di(N,2);
	}
	/*workbins will hold the total number of operations each thread will carry out to do
	matrix arithmetic:*/
	Eigen::VectorXd workbins(nThreadz);
	workbins.setZero(nThreadz);
	for(i=0; i<nThreadz; i++){
		for(j=0; j<rowbins[i].size(); j++){
			workbins[i] += diagcost + (rowbins[i](j))*offdiagcost;
		}
	}
	double rowslowest = workbins.maxCoeff();

	//Stuff for assessing slowest thread under cell-binning scheme:
	if(didUserProvideYhat){
		inicost = R_pow_di(N,3);
		diagcost = R_pow_di(N,3) + 3*R_pow_di(N,2) + 3*N;
	}
	else{
		inicost = (infoMatType==0) ? R_pow_di(N,2) : R_pow_di(N,3);
		/*^^^When the thread starts its work, and whenever it moves to a new row of
		 the AIM, it computes ytPdV_dtheta.*/
		/*Thread computes a gradient element whenever it computes a diagonal element
		 of the AIM:*/
		diagcost = (infoMatType==0) ? 2.5*R_pow_di(N,2) + 1.5*N : 2*R_pow_di(N,2) + N;
	}
	workbins.setConstant(nThreadz, inicost);
	int r=0, c=0;
	for(i=0; i<nThreadz; i++){
		for(j=0; j<AIMelembins[i].size(); j++){
			if(r==c){ //If we're at a diagonal element.
				//Add initial cost if we just moved to a new row:
				if(j>0){workbins(i) += inicost;}
				workbins(i) += diagcost;
			}
			else{workbins(i) += offdiagcost;}
			c++;
			/*If we're at the end of the row, then we move to the next one,
			and start at its diagonal element:*/
			if(c>=numExplicitFreePar){
				r++;
				c = r;
			}
		}
	}
	double cellslowest = workbins.maxCoeff();
	
	if(didUserSpecifyParallelDerivScheme){
		parallelDerivScheme = parallelDerivSchemeFromFrontend;
		return;
	}

	parallelDerivScheme = (rowslowest<=cellslowest) ? 2 : 3;
	return;
}


double omxGREMLFitState::pullAugVal(int thing, int row, int col){
	double val=0.0;
	if(row==-1 || col==-1){
		return(val);
	}
	switch(thing){
	case 0:
		if(aug){val = aug->data[0];}
		break;
	case 1:
		if(augGrad){val = augGrad->data[row+col];} //<--Remember that at least one of 'row' and 'col' should be 0.
		break;
	case 2:
		if(augHess){val = omxMatrixElement(augHess,row,col);}
		break;
	}
	return(val);
}


void omxGREMLFitState::recomputeAug(int thing, FitContext *fc){
	switch(thing){
	case 0:
		if(aug){omxRecompute(aug, fc);}
		break;
	case 1:
		if(augGrad){omxRecompute(augGrad, fc);}
		break;
	case 2:
		if(augHess){omxRecompute(augHess, fc);}
		break;
	}
}


void omxGREMLFitState::dVupdate(FitContext *fc){
	for(int i=0; i < numExplicitFreePar; i++){
		if(didUserGivedyhat[i] && dyhat[i]){
			if(OMX_DEBUG){
				mxLog("dyhat %d has matrix number? %s", i, dyhat[i]->hasMatrixNumber ? "True." : "False." );
				mxLog("dyhat %d is clean? %s", i, omxMatrixIsClean(dyhat[i]) ? "True." : "False." );
			}
			//Recompute if needs update and if NOT a parameter-independent algebra:
			if( omxNeedsUpdate(dyhat[i]) && !(indyAlg2[i]) ){
				if(OMX_DEBUG){
					mxLog("Recomputing dyhat %d, %s %s", i, dyhat[i]->getType(), dyhat[i]->name());
				}
				omxRecompute(dyhat[i], fc);
			}
		}
		if(didUserGivedV[i] && dV[i]){
			if(OMX_DEBUG){
				mxLog("dV %d has matrix number? %s", i, dV[i]->hasMatrixNumber ? "True." : "False." );
				mxLog("dV %d is clean? %s", i, omxMatrixIsClean(dV[i]) ? "True." : "False." );
			}
			//Recompute if needs update and if NOT a parameter-independent algebra:
			if( omxNeedsUpdate(dV[i]) && !(indyAlg[i]) ){
				if(OMX_DEBUG){
					mxLog("Recomputing dV %d, %s %s", i, dV[i]->getType(), dV[i]->name());
				}
				omxRecompute(dV[i], fc);
			}
		}
	}
}


void omxGREMLFitState::dVupdate_final(){
	for(int i=0; i < numExplicitFreePar; i++){
		if(didUserGivedyhat[i] && dyhat[i]){
			if(indyAlg2[i]){
				if(OMX_DEBUG){
					mxLog("dyhat %d has matrix number? %s", i, dyhat[i]->hasMatrixNumber ? "True." : "False." );
					mxLog("dyhat %d is clean? %s", i, omxMatrixIsClean(dyhat[i]) ? "True." : "False." );
				}
				if( omxNeedsUpdate(dyhat[i]) ){
					if(OMX_DEBUG){
						mxLog("Recomputing dyhat %d, %s %s", i, dyhat[i]->getType(), dyhat[i]->name());
					}
					omxRecompute(dyhat[i], NULL);
				}
			}
		}
		if(didUserGivedV[i] && dV[i]){
			if(indyAlg[i]){
				if(OMX_DEBUG){
					mxLog("dV %d has matrix number? %s", i, dV[i]->hasMatrixNumber ? "True." : "False." );
					mxLog("dV %d is clean? %s", i, omxMatrixIsClean(dV[i]) ? "True." : "False." );
				}
				if( omxNeedsUpdate(dV[i]) ){
					if(OMX_DEBUG){
						mxLog("Recomputing dV %d, %s %s", i, dV[i]->getType(), dV[i]->name());
					}
					omxRecompute(dV[i], NULL);
				}
			}
		}
	}
}




struct GRMFIMLFitState : omxFitFunction{
	int verbose;
	omxMatrix *y, *X, *invcov, *means;
	bool didUserProvideYhat;

	virtual void init() override;
	virtual void compute2(int want, FitContext *fc) override;
	virtual void populateAttr(SEXP algebra) override;
};

omxFitFunction *GRMFIMLFitInit()
{ return new GRMFIMLFitState; }

void GRMFIMLFitState::init()
{
	auto *oo = this;
	oo->openmpUser = false;
	oo->units = FIT_UNITS_MINUS2LL;

	//If user has provided rowwiseParallel=FALSE to frontend mxFitFunctionML(), then assume they want parallelized
	//numeric derivatives:
	ProtectedSEXP RrowwiseParallel(R_do_slot(rObj, Rf_install("rowwiseParallel")));
	oo->canDuplicate = !Rf_asLogical(RrowwiseParallel);

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose);
	
	omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
	didUserProvideYhat = oge->didUserProvideYhat;

	oo->y = omxGetExpectationComponent(expectation, "y");
	oo->invcov = omxGetExpectationComponent(expectation, "invcov");
	if(didUserProvideYhat){
		oo->means = omxGetExpectationComponent(expectation, "means");
	}
	else{
		oo->X = omxGetExpectationComponent(expectation, "X");
	}

}

void GRMFIMLFitState::compute2(int want, FitContext* fc){
	auto *oo = this;
	if(want & FF_COMPUTE_GRADIENT){invalidateGradient(fc);}
	const double NATLOG_2PI = 1.837877066409345483560659472811;	//<--log(2*pi)
	const double Scale = fabs(Global->llScale);
	omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
	Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(y), y->cols, 1);
	Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(invcov), invcov->rows, invcov->cols);
	//Eigen::Map< Eigen::MatrixXd > Yhat(omxMatrixDataColumnMajor(means), means->rows, 1);
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;
	if(want & FF_COMPUTE_FIT){
		omxExpectationCompute(fc, expectation, NULL);
		//Check for PD covariance matrix:
		if(oge->cholV_fail){
			oo->matrix->data[0] = NA_REAL;
			if (fc) fc->recordIterationError("expected covariance matrix is non-positive-definite");
			return;
		}
		/*Even if the quadratic form in X doesn't get used, if it's non-PD, then the yhats won't have been recalculated
		 at the current parameter values:*/
		if(oge->cholquadX_fail){
			oo->matrix->data[0] = NA_REAL;
			if (fc) fc->recordIterationError("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
			return;
		}
		if(didUserProvideYhat){
			Eigen::Map< Eigen::MatrixXd > Yhat(omxMatrixDataColumnMajor(means), means->rows, 1);
			Eigen::MatrixXd resids(means->rows,1);
			resids = Eigy - Yhat;
			//Compute ML fit:
			oo->matrix->data[0] = Scale*0.5*( (((double)oo->y->cols)*NATLOG_2PI) + oge->logdetV_om->data[0] +
				(resids.transpose() * Vinv.selfadjointView<Eigen::Lower>() * resids)(0,0) );
		}
		else{
			EigenMatrixAdaptor EigX(X);
			double ytPy = (( Eigy.transpose() * Vinv.selfadjointView<Eigen::Lower>() * Eigy ) -
				( Eigy.transpose() * oge->XtVinv.transpose() * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv * Eigy ))(0,0);
			if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
			oo->matrix->data[0] = Scale*0.5*( (((double)oo->y->cols) * NATLOG_2PI) + oge->logdetV_om->data[0] + ytPy);
		}
	}
	return;
}

void GRMFIMLFitState::populateAttr(SEXP algebra){
	//Not really anything to pass to the frontend fitfunction object.
	return;
}
