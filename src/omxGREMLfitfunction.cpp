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
	omxMatrix *y, *X, *cov, *invcov, *means, *origVdim_om;
	std::vector< omxMatrix* > dV;
	std::vector< const char* > dVnames;
	std::vector<int> indyAlg; //will keep track of which algebras don't get marked dirty after dropping cases
	std::vector<int> origdVdim;
	std::vector<int> didUserGivedV;
	void dVupdate(FitContext *fc);
	void dVupdate_final();
	int dVlength, usingGREMLExpectation, parallelDerivScheme, numExplicitFreePar, derivType, oldWantHess, infoMatType;
	double nll, REMLcorrection;
	Eigen::VectorXd gradient;
	Eigen::MatrixXd infoMat; //the Average Information matrix
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
	template <typename T1, typename T2>
	void gradientAndAIM1(
			int _nThreadz, int Eigyrows, FitContext *_fc, int _want, HessianBlock *_hb, omxGREMLExpectation *_oge, Eigen::MatrixBase<T1> &_P,
			double _Scale, Eigen::MatrixBase<T2> &_Py);
	template <typename T1, typename T2, typename T3>
	void gradientAndEIM1(
			int _nThreadz, int Eigyrows, FitContext *_fc, int _want, HessianBlock *_hb, omxGREMLExpectation *_oge, Eigen::MatrixBase<T1> &_P,
			double _Scale, Eigen::MatrixBase<T2> &_Py, Eigen::MatrixBase<T3> &_Eigy);

	template <typename T1, typename T2>
	void crude_numeric_dV(
		FitContext *_fc, Eigen::MatrixBase<T1> &_curEst, Eigen::MatrixBase<T2> &dV_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId);

	omxGREMLFitState() : jg(1, 1) {};
	virtual void init();
	virtual void compute(int want, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
};

struct GREMLSense {
	omxGREMLFitState *fs;
	int parNum, numFree=1, numcases2drop;
	std::vector< int > dropcase;
	Eigen::ArrayXd ref;
	Eigen::MatrixXd result; //<--necessary?
	FitContext *fc;

	GREMLSense(omxGREMLFitState *_fs, FitContext *_fc, int _numcases2drop,  std::vector< int > &_dropcase)
		: fs(_fs), fc(_fc){
		numcases2drop = _numcases2drop;
		dropcase.resize(_dropcase.size());
		for(int dci=0; dci < int(dropcase.size()); dci++){
			dropcase[dci] = _dropcase[dci];
		}
		Eigen::MatrixXd EigV(fs->cov->rows, fs->cov->cols);
		if(numcases2drop && fs->cov->rows > fs->y->cols){
			dropCasesAndEigenize(fs->cov, EigV, numcases2drop, dropcase, 1, fs->cov->rows);
		}
		else{
			EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(fs->cov), fs->cov->rows, fs->cov->cols);
		}
		/*We want to copy the elements of V at the current start values,
		 without copying V's algebra bindings or changing ref's type:*/
		int nrev = EigV.rows();
		ref.resize(nrev*(nrev+1)/2);
		for(int c=0; c < nrev; c++){
			for(int r=c; r < nrev; r++){
				ref( (c*nrev) + r ) = EigV(r,c);
			}
		}
	}

	template <typename T1>
	void operator()(double *myPars, int thrId, Eigen::ArrayBase<T1> &result1) const {
		FitContext *fc2 = thrId >= 0? fc->childList[thrId] : fc;
		fc2->setParamFromOptimizer(parNum, myPars[0]);
		std::vector< int > dropcase2 = dropcase;
		omxState *st = fc2->state;
		omxMatrix *mat = st->lookupDuplicate(fs->cov);
		omxRecompute(mat, fc2);
		Eigen::MatrixXd em(mat->rows, mat->cols);
		if(numcases2drop && mat->rows > fs->y->cols ){
			dropCasesAndEigenize(mat, em, numcases2drop, dropcase2, 1, mat->rows);
		}
		else{
			em = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols);
		}
		int nrem = em.rows();
		result1.resize(nrem*(nrem+1)/2, 1);
		/*We want to copy the elements of mat without copying its algebra bindings or changing result1's type:*/
		for(int c=0; c < nrem; c++){
			for(int r=c; r < nrem; r++){
				result1((c*nrem)+r, 0) = em(r,c);
			}
		}
	}
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
  newObj->usingGREMLExpectation = (strcmp(expectation->name, "MxExpectationGREML")==0 ? 1 : 0);
  if(!newObj->usingGREMLExpectation){
    //Maybe someday GREML fitfunction could be made compatible with another expectation, but not at present:
    mxThrow("GREML fitfunction is currently only compatible with GREML expectation");
  }
  else{
    omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
    oge->alwaysComputeMeans = 0;
  }

  newObj->y = omxGetExpectationComponent(expectation, "y");
  newObj->cov = omxGetExpectationComponent(expectation, "cov");
  newObj->invcov = omxGetExpectationComponent(expectation, "invcov");
  newObj->X = omxGetExpectationComponent(expectation, "X");
  newObj->means = omxGetExpectationComponent(expectation, "means");
  newObj->origVdim_om = omxGetExpectationComponent(expectation, "origVdim_om");
  newObj->nll = 0;
  newObj->REMLcorrection = 0;
  newObj->parallelDerivScheme = 0;
  newObj->varGroup = NULL;
  newObj->augGrad = NULL;
  newObj->augHess = NULL;
  newObj->dVlength = 0;
  newObj->oldWantHess = 0;

	//autoDerivType:
  {
  	ProtectedSEXP adt(R_do_slot(rObj, Rf_install("autoDerivType")));
  	if(strEQ(CHAR(Rf_asChar(adt)),"semiAnalyt")){derivType = 1;}
  	else if(strEQ(CHAR(Rf_asChar(adt)),"numeric")){derivType = 0;}
  	else{mxThrow("unrecognized character string provided for GREML fitfunction 'autoDerivType'");}
  }
  
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

  //Derivatives of V:
  if (R_has_slot(rObj, Rf_install("dV"))) {
  	ProtectedSEXP RdV(R_do_slot(rObj, Rf_install("dV")));
  	ProtectedSEXP RdVnames(R_do_slot(rObj, Rf_install("dVnames")));
  	newObj->dVlength = Rf_length(RdV);
  	newObj->dV.resize(newObj->dVlength);
  	newObj->indyAlg.resize(newObj->dVlength);
  	newObj->dVnames.resize(newObj->dVlength);
  	newObj->origdVdim.resize(newObj->dVlength);
  	if(newObj->dVlength){
  		if(!newObj->usingGREMLExpectation){
  			//Probably best not to allow use of dV if we aren't sure means will be calculated GREML-GLS way:
  			mxThrow("derivatives of 'V' matrix in GREML fitfunction only compatible with GREML expectation");
  		}
  		if(OMX_DEBUG) { mxLog("Processing derivatives of V."); }
  		int* dVint = INTEGER(RdV);
  		for(int i=0; i < newObj->dVlength; i++){
  			newObj->dV[i] = omxMatrixLookupFromStateByNumber(dVint[i], currentState);
  			SEXP elem;
  			{ScopedProtect p3(elem, STRING_ELT(RdVnames, i));
  				newObj->dVnames[i] = CHAR(elem);}
  		}
  	}
  }

  if(derivType==1 && !newObj->usingGREMLExpectation){
  	mxThrow("semi-analytic derivatives only compatible with GREML expectation");
  }

  if(newObj->dVlength || derivType==1){
    oo->hessianAvailable = true;
    newObj->rowbins.resize(Global->numThreads);
    newObj->AIMelembins.resize(Global->numThreads);
    for(int i=0; i < newObj->dVlength; i++){
    	/*Each dV must either (1) match the dimensions of V, OR (2) match the length of y if that is less than the
    	dimension of V (implying downsizing due to missing observations):*/
      if( ((newObj->dV[i]->rows == newObj->cov->rows)&&(newObj->dV[i]->cols == newObj->cov->cols)) ||
          ((newObj->y->cols < newObj->cov->rows)&&(newObj->dV[i]->rows == newObj->y->cols)&&
          	(newObj->dV[i]->cols == newObj->y->cols)) ){
      	newObj->origdVdim[i] = newObj->dV[i]->rows;
      }
      else{
        mxThrow("all derivatives of V must have the same dimensions as V");
  }}}

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
  			else{
  				mxThrow("if using semi-analytic derivatives and 'aug' has nonzero length, then 'augGrad' must as well");
  			}
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

void omxGREMLFitState::compute(int want, FitContext *fc)
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
 	int i;
 	Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(gff->y), gff->y->cols, 1);
 	Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(gff->invcov), gff->invcov->rows, gff->invcov->cols);
 	EigenMatrixAdaptor EigX(gff->X);
 	Eigen::MatrixXd P, Py;
 	double logdetV=0, logdetquadX=0, ytPy=0;

 	if(want & (FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 		if(gff->usingGREMLExpectation){
 			omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);

 			//Check that factorizations of V and the quadratic form in X succeeded:
 			if(oge->cholV_fail_om->data[0]){
 				oo->matrix->data[0] = NA_REAL;
 				if (fc) fc->recordIterationError("expected covariance matrix is non-positive-definite");
 				return;
 			}
 			if(oge->cholquadX_fail){
 				oo->matrix->data[0] = NA_REAL;
 				if (fc) fc->recordIterationError("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
 				return;
 			}

 			//Log determinant of V:
 			logdetV = oge->logdetV_om->data[0];

 			//Log determinant of quadX:
 			for(i=0; i < gff->X->cols; i++){
 				logdetquadX += log(oge->cholquadX_vectorD[i]);
 			}
 			logdetquadX *= 2;
 			gff->REMLcorrection = Scale*0.5*logdetquadX;

 			/*Finish computing fit (negative loglikelihood) if wanted.  P and Py will be needed later if analytic derivatives in use;
 			otherwise, extraneous calculations can be avoided:*/
 			if(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 				P.setZero(gff->invcov->rows, gff->invcov->cols);
 				P.triangularView<Eigen::Lower>() = (Vinv.selfadjointView<Eigen::Lower>() * //P = Vinv * (I-Hatmat)
 					(Eigen::MatrixXd::Identity(Vinv.rows(), Vinv.cols()) -
 					(EigX * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv))).triangularView<Eigen::Lower>();
 				Py = P.selfadjointView<Eigen::Lower>() * Eigy;
 				P.triangularView<Eigen::Upper>() = P.triangularView<Eigen::Lower>().transpose();
 				if(want & FF_COMPUTE_FIT){
 					ytPy = (Eigy.transpose() * Py)(0,0);
 					if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
 					oo->matrix->data[0] = gff->REMLcorrection +
 						Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + ytPy) + Scale*gff->pullAugVal(0L,0,0);
 					gff->nll = oo->matrix->data[0];
 					if(OMX_DEBUG){mxLog("augmentation is %3.3f",gff->pullAugVal(0L,0,0));}
 				}
 			}
 			/*ytPy can be calculated so that rate-limiting step is O(2kN^2), where k is the number of covariates,
 			and N is the dimension of Vinv (and typically N>>k):*/
 			else{
 				ytPy = (( Eigy.transpose() * Vinv.selfadjointView<Eigen::Lower>() * Eigy ) -
 					( Eigy.transpose() * oge->XtVinv.transpose() * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv * Eigy ))(0,0);
 				if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
 				oo->matrix->data[0] = gff->REMLcorrection +
 					Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + ytPy) + Scale*gff->pullAugVal(0L,0,0);
 				gff->nll = oo->matrix->data[0];
 				if(OMX_DEBUG){mxLog("augmentation is %3.3f",gff->pullAugVal(0L,0,0));}
 				return; //<--Since only fit value is wanted.
 			}
 		}
 		else{ //If not using GREML expectation, deal with means and cov in a general way to compute fit...
 			//Declare locals:
 			EigenMatrixAdaptor yhat(gff->means);
 			EigenMatrixAdaptor EigV(gff->cov);
 			double logdetV=0, logdetquadX=0;
 			Eigen::MatrixXd Vinv, quadX;
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
 			for(i=0; i < gff->X->rows; i++){
 				logdetV += log(cholV_vectorD[i]);
 			}
 			logdetV *= 2;

 			Vinv = cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse

 			quadX = EigX.transpose() * Vinv * EigX; //<--Quadratic form in X

 			cholquadX.compute(quadX); //<--Cholesky factorization of quadX
 			if(cholquadX.info() != Eigen::Success){
 				omxRaiseErrorf("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
 				oo->matrix->data[0] = NA_REAL;
 				return;
 			}
 			cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal();
 			for(i=0; i < gff->X->cols; i++){
 				logdetquadX += log(cholquadX_vectorD[i]);
 			}
 			logdetquadX *= 2;
 			gff->REMLcorrection = Scale*0.5*logdetquadX;

 			//Finish computing fit:
 			oo->matrix->data[0] = gff->REMLcorrection + Scale*0.5*( ((double)gff->y->rows * NATLOG_2PI) + logdetV +
 				( Eigy.transpose() * Vinv * (Eigy - yhat) )(0,0));
 			gff->nll = oo->matrix->data[0];
 			return;
 		}
 	}

 	if(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
 		//This part requires GREML expectation:
 		omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);

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
 				omxRaiseErrorf("GREML fitfunction cannot compute information matrix without analytic derivatives of V with respect to EVERY free parameter");
 			}
 			hb->vars.resize(gff->numExplicitFreePar);
 			hb->mat.resize(gff->numExplicitFreePar, gff->numExplicitFreePar);
 			gff->recomputeAug(2, fc);
 			wantHess = 1;
 		}
 		
 		if(oldWantHess != wantHess){
 			parallelDerivScheme = 0;
 			oldWantHess = wantHess;
 		}

 		if(gff->parallelDerivScheme==0){gff->planParallelDerivs(nThreadz,wantHess,gff->cov->rows);}

		/*If datapoints need to be dropped due to missingness, and we want the analytic Hessian,
 		then we need to resize any derivatives of V that come from front-end MxAlgebras
		 ahead of time,	in order to assure thread-safety of the parallelized code for
 		evaluating the gradient and Hessian (AIM):*/
		if(oge->numcases2drop && wantHess){
			for(i=0; i < numExplicitFreePar; i++){
				if(didUserGivedV[i] &&	gff->dV[i]->rows > Eigy.rows()){
					dropCasesFromAlgdV(gff->dV[i], oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[i]);
				}
			}
		}
		
 		//Begin parallelized evaluation of fitfunction derivatives:
 		switch(gff->parallelDerivScheme){
 		case 2: //bin by row
#pragma omp parallel num_threads(nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, r=0, c=0;
		double tr=0;
		Eigen::MatrixXd ytPdV_dtheta1;
		Eigen::MatrixXd dV_dtheta1(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter hrn.
		Eigen::MatrixXd dV_dtheta2(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter hcn.
		int threadID = omx_absolute_thread_num();
		int istart = 0;
		int iend = gff->rowbins[threadID].size();
		for(i=istart; i < iend; i++){
			tr=0;
			hrn = gff->rowbins[threadID](i); //Current row number of the AIM.
			if(gff->gradMap[hrn] < 0){continue;} //Check for negative parameter number.
			a1 = gff->dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
			if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){hb->vars[hrn] = gff->gradMap[hrn];}
			if( oge->numcases2drop && (gff->dV[hrn]->rows > Eigy.rows()) ){
				dropCasesAndEigenize(gff->dV[hrn], dV_dtheta1, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[hrn]);
			}
			else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[hrn]), gff->dV[hrn]->rows, gff->dV[hrn]->cols);}
			ytPdV_dtheta1 = Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
			for(hcn=hrn; hcn < gff->dVlength; hcn++){
				if(hcn==hrn){
					/*Need trace of P*dV_dtheta for gradient element...
					 Frustratingly, the selfadjointView has no row or column accessor function among its members.
					 But the trace of a product of two square symmetric matrices is the sum of the elements of
					 their elementwise product.*/
					for(c=0; c < gff->cov->rows; c++){
						for(r=c; r < gff->cov->rows; r++){
							tr += (r==c) ? P(r,c)*dV_dtheta1(r,c) : 2*P(r,c)*dV_dtheta1(r,c);
						}
					}
					gff->gradient(hrn) = Scale*0.5*(tr - (ytPdV_dtheta1 * Py)(0,0)) +
						Scale*gff->pullAugVal(1,a1,0);
					fc->gradZ(hrn) += gff->gradient(hrn);
					if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						gff->infoMat(hrn,hrn) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
							Scale*gff->pullAugVal(2,a1,a1);
					}
				}
				//I think it can be assumed at this point that the Hessian is wanted?:
				else{if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
					if(gff->gradMap[hcn] < 0){continue;}
					a2 = gff->dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
					if( oge->numcases2drop && (gff->dV[hcn]->rows > Eigy.rows()) ){
						dropCasesAndEigenize(gff->dV[hcn], dV_dtheta2, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[hcn]);
					}
					else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[hcn]), gff->dV[hcn]->rows, gff->dV[hcn]->cols);}
					gff->infoMat(hrn,hcn) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() *
						dV_dtheta2.selfadjointView<Eigen::Lower>() * Py)(0,0) + Scale*gff->pullAugVal(2,a1,a2);
					gff->infoMat(hcn,hrn) = gff->infoMat(hrn,hcn);
				}}}}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
 			break;
 		case 3: //bin by cell
#pragma omp parallel num_threads(nThreadz)
{
	try{
		int i=0, hrn=0, hcn=0, a1=0, a2=0, r=0, c=0, inielem=0;
		double tr=0;
		Eigen::MatrixXd ytPdV_dtheta1;
		Eigen::MatrixXd dV_dtheta1(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter hrn.
		Eigen::MatrixXd dV_dtheta2(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter hcn.
		int threadID = omx_absolute_thread_num();
		int iend = gff->AIMelembins[threadID].size();
		if(iend){inielem = gff->AIMelembins[threadID](0);}
		while(inielem > 0){
			hcn++;
			inielem--;
			if(hcn == gff->dVlength){
				hrn++;
				hcn=hrn;
			}
		}
		while(i < iend){
			if(gff->gradMap[hrn] < 0){continue;} //Check for negative parameter number.
			tr=0;
			a1 = gff->dAugMap[hrn]; //<--Index of augmentation derivatives to use for parameter hrn.
			if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){hb->vars[hrn] = gff->gradMap[hrn];}
			if(hrn==hcn || i==0){
				if( oge->numcases2drop && (gff->dV[hrn]->rows > Eigy.rows()) ){
					dropCasesAndEigenize(gff->dV[hrn], dV_dtheta1, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[hrn]);
				}
				else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[hrn]), gff->dV[hrn]->rows, gff->dV[hrn]->cols);}
				ytPdV_dtheta1 = Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
			}
			if(hrn==hcn){
				for(c=0; c < gff->cov->rows; c++){
					for(r=c; r < gff->cov->rows; r++){
						tr += (r==c) ? P(r,c)*dV_dtheta1(r,c) : 2*P(r,c)*dV_dtheta1(r,c);
					}
				}
				gff->gradient(hrn) = Scale*0.5*(tr - (ytPdV_dtheta1 * Py)(0,0)) +
					Scale*gff->pullAugVal(1,a1,0);
				fc->gradZ(hrn) += gff->gradient(hrn);
				if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
					gff->infoMat(hrn,hrn) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
						Scale*gff->pullAugVal(2,a1,a1);
				}
			}
			else{if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
				if(gff->gradMap[hcn] < 0){continue;}
				a2 = gff->dAugMap[hcn]; //<--Index of augmentation derivatives to use for parameter hcn.
				if( oge->numcases2drop && (gff->dV[hcn]->rows > Eigy.rows()) ){
					dropCasesAndEigenize(gff->dV[hcn], dV_dtheta2, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[hcn]);
				}
				else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[hcn]), gff->dV[hcn]->rows, gff->dV[hcn]->cols);}
				gff->infoMat(hrn,hcn) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() *
					dV_dtheta2.selfadjointView<Eigen::Lower>() * Py)(0,0) + Scale*gff->pullAugVal(2,a1,a2);
				gff->infoMat(hcn,hrn) = gff->infoMat(hrn,hcn);
			}}
			hcn++;
			i++;
			if(hcn == gff->dVlength){
				hrn++;
				hcn=hrn;
			}}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
 			break;
 		default: //bin naively (which is perfectly adequate for gradient-only, or for a single thread)
 			if(infoMatType==1){
 				gradientAndEIM1(nThreadz, Eigy.rows(), fc, want, hb, oge, P, Scale, Py, Eigy);
 			}
 			else{
#pragma omp parallel num_threads(nThreadz)
{
	try{
		int i=0, j=0, t1=0, t2=0, a1=0, a2=0, r=0, c=0;
		double tr=0;
		Eigen::MatrixXd ytPdV_dtheta1;
		Eigen::MatrixXd dV_dtheta1(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter i.
		//GREMLSense sense(this, fc, oge->numcases2drop, oge->dropcase);
		Eigen::VectorXd curEst(numExplicitFreePar);
		fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = threadID * numExplicitFreePar / nThreadz;
		int iend = (threadID+1) * numExplicitFreePar / nThreadz;
		if(threadID == nThreadz-1){iend = numExplicitFreePar;}
		for(i=istart; i < iend; i++){
			tr=0;
			t1 = gff->gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			if(didUserGivedV[t1] || derivType==1){
				a1 = gff->dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
				if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){hb->vars[i] = t1;}
				if(didUserGivedV[t1]){
					if( oge->numcases2drop && (gff->dV[i]->rows > Eigy.rows()) ){
						dropCasesAndEigenize(gff->dV[i], dV_dtheta1, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[i]);
					}
					else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[i]), gff->dV[i]->rows, gff->dV[i]->cols);}
				}
				else{
					/*
					sense.parNum = t1;
					curEst1p[0] = curEst[t1];
					Eigen::MatrixXd tmpresult(sense.ref.size(), 1);
					jg(sense, sense.ref, curEst1p, false, tmpresult);
					int dVdim = int(sqrt(tmpresult.rows() * 2.0));
					dV_dtheta1.resize(dVdim,dVdim);
					for(int cc=0; cc < dVdim; cc++){
						for(int rr=cc; rr < dVdim; rr++){
							dV_dtheta1(rr,cc) = tmpresult((cc*dVdim)+rr, 0);
						}
					}
					*/
					crude_numeric_dV(fc, curEst, dV_dtheta1, t1, oge, (nThreadz>1 ? threadID : -1));
				}
				ytPdV_dtheta1 = Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
				for(j=i; j < numExplicitFreePar; j++){
					if(j==i){
						/*Need trace of P*dV_dtheta for gradient element...
						 Frustratingly, the selfadjointView has no row or column accessor function among its members.
						 But the trace of a product of two square symmetric matrices is the sum of the elements of
						 their elementwise product.*/
						for(c=0; c < gff->cov->rows; c++){
							for(r=c; r < gff->cov->rows; r++){
								tr += (r==c) ? P(r,c)*dV_dtheta1(r,c) : 2*P(r,c)*dV_dtheta1(r,c);
							}
						}
						gff->gradient(t1) = Scale*0.5*(tr - (ytPdV_dtheta1 * Py)(0,0)) +
							Scale*gff->pullAugVal(1,a1,0);
						fc->gradZ(t1) += gff->gradient(t1);
						if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							gff->infoMat(t1,t1) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) +
								Scale*gff->pullAugVal(2,a1,a1);
						}
					}
					else{
						if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
							Eigen::MatrixXd dV_dtheta2(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter j.
							t2 = gff->gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							a2 = gff->dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
							if(didUserGivedV[t2]){
								if( oge->numcases2drop && (gff->dV[j]->rows > Eigy.rows()) ){
									dropCasesAndEigenize(gff->dV[j], dV_dtheta2, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[j]);
								}
								else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[j]), gff->dV[j]->rows, gff->dV[j]->cols);}
							}
							else{
								/*
								sense.parNum = t2;
								curEst1p[0] = curEst[t2];
								Eigen::MatrixXd tmpresult(sense.ref.size(), 1);
								jg(sense, sense.ref, curEst1p, false, tmpresult);
								int dVdim = int(sqrt(tmpresult.rows() * 2.0));
								dV_dtheta2.resize(dVdim,dVdim);
								for(int cc=0; cc < dVdim; cc++){
									for(int rr=cc; rr < dVdim; rr++){
										dV_dtheta2(rr,cc) = tmpresult((cc*dVdim)+rr, 0);
									}
								}
								*/
								crude_numeric_dV(fc, curEst, dV_dtheta2, t2, oge, (nThreadz>1 ? threadID : -1));
							}
							gff->infoMat(t1,t2) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() *
								dV_dtheta2.selfadjointView<Eigen::Lower>() * Py)(0,0) + Scale*gff->pullAugVal(2,a1,a2);
							gff->infoMat(t2,t1) = gff->infoMat(t1,t2);
						}
					}
				}
			}
			else{
				gff->gradient(t1) = NA_REAL;
				fc->gradZ(t1) = NA_REAL;
			}
		}
	} catch (const std::exception& e) {
		omxRaiseErrorf("%s", e.what());
	} catch (...) {
		omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
	}
}
} 		
 		break;
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
 	}
 	return;
 }


//Foundation for separate functions that compute fitfunction derivatives from dV:
template <typename T1, typename T2>
void omxGREMLFitState::gradientAndAIM1(
		int _nThreadz, int Eigyrows, FitContext *_fc, int _want, HessianBlock *_hb, omxGREMLExpectation *_oge, Eigen::MatrixBase<T1> &_P,
    double _Scale, Eigen::MatrixBase<T2> &_Py){
	return;
}
template <typename T1, typename T2, typename T3>
void omxGREMLFitState::gradientAndEIM1(
		int _nThreadz, int Eigyrows, FitContext *_fc, int _want, HessianBlock *_hb, omxGREMLExpectation *_oge, Eigen::MatrixBase<T1> &_P,
		double _Scale, Eigen::MatrixBase<T2> &_Py, Eigen::MatrixBase<T3> &_Eigy){
#pragma omp parallel num_threads(_nThreadz)
{
	try{
		int i=0, t1=0, a1=0, r=0, c=0;
		double tr1=0;
		Eigen::MatrixXd dV_dtheta1(Eigyrows, Eigyrows); //<--Derivative of V w/r/t parameter i.
		//GREMLSense sense(this, _fc, _oge->numcases2drop, _oge->dropcase);
		Eigen::VectorXd curEst(numExplicitFreePar);
		_fc->copyEstToOptimizer(curEst);
		Eigen::VectorXd curEst1p(1);
		int threadID = omx_absolute_thread_num();
		int istart = threadID * numExplicitFreePar / _nThreadz;
		int iend = (threadID+1) * numExplicitFreePar / _nThreadz;
		if(threadID == _nThreadz-1){iend = numExplicitFreePar;}
		for(i=istart; i < iend; i++){
			tr1=0;
			t1 = gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			//If the Hessian isn't wanted, then we compute the gradient elements the fast way.
			//If the Hessian is wanted, then we compute the gradient in a slower way, but 
			//which saves time when computing the expected information matrix:
			if(_want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
				Eigen::MatrixXd dV_dtheta2(Eigyrows, Eigyrows), dV_dtheta1P, dV_dtheta2P;
				int j=0, t2=0, a2=0;
				double tr2=0;
				_hb->vars[i] = t1;
				if(didUserGivedV[t1] || derivType==1){
					a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
					if(didUserGivedV[t1]){
						if( _oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
							dropCasesAndEigenize(dV[i], dV_dtheta1, _oge->numcases2drop, _oge->dropcase, 1, origdVdim[i]);
						}
						else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(dV[i]), dV[i]->rows, dV[i]->cols);}
					}
					else{
						/*
						 sense.parNum = t1;
						 curEst1p[0] = curEst[t1];
						 Eigen::MatrixXd tmpresult(sense.ref.size(), 1);
						 jg(sense, sense.ref, curEst1p, false, tmpresult);
						 int dVdim = int(sqrt(tmpresult.rows() * 2.0));
						 dV_dtheta1.resize(dVdim,dVdim);
						 for(int cc=0; cc < dVdim; cc++){
						 for(int rr=cc; rr < dVdim; rr++){
						 dV_dtheta1(rr,cc) = tmpresult((cc*dVdim)+rr, 0);
						 }
						 }
						 */
						crude_numeric_dV(_fc, curEst, dV_dtheta1, t1, _oge, (_nThreadz>1 ? threadID : -1));
					}
					//Eigen does not define a multiplication operator for selfadjointView times selfadjointView:
					dV_dtheta1P = dV_dtheta1.selfadjointView<Eigen::Lower>() * _P;
					tr1 = dV_dtheta1P.trace();
					for(j=i; j < numExplicitFreePar; j++){
						tr2=0;
						if(j==i){
							gradient(t1) = _Scale*0.5*(tr1 - (_Py.transpose() * dV_dtheta1P * _Eigy)(0,0)) +
								_Scale*pullAugVal(1,a1,0);
							_fc->gradZ(t1) += gradient(t1);
							for(c=0; c < cov->rows; c++){
								for(r=0; r < cov->rows; r++){
									tr2 += dV_dtheta1P(r,c) * dV_dtheta1P(r,c);
								}
							}
							infoMat(t1,t1) = _Scale*0.5*tr2 + _Scale*pullAugVal(2,a1,a1);
						}
						else{
							t2 = gradMap[j]; //<--Parameter number for parameter j.
							if(t2 < 0){continue;}
							a2 = dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
							if(didUserGivedV[t2]){
								if( _oge->numcases2drop && (dV[j]->rows > Eigyrows) ){
									dropCasesAndEigenize(dV[j], dV_dtheta2, _oge->numcases2drop, _oge->dropcase, 1, origdVdim[j]);
								}
								else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(dV[j]), dV[j]->rows, dV[j]->cols);}
							}
							else{
								/*
								 sense.parNum = t2;
								 curEst1p[0] = curEst[t2];
								 Eigen::MatrixXd tmpresult(sense.ref.size(), 1);
								 jg(sense, sense.ref, curEst1p, false, tmpresult);
								 int dVdim = int(sqrt(tmpresult.rows() * 2.0));
								 dV_dtheta2.resize(dVdim,dVdim);
								 for(int cc=0; cc < dVdim; cc++){
								 for(int rr=cc; rr < dVdim; rr++){
								 dV_dtheta2(rr,cc) = tmpresult((cc*dVdim)+rr, 0);
								 }
								 }
								 */
								crude_numeric_dV(_fc, curEst, dV_dtheta2, t2, _oge, (_nThreadz>1 ? threadID : -1));
							}
							dV_dtheta2P = dV_dtheta2.selfadjointView<Eigen::Lower>() * _P;//.template selfadjointView<Eigen::Lower>();
							for(c=0; c < cov->rows; c++){
								for(r=0; r < cov->rows; r++){
									tr2 += dV_dtheta1P(r,c) * dV_dtheta2P(r,c);
								}
							}
							infoMat(t1,t2) = _Scale*0.5*tr2 + _Scale*pullAugVal(2,a1,a2);
							infoMat(t2,t1) = infoMat(t1,t2);
						}
					}
				}
				else{
					gradient(t1) = NA_REAL;
					_fc->gradZ(t1) = NA_REAL;
				}
			}
			
			else{
				Eigen::MatrixXd ytPdV_dtheta1;
				if(didUserGivedV[t1] || derivType==1){
					a1 = dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
					if(didUserGivedV[t1]){
						if( _oge->numcases2drop && (dV[i]->rows > Eigyrows) ){
							dropCasesAndEigenize(dV[i], dV_dtheta1, _oge->numcases2drop, _oge->dropcase, 1, origdVdim[i]);
						}
						else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(dV[i]), dV[i]->rows, dV[i]->cols);}
					}
					else{
						/*
						 sense.parNum = t1;
						 curEst1p[0] = curEst[t1];
						 Eigen::MatrixXd tmpresult(sense.ref.size(), 1);
						 jg(sense, sense.ref, curEst1p, false, tmpresult);
						 int dVdim = int(sqrt(tmpresult.rows() * 2.0));
						 dV_dtheta1.resize(dVdim,dVdim);
						 for(int cc=0; cc < dVdim; cc++){
						 for(int rr=cc; rr < dVdim; rr++){
						 dV_dtheta1(rr,cc) = tmpresult((cc*dVdim)+rr, 0);
						 }
						 }
						 */
						crude_numeric_dV(_fc, curEst, dV_dtheta1, t1, _oge, (_nThreadz>1 ? threadID : -1));
					}
					ytPdV_dtheta1 = _Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
					/*Need trace of P*dV_dtheta for gradient element...
					 Frustratingly, the selfadjointView has no row or column accessor function among its members.
					 But the trace of a product of two square symmetric matrices is the sum of the elements of
					 their elementwise product.*/
					for(c=0; c < cov->rows; c++){
						for(r=c; r < cov->rows; r++){
							tr1 += (r==c) ? _P(r,c)*dV_dtheta1(r,c) : 2*_P(r,c)*dV_dtheta1(r,c);
						}
					}
					gradient(t1) = _Scale*0.5*(tr1 - (ytPdV_dtheta1 * _Py)(0,0)) +
						_Scale*pullAugVal(1,a1,0);
					_fc->gradZ(t1) += gradient(t1);
				}
				else{
					gradient(t1) = NA_REAL;
					_fc->gradZ(t1) = NA_REAL;
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


template <typename T1, typename T2>
void omxGREMLFitState::crude_numeric_dV(
		FitContext *_fc, Eigen::MatrixBase<T1> &_curEst, Eigen::MatrixBase<T2> &dV_dtheta, int Parnum, omxGREMLExpectation *ge, int thrId)
{
	int c, r;
	FitContext *fc2 = _fc;
	if(thrId >= 0){
		fc2 = _fc->childList[thrId];
	}
	omxState *st = fc2->state;

	//Store current elements of V:
	Eigen::MatrixXd EigV(cov->rows, cov->cols);
	if(ge->numcases2drop && cov->rows > y->cols){
		dropCasesAndEigenize(cov, EigV, ge->numcases2drop, ge->dropcase, 1, cov->rows);
	}
	else{
		EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(cov), cov->rows, cov->cols);
	}
	Eigen::MatrixXd Vcurr(EigV.rows(), EigV.cols());
	for(c=0; c < EigV.cols(); c++){
		for(r=c; r < EigV.rows(); r++){
			Vcurr(r,c) = EigV(r,c);
		}
	}

	//Shift parameter of interest and compute V at new point in parameter space:
	_curEst[Parnum] += 1e-4;
	fc2->setEstFromOptimizer(_curEst);
	omxMatrix *mat = st->lookupDuplicate(cov);
	omxRecompute(mat, fc2);

	if( ge->numcases2drop && mat->rows > y->cols ){
		dropCasesAndEigenize(mat, dV_dtheta, ge->numcases2drop, ge->dropcase, 1, mat->rows);
	}
	else{
		dV_dtheta = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols);
	}

	//Does this actually save memory...?:
	dV_dtheta.template triangularView<Eigen::Lower>() = (dV_dtheta - Vcurr)/1e-4;

	//Put things back the way they were:
	_curEst[Parnum] -= 1e-4;
	fc2->setEstFromOptimizer(_curEst);
	if(thrId<0){
		omxRecompute(cov, fc2);
	}
	else{
		EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(cov), cov->rows, cov->cols);
		for(c=0; c < cov->rows; c++){
			for(r=c; r < cov->rows; r++){
				EigV(r,c) = Vcurr(r,c);
			}
		}
	}
	return;
}



void omxGREMLFitState::populateAttr(SEXP algebra)
{
	auto *oo = this;
	omxGREMLFitState *gff = (omxGREMLFitState*)this;
	gff->dVupdate_final();
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }

  SEXP nval, mlfitval;
  int userSuppliedDataNumObs = (int)(( (omxGREMLExpectation*)(oo->expectation) )->data2->numObs);

  //Tell the frontend fitfunction counterpart how many observations there are...:
  {
  ScopedProtect p1(nval, Rf_allocVector(INTSXP, 1));
  INTEGER(nval)[0] = 1L - userSuppliedDataNumObs;
  R_do_slot_assign(rObj, Rf_install("numObs"), nval);
  /*^^^^The end result is that number of observations will be reported as 1 in summary()...
  which is always correct with GREML.  This is a bit of a hack, since it is sneaking this
  negative numObs into the pre-backend fitfunction that summary() looks at...*/
	}

	{
	ScopedProtect p1(mlfitval, Rf_allocVector(REALSXP, 1));
	REAL(mlfitval)[0] = gff->nll - gff->REMLcorrection;
	Rf_setAttrib(algebra, Rf_install("MLfit"), mlfitval);
	}
}


void omxGREMLFitState::buildParamMap(FreeVarGroup *newVarGroup)
{
	if(OMX_DEBUG) { mxLog("Building parameter map for GREML fitfunction."); }
	varGroup = newVarGroup;
	numExplicitFreePar = int(varGroup->vars.size());
	gradient.setZero(numExplicitFreePar);
	infoMat.setZero(numExplicitFreePar,numExplicitFreePar);
	didUserGivedV.resize(numExplicitFreePar);
	didUserGivedV.assign(size_t(numExplicitFreePar),0);
	gradMap.resize(numExplicitFreePar);
	dAugMap.resize(numExplicitFreePar);
	dAugMap.assign(size_t(numExplicitFreePar),-1);
	indyAlg.resize(numExplicitFreePar);
	indyAlg.assign(size_t(numExplicitFreePar),0);
	int gx=0;
	if(dVlength){
		if(dVlength > numExplicitFreePar){
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
					//Remember that didUserGivedV was set to all zeroes just above, and dAugMap to all -1s.
					++gx;
					break;
				}
				if (strEQ(dVnames_temp[nx], varGroup->vars[vx]->name)) {
					gradMap[gx] = vx;
					dV[gx] = dV_temp[nx];
					dVnames[gx] = dVnames_temp[nx]; //<--Probably not strictly necessary...
					origdVdim[gx] = origdVdim_temp[nx];
					dAugMap[gx] = nx;
					indyAlg[gx] = ( dV_temp[nx]->algebra && !(dV_temp[nx]->dependsOnParameters()) ) ? 1 : 0;
					didUserGivedV[gx] = 1;
					++gx;
					break;
				}
			}
		}
		if (gx != numExplicitFreePar) mxThrow("Problem in dVnames mapping");
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
			}}
		}
	}
	else{
		for(gx=0; gx < numExplicitFreePar; ++gx){
			gradMap[gx] = gx;
		}
	}
}


void omxGREMLFitState::planParallelDerivs(int nThreadz, int wantHess, int Vrows){
	//Note: AIM = Average Information Matrix (Hessian)
	if(wantHess==0 || nThreadz<2 || numExplicitFreePar==1 || dVlength < numExplicitFreePar || infoMatType==1){
		//^^^TODO: relax the last two conditions.
		parallelDerivScheme = 1; //Divvy up parameters the old, naive way.
		return;
	}

	/*Under the AIM-row-binning scheme, each thread will calculate the gradient element and the row
	of the AIM (starting with its diagonal element) for each of its parameters.*/

	//Stuff for row binning:
	int i, j, minbin;
	std::vector<int> rownums(dVlength,0);
	Eigen::VectorXi rowbinsums(nThreadz);
	rowbinsums.setZero(nThreadz);
	for(i=dVlength; i>0; i--){rownums[dVlength-i] = i;}

	//Greedy partitioning algorithm to bin rows of the AIM:
	for(i=0; i<dVlength; i++){
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
	int numcells = dVlength*(dVlength+1)/2;
	int targ = int(trunc(numcells/nThreadz));
	targ = (targ<1) ? 1 : targ;
	int remainder = numcells % nThreadz;
	int jlim, cellnum=0;

	//Bin AIM elements:
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

	//Stuff for assessing slowest thread under row-binning scheme:
	double N = double(Vrows);
	/*The computational cost of computing a diagonal element includes the upfront cost of
	computing ytPdV_dtheta, and the cost of computing the gradient element.
	2*N^2 for ytPdV_dtheta
	1.5*N^2 + 0.5*N to efficiently calculate trace of PdV_dtheta
	2*N to finish gradient element
	(2*N^2) + 2*N for diagonal element:*/
	double diagcost = 5.5*R_pow_di(N,2) + 4.5*N;
	double offdiagcost = 4*R_pow_di(N,2) + 2*N;
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
	double inicost = 2*R_pow_di(N,2);
	/*^^^When the thread starts its work, and whenever it moves to a new row of
	the AIM, it computes ytPdV_dtheta.*/
	/*Thread computes a gradient element whenever it computes a diagonal element
	of the AIM:*/
	diagcost = 3.5*R_pow_di(N,2) + 4.5*N;
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

	parallelDerivScheme = (rowslowest<=cellslowest) ? 2 : 3;
	//The bin-by-row code assumes nThreadz <= numExplicitFreePar:
	if(nThreadz > numExplicitFreePar){parallelDerivScheme = 3;}
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
		if(!didUserGivedV[i]){
			continue;
		}
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


void omxGREMLFitState::dVupdate_final(){
	for(int i=0; i < numExplicitFreePar; i++){
		if(!didUserGivedV[i]){
			continue;
		}
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



struct GRMFIMLFitState : omxFitFunction{
	int verbose;
	omxMatrix *y, *invcov, *means;

	virtual void init();
	virtual void compute(int want, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
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

	oo->y = omxGetExpectationComponent(expectation, "y");
	oo->invcov = omxGetExpectationComponent(expectation, "invcov");
	oo->means = omxGetExpectationComponent(expectation, "means");

}

void GRMFIMLFitState::compute(int want, FitContext* fc){
	auto *oo = this;
	if(want & FF_COMPUTE_GRADIENT){invalidateGradient(fc);}
	const double NATLOG_2PI = 1.837877066409345483560659472811;	//<--log(2*pi)
	const double Scale = fabs(Global->llScale);
	omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation);
	Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(y), y->cols, 1);
	Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(invcov), invcov->rows, invcov->cols);
	Eigen::Map< Eigen::MatrixXd > Yhat(omxMatrixDataColumnMajor(means), means->rows, 1);
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;
	if(want & FF_COMPUTE_FIT){
		omxExpectationCompute(fc, expectation, NULL);
		//Check for PD covariance matrix:
		if(oge->cholV_fail_om->data[0]){
			oo->matrix->data[0] = NA_REAL;
			if (fc) fc->recordIterationError("expected covariance matrix is non-positive-definite");
			return;
		}
		/*This function doesn't use the quadratic form in X, but if it's non-PD, then the yhats won't have been recalculated
		  at the current parameter values:*/
		if(oge->cholquadX_fail){
			oo->matrix->data[0] = NA_REAL;
			if (fc) fc->recordIterationError("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
			return;
		}
		Eigen::MatrixXd resids(means->rows,1);
		resids = Eigy - Yhat;
		//Compute ML fit:
		oo->matrix->data[0] = Scale*0.5*( (((double)oo->y->cols)*NATLOG_2PI) + oge->logdetV_om->data[0] +
			(resids.transpose() * Vinv.selfadjointView<Eigen::Lower>() * resids)(0,0) );
	}
	return;
}

void GRMFIMLFitState::populateAttr(SEXP algebra){
	//Not really anything to pass to the frontend fitfunction object.
	return;
}
