 /*
 *  Copyright 2007-2016 The OpenMx Project
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
#include "omxGREMLfitfunction.h"
#include "omxGREMLExpectation.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Rmath.h>


 void omxInitGREMLFitFunction(omxFitFunction *oo){
  
  if(OMX_DEBUG) { mxLog("Initializing GREML fitfunction."); }
  SEXP rObj = oo->rObj;
  SEXP dV, dVnames, aug, augGrad, augHess;
  int i=0;
  
  oo->units = FIT_UNITS_MINUS2LL;
  oo->computeFun = omxCallGREMLFitFunction;
  oo->ciFun = loglikelihoodCIFun;
  oo->destructFun = omxDestroyGREMLFitFunction;
  oo->populateAttrFun = omxPopulateGREMLAttributes;
  
  omxGREMLFitState *newObj = new omxGREMLFitState;
  oo->argStruct = (void*)newObj;
  omxExpectation* expectation = oo->expectation;
  omxState* currentState = expectation->currentState;
  newObj->usingGREMLExpectation = (strcmp(expectation->expType, "MxExpectationGREML")==0 ? 1 : 0);
  if(!newObj->usingGREMLExpectation){
    //Maybe someday GREML fitfunction could be made compatible with another expectation, but not at present:
    Rf_error("GREML fitfunction is currently only compatible with GREML expectation");
  }
  else{
    omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
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
  
  //Augmentation:
  {
  ScopedProtect p1(aug, R_do_slot(rObj, Rf_install("aug")));
  if(Rf_length(aug)){
  	int* augint = INTEGER(aug);
  	newObj->aug = omxMatrixLookupFromStateByNumber(augint[0], currentState);
  }
  else{newObj->aug = NULL;}
  }
  
  //Derivatives of V:
  {
  ScopedProtect p1(dV, R_do_slot(rObj, Rf_install("dV")));
  ScopedProtect p2(dVnames, R_do_slot(rObj, Rf_install("dVnames")));
  newObj->dVlength = Rf_length(dV);  
  newObj->dV.resize(newObj->dVlength);
  newObj->indyAlg.resize(newObj->dVlength);
  newObj->dVnames.resize(newObj->dVlength);
  newObj->origdVdim.resize(newObj->dVlength);
	if(newObj->dVlength){
    if(!newObj->usingGREMLExpectation){
      //Probably best not to allow use of dV if we aren't sure means will be calculated GREML-GLS way:
      Rf_error("derivatives of 'V' matrix in GREML fitfunction only compatible with GREML expectation");
    }
    if(OMX_DEBUG) { mxLog("Processing derivatives of V."); }
		int* dVint = INTEGER(dV);
    for(i=0; i < newObj->dVlength; i++){
      newObj->dV[i] = omxMatrixLookupFromStateByNumber(dVint[i], currentState);
      SEXP elem;
      {ScopedProtect p3(elem, STRING_ELT(dVnames, i));
			newObj->dVnames[i] = CHAR(elem);}
	}}
  }
  
  if(newObj->dVlength){
    oo->gradientAvailable = true;
    newObj->gradient.setZero(newObj->dVlength,1);
    oo->hessianAvailable = true;
    newObj->avgInfo.setZero(newObj->dVlength,newObj->dVlength);
    newObj->rowbins.resize(Global->numThreads);
    newObj->AIMelembins.resize(Global->numThreads);
    for(i=0; i < newObj->dVlength; i++){
    	/*Each dV must either (1) match the dimensions of V, OR (2) match the length of y if that is less than the 
    	dimension of V (implying downsizing due to missing observations):*/
      if( ((newObj->dV[i]->rows == newObj->cov->rows)&&(newObj->dV[i]->cols == newObj->cov->cols)) ||
          ((newObj->y->cols < newObj->cov->rows)&&(newObj->dV[i]->rows == newObj->y->cols)&&
          	(newObj->dV[i]->cols == newObj->y->cols)) ){
      	newObj->origdVdim[i] = newObj->dV[i]->rows;
      }
      else{
        Rf_error("all derivatives of V must have the same dimensions as V");
  }}}
  
  //Augmentation derivatives:
	if(newObj->dVlength && newObj->aug){
	//^^^Ignore derivatives of aug unless aug itself and objective derivatives are supplied.	
		ScopedProtect p1(augGrad, R_do_slot(rObj, Rf_install("augGrad")));
		ScopedProtect p2(augHess, R_do_slot(rObj, Rf_install("augHess")));
		if(!Rf_length(augGrad)){
			if(Rf_length(augHess)){
				Rf_error("if argument 'augHess' has nonzero length, then argument 'augGrad' must as well");
			}
			else{
				Rf_error("if arguments 'dV' and 'aug' have nonzero length, then 'augGrad' must as well");
		}}
		else{
			int* augGradint = INTEGER(augGrad);
			newObj->augGrad = omxMatrixLookupFromStateByNumber(augGradint[0], currentState);
			if(Rf_length(augHess)){
				int* augHessint = INTEGER(augHess);
				newObj->augHess = omxMatrixLookupFromStateByNumber(augHessint[0], currentState);
			}
			else{oo->hessianAvailable = false;}
		}
	}
}




void omxCallGREMLFitFunction(omxFitFunction *oo, int want, FitContext *fc){
  if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
  
  //Recompute Expectation:
  omxExpectation* expectation = oo->expectation;
  omxExpectationCompute(fc, expectation, NULL);
    
  omxGREMLFitState *gff = (omxGREMLFitState*)oo->argStruct; //<--Cast generic omxFitFunction to omxGREMLFitState
  
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
  P.setZero(gff->invcov->rows, gff->invcov->cols);
  double logdetV=0, logdetquadX=0, ytPy=0;
  
  if(want & (FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
    if(gff->usingGREMLExpectation){
      omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
      
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
      
      //Finish computing fit (negative loglikelihood):
      P.triangularView<Eigen::Lower>() = (Vinv.selfadjointView<Eigen::Lower>() * //P = Vinv * (I-Hatmat)
        (Eigen::MatrixXd::Identity(Vinv.rows(), Vinv.cols()) - 
          (EigX * oge->quadXinv.selfadjointView<Eigen::Lower>() * oge->XtVinv))).triangularView<Eigen::Lower>();
      Py = P.selfadjointView<Eigen::Lower>() * Eigy;
      ytPy = (Eigy.transpose() * Py)(0,0);
      if(OMX_DEBUG) {mxLog("ytPy is %3.3f",ytPy);}
      oo->matrix->data[0] = gff->REMLcorrection + 
      	Scale*0.5*( (((double)gff->y->cols) * NATLOG_2PI) + logdetV + ytPy) + Scale*gff->pullAugVal(0L,0,0);
      gff->nll = oo->matrix->data[0];
      if(OMX_DEBUG){mxLog("augmentation is %3.3f",gff->pullAugVal(0L,0,0));}
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
    omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
  	
  	//Recompute derivatives:
  	gff->dVupdate(fc);
  	gff->recomputeAug(1, fc);
    
    //Declare local variables for this scope:
    int nThreadz = Global->numThreads;
    int wantHess = 0;
    
    fc->grad.resize(gff->dVlength); //<--Resize gradient in FitContext
    
    //Set up new HessianBlock:
    HessianBlock *hb = new HessianBlock;
    if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
      hb->vars.resize(gff->dVlength);
      hb->mat.resize(gff->dVlength, gff->dVlength);
      gff->recomputeAug(2, fc);
      wantHess = 1;
    }
    
    if(gff->parallelDerivScheme==0){gff->planParallelDerivs(nThreadz,wantHess,gff->cov->rows);}
    
    //Begin looping thru free parameters:
#pragma omp parallel num_threads(nThreadz)
{
		int i=0, j=0, t1=0, t2=0, a1=0, a2=0, r=0, c=0;
		double tr=0;
		Eigen::MatrixXd ytPdV_dtheta1;
		//Eigen::VectorXd diagPdV_dtheta1;
		Eigen::MatrixXd dV_dtheta1(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter i.
		Eigen::MatrixXd dV_dtheta2(Eigy.rows(), Eigy.rows()); //<--Derivative of V w/r/t parameter j.
		int threadID = omx_absolute_thread_num();
		int istart = threadID * gff->dVlength / nThreadz;
		int iend = (threadID+1) * gff->dVlength / nThreadz;
		if(threadID == nThreadz-1){iend = gff->dVlength;}
		for(i=istart; i < iend; i++){
			tr=0;
			t1 = gff->gradMap[i]; //<--Parameter number for parameter i.
			if(t1 < 0){continue;}
			a1 = gff->dAugMap[i]; //<--Index of augmentation derivatives to use for parameter i.
			if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){hb->vars[i] = t1;}
			if( oge->numcases2drop && (gff->dV[i]->rows > Eigy.rows()) ){
				dropCasesAndEigenize(gff->dV[i], dV_dtheta1, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[i]);
			}
			else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[i]), gff->dV[i]->rows, gff->dV[i]->cols);}
			ytPdV_dtheta1 = Py.transpose() * dV_dtheta1.selfadjointView<Eigen::Lower>();
			for(j=i; j < gff->dVlength; j++){
				if(j==i){
					/*Need trace of P*dV_dtheta for gradient element...
					Frustratingly, the selfadjointView has no row or column accessor function among its members.
					But the trace of a product of two square symmetric matrices is the sum of the elements of
					their elementwise product.*/
					//diagPdV_dtheta1(k) = (P.selfadjointView<Eigen::Lower>()).row(k) * (dV_dtheta1.selfadjointView<Eigen::Lower>()).col(k);
					for(c=0; c < gff->cov->rows; c++){
						for(r=c; r < gff->cov->rows; r++){
							tr += (r==c) ? P(r,c)*dV_dtheta1(r,c) : 2*P(r,c)*dV_dtheta1(r,c);
						}
					}
					gff->gradient(t1) = Scale*0.5*(tr - (ytPdV_dtheta1 * Py)(0,0)) + 
						Scale*gff->pullAugVal(1,a1,0);
					fc->grad(t1) += gff->gradient(t1);
					if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
						gff->avgInfo(t1,t1) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() * ytPdV_dtheta1.transpose())(0,0) + 
							Scale*gff->pullAugVal(2,a1,a1);
					}
				}
				else{if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
					t2 = gff->gradMap[j]; //<--Parameter number for parameter j.
					if(t2 < 0){continue;}
					a2 = gff->dAugMap[j]; //<--Index of augmentation derivatives to use for parameter j.
					if( oge->numcases2drop && (gff->dV[j]->rows > Eigy.rows()) ){
						dropCasesAndEigenize(gff->dV[j], dV_dtheta2, oge->numcases2drop, oge->dropcase, 1, gff->origdVdim[j]);
					}
					else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[j]), gff->dV[j]->rows, gff->dV[j]->cols);}
					gff->avgInfo(t1,t2) = Scale*0.5*(ytPdV_dtheta1 * P.selfadjointView<Eigen::Lower>() * 
						dV_dtheta2.selfadjointView<Eigen::Lower>() * Py)(0,0) + Scale*gff->pullAugVal(2,a1,a2);
					gff->avgInfo(t2,t1) = gff->avgInfo(t1,t2);
				}}}}
}
    //Assign upper triangle elements of avgInfo to the HessianBlock:
    if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
      for (size_t d1=0, h1=0; h1 < gff->dV.size(); ++h1) {
		    for (size_t d2=0, h2=0; h2 <= h1; ++h2) {
				  	hb->mat(d2,d1) = gff->avgInfo(h2,h1);
				    ++d2;
        }
			  ++d1;	
		  }
		  fc->queue(hb);
  }}
  return;
}



void omxDestroyGREMLFitFunction(omxFitFunction *oo){
  if(OMX_DEBUG) {mxLog("Freeing GREML FitFunction.");}
    if(oo->argStruct == NULL) return;
    omxGREMLFitState* owo = ((omxGREMLFitState*)oo->argStruct);
    delete owo;
}


static void omxPopulateGREMLAttributes(omxFitFunction *oo, SEXP algebra){
	omxGREMLFitState *gff = (omxGREMLFitState*)oo->argStruct;
	gff->dVupdate_final();
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }
  SEXP rObj = oo->rObj;
  SEXP nval, mlfitval;
  int userSuppliedDataNumObs = (int)(( (omxGREMLExpectation*)(oo->expectation->argStruct) )->data2->numObs);
  
  //Tell the frontend fitfunction counterpart how many observations there are...:
  {
  //ScopedProtect p1(nval, R_do_slot(rObj, Rf_install("numObs")));
  ScopedProtect p1(nval, Rf_allocVector(INTSXP, 1));
  INTEGER(nval)[0] = 1L - userSuppliedDataNumObs;
  R_do_slot_assign(rObj, Rf_install("numObs"), nval);
  /*^^^^The end result is that number of observations will be reported as 1 in summary()...
  which is always correct with GREML.  This is a bit of a hack, since it is sneaking this
  negative numObs into the pre-backend fitfunction that summary() looks at...*/
	}
	
	{
	//ScopedProtect p1(mlfitval, R_do_slot(rObj, Rf_install("MLfit")));
	ScopedProtect p1(mlfitval, Rf_allocVector(REALSXP, 1));
	REAL(mlfitval)[0] = gff->nll - gff->REMLcorrection;
	Rf_setAttrib(algebra, Rf_install("MLfit"), mlfitval);
	}
}


void omxGREMLFitState::buildParamMap(FreeVarGroup *newVarGroup)
{
	if(OMX_DEBUG) { mxLog("Building parameter map for GREML fitfunction."); }
	varGroup = newVarGroup;
	if(dVlength){
		/*The pointers to the derivatives of V, their names, and their original dimensions get temporariliy 
		copied here:*/
		std::vector< omxMatrix* > dV_temp = dV;
		std::vector< const char* > dVnames_temp = dVnames;
		std::vector<int> origdVdim_temp = origdVdim;
		gradMap.resize(dVlength);
		dAugMap.resize(dVlength);
		int gx=0;
		/*If there are no problems, then every time vx gets incremented, it should become equal to the current 
		value of gx*/
		for (int vx=0; vx < int(varGroup->vars.size()); ++vx) {
			for (int nx=0; nx < dVlength; ++nx) {
				if (strEQ(dVnames_temp[nx], varGroup->vars[vx]->name)) {
					gradMap[gx] = vx;
					dV[gx] = dV_temp[nx];
					dVnames[gx] = dVnames_temp[nx]; //<--Probably not strictly necessary...
					origdVdim[gx] = origdVdim_temp[nx];
					dAugMap[gx] = nx;
					indyAlg[gx] = ( dV_temp[nx]->algebra && !(dV_temp[nx]->dependsOnParameters()) ) ? 1 : 0;
					++gx;
					break;
				}
			}
		}/*By the end of the loop, the member objects of the omxGREMLFitState (dV, dVnames, etc.) should have their
		elements arranged to match the order in which the free parameters appear in the freeVarGroup*/
		if (gx != dVlength) Rf_error("Problem in dVnames mapping"); //possibly, argument 'dV' has elements not named with free parameter labels
		if( gx < int(varGroup->vars.size()) ){Rf_error("At least one free parameter has no corresponding element in 'dV'");}
		
		if(augGrad){
			int ngradelem = std::max(augGrad->rows, augGrad->cols);
			if(ngradelem != dVlength){
				Rf_error("matrix referenced by 'augGrad' must have same number of elements as argument 'dV'");
			}
			if(augHess){
				if (augHess->rows != augHess->cols) {
					Rf_error("matrix referenced by 'augHess' must be square (instead of %dx%d)",
              augHess->rows, augHess->cols);
				}
				if(augHess->rows != ngradelem){
					Rf_error("Augmentation derivatives non-conformable (gradient is size %d and Hessian is %dx%d)",
              ngradelem, augHess->rows, augHess->cols);
			}}
		}
	}
}


void omxGREMLFitState::planParallelDerivs(int nThreadz, int wantHess, int Vrows){
	//Note: AIM = Average Information Matrix (Hessian)
	if(wantHess==0 || nThreadz<2){
		parallelDerivScheme = 1; //Divvy up parameters the old, naive way.
		return;
	}

	/*Under the AIM-row-binning scheme, each thread will calculate the gradient element and the row
	of the AIM (starting with its diagonal element) for each of its parameters.*/
	
	//Stuff for row binning:
	int rownums[dVlength], i, j, minbin;
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
	2*N^2 + N to efficiently calculate trace of PdV_dtheta
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
			if(c>=dVlength){
				r++;
				c = r;
			}
		}
	}
	double cellslowest = workbins.maxCoeff();
	
	parallelDerivScheme = (rowslowest<cellslowest) ? 2 : 3;
	return;
}
 

double omxGREMLFitState::pullAugVal(int thing, int row, int col){
	double val=0;
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
	for(int i=0; i < dVlength; i++){
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
	for(int i=0; i < dVlength; i++){
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


