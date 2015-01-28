 /*
 *  Copyright 2007-2015 The OpenMx Project
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
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Dense"

struct omxGREMLFitState { 
  //TODO: Some of these members might be redundant with what's stored in the FitContext, 
  //and could therefore be cut
  omxMatrix* y;
  omxMatrix* X;
  omxMatrix* V;
  std::vector< omxMatrix* > dV;
  std::vector< const char* > dVnames;
  int dVlength;
  int* do_fixeff;
  double nll;
  Eigen::MatrixXd XtVinv;
  Eigen::MatrixXd quadXinv;
  Eigen::MatrixXd P;
  Eigen::MatrixXd Py;
  Eigen::VectorXd gradient;
  Eigen::MatrixXd avgInfo; //the Average Information matrix
  FreeVarGroup *varGroup;
	std::vector<int> gradMap;
  void buildParamMap(FreeVarGroup *newVarGroup);
}; 

void omxCallGREMLFitFunction(omxFitFunction *oo, int want, FitContext *fc){
  if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
  
  //Recompute Expectation:
  omxExpectation* expectation = oo->expectation;
  omxExpectationCompute(expectation, NULL);
    
  omxGREMLFitState *gff = (omxGREMLFitState*)oo->argStruct; //<--Cast generic omxFitFunction to omxGREMLFitState
  
  //Ensure that the pointer in the GREML fitfunction is directed at the right FreeVarGroup (is this necessary?):
  if(fc && gff->varGroup != fc->varGroup){
    gff->buildParamMap(fc->varGroup);
	}
  
  //Declare local variables used in more than one scope in this function:
  const double Scale = fabs(Global->llScale); //<--absolute value of loglikelihood scale
  int i;
  EigenMatrixAdaptor Eigy(gff->y);
  
  if(want & (FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
    
    //Declare local variables for this scope:
    double logdetV=0, logdetquadX=0;
    Eigen::MatrixXd Vinv, quadX;
    EigenMatrixAdaptor EigX(gff->X);
    EigenMatrixAdaptor EigV(gff->V);
    Eigen::LDLT< Eigen::MatrixXd > rbstcholV(gff->y->rows);
    Eigen::LDLT< Eigen::MatrixXd > rbstcholquadX(gff->X->cols);
    
    rbstcholV.compute(EigV); //<--Robust Cholesky factorization of V
    //Check that factorization succeeded and that V is positive definite:
    if(rbstcholV.info() != Eigen::Success){
      omxRaiseErrorf("Cholesky factorization failed due to unknown numerical error (is the expected covariance matrix asymmetric?)");
      oo->matrix->data[0] = NA_REAL;
      return;
    }
    if(!rbstcholV.isPositive()){
      oo->matrix->data[0] = NA_REAL;
      if (fc) fc->recordIterationError("Expected covariance matrix is non-positive-definite");
      return;
    }
    //Log determinant of V:
    for(i=0; i < gff->y->rows; i++){
      logdetV += log(rbstcholV.vectorD()[i]);
    }
    Vinv = rbstcholV.solve(Eigen::MatrixXd::Identity(gff->V->rows, gff->V->cols)); //<-- V inverse
    
    gff->XtVinv = EigX.transpose() * Vinv;
    quadX = gff->XtVinv * EigX; //<--Quadratic form in X
    
    //Do for XtVinvX as was done for V:
    rbstcholquadX.compute(quadX);
    if(rbstcholquadX.info() != Eigen::Success){
      omxRaiseErrorf("Cholesky factorization failed due to unknown numerical error");
      oo->matrix->data[0] = NA_REAL;
      return;
    }
    if(!rbstcholquadX.isPositive()){
      oo->matrix->data[0] = NA_REAL;
      if (fc) fc->recordIterationError("Cholesky factorization failed; possibly, the matrix of covariates is rank-deficient");
      return;
    }
    for(i=0; i < gff->X->cols; i++){
      logdetquadX += log(rbstcholquadX.vectorD()[i]);
    }
    gff->quadXinv = rbstcholquadX.solve(Eigen::MatrixXd::Identity(gff->X->cols, gff->X->cols));
    
    //Finish computing fit (negative loglikelihood):
    gff->P = Vinv - (gff->XtVinv.transpose() * gff->quadXinv * gff->XtVinv);
    gff->Py = gff->P * Eigy;
    oo->matrix->data[0] = Scale*0.5*(logdetV + logdetquadX + (Eigy.transpose() * gff->Py)(0,0));
    gff->nll = oo->matrix->data[0]; 
  }
  
  if(want & (FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
    //Declare local variables for this scope:
    int j=0, t1=0, t2=0;
    Eigen::MatrixXd PdV_dtheta1, PdV_dtheta2;
    
    fc->grad.resize(gff->dVlength); //<--Resize gradient in FitContext
    
    //Set up new HessianBlock:
    HessianBlock *hb = new HessianBlock;
    if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
      hb->vars.resize(gff->dVlength);
      hb->mat.resize(gff->dVlength, gff->dVlength);
    }
    
    //Begin looping thru free parameters:
    for(i=0; i < gff->dVlength ; i++){
      t1 = gff->gradMap[i]; //<--Parameter number for parameter i.
      if(t1 < 0){continue;}
      if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){hb->vars[i] = t1;}
      EigenMatrixAdaptor dV_dtheta1(gff->dV[i]); //<--Derivative of V w/r/t parameter i.
      PdV_dtheta1 = gff->P * dV_dtheta1;
      for(j=i; j < gff->dVlength; j++){
        if(j==i){
          gff->gradient(t1) = Scale*0.5*(PdV_dtheta1.trace() - (Eigy.transpose() * PdV_dtheta1 * gff->Py)(0,0));
          fc->grad(t1) += gff->gradient(t1);
          if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
            gff->avgInfo(t1,t1) = Scale*0.5*(Eigy.transpose() * PdV_dtheta1 * PdV_dtheta1 * gff->Py)(0,0);
          }
        }
        else{if(want & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
          t2 = gff->gradMap[j]; //<--Parameter number for parameter j.
          if(t2 < 0){continue;}
          EigenMatrixAdaptor dV_dtheta2(gff->dV[j]); //<--Derivative of V w/r/t parameter j.
          gff->avgInfo(t1,t2) = Scale*0.5*(Eigy.transpose() * PdV_dtheta1 * gff->P * dV_dtheta2 * gff->Py)(0,0);
          gff->avgInfo(t2,t1) = gff->avgInfo(t1,t2);
    }}}}
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



void omxInitGREMLFitFunction(omxFitFunction *oo){
  
  if(OMX_DEBUG) { mxLog("Initializing GREML fit function."); }
  oo->computeFun = omxCallGREMLFitFunction;
  oo->destructFun = omxDestroyGREMLFitFunction;
  oo->populateAttrFun = omxPopulateGREMLAttributes;
  //omxGREMLFitState *newObj = (omxGREMLFitState*) R_alloc(1, sizeof(omxGREMLFitState));
  omxGREMLFitState *newObj = new omxGREMLFitState;
  oo->argStruct = (void*)newObj;
  omxExpectation* expectation = oo->expectation;
  newObj->y = omxGetExpectationComponent(expectation, oo, "y");
  newObj->V = omxGetExpectationComponent(expectation, oo, "V");
  newObj->X = omxGetExpectationComponent(expectation, oo, "X");
  newObj->nll = 0;
  newObj->varGroup = NULL;
  
  omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
  newObj->do_fixeff = oge->do_fixeff;
  newObj->dV = oge->dV;
  newObj->dVnames = oge->dVnames;
  newObj->dVlength = oge->dVlength;
  if(newObj->dVlength){
    oo->gradientAvailable = true;
    newObj->gradient.setZero(newObj->dVlength,1);
    oo->hessianAvailable = true;
    newObj->avgInfo.setZero(newObj->dVlength,newObj->dVlength);
  }
  //omxRaiseErrorf("Best to stop here for now");
}


void omxDestroyGREMLFitFunction(omxFitFunction *oo){
  if(OMX_DEBUG) {mxLog("Freeing GREML FitFunction.");}
    if(oo->argStruct == NULL) return;
    omxGREMLFitState* owo = ((omxGREMLFitState*)oo->argStruct);
    delete owo;
}


static void omxPopulateGREMLAttributes(omxFitFunction *oo, SEXP algebra){
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }
  omxGREMLFitState *gff = ((omxGREMLFitState*)oo->argStruct);
  if(gff->do_fixeff){
    EigenMatrixAdaptor Eigy(gff->y);
	Eigen::MatrixXd GREML_b = gff->quadXinv * gff->XtVinv * Eigy;
    SEXP b_ext, bcov_ext;
    Rf_protect(b_ext = Rf_allocMatrix(REALSXP, GREML_b.rows(), 1));
    for(int row = 0; row < GREML_b.rows(); row++){
  			REAL(b_ext)[0 * GREML_b.rows() + row] = GREML_b(row,0);
    }
  
  	Rf_protect(bcov_ext = Rf_allocMatrix(REALSXP, gff->quadXinv.rows(), gff->quadXinv.cols()));
    for(int row = 0; row < gff->quadXinv.rows(); row++){
		  for(int col = 0; col < gff->quadXinv.cols(); col++){
			  REAL(bcov_ext)[col * gff->quadXinv.rows() + row] = gff->quadXinv(row,col);
		}}
  
  	Rf_setAttrib(algebra, Rf_install("b"), b_ext);
  	Rf_setAttrib(algebra, Rf_install("bcov"), bcov_ext);
  }
}

//Alternate way to do fixed effects using QR solve:
/*  }
  if(want & (FF_COMPUTE_FIXEDEFFECTS)){
    Eigen::MatrixXd S, Sinv, SinvX, Sinvy, quadX;
    EigenMatrixAdaptor Eigy = EigenMatrixAdaptor(gff->y);
    EigenMatrixAdaptor EigX = EigenMatrixAdaptor(gff->X);
    EigenMatrixAdaptor EigV = EigenMatrixAdaptor(gff->V);
    Eigen::LLT< Eigen::MatrixXd > cholV(gff->y->rows);
    Eigen::LLT< Eigen::MatrixXd > cholquadX(gff->X->cols);
    
    cholV.compute(EigV);
    if(cholV.info() != Eigen::Success){
      omxRaiseErrorf("Cholesky factorization failed due to unknown numerical error (is the expected covariance matrix asymmetric?)");
      return;
    }
    S = cholV.matrixL();
    Sinv = S.inverse();
    SinvX = Sinv * EigX;
    Sinvy = Sinv * Eigy;
    fc->GREML_b = SinvX.colPivHouseholderQr().solve(Sinvy);
    quadX = EigX.transpose() * Sinv * Sinv.transpose() * EigX;
    cholquadX.compute(quadX);
    fc->GREML_bcov = cholquadX.solve(Eigen::MatrixXd::Identity(gff->X->cols, gff->X->cols));
    return;    
  } */


void omxGREMLFitState::buildParamMap(FreeVarGroup *newVarGroup)
{
  varGroup = newVarGroup;
	gradMap.resize(dV.size());
	for (size_t nx=0; nx < dV.size(); ++nx) {
		int to = varGroup->lookupVar(dVnames[nx]);
		gradMap[nx] = to;
	}
}
