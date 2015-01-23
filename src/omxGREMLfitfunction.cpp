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
  omxMatrix* y;
  omxMatrix* X;
  omxMatrix* V;
  int* do_fixeff;
  Eigen::MatrixXd XtVinv;
  Eigen::MatrixXd quadXinv;
  Eigen::MatrixXd P;
  Eigen::MatrixXd ytP;
}; 

void omxCallGREMLFitFunction(omxFitFunction *oo, int want, FitContext *fc){
  if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
  
  //Recompute Expectation:
  omxExpectation* expectation = oo->expectation;
  omxExpectationCompute(expectation, NULL);
  
  omxGREMLFitState *gff = (omxGREMLFitState*)oo->argStruct; //<--Cast generic omxFitFunction to omxGREMLFitState
  
  //if(want & (FF_COMPUTE_FIT)){
    //Declare local variables:
    int i;
    double logdetV=0, logdetquadX=0, nll=0;
    Eigen::MatrixXd Vinv, quadX;
    EigenMatrixAdaptor Eigy(gff->y);
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
    quadX = gff->XtVinv * EigX;
    
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
    
        //Finish computing fit:
    gff->P = Vinv - (gff->XtVinv.transpose() * gff->quadXinv * gff->XtVinv);
    gff->ytP = Eigy.transpose() * gff->P;
    nll = 0.5*(logdetV + logdetquadX + (gff->ytP * Eigy)(0,0));
    oo->matrix->data[0] = nll;
    return;
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
  omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
  newObj->do_fixeff = oge->do_fixeff;
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
