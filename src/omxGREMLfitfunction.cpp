 /*
 *  Copyright 2007-2014 The OpenMx Project
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
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Dense"

struct omxGREMLFitFunction {
  omxMatrix* y;
  omxMatrix* X;
  omxMatrix* V;
  /*
  Eigen::MatrixXd Eigy;
  Eigen::MatrixXd EigV;
  Eigen::MatrixXd EigX;
  Eigen::MatrixXd P;
  Eigen::MatrixXd ytP;  */
}; 

void omxCallGREMLFitFunction(omxFitFunction *oo, int want, FitContext *fc){
  if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
  
  omxExpectation* expectation = oo->expectation;
  omxExpectationCompute(expectation, NULL);
  
  omxGREMLFitFunction *gff = (omxGREMLFitFunction*)oo->argStruct;
  int i;
  double logdetV=0, logdetquadX=0, nll=0;
  Eigen::MatrixXd Vinv, XtVinv, quadX, quadXinv, P, ytP;
  EigenMatrixAdaptor Eigy = EigenMatrixAdaptor(gff->y);
  EigenMatrixAdaptor EigX = EigenMatrixAdaptor(gff->X);
  EigenMatrixAdaptor EigV = EigenMatrixAdaptor(gff->V);
  Eigen::LDLT< Eigen::MatrixXd > rbstcholV(gff->y->rows);
  Eigen::LDLT< Eigen::MatrixXd > rbstcholquadX(gff->X->cols);
  
  rbstcholV.compute(EigV);
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
  for(i=0; i < gff->y->rows; i++){
    logdetV += log(rbstcholV.vectorD()[i]);
  }
  Vinv = rbstcholV.solve(Eigen::MatrixXd::Identity(gff->V->rows, gff->V->cols));
  
  XtVinv = EigX.transpose() * Vinv;
  quadX = XtVinv * EigX;
  
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
  quadXinv = rbstcholquadX.solve(Eigen::MatrixXd::Identity(gff->X->cols, gff->X->cols));
  
  P = Vinv - (XtVinv.transpose() * quadXinv * XtVinv);
  ytP = Eigy.transpose() * P;
  nll = 0.5*(logdetV + logdetquadX + (ytP * Eigy)(0,0));
  oo->matrix->data[0] = nll;
}


void omxInitGREMLFitFunction(omxFitFunction *oo){
  
  if(OMX_DEBUG) { mxLog("Initializing GREML fit function."); }
  oo->computeFun = omxCallGREMLFitFunction;
  oo->destructFun = omxDestroyGREMLFitFunction;
  oo->populateAttrFun = omxPopulateGREMLAttributes;
  omxGREMLFitFunction *newObj = (omxGREMLFitFunction*) R_alloc(1, sizeof(omxGREMLFitFunction));
  oo->argStruct = (void*)newObj;
  omxExpectation* expectation = oo->expectation;
  newObj->y = omxGetExpectationComponent(expectation, oo, "y");
  newObj->V = omxGetExpectationComponent(expectation, oo, "V");
  newObj->X = omxGetExpectationComponent(expectation, oo, "X");
  //newObj->Eigy = Eigen::Map< Eigen::MatrixXd >(newObj->y->data, newObj->y->rows, 1);
  //newObj->EigX = Eigen::Map< Eigen::MatrixXd >(newObj->X->data, newObj->X->rows, newObj->X->cols);
}


void omxDestroyGREMLFitFunction(omxFitFunction *oo){
  if(OMX_DEBUG) {mxLog("Freeing GREML FitFunction.");}
    if(oo->argStruct == NULL) return;
    omxGREMLFitFunction* owo = ((omxGREMLFitFunction*)oo->argStruct);
    delete owo;
}


static void omxPopulateGREMLAttributes(omxFitFunction *oo, SEXP algebra){
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }
}
