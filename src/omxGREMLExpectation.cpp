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
 
#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxDefines.h"
#include "omxGREMLExpectation.h"
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Dense"
 
void omxInitGREMLExpectation(omxExpectation* ox){
  
  SEXP rObj = ox->rObj;
  SEXP Rmtx, casesToDrop; //dV, dVnames
  int i=0;
  omxState* currentState = ox->currentState;
  
  if(OMX_DEBUG) { mxLog("Initializing GREML expectation."); }
  
  //omxGREMLExpectation *oge = (omxGREMLExpectation*) R_alloc(1, sizeof(omxGREMLExpectation));
  omxGREMLExpectation *oge = new omxGREMLExpectation;
  
  /* Set Expectation Calls and Structures */
  ox->computeFun = omxComputeGREMLExpectation;
	ox->destructFun = omxDestroyGREMLExpectation;
	ox->componentFun = omxGetGREMLExpectationComponent;
	ox->populateAttrFun = omxPopulateGREMLAttributes;
	ox->argStruct = (void*) oge;
  
    /* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Processing V."); }
	oge->cov = omxNewMatrixFromSlot(rObj, currentState, "V");
  if( oge->cov->rows != oge->cov->cols ){
    Rf_error("V matrix is not square");
  }
	if(OMX_DEBUG) { mxLog("Processing X."); }
  {ScopedProtect p1(Rmtx, R_do_slot(rObj, Rf_install("X")));
	oge->X = omxNewMatrixFromRPrimitive(Rmtx, currentState, 0, 0);
  }
  {ScopedProtect p1(Rmtx, R_do_slot(rObj, Rf_install("y")));
	oge->y = omxNewMatrixFromRPrimitive(Rmtx, currentState, 0, 0);
  }
  oge->means = omxInitMatrix(oge->y->rows, 1, 1, currentState);
  for(i=0; i < oge->y->rows; i++){oge->means->data[i] = 0;}

  //Deal with missing data:
  int* casesToDrop_intptr;
  oge->numcases2drop = 0;
  {
  ScopedProtect p1(casesToDrop, R_do_slot(rObj, Rf_install("casesToDrop")));
  if(Rf_length(casesToDrop)){
    if(OMX_DEBUG) { mxLog("Preparing GREML expectation to handle missing data."); }
    oge->numcases2drop = Rf_length(casesToDrop);
    casesToDrop_intptr = INTEGER(casesToDrop);
    oge->dropcase.assign(oge->cov->rows,0);
    for(i=0; i < Rf_length(casesToDrop); i++){
      if(casesToDrop_intptr[i] >= oge->cov->rows){
        Rf_warning("casesToDrop vector in GREML expectation contains indices greater than the number of observations");
        oge->numcases2drop--; 
      }
      //Need to subtract 1 from the index because R begins array indexing with 1, not 0:
      else{oge->dropcase[casesToDrop_intptr[i]-1] = 1;}
  }}
  }
  if(oge->y->rows != oge->cov->rows - oge->numcases2drop){
    Rf_error("y and V matrices do not have equal numbers of rows");
  }
  
  //Initially compute everything involved in computing means:
  oge->alwaysComputeMeans = 1;
  oge->cholV_fail = 0; 
  oge->cholquadX_fail = 0;
  EigenMatrixAdaptor EigX(oge->X);
  EigenMatrixAdaptor Eigy(oge->y);
  Eigen::Map< Eigen::MatrixXd > yhat(omxMatrixDataColumnMajor(oge->means), oge->means->rows, oge->means->cols);
  Eigen::MatrixXd EigV(oge->y->rows, oge->y->rows);
  Eigen::MatrixXd quadX;
  Eigen::LLT< Eigen::MatrixXd > cholV(oge->y->rows);
  Eigen::LLT< Eigen::MatrixXd > cholquadX(oge->X->cols);
  if( oge->numcases2drop ){
    dropCasesAndEigenize(oge->cov, EigV, oge->numcases2drop, oge->dropcase);
  }
  else{EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(oge->cov), oge->cov->rows, oge->cov->cols);}
  cholV.compute(EigV.selfadjointView<Eigen::Lower>());
  if(cholV.info() != Eigen::Success){
    Rf_error("Expected covariance matrix is non-positive-definite at initial values");
  }
  oge->cholV_vectorD = (( Eigen::MatrixXd )(cholV.matrixL())).diagonal();
  oge->Vinv = cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse
  oge->XtVinv = EigX.transpose() * oge->Vinv;
  quadX = oge->XtVinv * EigX;
  cholquadX.compute(quadX.selfadjointView<Eigen::Lower>());
  if(cholquadX.info() != Eigen::Success){
    Rf_error("Cholesky factorization failed at initial values; possibly, the matrix of covariates is rank-deficient");
  }
  oge->cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal();
  oge->quadXinv = cholquadX.solve(Eigen::MatrixXd::Identity(oge->X->cols, oge->X->cols));
  oge->b = oge->quadXinv * oge->XtVinv * Eigy;
  yhat = EigX * oge->b;
  //Rf_error("Best to stop here");
  
  
/*  {ScopedProtect p1(dV, R_do_slot(rObj, Rf_install("dV")));
	ScopedProtect p2(dVnames, R_do_slot(rObj, Rf_install("dVnames")));
  oge->dVlength = Rf_length(dV);  
  oge->dV.resize(oge->dVlength);
  oge->dVnames.resize(oge->dVlength);
	if(oge->dVlength){
    if(OMX_DEBUG) { mxLog("Processing derivatives of V."); }
		int* dVint = INTEGER(dV);
    for(i=0; i < oge->dVlength; i++){
      oge->dV[i] = omxMatrixLookupFromState1(dVint[i], currentState);
      SEXP elem;
      {ScopedProtect p3(elem, STRING_ELT(dVnames, i));
			oge->dVnames[i] = CHAR(elem);}
	}}
  } */
}


void omxComputeGREMLExpectation(omxExpectation* ox, const char *, const char *) {
  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
	omxRecompute(oge->cov, NULL);
  oge->cholV_fail = 0;
  oge->cholquadX_fail = 0;
  
  EigenMatrixAdaptor EigX(oge->X);
  EigenMatrixAdaptor Eigy(oge->y);
  Eigen::Map< Eigen::MatrixXd > yhat(omxMatrixDataColumnMajor(oge->means), oge->means->rows, oge->means->cols);
  Eigen::MatrixXd EigV(oge->y->rows, oge->y->rows);
  Eigen::MatrixXd quadX;
  Eigen::LLT< Eigen::MatrixXd > cholV(oge->y->rows);
  Eigen::LLT< Eigen::MatrixXd > cholquadX(oge->X->cols);
  if( oge->numcases2drop ){
    dropCasesAndEigenize(oge->cov, EigV, oge->numcases2drop, oge->dropcase);
  }
  else{EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(oge->cov), oge->cov->rows, oge->cov->cols);}
  cholV.compute(EigV);
  if(cholV.info() != Eigen::Success){
    oge->cholV_fail = 1;
    return;
  }
  oge->cholV_vectorD = (( Eigen::MatrixXd )(cholV.matrixL())).diagonal();
  oge->Vinv = cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse
  oge->XtVinv = EigX.transpose() * oge->Vinv;
  quadX = oge->XtVinv * EigX;
  cholquadX.compute(quadX);
  if(cholquadX.info() != Eigen::Success){ 
    oge->cholquadX_fail = 1;
    return;
  }
  oge->cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal();
  oge->quadXinv = cholquadX.solve(Eigen::MatrixXd::Identity(oge->X->cols, oge->X->cols));
  if(oge->alwaysComputeMeans){
    oge->b = oge->quadXinv * oge->XtVinv * Eigy;
    yhat = EigX * oge->b;
  }
  
/*  if(oge->dVlength){
    for(int i=0; i < oge->dVlength; i++){
      omxRecompute(oge->dV[i], NULL);
  }} */
}


void omxDestroyGREMLExpectation(omxExpectation* ox) {
	if(OMX_DEBUG) { mxLog("Destroying GREML Expectation."); }
  omxGREMLExpectation* argStruct = (omxGREMLExpectation*)(ox->argStruct);
  omxFreeMatrix(argStruct->means);
}


void omxPopulateGREMLAttributes(omxExpectation *ox, SEXP algebra) {
  if(OMX_DEBUG) { mxLog("Populating GREML expectation attributes."); }

  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
  
  Rf_setAttrib(algebra, Rf_install("numStats"), Rf_ScalarReal(oge->y->rows));
  Rf_setAttrib(algebra, Rf_install("numFixEff"), Rf_ScalarInteger(oge->X->cols));
}


omxMatrix* omxGetGREMLExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component){
/* Return appropriate parts of Expectation to the Fit Function */
  if(OMX_DEBUG) { mxLog("GREML expectation: %s requested--", component); }

	omxGREMLExpectation* oge = (omxGREMLExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = oge->cov;
	} else if(strEQ("X", component)) {
		retval = oge->X;
	} else if(strEQ("y", component)) {
		retval = oge->y;
	}
  else if(strEQ("means", component)) {
  	retval = oge->means;
  }
	if (retval) omxRecompute(retval, NULL);
	
	return retval;
}


void dropCasesAndEigenize(omxMatrix* om, Eigen::MatrixXd &em, int num2drop, std::vector< int > todrop){
/*TODO: Assuming this function is only to be used with symmetric matrices, rewrite it to ignore nonunique
matrix elements*/
  
  if(OMX_DEBUG) { mxLog("Trimming out cases with missing data..."); }
  
  if(num2drop < 1){ return; }
  
  omxEnsureColumnMajor(om);
  
  if(om->algebra == NULL){ //i.e., if omxMatrix is from a frontend MxMatrix
  
    em.setZero(om->rows - num2drop, om->cols - num2drop);
  
    int nextCol = 0;
    int nextRow = 0;
    
    for(int j = 0; j < om->cols; j++) {
  	  if(todrop[j]) continue;
  		nextRow = 0;
  		for(int k = 0; k < om->rows; k++) {
  			if(todrop[k]) continue;
  			em(nextRow,nextCol) = omxAliasedMatrixElement(om, k, j);
  			nextRow++;
  		}
  		nextCol++;
  	}
  }
  else{ /*If the omxMatrix is from an algebra, then copying is not necessary; it can be resized directly
  and and Eigen-mapped, since the algebra will be recalculated back to its original dimensions anyhow.*/
    if(om->originalRows == 0 || om->originalCols == 0) Rf_error("Not allocated");
    if (om->rows != om->originalRows || om->cols != om->originalCols) {
      // Feasible, but the code is currently not robust to this case
      Rf_error("Can only omxRemoveRowsAndColumns once");
    }
    
    int oldRows = om->originalRows;
    int oldCols = om->originalCols;
    
    int nextCol = 0;
    int nextRow = 0;
    
    om->rows = oldRows - num2drop;
    om->cols = oldCols - num2drop;
    
    for(int j = 0; j < oldCols; j++) {
      if(todrop[j]) continue;
      nextRow = 0;
      for(int k = 0; k < oldRows; k++) {
        if(todrop[k]) continue;
        omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j));
        nextRow++;
      }
      nextCol++;
    }
    em = Eigen::Map< Eigen::MatrixXd >(om->data, om->rows, om->cols);
    omxMarkDirty(om);
  }
  if(OMX_DEBUG) { mxLog("Finished trimming out cases with missing data..."); }
}