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
  std::vector< int > dropcase;
  int numcases2drop;
  std::vector< omxMatrix* > dV;
  std::vector< const char* > dVnames;
  int dVlength;
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
  
  //Ensure that the pointer in the GREML fitfunction is directed at the right FreeVarGroup
  //(not necessary for most compute plans):
  if(fc && gff->varGroup != fc->varGroup){
    gff->buildParamMap(fc->varGroup);
	}
  
  //Declare local variables used in more than one scope in this function:
  const double Scale = fabs(Global->llScale); //<--absolute value of loglikelihood scale
  int i;
  EigenMatrixAdaptor Eigy(gff->y);
  
  //Trim out cases with missing data from V, if necessary:
  Eigen::MatrixXd EigV(gff->y->rows, gff->y->rows);
  if( gff->numcases2drop ){
    dropCasesAndEigenize(gff->V, EigV, gff->numcases2drop, gff->dropcase);
  }
  else{EigV = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->V), gff->V->rows, gff->V->cols);}
  
  if(want & (FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)){
    
    //Declare local variables for this scope:
    double logdetV=0, logdetquadX=0;
    Eigen::MatrixXd Vinv, quadX;
    EigenMatrixAdaptor EigX(gff->X);
    //EigenMatrixAdaptor EigV(gff->V);
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
    Vinv = rbstcholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse
    
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
    Eigen::MatrixXd PdV_dtheta1, PdV_dtheta2;//, dV_dtheta1, dV_dtheta2;
    Eigen::MatrixXd dV_dtheta1(gff->y->rows, gff->y->rows); //<--Derivative of V w/r/t parameter i.
    Eigen::MatrixXd dV_dtheta2(gff->y->rows, gff->y->rows); //<--Derivative of V w/r/t parameter j.
    
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
      if( gff->numcases2drop ){
        dropCasesAndEigenize(gff->dV[i], dV_dtheta1, gff->numcases2drop, gff->dropcase);
      }
      else{dV_dtheta1 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[i]), gff->dV[i]->rows, gff->dV[i]->cols);}
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
          if( gff->numcases2drop ){
            dropCasesAndEigenize(gff->dV[j], dV_dtheta2, gff->numcases2drop, gff->dropcase);
          }
          else{dV_dtheta2 = Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(gff->dV[j]), gff->dV[j]->rows, gff->dV[j]->cols);}
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
  
  if(OMX_DEBUG) { mxLog("Initializing GREML fitfunction."); }
  oo->units = FIT_UNITS_UNKNOWN;  // should be FIT_UNITS_MINUS2LL ? TODO
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
  
  if( newObj->X->rows != newObj->y->rows ){
    Rf_error("X and y matrices do not have equal numbers of rows");
  }
   if( newObj->V->rows != newObj->V->cols ){
    Rf_error("V matrix is not square");
  }
  
  //Deal with missing data:
  SEXP rObj = oo->rObj;
  SEXP casesToDrop, do_drop;
  int* casesToDrop_intptr;
  int i=0;
  newObj->numcases2drop = 0;
  {
  ScopedProtect p1(casesToDrop, R_do_slot(rObj, Rf_install("casesToDrop")));
  ScopedProtect p2(do_drop, R_do_slot(rObj, Rf_install("dropNAfromV")));
  if(Rf_length(casesToDrop) && Rf_asInteger(do_drop)){
    if(OMX_DEBUG) { mxLog("Preparing GREML fitfunction to handle missing data."); }
    newObj->numcases2drop = Rf_length(casesToDrop);
    casesToDrop_intptr = INTEGER(casesToDrop);
    newObj->dropcase.assign(newObj->V->rows,0);
    for(i=0; i < Rf_length(casesToDrop); i++){
      if(casesToDrop_intptr[i] >= newObj->V->rows){
        Rf_warning("casesToDrop vector in GREML fitfunction contains indices greater than the number of observations");
        newObj->numcases2drop--; 
      }
      //Need to subtract 1 from the index because R begins array indexing with 1, not 0:
      else{newObj->dropcase[casesToDrop_intptr[i]-1] = 1;}
  }}
  }
  
  if(newObj->y->rows != newObj->V->rows - newObj->numcases2drop){
    Rf_error("y and V matrices do not have equal numbers of rows");
  }
  
  //Tell the frontend fitfunction counterpart how many observations there are:
  SEXP nval;
  ScopedProtect p3(nval, R_do_slot(rObj, Rf_install("numObs")));
  int* numobs = INTEGER(nval);
  numobs[0] = newObj->y->rows;
  
  omxGREMLExpectation* oge = (omxGREMLExpectation*)(expectation->argStruct);
  if(OMX_DEBUG) { mxLog("Beginning last steps in initializing GREML fitfunction."); }
  newObj->dV = oge->dV;
  newObj->dVnames = oge->dVnames;
  newObj->dVlength = oge->dVlength;
  if(newObj->dVlength){
    oo->gradientAvailable = true;
    newObj->gradient.setZero(newObj->dVlength,1);
    oo->hessianAvailable = true;
    newObj->avgInfo.setZero(newObj->dVlength,newObj->dVlength);
    for(i=0; i < newObj->dVlength; i++){
      if( (newObj->dV[i]->rows != newObj->V->rows) || (newObj->dV[i]->cols != newObj->V->cols) ){
        Rf_error("all derivatives of V must have the same dimensions as V");
}}}}


void omxDestroyGREMLFitFunction(omxFitFunction *oo){
  if(OMX_DEBUG) {mxLog("Freeing GREML FitFunction.");}
    if(oo->argStruct == NULL) return;
    omxGREMLFitState* owo = ((omxGREMLFitState*)oo->argStruct);
    delete owo;
}


static void omxPopulateGREMLAttributes(omxFitFunction *oo, SEXP algebra){
  if(OMX_DEBUG) { mxLog("Populating GREML Attributes."); }

  omxGREMLFitState *gff = ((omxGREMLFitState*)oo->argStruct);
  
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
  if(OMX_DEBUG) { mxLog("Building parameter map for GREML fitfunction."); }
  varGroup = newVarGroup;
	gradMap.resize(dV.size());
	for (size_t nx=0; nx < dV.size(); ++nx) {
		int to = varGroup->lookupVar(dVnames[nx]);
		gradMap[nx] = to;
	}
}


static double omxAliasedMatrixElement(omxMatrix *om, int row, int col)
{
  int index = 0;
	if(row >= om->originalRows || col >= om->originalCols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.", 
			row + 1, col + 1, om->originalRows, om->originalCols);
		Rf_error(errstr);
		free(errstr);  // TODO not reached
        return (NA_REAL);
	}
	if(om->colMajor) {
		index = col * om->originalRows + row;
	} else {
		index = row * om->originalCols + col;
	}
	return om->data[index];
}


void dropCasesAndEigenize(omxMatrix* om, Eigen::MatrixXd &em, int num2drop, std::vector< int > todrop){
  
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
