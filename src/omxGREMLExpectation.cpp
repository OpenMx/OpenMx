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

#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxDefines.h"
#include "omxGREMLExpectation.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include "EnableWarnings.h"

omxExpectation *omxInitGREMLExpectation(omxState *st, int num)
{ return new omxGREMLExpectation(st, num); }

void omxGREMLExpectation::connectToData()
{
  setConnectedToData(true);

	auto dc = getDataColumns();
	for (int cx=0; cx < int(dc.size()); ++cx) {
		int var = dc[cx];
		data->assertColumnIsData(var, OMXDATA_REAL);
	}
}

void omxGREMLExpectation::init()
{
	loadDataColFromR();

  SEXP Rmtx, RdoREML, Ryhat, casesToDrop, RyXcolnames;
  int i=0;

  if(OMX_DEBUG) { mxLog("Initializing GREML expectation."); }

  omxGREMLExpectation *oge = this;

    /* Set up expectation structures */
  //y:
  if(OMX_DEBUG) { mxLog("Processing y."); }
  oge->y = new omxData();
  {ScopedProtect p1(Rmtx, R_do_slot(rObj, Rf_install("y")));
	  oge->y->newDataStatic(currentState, Rmtx);
  }
  //V:
	if(OMX_DEBUG) { mxLog("Processing V."); }
	oge->cov = omxNewMatrixFromSlot(rObj, currentState, "V");
  if( oge->cov->rows != oge->cov->cols ){
    mxThrow("'V' matrix is not square");
  }
  //X:
	if(OMX_DEBUG) { mxLog("Processing X."); }
  {ScopedProtect p1(Rmtx, R_do_slot(rObj, Rf_install("X")));
	oge->X = omxNewMatrixFromRPrimitive(Rmtx, currentState, 0, 0);
  }
	oge->didUserProvideYhat = false;
	//REML:
  {
  	ScopedProtect p1(RdoREML, R_do_slot(rObj, Rf_install("REML")));
  	oge->doREML = Rf_asLogical(RdoREML);
  }
  //Eigy (local) will have however many rows and 1 column:
  Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(oge->y->dataMat), oge->y->dataMat->cols, 1);
  if(oge->X->rows != Eigy.rows()){mxThrow("'X' and 'y' matrices have different numbers of rows");}
  //yhat, provided by user:
  if(!oge->doREML){
  	ScopedProtect p1(Ryhat, R_do_slot(rObj, Rf_install("yhat")));
  	if(Rf_length(Ryhat)){
  		oge->didUserProvideYhat = true;
  		oge->yhatFromUser = omxNewMatrixFromSlot(rObj, currentState, "yhat");
  		oge->residual.setZero(Eigy.rows(),1);
  	}
  }
  //means:
  oge->means = omxInitMatrix(Eigy.rows(), 1, 1, currentState);
  //logdetV_om:
  oge->logdetV_om = omxInitMatrix(1, 1, 1, currentState);
  oge->logdetV_om->data[0] = 0;
  //cholV_fail:
  oge->cholV_fail = false;
  //quadXinv:
  oge->quadXinv.setZero(oge->X->cols, oge->X->cols);
  //original dimensions of V:
  oge->origVdim = oge->cov->rows;
  
  /*if(oge->didUserProvideYhat){
  	mxThrow("the first element of the yhat vector is %f",oge->yhatFromUser->data[0]);
  }*/
		


  //Deal with missing data:
  int* casesToDrop_intptr;
  oge->numcases2drop = 0;
  {
  ScopedProtect p1(casesToDrop, R_do_slot(rObj, Rf_install("casesToDrop")));
  if(Rf_length(casesToDrop)){
    if(OMX_DEBUG) { mxLog("Preparing GREML expectation to handle missing data."); }
    oge->numcases2drop = Rf_length(casesToDrop);
    casesToDrop_intptr = INTEGER(casesToDrop);
    oge->dropcase.assign(oge->cov->rows,false);
    for(i=0; i < Rf_length(casesToDrop); i++){
      if(casesToDrop_intptr[i] > oge->cov->rows){
        Rf_warning("casesToDrop vector in GREML expectation contains indices greater than the number of datapoints");
        oge->numcases2drop--;
      }
      //Need to subtract 1 from the index because R begins array indexing with 1, not 0:
      else{oge->dropcase[casesToDrop_intptr[i]-1] = true;}
  }}
  }
  if(Eigy.rows() != oge->cov->rows - oge->numcases2drop){
    mxThrow("y and V matrices do not have equal numbers of rows");
  }
  if(oge->yhatFromUser){
  	if(Eigy.rows() != oge->yhatFromUser->rows - oge->numcases2drop){
  		mxThrow("'y' and 'yhat' vectors have different numbers of rows");
  	}
  }
  

  //column names of y and X:
  {
  	ScopedProtect p1(RyXcolnames, R_do_slot(rObj, Rf_install("yXcolnames")));
  	if(Rf_length(RyXcolnames)){
  		oge->yXcolnames.resize(Rf_length(RyXcolnames));
  		for(i=0; i < Rf_length(RyXcolnames); i++){
  			SEXP elem;
  			{ScopedProtect p2(elem, STRING_ELT(RyXcolnames, i));
  				oge->yXcolnames[i] = CHAR(elem);}
  		}
  	}
  }

  oge->cholquadX_fail = false;
  EigenMatrixAdaptor EigX(oge->X);
  double *ptrToMatrix;
  Eigen::LLT< Eigen::MatrixXd > cholV(Eigy.rows());
  if( oge->numcases2drop && (oge->cov->rows > Eigy.rows()) ){
  	dropCasesAndEigenizeSquareMatrix(oge->cov, EigV_filtered, ptrToMatrix, oge->numcases2drop, oge->dropcase, true, oge->origVdim, false);
  }
  else{
  	ptrToMatrix = omxMatrixDataColumnMajor(oge->cov);
  }
  Eigen::Map< Eigen::MatrixXd > EigV( ptrToMatrix, Eigy.rows(), Eigy.rows() );
  //invcov:
  oge->invcov = omxInitMatrix(EigV.rows(), EigV.cols(), 1, currentState);
  Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(oge->invcov), EigV.rows(), EigV.cols());
  cholV.compute(EigV.selfadjointView<Eigen::Lower>());
  if(cholV.info() != Eigen::Success){
  	mxThrow("Expected covariance matrix is non-positive-definite at initial values");
  }
  oge->cholV_vectorD = (( Eigen::MatrixXd )(cholV.matrixL())).diagonal();
  for(i=0; i < oge->X->rows; i++){
  	oge->logdetV_om->data[0] += log(oge->cholV_vectorD[i]);
  }
  oge->logdetV_om->data[0] *= 2;
  Vinv = cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )); //<-- V inverse
  
  if(!oge->didUserProvideYhat){
  	Eigen::MatrixXd quadX(oge->X->cols, oge->X->cols);
  	//Apparently you need to initialize a matrix's elements before you try to write to its lower triangle:
  	quadX.setZero(oge->X->cols, oge->X->cols);
  	Eigen::LLT< Eigen::MatrixXd > cholquadX(oge->X->cols);
  	
  	oge->XtVinv = EigX.transpose() * Vinv;
  	quadX.triangularView<Eigen::Lower>() = oge->XtVinv * EigX;
  	cholquadX.compute(quadX.selfadjointView<Eigen::Lower>());
  	if(cholquadX.info() != Eigen::Success){
  		mxThrow("Cholesky factorization failed at initial values; possibly, the matrix of covariates is rank-deficient");
  	}
  	if(doREML){ oge->cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal(); }
  	oge->quadXinv = ( cholquadX.solve(Eigen::MatrixXd::Identity(oge->X->cols, oge->X->cols)) ).triangularView<Eigen::Lower>();
  }
  
  if(oge->didUserProvideYhat){
  	double *ptrToVector;
  	//Eigen::Map< Eigen::MatrixXd > Eigyhat(omxMatrixDataColumnMajor(oge->yhatFromUser->dataMat), oge->yhatFromUser->dataMat->cols, 1);
  	if( oge->numcases2drop && (oge->yhatFromUser->rows > Eigy.rows()) ){
  		dropCasesAndEigenizeColumnVector(oge->yhatFromUser, EigyhatFromUser_filtered, ptrToVector, oge->numcases2drop, oge->dropcase, false, oge->origVdim, false);
  	}
  	else{
  		ptrToVector = omxMatrixDataColumnMajor(oge->yhatFromUser);
  	}
  	Eigen::Map< Eigen::MatrixXd > Eigyhat( ptrToVector, Eigy.rows(), 1 );
  	residual = Eigy - Eigyhat;
  }

  /*Prepare y as the data that the FIML fitfunction will use:*/
  data2 = data;  // for safekeeping
  data = oge->y;
  if (data2->hasDefinitionVariables()) {
	  mxThrow("definition variables are incompatible (and unnecessary) with GREML expectation");
  }
}


void omxGREMLExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	omxGREMLExpectation* oge = this;
	omxRecompute(oge->cov, fc);
  int i=0;
  oge->cholV_fail = false;
  oge->cholquadX_fail = false;
  oge->logdetV_om->data[0] = 0;

  EigenMatrixAdaptor EigX(oge->X);
  Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(oge->y->dataMat), oge->y->dataMat->cols, 1);
  double *ptrToMatrix;
  Eigen::Map< Eigen::MatrixXd > Vinv(omxMatrixDataColumnMajor(oge->invcov), oge->invcov->rows, oge->invcov->cols);
  Eigen::LLT< Eigen::MatrixXd > cholV(oge->y->dataMat->rows);
  if( oge->numcases2drop && (oge->cov->rows > Eigy.rows()) ){
  	dropCasesAndEigenizeSquareMatrix(oge->cov, EigV_filtered, ptrToMatrix, oge->numcases2drop, oge->dropcase, true, oge->origVdim, false);
  }
  else{
  	ptrToMatrix = omxMatrixDataColumnMajor(oge->cov);
  }
  Eigen::Map< Eigen::MatrixXd > EigV( ptrToMatrix, Eigy.rows(), Eigy.rows() );
  cholV.compute(EigV.selfadjointView<Eigen::Lower>());
  if(cholV.info() != Eigen::Success){
  	oge->cholV_fail = true;
  	return;
  }
  oge->cholV_vectorD = (( Eigen::MatrixXd )(cholV.matrixL())).diagonal();
  for(i=0; i < oge->X->rows; i++){
  	oge->logdetV_om->data[0] += log(oge->cholV_vectorD[i]);
  }
  oge->logdetV_om->data[0] *= 2;
  //V inverse:
  Vinv.triangularView<Eigen::Lower>() = ( cholV.solve(Eigen::MatrixXd::Identity( EigV.rows(), EigV.cols() )) ).triangularView<Eigen::Lower>();
  
  if(!oge->didUserProvideYhat){
  	Eigen::MatrixXd quadX(oge->X->cols, oge->X->cols);
  	quadX.setZero(oge->X->cols, oge->X->cols);
  	Eigen::LLT< Eigen::MatrixXd > cholquadX(oge->X->cols);
  	
  	oge->XtVinv = EigX.transpose() * Vinv.selfadjointView<Eigen::Lower>();
  	quadX.triangularView<Eigen::Lower>() = oge->XtVinv * EigX;
  	cholquadX.compute(quadX.selfadjointView<Eigen::Lower>());
  	if(cholquadX.info() != Eigen::Success){
  		oge->cholquadX_fail = true;
  		return;
  	}
  	if(doREML){ oge->cholquadX_vectorD = (( Eigen::MatrixXd )(cholquadX.matrixL())).diagonal(); }
  	oge->quadXinv = ( cholquadX.solve(Eigen::MatrixXd::Identity(oge->X->cols, oge->X->cols)) ).triangularView<Eigen::Lower>();
  }
  
  if(oge->didUserProvideYhat){
  	omxRecompute(oge->yhatFromUser, fc);
  	double *ptrToVector;
  	//Eigen::Map< Eigen::MatrixXd > Eigyhat(omxMatrixDataColumnMajor(oge->yhatFromUser->dataMat), oge->yhatFromUser->dataMat->cols, 1);
  	if( oge->numcases2drop && (oge->yhatFromUser->rows > Eigy.rows()) ){
  		dropCasesAndEigenizeColumnVector(oge->yhatFromUser, EigyhatFromUser_filtered, ptrToVector, oge->numcases2drop, oge->dropcase, false, oge->yhatFromUser->rows, false);
  	}
  	else{
  		ptrToVector = omxMatrixDataColumnMajor(oge->yhatFromUser);
  	}
  	Eigen::Map< Eigen::MatrixXd > Eigyhat( ptrToVector, Eigy.rows(), 1 );
  	residual = Eigy - Eigyhat;
  }

	super::compute(fc, what, how);
}


omxGREMLExpectation::~omxGREMLExpectation()
{
	if(OMX_DEBUG) { mxLog("Destroying GREML Expectation."); }
	omxGREMLExpectation* argStruct = this;
  omxFreeMatrix(argStruct->means);
  omxFreeMatrix(argStruct->invcov);
  omxFreeMatrix(argStruct->logdetV_om);
  omxFreeMatrix(argStruct->y->dataMat);
  delete argStruct->y;
  omxFreeMatrix(argStruct->X);
  //omxFreeMatrix(argStruct->yhatFromUser);
}


/*Possible TODO: it will require some additional computation, but it is probably best to calculate the final
regression coefficients using QR, which is more numerically stable*/
void omxGREMLExpectation::populateAttr(SEXP algebra) {
  if(OMX_DEBUG) { mxLog("Populating GREML expectation attributes."); }

  omxGREMLExpectation* oge = this;
  
  SEXP RyXcolnames, Rresid, Rfitted;
  Eigen::Map< Eigen::MatrixXd > Eigy(omxMatrixDataColumnMajor(oge->y->dataMat), oge->y->dataMat->cols, 1);

  {
  ProtectedSEXP RnumStat(Rf_ScalarReal(oge->y->dataMat->cols));
  Rf_setAttrib(algebra, Rf_install("numStats"), RnumStat);
  }
  
  if(!didUserProvideYhat){
  	{
  		ProtectedSEXP RnumFixEff(Rf_ScalarInteger(oge->X->cols));
  		Rf_setAttrib(algebra, Rf_install("numFixEff"), RnumFixEff);
  	}
  	
  	SEXP b_ext, bcov_ext; 
  	oge->quadXinv = oge->quadXinv.selfadjointView<Eigen::Lower>();
  	Eigen::MatrixXd GREML_b = oge->quadXinv * (oge->XtVinv * Eigy);
  	
  	{
  		ScopedProtect p1(b_ext, Rf_allocMatrix(REALSXP, GREML_b.rows(), 1));
  		for(int row = 0; row < GREML_b.rows(); row++){
  			REAL(b_ext)[0 * GREML_b.rows() + row] = GREML_b(row,0);
  		}
  		Rf_setAttrib(algebra, Rf_install("b"), b_ext);
  	}
  	
  	{
  		ScopedProtect p1(bcov_ext, Rf_allocMatrix(REALSXP, oge->quadXinv.rows(),
                                              oge->quadXinv.cols()));
  		for(int row = 0; row < oge->quadXinv.rows(); row++){
  			for(int col = 0; col < oge->quadXinv.cols(); col++){
  				REAL(bcov_ext)[col * oge->quadXinv.rows() + row] = oge->quadXinv(row,col);
  			}}
  		Rf_setAttrib(algebra, Rf_install("bcov"), bcov_ext);
  	}
  	
  	EigenMatrixAdaptor EigX(oge->X);
  	Eigen::MatrixXd fitted = EigX * GREML_b;
  	Eigen::MatrixXd resid = Eigy - fitted;
  	
  	{
  		ScopedProtect p1(Rfitted, Rf_allocMatrix(REALSXP, Eigy.rows(), 1));
  		ScopedProtect p2(Rresid, Rf_allocMatrix(REALSXP, Eigy.rows(), 1));
  		for(int row=0; row < Eigy.rows(); row++){
  			REAL(Rfitted)[row] = fitted(row,0);
  			REAL(Rresid)[row] = resid(row,0);
  		}
  		Rf_setAttrib(algebra, Rf_install("fitted.values"), Rfitted);
  		Rf_setAttrib(algebra, Rf_install("residuals"), Rresid);
  	}
  }
  else{
  	Eigen::Map< Eigen::MatrixXd > Eigyhat(omxMatrixDataColumnMajor(yhatFromUser), Eigy.rows(), 1);
  	{
  		ScopedProtect p1(Rfitted, Rf_allocMatrix(REALSXP, Eigy.rows(), 1));
  		ScopedProtect p2(Rresid, Rf_allocMatrix(REALSXP, Eigy.rows(), 1));
  		for(int row=0; row < Eigy.rows(); row++){
  			REAL(Rfitted)[row] = Eigyhat(row,0);
  			REAL(Rresid)[row] = residual(row,0);
  		}
  		Rf_setAttrib(algebra, Rf_install("fitted.values"), Rfitted);
  		Rf_setAttrib(algebra, Rf_install("residuals"), Rresid);
  	}
  }

  //yXcolnames:
  {
  	if(oge->yXcolnames.size()){
  		ScopedProtect p1(RyXcolnames, Rf_allocVector(STRSXP, oge->yXcolnames.size()));
  		for(int i=0; i < (int)(oge->yXcolnames.size()); i++){
  			SET_STRING_ELT(RyXcolnames, i, Rf_mkChar(oge->yXcolnames[i]));
  		}
  		Rf_setAttrib(algebra, Rf_install("yXcolnames"), RyXcolnames);
  	}
  }
  

}

omxMatrix *omxGREMLExpectation::getComponent(const char* component){
/* Return appropriate parts of Expectation to the Fit Function */
  if(OMX_DEBUG) { mxLog("GREML expectation: %s requested--", component); }

  omxGREMLExpectation* oge = this;
	omxMatrix* retval = NULL;


  if(strEQ("y", component)) {
    retval = oge->y->dataMat;
	}
  else if(strEQ("invcov", component)) {
    retval = oge->invcov;
  }
  else if(strEQ("means", component)) {
  	retval = oge->means;
  }
  else if(strEQ("logdetV_om", component)){
    retval = oge->logdetV_om;
  }
  else if(strEQ("cov", component)) {
		retval = oge->cov;
	}
  else if(strEQ("X", component)) {
		retval = oge->X;
	}
  else if(strEQ("yhatFromUser", component)) {
  	retval = oge->yhatFromUser;
  }

	if (retval) omxRecompute(retval, NULL);

	return retval;
}



double omxAliasedMatrixElement(omxMatrix *om, int row, int col, int origRows, int origCols)
{
  int index = 0;
  if(row >= origRows || col >= origCols){
		mxThrow("Requested improper value (%d, %d) from (%d x %d) matrix %s",
           row + 1, col + 1, origRows, origCols, om->name());
		return (NA_REAL);
	}
	index = col * origRows + row; //<--om should always be column-major by this point.
	return om->data[index];
}



void dropCasesFromAlgdV(omxMatrix* om, int num2drop, std::vector< bool > &todrop, int symmetric, int origDim){

	if(OMX_DEBUG) { mxLog("Trimming out cases with missing data..."); }

	if(num2drop < 1 || om->algebra == NULL){ return; }

	omxEnsureColumnMajor(om);

	if(origDim==0){mxThrow("Memory not allocated for algebra %s at downsize time",
    om->name());}
	if(om->rows != origDim || om->cols != origDim){
		//Not sure if there are cases where this should be allowed
		mxThrow("More than one attempt made to downsize algebra %s", om->name());
		//return;
	}

	int nextCol = 0;
	int nextRow = 0;

	om->rows = origDim - num2drop;
	om->cols = origDim - num2drop;

	for(int j = 0; j < origDim; j++){ //<--j indexes columns
		if(todrop[j]) continue;
		nextRow = (symmetric ? nextCol : 0);
		for(int k = (symmetric ? j : 0); k < origDim; k++){ //<--k indexes rows
			if(todrop[k]) continue;
			omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j, origDim, origDim));
			nextRow++;
		}
		nextCol++;
	}
	omxMarkDirty(om); //<--Need to mark it dirty so that it eventually gets recalculated back to original dimensions.
	//^^^Algebras that do not depend upon free parameters, and upon which V does not depend, will not be
	//recalculated back to full size until optimization is complete (the GREML fitfunction is smart about that).
	if(OMX_DEBUG) { mxLog("Finished trimming out cases with missing data..."); }
}
