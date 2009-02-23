#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_FIML_OBJECTIVE_
#define _OMX_FIML_OBJECTIVE_ TRUE

typedef struct {

	omxMatrix* cov;
	omxMatrix* means;
	omxMatrix* data;
	omxMatrix* smallRow;
	omxMatrix* smallCov;
	omxMatrix* RCX;

} omxFIMLObjective;

void omxDestroyFIMLObjective(omxObjective *oo) {

}

void omxCallFIMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n");}
	// Requires: Data, means, covariances.

	SEXP matrixDims;
	int *dimList;
	double sum;
	char u = 'U';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	int onei = 1;
	int mainDist = 0;
	double Q = 0.0;
	double logDet = 0;
	int nextRow, nextCol, numCols, numRemoves;

	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *dataRows;

	cov 		= ((omxFIMLObjective*)oo->argStruct)->cov;			// Locals, for readability.  Compiler should cut through this.
	means		= ((omxFIMLObjective*)oo->argStruct)->means;
	smallRow 	= ((omxFIMLObjective*)oo->argStruct)->smallRow;
	smallCov 	= ((omxFIMLObjective*)oo->argStruct)->smallCov;
	RCX 		= ((omxFIMLObjective*)oo->argStruct)->RCX;
	dataRows	= ((omxFIMLObjective*)oo->argStruct)->data;
	
	omxRecomputeMatrix(cov);							// We assume data won't need to be recomputed.

	int toRemove[cov->cols];

	sum = 0.0;

	for(int row = 0; row < dataRows->rows; row++) {
		logDet = 0.0;
		Q = 0.0;

		/** HANDLE MISSINGNESS HERE **/
		// Note:  This really aught to be done using a matrix multiply.  Why isn't it?
		numCols = 0;
		numRemoves = 0;

		omxResizeMatrix(smallRow, 1, cov->cols, FALSE); // Reset Row size //TODO: Test to see if aliasing is faster.

		// Determine how many rows/cols to remove.
		for(int j = 0; j < dataRows->cols; j++) {
			if(ISNA(omxMatrixElement(dataRows, row, j))) {
				numRemoves++;
				toRemove[j] = 1;
			} else {
				toRemove[j] = 0;

			}
		}
		
		omxResizeMatrix(smallRow, 1, cov->cols - numRemoves, TRUE);	// Subsample this Row
		
		for(int j = 0; j < dataRows->cols; j++) {
			if(!toRemove[j])
				omxSetMatrixElement(smallRow, 0, numCols++, omxMatrixElement(dataRows, row, j));
		}
		if(numCols==0) continue;

		omxResetAliasedMatrix(smallCov);						// Subsample covariance matrix
		omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);


		/* The Calculation */
		F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) error("Covariance Matrix is not positive-definite.");
		for(int diag = 0; diag < (smallCov->rows); diag++) {
			logDet += log(fabs(smallCov->data[diag + (diag * smallCov->rows)]));
		}
		logDet *= 2.0;
		F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) error("Cannot invert covariance matrix.");
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		for(int col = 0; col < smallRow->cols; col++) {
			Q += RCX->data[col] * smallRow->data[col];
		}
		sum += logDet + Q + (log(2 * M_PI) * smallRow->cols);
	
	}
	
	oo->myMatrix->data[0] = sum;

}

unsigned short int omxNeedsUpdateFIMLObjective(omxObjective* oo) {
	return omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->cov)
		|| omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->means);
}

void omxInitFIMLObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	SEXP nextMatrix;
	int index;
	omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(sizeof(omxFIMLObjective), 1);
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	newObj->means = omxNewMatrixFromMxMatrixPtr(nextMatrix);
	UNPROTECT(1);
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	newObj->cov = omxNewMatrixFromMxMatrixPtr(nextMatrix);
	UNPROTECT(1);
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));   // TODO: Need better way to process data elements.
	index = round(REAL(nextMatrix)[0]);
	PROTECT(nextMatrix = VECTOR_ELT(dataList, index));
	newObj->data = omxNewMatrixFromMxMatrix(nextMatrix);
	UNPROTECT(2);
	
	/* Temporary storage for calculation */
	newObj->smallRow = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE);
	newObj->smallCov = omxInitMatrix(NULL, newObj->cov->rows, newObj->cov->cols, TRUE);
	newObj->RCX = omxInitMatrix(NULL, 1, newObj->data->cols, TRUE);
	
	omxAliasMatrix(newObj->smallCov, newObj->cov);					// Will keep its aliased state from here on.
	
	oo->needsUpdateFun = omxNeedsUpdateFIMLObjective;
	oo->argStruct = (void*) newObj;
	
}


#endif /* _OMX_FIML_OBJECTIVE_ */