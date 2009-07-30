/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <R.h> 
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_FIML_OBJECTIVE_
#define _OMX_FIML_OBJECTIVE_ TRUE

extern omxMatrix** matrixList;

typedef struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
	int numLocations;		// Num locations
	double** location;		// And where it goes
	omxMatrix** matrices;	// Matrix numbers for dirtying

} omxDefinitionVar;

typedef struct omxFIMLObjective {

	omxMatrix* cov;
	omxMatrix* means;
	omxMatrix* data;
	omxMatrix* dataRow;
	omxMatrix* smallRow;
	omxMatrix* smallCov;
	omxMatrix* RCX;
//	double* zeros;
	omxDefinitionVar* defVars;
	int numDefs;

} omxFIMLObjective;

void handleDefinitionVarList(omxMatrix* data, int row, omxDefinitionVar* defVars, int numDefs) {

	if(OMX_DEBUG) { Rprintf("Processing Definition Vars.\n"); }

	/* Fill in Definition Var Estimates */
	for(int k = 0; k < numDefs; k++) {
		for(int l = 0; l < defVars[k].numLocations; l++) {
			if(OMX_DEBUG) {
				Rprintf("Populating column %d (value %3.2f) into matrix %d.\n", defVars[k].column, omxMatrixElement(data, row, defVars[k].column), defVars[k].matrices[l]);
			}
			*(defVars[k].location[l]) = omxMatrixElement(data, row, defVars[k].column);
			omxMarkDirty(defVars[k].matrices[l]);
			if(ISNA(omxMatrixElement(data, row, defVars[k].column))) {
				error("Error NYI: Missing Definition Vars Not Yet Implemented.");
			}
		}
	}
}

void omxDestroyFIMLObjective(omxObjective *oo) {

}

void omxCallFIMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n");}
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

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
	int numDefs;
	int nextRow, nextCol, numCols, numRemoves;

	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *data, *dataRow;
	omxDefinitionVar* defVars;

	cov 		= ((omxFIMLObjective*)oo->argStruct)->cov;		// Locals, for readability.  Should compile out.
	means		= ((omxFIMLObjective*)oo->argStruct)->means;
	smallRow 	= ((omxFIMLObjective*)oo->argStruct)->smallRow;
	smallCov 	= ((omxFIMLObjective*)oo->argStruct)->smallCov;
	RCX 		= ((omxFIMLObjective*)oo->argStruct)->RCX;
	data		= ((omxFIMLObjective*)oo->argStruct)->data;
	dataRow		= ((omxFIMLObjective*)oo->argStruct)->dataRow;
	defVars		= ((omxFIMLObjective*)oo->argStruct)->defVars;
	numDefs		= ((omxFIMLObjective*)oo->argStruct)->numDefs;

	if(numDefs == 0) {
		omxRecompute(cov);			// Only recompute this here if there are no definition vars
		omxRecompute(means);
	}
	
	if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }

	int toRemove[cov->cols];
	int ipiv[cov->rows];
	int lwork = 2*cov->rows;
	double work[lwork];
	int zeros[cov->cols];

	sum = 0.0;

	for(int row = 0; row < data->rows; row++) {
		oo->matrix->currentState->currentRow = row;		// Set to a new row.
		logDet = 1.0;
		Q = 0.0;

		// Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
		numCols = 0;
		numRemoves = 0;

		omxResetAliasedMatrix(smallRow); 										// Reset Row size
		
		// Handle Definition Variables.
		if(numDefs != 0) {
			handleDefinitionVarList(data, row, defVars, numDefs);
			omxStateNextRow(oo->matrix->currentState);							// Advance Row
			omxRecompute(cov);
			omxRecompute(means);
		}
		
		// Determine how many rows/cols to remove.
		for(int j = 0; j < dataRow->cols; j++) {
			int value = (int) omxMatrixElement(dataRow, 0, j);
			double dataValue = omxMatrixElement(data, row, j);
			if(isnan(dataValue) || dataValue == NA_REAL) {
				numRemoves++;
				toRemove[j] = 1;
			} else {
				if(means != NULL) {
					omxSetMatrixElement(smallRow, 0, j, (dataValue -  omxVectorElement(means, j)));
				} else {
					omxSetMatrixElement(smallRow, 0, j, dataValue);
				}
				numCols++;
				toRemove[j] = 0;
			}
			zeros[j] = 0;
		}
		
		if(dataRow->cols <= numRemoves) continue;
		omxRemoveRowsAndColumns(smallRow, 0, numRemoves, zeros, toRemove); 	// Reduce it.
		
		omxResetAliasedMatrix(smallCov);						// Subsample covariance matrix
		omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);
		
//		if(OMX_DEBUG) { omxPrint(smallCov, "Local Covariance Matrix"); }

		/* The Calculation */
		F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			char errStr[250];
			sprintf(errStr, "Backing out of parameter space region where covariance matrix is not positive-definite.");
			omxRaiseError(oo->matrix->currentState, -1, errStr);
			return;
		}
		for(int diag = 0; diag < (smallCov->rows); diag++) {
			logDet *= omxMatrixElement(smallCov, diag, diag);
		}
		logDet = log(logDet * logDet);
		
		F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			char errstr[250];
			sprintf(errstr, "Cannot invert covariance matrix. Error %d.", info);
			omxRaiseError(oo->matrix->currentState, -1, errstr);
			return;
		}
//		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);

		sum += logDet + Q + (log(2 * M_PI) * smallRow->cols);
		if(OMX_DEBUG) {Rprintf("Change in Total Likelihood is %3.3f + %3.3f + %3.3f = %3.3f, total Likelihood is %3.3f\n", logDet, Q, (log(2 * M_PI) * smallRow->cols), logDet + Q + (log(2 * M_PI) * smallRow->cols), sum);}
	}

	if(OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
	oo->matrix->data[0] = sum;

}

unsigned short int omxNeedsUpdateFIMLObjective(omxObjective* oo) {
	return omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->cov)
		|| omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->means);
}

void omxInitFIMLObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {

	if(OMX_DEBUG) { Rprintf("Initializing FIML objective function.\n"); }

	SEXP nextMatrix, itemList, nextItem, dataSource, columnSource;
	int nextDef, index, data, column;
	omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(1, sizeof(omxFIMLObjective));

	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	newObj->means = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->means == NULL) { error("No means in FIML evaluation.");}
	UNPROTECT(1);

	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	newObj->cov = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing data source.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data"))); // TODO: Need better way to process data elements.
	index = (int) REAL(nextMatrix)[0];
	PROTECT(nextMatrix = VECTOR_ELT(dataList, index));
	PROTECT(nextMatrix = GET_SLOT(nextMatrix, install("observed")));
	newObj->data = omxNewMatrixFromMxMatrix(nextMatrix, oo->matrix->currentState);
	UNPROTECT(3);
	
	if(OMX_DEBUG) {Rprintf("Accessing variable mapping structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("dataRow")));
	newObj->dataRow = omxNewMatrixFromMxMatrix(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);
	
	if(OMX_DEBUG) {Rprintf("Accessing definition variables structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
	newObj->numDefs = length(nextMatrix);
	if(OMX_DEBUG) {Rprintf("Number of definition variables is %d.\n", newObj->numDefs); }
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
		PROTECT(itemList = VECTOR_ELT(nextMatrix, nextDef));
		PROTECT(dataSource = VECTOR_ELT(itemList, 0));
		if(OMX_DEBUG) {Rprintf("Data source number is %d.\n", (int) REAL(dataSource)[0]); }
		newObj->defVars[nextDef].data = (int) REAL(dataSource)[0];
		PROTECT(columnSource = VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG) {Rprintf("Data column number is %d.\n", (int) REAL(columnSource)[0]); }
		newObj->defVars[nextDef].column = (int) REAL(columnSource)[0];
		UNPROTECT(2); // unprotect dataSource and columnSource
		newObj->defVars[nextDef].numLocations = length(itemList) - 2;
		newObj->defVars[nextDef].location = (double **) R_alloc(length(itemList) - 2, sizeof(double*));
		newObj->defVars[nextDef].matrices = (omxMatrix **) R_alloc(length(itemList) - 2, sizeof(omxMatrix*));
		for(index = 2; index < length(itemList); index++) {
			PROTECT(nextItem = VECTOR_ELT(itemList, index));
			newObj->defVars[nextDef].location[index-2] = omxLocationOfMatrixElement(
				oo->matrix->currentState->matrixList[(int) REAL(nextItem)[0]],
				(int) REAL(nextItem)[1], (int) REAL(nextItem)[2]);
			newObj->defVars[nextDef].matrices[index-2] = oo->matrix->currentState->matrixList[(int) REAL(nextItem)[0]];
			UNPROTECT(1); // unprotect nextItem
		}
		UNPROTECT(1); // unprotect itemList
	}
	UNPROTECT(1); // unprotect nextMatrix

	/* Temporary storage for calculation */
	newObj->smallRow = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);
	newObj->smallCov = omxInitMatrix(NULL, newObj->cov->rows, newObj->cov->cols, TRUE, oo->matrix->currentState);
	newObj->RCX = omxInitMatrix(NULL, 1, newObj->data->cols, TRUE, oo->matrix->currentState);
//	newObj->zeros = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);

	omxAliasMatrix(newObj->smallCov, newObj->cov);					// Will keep its aliased state from here on.
	// We can alias smallrow to itself because: 1) it has local data 2) we populate it explicitly in FIML.
	// TODO: do this more cleanly.
	omxAliasMatrix(newObj->smallRow, newObj->smallRow);
	
	oo->objectiveFun = omxCallFIMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateFIMLObjective;
	oo->destructFun = omxDestroyFIMLObjective;
	oo->repopulateFun = NULL;

	oo->argStruct = (void*) newObj;
	if(OMX_DEBUG) {Rprintf("Finished import of FIML objective function.\n"); }
}

#endif /* _OMX_FIML_OBJECTIVE_ */
