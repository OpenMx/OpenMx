/*
 *  Copyright 2007-2010 The OpenMx Project
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
#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"

#ifndef _OMX_ROW_OBJECTIVE_
#define _OMX_ROW_OBJECTIVE_ TRUE

// TODO: Migrate omxDefinitionVar struct to a more central location.
 typedef struct omxDefinitionVar {		 	// Definition Var

 	int data, column;		// Where it comes from
 	omxData* source;		// Data source
 	int numLocations;		// Num locations
 	double** location;		// And where it goes
 	omxMatrix** matrices;	// Matrix numbers for dirtying

 } omxDefinitionVar;

 extern int handleDefinitionVarList(omxData* data, int row, omxDefinitionVar* defVars, double* oldDefs, int numDefs);

typedef struct omxRowObjective {

	/* Parts of the R  MxRowObjective Object */
	omxMatrix* rowAlgebra;		// Row-by-row algebra
	omxMatrix* rowResults;		// Aggregation of row algebra results
	omxMatrix* reduceAlgebra;	// Algebra performed after row-by-row computation
    omxMatrix* filteredDataRow; // Data row minus NAs
    omxMatrix* existenceVector; // Set of NAs
    omxMatrix* dataColumns;		// The order of columns in the data matrix

    /* Contiguous data note for contiguity speedup */
	omxContiguousData contiguous;		// Are the dataColumns contiguous within the data set

	/* Structures determined from info in the MxRowObjective Object*/
    double* oldDefs;            // The previous defVar vector.  To avoid recalculations where possible.         // NYI: This element currently unused.
	omxDefinitionVar* defVars;	// A list of definition variables
	int numDefs;				// The length of the defVars list
	omxMatrix* dataRow;         // One row of data, kept for aliasing only
	omxData*   data;			// The data
    

} omxRowObjective;

void omxDestroyRowObjective(omxObjective *oo) {

}

omxRListElement* omxSetFinalReturnsRowObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 0;
	omxRListElement* retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));

	retVal[0].numValues = 0;

	return retVal;
}


void omxCopyMatrixToRow(omxMatrix* source, int row, omxMatrix* target) {
	
	int i;
	for(i = 0; i < source->cols; i++) {
		omxSetMatrixElement(target, row, i, omxMatrixElement(source, 0, i));
	}

}

void omxCallRowObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
    if(OMX_DEBUG) { Rprintf("Beginning Row Evaluation.\n");}
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

	int numDefs;

    omxMatrix *rowAlgebra, *rowResults, *reduceAlgebra;
    omxMatrix *filteredDataRow, *dataRow, *existenceVector;
    omxMatrix *dataColumns;
	omxDefinitionVar* defVars;
	omxData *data;
	int isContiguous, contiguousStart, contiguousLength;
    double* oldDefs;
    int numCols, numRemoves;

    omxRowObjective* oro = ((omxRowObjective*)oo->argStruct);

	rowAlgebra	    = oro->rowAlgebra;	 // Locals, for readability.  Should compile out.
	rowResults	    = oro->rowResults;
	reduceAlgebra   = oro->reduceAlgebra;
	data		    = oro->data;
	defVars		    = oro->defVars;
	numDefs		    = oro->numDefs;
    oldDefs         = oro->oldDefs;
    dataColumns     = oro->dataColumns;
    dataRow         = oro->dataRow;
    filteredDataRow = oro->filteredDataRow;
    existenceVector = oro->existenceVector;
    
    isContiguous    = oro->contiguous.isContiguous;
	contiguousStart = oro->contiguous.start;
	contiguousLength = oro->contiguous.length;
	
	
	if(rowResults->cols != rowAlgebra->cols || rowResults->rows != data->rows) {
		if(OMX_DEBUG_ROWS) { Rprintf("Resizing rowResults from %dx%d to %dx%d.\n", rowResults->rows, rowResults->cols, data->rows, rowAlgebra->cols); }
		omxResizeMatrix(rowResults, data->rows, rowAlgebra->cols, FALSE);
	}
		
    numCols = dataColumns->cols;
	int toRemove[dataColumns->cols];
	int zeros[dataColumns->cols];
	memset(zeros, 0, sizeof(int) * dataColumns->cols);  // Shouldn't be required.

	for(int row = 0; row < data->rows; row++) {

		// Handle Definition Variables.
        if(OMX_DEBUG_ROWS) { Rprintf("numDefs is %d", numDefs);}
		if(numDefs != 0) {		// With no defs, just copy repeatedly to the rowResults matrix.
			handleDefinitionVarList(data, row, defVars, oldDefs, numDefs);
		}

		omxStateNextRow(oo->matrix->currentState);						// Advance row
		
        // Populate data row
		numRemoves = 0;
	
		if (isContiguous) {
			omxContiguousDataRow(data, row, contiguousStart, contiguousLength, dataRow);
		} else {
			omxDataRow(data, row, dataColumns, dataRow);	// Populate data row
		}
		
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxMatrixElement(dataRow, 0, j);
			if(isnan(dataValue) || dataValue == NA_REAL) {
				numRemoves++;
				toRemove[j] = 1;
                omxSetMatrixElement(existenceVector, 0, j, 0);
			} else {
			    toRemove[j] = 0;
                omxSetMatrixElement(existenceVector, 0, j, 1);
			}
		}		
		// TODO: Determine if this is the correct response.
		
		if(numRemoves == numCols) {
		    char *errstr = calloc(250, sizeof(char));
			sprintf(errstr, "Row %d completely missing.  omxRowObjective cannot have completely missing rows.", omxDataIndex(data, row));
			omxRaiseError(oo->matrix->currentState, -1, errstr);
			free(errstr);
			return;
		}

		omxResetAliasedMatrix(filteredDataRow); 			// Reset the row
		omxRemoveRowsAndColumns(filteredDataRow, 0, numRemoves, zeros, toRemove);

		omxRecompute(rowAlgebra);							// Compute this row

		omxCopyMatrixToRow(rowAlgebra, omxDataIndex(data, row), rowResults);
	}

	omxRecompute(reduceAlgebra);

	omxCopyMatrix(oo->matrix, reduceAlgebra);

}

unsigned short int omxNeedsUpdateRowObjective(omxObjective* oo) {
	return omxMatrixNeedsUpdate(((omxRowObjective*)oo->argStruct)->rowAlgebra)
		|| omxMatrixNeedsUpdate(((omxRowObjective*)oo->argStruct)->reduceAlgebra);
}

void omxInitRowObjective(omxObjective* oo, SEXP rObj) {

	if(OMX_DEBUG) { Rprintf("Initializing Row/Reduce objective function.\n"); }

	SEXP nextMatrix, itemList, nextItem, dataSource, columnSource;
	int nextDef, index;
	omxRowObjective *newObj = (omxRowObjective*) R_alloc(1, sizeof(omxRowObjective));

	if(OMX_DEBUG) {Rprintf("Accessing data source.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));
	newObj->data = omxNewDataFromMxDataPtr(nextMatrix, oo->matrix->currentState);
	if(newObj->data == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No data provided to omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	UNPROTECT(1); // nextMatrix

	PROTECT(nextMatrix = GET_SLOT(rObj, install("rowAlgebra")));
	newObj->rowAlgebra = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->rowAlgebra == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No row-wise algebra in omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	UNPROTECT(1);// nextMatrix

	PROTECT(nextMatrix = GET_SLOT(rObj, install("filteredDataRow")));
	newObj->filteredDataRow = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->filteredDataRow == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	// Create the original data row from which to filter.
    newObj->dataRow = omxInitMatrix(NULL, newObj->filteredDataRow->rows, newObj->filteredDataRow->cols, TRUE, oo->matrix->currentState);
    omxAliasMatrix(newObj->filteredDataRow, newObj->dataRow);
	UNPROTECT(1);// nextMatrix

	PROTECT(nextMatrix = GET_SLOT(rObj, install("existenceVector")));
	newObj->existenceVector = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
    // Do we allow NULL existence?  (Whoa, man. That's, like, deep, or something.)
	if(newObj->existenceVector == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No existance matrix in omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	UNPROTECT(1);// nextMatrix


	PROTECT(nextMatrix = GET_SLOT(rObj, install("rowResults")));
	newObj->rowResults = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->rowResults == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	UNPROTECT(1);// nextMatrix

	PROTECT(nextMatrix = GET_SLOT(rObj, install("reduceAlgebra")));
	newObj->reduceAlgebra = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->reduceAlgebra == NULL) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "No row reduction algebra in omxRowObjective.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	UNPROTECT(1);// nextMatrix
	
	if(OMX_DEBUG) {Rprintf("Accessing variable mapping structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("dataColumns")));
	newObj->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, oo->matrix->currentState, 0, 0);
	if(OMX_DEBUG) { omxPrint(newObj->dataColumns, "Variable mapping"); }
	UNPROTECT(1);
	

	if(OMX_DEBUG) {Rprintf("Accessing definition variables structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
	newObj->numDefs = length(nextMatrix);
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars
	if(OMX_DEBUG) {Rprintf("Number of definition variables is %d.\n", newObj->numDefs); }
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
		PROTECT(itemList = VECTOR_ELT(nextMatrix, nextDef));
		PROTECT(dataSource = VECTOR_ELT(itemList, 0));
		if(OMX_DEBUG) {Rprintf("Data source number is %d.\n", INTEGER(dataSource)[0]); }
		newObj->defVars[nextDef].data = INTEGER(dataSource)[0];
		newObj->defVars[nextDef].source = oo->matrix->currentState->dataList[INTEGER(dataSource)[0]];
		PROTECT(columnSource = VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG) {Rprintf("Data column number is %d.\n", INTEGER(columnSource)[0]); }
		newObj->defVars[nextDef].column = INTEGER(columnSource)[0];
		UNPROTECT(2); // unprotect dataSource and columnSource
		newObj->defVars[nextDef].numLocations = length(itemList) - 2;
		newObj->defVars[nextDef].location = (double **) R_alloc(length(itemList) - 2, sizeof(double*));
		newObj->defVars[nextDef].matrices = (omxMatrix **) R_alloc(length(itemList) - 2, sizeof(omxMatrix*));
		for(index = 2; index < length(itemList); index++) {
			PROTECT(nextItem = VECTOR_ELT(itemList, index));
			newObj->defVars[nextDef].location[index-2] = omxLocationOfMatrixElement(
				oo->matrix->currentState->matrixList[INTEGER(nextItem)[0]],
				INTEGER(nextItem)[1], INTEGER(nextItem)[2]);
			newObj->defVars[nextDef].matrices[index-2] = oo->matrix->currentState->matrixList[INTEGER(nextItem)[0]];
			UNPROTECT(1); // unprotect nextItem
		}
		UNPROTECT(1); // unprotect itemList
	}
	UNPROTECT(1); // unprotect nextMatrix
	
	/* Set up data columns */
	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);

	oo->objectiveFun = omxCallRowObjective;
	oo->needsUpdateFun = omxNeedsUpdateRowObjective;
	oo->setFinalReturns = omxSetFinalReturnsRowObjective;
	oo->destructFun = omxDestroyRowObjective;
	oo->repopulateFun = NULL;

	oo->argStruct = (void*) newObj;
}

#endif /* _OMX_ROW_OBJECTIVE_ */
