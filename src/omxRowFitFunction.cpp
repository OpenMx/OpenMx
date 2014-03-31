/*
 *  Copyright 2007-2014 The OpenMx Project
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

#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxRowFitFunction.h"
#include "omxFIMLFitFunction.h"

void omxDestroyRowFitFunction(omxFitFunction *oo) {

	omxRowFitFunction* argStruct = (omxRowFitFunction*)(oo->argStruct);

	omxFreeMatrix(argStruct->dataRow);
}

void omxCopyMatrixToRow(omxMatrix* source, int row, omxMatrix* target) {
	
	int i;
	for(i = 0; i < source->cols; i++) {
		omxSetMatrixElement(target, row, i, omxMatrixElement(source, 0, i));
	}

}

void markDataRowDependencies(omxState* os, omxRowFitFunction* orff) {

	int numDeps = orff->numDataRowDeps;
	int *deps = orff->dataRowDeps;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(os->matrixList[~value]);
		} else {
			omxMarkDirty(os->algebraList[value]);
		}
	}

}

void omxRowFitFunctionSingleIteration(omxFitFunction *localobj, omxFitFunction *sharedobj, int rowbegin, int rowcount) {

    omxRowFitFunction* oro = ((omxRowFitFunction*) localobj->argStruct);
    omxRowFitFunction* shared_oro = ((omxRowFitFunction*) sharedobj->argStruct);

	int numDefs;

    omxMatrix *rowAlgebra, *rowResults;
    omxMatrix *filteredDataRow, *dataRow, *existenceVector;
    omxMatrix *dataColumns;
	omxDefinitionVar* defVars;
	omxData *data;
	int isContiguous, contiguousStart, contiguousLength;
    double* oldDefs;
    int numCols, numRemoves;

	rowAlgebra	    = oro->rowAlgebra;
	rowResults	    = shared_oro->rowResults;
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

	numCols = dataColumns->cols;
	int *toRemove = (int*) malloc(sizeof(int) * dataColumns->cols);
	int *zeros = (int*) calloc(dataColumns->cols, sizeof(int));

    resetDefinitionVariables(oldDefs, numDefs);

	for(int row = rowbegin; row < data->rows && (row - rowbegin) < rowcount; row++) {

		// Handle Definition Variables.
        if(OMX_DEBUG_ROWS(row)) { mxLog("numDefs is %d", numDefs);}
		if(numDefs != 0) {		// With defs, just copy repeatedly to the rowResults matrix.
			handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs);
		}

		omxStateNextRow(localobj->matrix->currentState);						// Advance row
		
        // Populate data row
		numRemoves = 0;
	
		if (isContiguous) {
			omxContiguousDataRow(data, row, contiguousStart, contiguousLength, dataRow);
		} else {
			omxDataRow(data, row, dataColumns, dataRow);	// Populate data row
		}

		markDataRowDependencies(localobj->matrix->currentState, oro);
		
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxVectorElement(dataRow, j);
			if(std::isnan(dataValue)) {
				numRemoves++;
				toRemove[j] = 1;
                omxSetVectorElement(existenceVector, j, 0);
			} else {
			    toRemove[j] = 0;
                omxSetVectorElement(existenceVector, j, 1);
			}
		}		
		// TODO: Determine if this is the correct response.
		
		if(numRemoves == numCols) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Row %d completely missing.  omxRowFitFunction cannot have completely missing rows.", omxDataIndex(data, row));
			omxRaiseError(localobj->matrix->currentState, -1, errstr);
			free(errstr);
			continue;
		}

		omxResetAliasedMatrix(filteredDataRow); 			// Reset the row
		omxRemoveRowsAndColumns(filteredDataRow, 0, numRemoves, zeros, toRemove);

		omxRecompute(rowAlgebra);							// Compute this row

		omxCopyMatrixToRow(rowAlgebra, omxDataIndex(data, row), rowResults);
	}
	free(toRemove);
	free(zeros);
}

static void omxCallRowFitFunction(omxFitFunction *oo, int want, FitContext *) {
	if (want & (FF_COMPUTE_PREOPTIMIZE)) return;

    if(OMX_DEBUG) { mxLog("Beginning Row Evaluation.");}
	// Requires: Data, means, covariances.

	omxMatrix* objMatrix  = oo->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState==globalState? Global->numChildren : 0;

    omxMatrix *reduceAlgebra;
	omxData *data;

    omxRowFitFunction* oro = ((omxRowFitFunction*) oo->argStruct);

	reduceAlgebra   = oro->reduceAlgebra;
	data		    = oro->data;

	/* Michael Spiegel, 7/31/12
	* The demo "RowFitFunctionSimpleExamples" will fail in the parallel 
	* Hessian calculation if the resizing operation is performed.
	*
	omxMatrix *rowAlgebra, *rowResults
	rowAlgebra	    = oro->rowAlgebra;
	rowResults	    = oro->rowResults;

	if(rowResults->cols != rowAlgebra->cols || rowResults->rows != data->rows) {
		if(OMX_DEBUG_ROWS(1)) { 
			mxLog("Resizing rowResults from %dx%d to %dx%d.", 
				rowResults->rows, rowResults->cols, 
				data->rows, rowAlgebra->cols); 
		}
		omxResizeMatrix(rowResults, data->rows, rowAlgebra->cols, FALSE);
	}
	*/
		
    int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				omxRowFitFunctionSingleIteration(childFit, oo, stride * i, data->rows - stride * i);
			} else {
				omxRowFitFunctionSingleIteration(childFit, oo, stride * i, stride);
			}
		}
	} else {
		omxRowFitFunctionSingleIteration(oo, oo, 0, data->rows);
	}

	omxRecompute(reduceAlgebra);

	omxCopyMatrix(oo->matrix, reduceAlgebra);

}

void omxInitRowFitFunction(omxFitFunction* oo) {

	if(OMX_DEBUG) { mxLog("Initializing Row/Reduce fit function."); }

	SEXP rObj = oo->rObj;
	SEXP nextMatrix, itemList, nextItem;
	int nextDef, index, numDeps;

	omxRowFitFunction *newObj = (omxRowFitFunction*) R_alloc(1, sizeof(omxRowFitFunction));

	if(OMX_DEBUG) {mxLog("Accessing data source."); }
	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("data")));
	newObj->data = omxDataLookupFromState(nextMatrix, oo->matrix->currentState);
	if(newObj->data == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No data provided to omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	Rf_unprotect(1); // nextMatrix

	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("rowAlgebra")));
	newObj->rowAlgebra = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	if(newObj->rowAlgebra == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row-wise algebra in omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	Rf_unprotect(1);// nextMatrix

	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("filteredDataRow")));
	newObj->filteredDataRow = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	if(newObj->filteredDataRow == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	// Create the original data row from which to filter.
    newObj->dataRow = omxInitMatrix(NULL, newObj->filteredDataRow->rows, newObj->filteredDataRow->cols, TRUE, oo->matrix->currentState);
    omxAliasMatrix(newObj->filteredDataRow, newObj->dataRow);
	Rf_unprotect(1);// nextMatrix

	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("existenceVector")));
	newObj->existenceVector = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
    // Do we allow NULL existence?  (Whoa, man. That's, like, deep, or something.)
	if(newObj->existenceVector == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No existance matrix in omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	Rf_unprotect(1);// nextMatrix


	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("rowResults")));
	newObj->rowResults = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	if(newObj->rowResults == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	Rf_unprotect(1);// nextMatrix

	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("reduceAlgebra")));
	newObj->reduceAlgebra = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	if(newObj->reduceAlgebra == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row reduction algebra in omxRowFitFunction.");
		omxRaiseError(oo->matrix->currentState, -1, errstr);
		free(errstr);
	}
	Rf_unprotect(1);// nextMatrix
	
	if(OMX_DEBUG) {mxLog("Accessing variable mapping structure."); }
	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("dataColumns")));
	newObj->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, oo->matrix->currentState, 0, 0);
	if(OMX_DEBUG) { omxPrint(newObj->dataColumns, "Variable mapping"); }
	Rf_unprotect(1);

	if(OMX_DEBUG) {mxLog("Accessing data row dependencies."); }
	Rf_protect(nextItem = R_do_slot(rObj, Rf_install("dataRowDeps")));
	numDeps = LENGTH(nextItem);
	newObj->numDataRowDeps = numDeps;
	newObj->dataRowDeps = (int*) R_alloc(numDeps, sizeof(int));
	for(int i = 0; i < numDeps; i++) {
		newObj->dataRowDeps[i] = INTEGER(nextItem)[i];
	}
	Rf_unprotect(1);

	if(OMX_DEBUG) {mxLog("Accessing definition variables structure."); }
	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("definitionVars")));
	newObj->numDefs = Rf_length(nextMatrix);
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars
	if(OMX_DEBUG) {mxLog("Number of definition variables is %d.", newObj->numDefs); }
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
		SEXP dataSource, columnSource, depsSource; 

		Rf_protect(itemList = VECTOR_ELT(nextMatrix, nextDef));
		Rf_protect(dataSource = VECTOR_ELT(itemList, 0));
		if(OMX_DEBUG) {mxLog("Data source number is %d.", INTEGER(dataSource)[0]); }
		newObj->defVars[nextDef].data = INTEGER(dataSource)[0];
		newObj->defVars[nextDef].source = oo->matrix->currentState->dataList[INTEGER(dataSource)[0]];
		Rf_protect(columnSource = VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG) {mxLog("Data column number is %d.", INTEGER(columnSource)[0]); }
		newObj->defVars[nextDef].column = INTEGER(columnSource)[0];
		Rf_protect(depsSource = VECTOR_ELT(itemList, 2));
		numDeps = LENGTH(depsSource);
		newObj->defVars[nextDef].numDeps = numDeps;
		newObj->defVars[nextDef].deps = (int*) R_alloc(numDeps, sizeof(int));
		for(int i = 0; i < numDeps; i++) {
			newObj->defVars[nextDef].deps[i] = INTEGER(depsSource)[i];
		}
		Rf_unprotect(3); // unprotect dataSource, columnSource, and depsSource

		newObj->defVars[nextDef].numLocations = Rf_length(itemList) - 3;
		newObj->defVars[nextDef].matrices = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));
		newObj->defVars[nextDef].rows = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));
		newObj->defVars[nextDef].cols = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));

		for(index = 3; index < Rf_length(itemList); index++) {
			Rf_protect(nextItem = VECTOR_ELT(itemList, index));
			newObj->defVars[nextDef].matrices[index-3] = INTEGER(nextItem)[0];
			newObj->defVars[nextDef].rows[index-3]     = INTEGER(nextItem)[1];
			newObj->defVars[nextDef].cols[index-3]     = INTEGER(nextItem)[2];
			Rf_unprotect(1); // unprotect nextItem
		}
		Rf_unprotect(1); // unprotect itemList
	}
	Rf_unprotect(1); // unprotect nextMatrix
	
	/* Set up data columns */
	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);

	oo->computeFun = omxCallRowFitFunction;
	oo->destructFun = omxDestroyRowFitFunction;
	oo->usesChildModels = TRUE;

	oo->argStruct = (void*) newObj;
}


