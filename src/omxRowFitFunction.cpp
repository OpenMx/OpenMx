/*
 *  Copyright 2007-2015 The OpenMx Project
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
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxRowFitFunction.h"
#include "omxFIMLFitFunction.h"

void omxDestroyRowFitFunction(omxFitFunction *oo) {

	omxRowFitFunction* argStruct = (omxRowFitFunction*)(oo->argStruct);

	omxFreeMatrix(argStruct->dataRow);
	omxFreeMatrix(argStruct->dataColumns);
	delete argStruct;
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

static void omxRowFitFunctionSingleIteration(omxFitFunction *localobj, omxFitFunction *sharedobj, int rowbegin, int rowcount,
					     FitContext *fc) {

    omxRowFitFunction* oro = ((omxRowFitFunction*) localobj->argStruct);
    omxRowFitFunction* shared_oro = ((omxRowFitFunction*) sharedobj->argStruct);

    omxMatrix *rowAlgebra, *rowResults;
    omxMatrix *filteredDataRow, *dataRow, *existenceVector;
    omxMatrix *dataColumns;
	omxData *data;
	int isContiguous, contiguousStart, contiguousLength;
    int numCols;

	rowAlgebra	    = oro->rowAlgebra;
	rowResults	    = shared_oro->rowResults;
	data		    = oro->data;
    dataColumns     = oro->dataColumns;
    dataRow         = oro->dataRow;
    filteredDataRow = oro->filteredDataRow;
    existenceVector = oro->existenceVector;
    
    isContiguous    = oro->contiguous.isContiguous;
	contiguousStart = oro->contiguous.start;
	contiguousLength = oro->contiguous.length;

	Eigen::VectorXd oldDefs;
	oldDefs.resize(data->defVars.size());
	oldDefs.setConstant(NA_REAL);

	numCols = dataColumns->cols;
	int *toRemove = (int*) malloc(sizeof(int) * dataColumns->cols);
	int *zeros = (int*) calloc(dataColumns->cols, sizeof(int));

	for(int row = rowbegin; row < data->rows && (row - rowbegin) < rowcount; row++) {
		mxLogSetCurrentRow(row);

		data->handleDefinitionVarList(localobj->matrix->currentState, row, oldDefs.data());

        // Populate data row
		if (isContiguous) {
			omxContiguousDataRow(data, row, contiguousStart, contiguousLength, dataRow);
		} else {
			omxDataRow(data, row, dataColumns, dataRow);	// Populate data row
		}

		markDataRowDependencies(localobj->matrix->currentState, oro);
		
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxVectorElement(dataRow, j);
			if(std::isnan(dataValue)) {
				toRemove[j] = 1;
                omxSetVectorElement(existenceVector, j, 0);
			} else {
			    toRemove[j] = 0;
                omxSetVectorElement(existenceVector, j, 1);
			}
		}		
		
		omxCopyMatrix(filteredDataRow, dataRow);
		omxRemoveRowsAndColumns(filteredDataRow, zeros, toRemove);

		omxRecompute(rowAlgebra, fc);

		omxCopyMatrixToRow(rowAlgebra, omxDataIndex(data, row), rowResults);
	}
	free(toRemove);
	free(zeros);
}

static void omxCallRowFitFunction(omxFitFunction *oo, int want, FitContext *fc) {
	if (want & (FF_COMPUTE_PREOPTIMIZE)) return;

    if(OMX_DEBUG) { mxLog("Beginning Row Evaluation.");}
	// Requires: Data, means, covariances.

	omxMatrix* objMatrix  = oo->matrix;
	int numChildren = fc->childList.size();

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
		omxResizeMatrix(rowResults, data->rows, rowAlgebra->cols);
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
			FitContext *kid = fc->childList[i];
			omxMatrix *childMatrix = kid->lookupDuplicate(objMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				omxRowFitFunctionSingleIteration(childFit, oo, stride * i, data->rows - stride * i, fc);
			} else {
				omxRowFitFunctionSingleIteration(childFit, oo, stride * i, stride, fc);
			}
		}
	} else {
		omxRowFitFunctionSingleIteration(oo, oo, 0, data->rows, fc);
	}

	omxRecompute(reduceAlgebra, fc);

	omxCopyMatrix(oo->matrix, reduceAlgebra);

}

void omxInitRowFitFunction(omxFitFunction* oo) {

	if(OMX_DEBUG) { mxLog("Initializing Row/Reduce fit function."); }

	SEXP rObj = oo->rObj;
	SEXP nextMatrix, nextItem;
	int numDeps;

	omxRowFitFunction *newObj = new omxRowFitFunction;

	if(OMX_DEBUG) {mxLog("Accessing data source."); }
	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("data")));
	newObj->data = omxDataLookupFromState(nextMatrix, oo->matrix->currentState);
	if(newObj->data == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No data provided to omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}
	}

	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("rowAlgebra")));
	newObj->rowAlgebra = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	if(newObj->rowAlgebra == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row-wise algebra in omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}
	}

	{
		ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("units")));
		oo->setUnitsFromName(CHAR(STRING_ELT(nextMatrix, 0)));
	}

	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("filteredDataRow")));
	newObj->filteredDataRow = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	}
	if(newObj->filteredDataRow == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}
	// Create the original data row from which to filter.
	newObj->dataRow = omxInitMatrix(newObj->filteredDataRow->rows,
					newObj->filteredDataRow->cols, TRUE, oo->matrix->currentState);
	omxCopyMatrix(newObj->filteredDataRow, newObj->dataRow);

	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("existenceVector")));
	newObj->existenceVector = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	}
    // Do we allow NULL existence?  (Whoa, man. That's, like, deep, or something.)
	if(newObj->existenceVector == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No existance matrix in omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}


	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("rowResults")));
	newObj->rowResults = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	}
	if(newObj->rowResults == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row results matrix in omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}

	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("reduceAlgebra")));
	newObj->reduceAlgebra = omxMatrixLookupFromState1(nextMatrix, oo->matrix->currentState);
	}
	if(newObj->reduceAlgebra == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "No row reduction algebra in omxRowFitFunction.");
		omxRaiseError(errstr);
		free(errstr);
	}
	
	if(OMX_DEBUG) {mxLog("Accessing variable mapping structure."); }
	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("dataColumns")));
	newObj->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, oo->matrix->currentState, 0, 0);
	}
	if(OMX_DEBUG) { omxPrint(newObj->dataColumns, "Variable mapping"); }

	if(OMX_DEBUG) {mxLog("Accessing data row dependencies."); }
	{ ScopedProtect p1(nextItem, R_do_slot(rObj, Rf_install("dataRowDeps")));
	numDeps = LENGTH(nextItem);
	newObj->numDataRowDeps = numDeps;
	newObj->dataRowDeps = (int*) R_alloc(numDeps, sizeof(int));
	for(int i = 0; i < numDeps; i++) {
		newObj->dataRowDeps[i] = INTEGER(nextItem)[i];
	}
	}

	/* Set up data columns */
	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);

	oo->computeFun = omxCallRowFitFunction;
	oo->destructFun = omxDestroyRowFitFunction;

	oo->argStruct = (void*) newObj;
}


