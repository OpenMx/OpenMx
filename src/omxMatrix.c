/*
 *  Copyright 2007-2012 The OpenMx Project
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

/***********************************************************
*
*  omxMatrix.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Contains code for the omxMatrix class
*   omxDataMatrices hold necessary information to simplify
* 	dealings between the OpenMX back end and BLAS.
*
**********************************************************/
#include "omxMatrix.h"

// forward declarations
omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber);

const char omxMatrixMajorityList[3] = "Tn";		// BLAS Column Majority.

void omxPrintMatrix(omxMatrix *source, char* header) {
	int j, k;

	Rprintf("%s: (%d x %d) [%s-major]\n", header, source->rows, source->cols, (source->colMajor?"col":"row"));
	if(OMX_DEBUG_MATRIX) {Rprintf("Matrix Printing is at %0x\n", source);}

	if(source->colMajor) {
		for(j = 0; j < source->rows; j++) {
			for(k = 0; k < source->cols; k++) {
				Rprintf("\t%3.6f", source->data[k*source->rows+j]);
			}
			Rprintf("\n");
		}
	} else {
		for(j = 0; j < source->cols; j++) {
			for(k = 0; k < source->rows; k++) {
				Rprintf("\t%3.6f", source->data[k*source->cols+j]);
			}
			Rprintf("\n");
		}
	}
}

omxMatrix* omxInitMatrix(omxMatrix* om, int nrows, int ncols, unsigned short isColMajor, omxState* os) {

	if(om == NULL) om = (omxMatrix*) R_alloc(1, sizeof(omxMatrix));
	if(OMX_DEBUG_MATRIX) { Rprintf("Initializing matrix 0x%0x to (%d, %d) with state at 0x%x.\n", om, nrows, ncols, os); }

	om->hasMatrixNumber = 0;
	om->rows = nrows;
	om->cols = ncols;
	om->colMajor = (isColMajor?1:0);

	om->originalRows = om->rows;
	om->originalCols = om->cols;
	om->originalColMajor=om->colMajor;

	if(om->rows == 0 || om->cols == 0) {
		om->data = NULL;
		om->localData = FALSE;
	} else {
		om->data = (double*) Calloc(nrows * ncols, double);
		om->localData = TRUE;
	}

	om->populateFrom = NULL;
	om->populateFromCol = NULL;
	om->populateFromRow = NULL;
	om->populateToCol = NULL;
	om->populateToRow = NULL;

	om->numPopulateLocations = 0;

	om->aliasedPtr = NULL;
	om->algebra = NULL;
	om->objective = NULL;

	om->currentState = os;
	om->lastCompute = -2;
	om->lastRow = -2;
	om->isTemporary = FALSE;
	om->name = NULL;

	omxMatrixCompute(om);

	return om;

}

omxMatrix* omxInitTemporaryMatrix(omxMatrix* om, int nrows, int ncols, unsigned short isColMajor, omxState* os) {

	if(om == NULL) {
		om = (omxMatrix*) Calloc(1, omxMatrix);
	}

	om = omxInitMatrix(om, nrows, ncols, isColMajor, os);
	om->isTemporary = TRUE;
	
	return(om);

}

void omxUpdateMatrixHelper(omxMatrix *dest, omxMatrix *orig) {

	dest->rows = orig->rows;
	dest->cols = orig->cols;
	dest->colMajor = orig->colMajor;
	dest->originalRows = dest->rows;
	dest->originalCols = dest->cols;
	dest->originalColMajor = dest->colMajor;
	dest->lastCompute = orig->lastCompute;
	dest->lastRow = orig->lastRow;
	dest->localData = orig->localData;

	if ((dest->data == NULL) && (dest->rows > 0) && (dest->cols > 0)) {
		dest->data = (double*) Calloc(dest->rows * dest->cols, double);
	}

	memcpy(dest->data, orig->data, dest->rows * dest->cols * sizeof(double));

	if (dest->aliasedPtr != NULL && orig->aliasedPtr != NULL) {
		omxUpdateMatrix(dest->aliasedPtr, orig->aliasedPtr);
	}

	omxMatrixCompute(dest);
}

void omxUpdateMatrix(omxMatrix *dest, omxMatrix *orig) {
	omxUpdateMatrixHelper(dest, orig);

	if (orig->algebra) {
		omxUpdateAlgebra(dest->algebra, orig->algebra);
	} else if (orig->objective) {
		omxUpdateObjectiveFunction(dest->objective, orig->objective);
	}

}

void omxCopyMatrix(omxMatrix *dest, omxMatrix *orig) {
	/* Copy a matrix.  NOTE: Matrix maintains its algebra bindings. */

	if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("omxCopyMatrix"); }

	int regenerateMemory = TRUE;
	int numPopLocs = orig->numPopulateLocations;

	if(dest->localData && (dest->originalRows == orig->rows && dest->originalCols == orig->cols)) {
		regenerateMemory = FALSE;				// If it's local data and the right size, we can keep memory.
	}

	dest->rows = orig->rows;
	dest->cols = orig->cols;
	dest->colMajor = orig->colMajor;
	dest->originalRows = dest->rows;
	dest->originalCols = dest->cols;
	dest->originalColMajor = dest->colMajor;
	dest->lastCompute = orig->lastCompute;
	dest->lastRow = orig->lastRow;

	dest->numPopulateLocations = numPopLocs;
	if (numPopLocs > 0) {
		dest->populateFrom = (int*)R_alloc(numPopLocs, sizeof(int));
		dest->populateFromRow = (int*)R_alloc(numPopLocs, sizeof(int));
		dest->populateFromCol = (int*)R_alloc(numPopLocs, sizeof(int));
		dest->populateToRow = (int*)R_alloc(numPopLocs, sizeof(int));
		dest->populateToCol = (int*)R_alloc(numPopLocs, sizeof(int));
		
		memcpy(dest->populateFrom, orig->populateFrom, numPopLocs * sizeof(int));
		memcpy(dest->populateFromRow, orig->populateFromRow, numPopLocs * sizeof(int));
		memcpy(dest->populateFromCol, orig->populateFromCol, numPopLocs * sizeof(int));
		memcpy(dest->populateToRow, orig->populateToRow, numPopLocs * sizeof(int));
		memcpy(dest->populateToCol, orig->populateToCol, numPopLocs * sizeof(int));
	}

	if(dest->rows == 0 || dest->cols == 0) {
		omxFreeMatrixData(dest);
		dest->data = NULL;
		dest->localData=FALSE;
	} else {
		if(regenerateMemory) {
			omxFreeMatrixData(dest);											// Free and regenerate memory
			dest->data = (double*) Calloc(dest->rows * dest->cols, double);
		}
		memcpy(dest->data, orig->data, dest->rows * dest->cols * sizeof(double));
		dest->localData = TRUE;
	}

	dest->aliasedPtr = NULL;

	omxMatrixCompute(dest);
}

void omxAliasMatrix(omxMatrix *dest, omxMatrix *src) {
	omxCopyMatrix(dest, src);
	dest->aliasedPtr = src;					// Alias now follows back matrix precisely.
}

void omxFreeMatrixData(omxMatrix * om) {

	if(om->localData && om->data != NULL) {
		if(OMX_DEBUG_MATRIX) { Rprintf("Freeing matrix at 0x%0x. Localdata = %d.\n", om->data, om->localData); }
		Free(om->data);
		om->data = NULL;
		om->localData = FALSE;
	}

}

void omxFreeAllMatrixData(omxMatrix *om) {
    
    if(om == NULL) return;

	if(OMX_DEBUG) { 
	    Rprintf("Freeing matrix at 0x%0x with data = 0x%x, algebra 0x%x, and objective 0x%x.\n", 
	        om, om->data, om->algebra); 
	}

	if(om->localData && om->data != NULL) {
		Free(om->data);
		om->data = NULL;
		om->localData = FALSE;
	}

	if(om->algebra != NULL) {
		omxFreeAlgebraArgs(om->algebra);
		om->algebra = NULL;
	}

	if(om->objective != NULL) {
		omxFreeObjectiveArgs(om->objective);
		om->objective = NULL;
	}
	
	if(om->isTemporary) {
		Free(om);
		om = NULL;
	}

}

void omxZeroByZeroMatrix(omxMatrix *om) {
	if (om->rows > 0 || om->cols > 0) {
		omxResizeMatrix(om, 0, 0, FALSE);
	}
}

omxMatrix* omxNewIdentityMatrix(int nrows, omxState* state) {
	omxMatrix* newMat = NULL;
	int l,k;

	newMat = omxInitMatrix(newMat, nrows, nrows, FALSE, state);
	for(k =0; k < newMat->rows; k++) {
		for(l = 0; l < newMat->cols; l++) {
			if(l == k) {
				omxSetMatrixElement(newMat, k, l, 1);
			} else {
				omxSetMatrixElement(newMat, k, l, 0);
			}
		}
	}
	return newMat;
}

omxMatrix* omxDuplicateMatrix(omxMatrix* src, omxState* newState) {
	omxMatrix* newMat;
    
	if(src == NULL) return NULL;
	newMat = omxInitMatrix(NULL, src->rows, src->cols, FALSE, newState);
	omxCopyMatrix(newMat, src);
	newMat->hasMatrixNumber = src->hasMatrixNumber;
	newMat->matrixNumber    = src->matrixNumber;
    
    return newMat;    
}

void omxResizeMatrix(omxMatrix *om, int nrows, int ncols, unsigned short keepMemory) {
	// Always Recompute() before you Resize().
	if(OMX_DEBUG_MATRIX) { 
		Rprintf("Resizing matrix from (%d, %d) to (%d, %d) (keepMemory: %d)", 
			om->rows, om->cols, 
			nrows, ncols, keepMemory);
	}
	if((keepMemory == FALSE) && (om->rows != nrows || om->cols != ncols)) {
		if(OMX_DEBUG_MATRIX) { Rprintf(" and regenerating memory to do it"); }
		omxFreeMatrixData(om);
		om->data = (double*) Calloc(nrows * ncols, double);
		om->localData = TRUE;
	} else if(om->originalRows * om->originalCols < nrows * ncols) {
		warning("Upsizing an existing matrix may cause undefined behavior.\n"); // TODO: Define this behavior?
	}

	if(OMX_DEBUG_MATRIX) { Rprintf(".\n"); }
	om->rows = nrows;
	om->cols = ncols;
	if(keepMemory == FALSE) {
		om->originalRows = om->rows;
		om->originalCols = om->cols;
	}

	omxMatrixCompute(om);
}

void omxResetAliasedMatrix(omxMatrix *om) {
	om->rows = om->originalRows;
	om->cols = om->originalCols;
	if(om->aliasedPtr != NULL) {
		memcpy(om->data, om->aliasedPtr->data, om->rows*om->cols*sizeof(double));
		om->colMajor = om->aliasedPtr->colMajor;
	}
	omxMatrixCompute(om);
}

void omxMatrixCompute(omxMatrix *om) {

	if(OMX_DEBUG_MATRIX) { Rprintf("Matrix compute: 0x%0x, 0x%0x, algebra: 0x%x.\n", om, om->currentState, om->algebra); }
	om->majority = &(omxMatrixMajorityList[(om->colMajor?1:0)]);
	om->minority = &(omxMatrixMajorityList[(om->colMajor?0:1)]);
	om->leading = (om->colMajor?om->rows:om->cols);
	om->lagging = (om->colMajor?om->cols:om->rows);

	for(int i = 0; i < om->numPopulateLocations; i++) {
		int index = om->populateFrom[i];
		omxMatrix* sourceMatrix;
		if (index < 0) {
			sourceMatrix = om->currentState->matrixList[~index];
		} else {
			sourceMatrix = om->currentState->algebraList[index];
		}
		if (sourceMatrix != NULL) {
			omxRecompute(sourceMatrix);				// Make sure it's up to date
			double value = omxMatrixElement(sourceMatrix, om->populateFromRow[i], om->populateFromCol[i]);
			omxSetMatrixElement(om, om->populateToRow[i], om->populateToCol[i], value);
			// And then fill in the details.  Use the Helper here in case of transposition/downsampling.
		}
	}

	om->isDirty = FALSE;
	om->lastCompute = om->currentState->computeCount;
	om->lastRow = om->currentState->currentRow;

}

double* omxLocationOfMatrixElement(omxMatrix *om, int row, int col) {
	int index = 0;
	if(om->colMajor) {
		index = col * om->rows + row;
	} else {
		index = row * om->cols + col;
	}
	return om->data + index;
}

void vectorElementError(int index, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	if ((numrow > 1) && (numcol > 1)) {
		sprintf(errstr, "Requested improper index (%d) from a malformed vector of dimensions (%d, %d).", 
			index, numrow, numcol);
	} else {
		int length = (numrow > 1) ? numrow : numcol;
		sprintf(errstr, "Requested improper index (%d) from vector of length (%d).", 
			index, length);
	}
	error(errstr);
	free(errstr);
}

void setMatrixError(omxMatrix *om, int row, int col, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	static const char *matrixString = "matrix";
	static const char *algebraString = "algebra";
	static const char *objectiveString = "objective";
	const char *typeString;
	if (om->algebra != NULL) {
		typeString = algebraString;
	} else if (om->objective != NULL) {
		typeString = objectiveString;
	} else {
		typeString = matrixString;
	}
	if (om->name == NULL) {
		sprintf(errstr, "Attempted to set row and column (%d, %d) in %s with dimensions %d x %d.", 
			row, col, typeString, numrow, numcol);
	} else {
		sprintf(errstr, "Attempted to set row and column (%d, %d) in %s \"%s\" with dimensions %d x %d.", 
			row, col, typeString, om->name, numrow, numcol);
	}
	error(errstr);
	free(errstr);
}

void matrixElementError(int row, int col, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.",
		row, col, numrow, numcol);
	error(errstr);
	free(errstr);
}

void setVectorError(int index, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	if ((numrow > 1) && (numcol > 1)) {
		sprintf(errstr, "Attempting to set improper index (%d) from a malformed vector of dimensions (%d, %d).", 
			index, numrow, numcol);
	} else {
		int length = (numrow > 1) ? numrow : numcol;
		sprintf(errstr, "Setting improper index (%d) from vector of length %d.", 
			index, length);
	}
	error(errstr);
	free(errstr);
}

double omxAliasedMatrixElement(omxMatrix *om, int row, int col) {
	int index = 0;
	if(row >= om->originalRows || col >= om->originalCols) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.", row+1, col+1, om->originalRows, om->originalCols);
		error(errstr);
		free(errstr);
        return (NA_REAL);
	}
	if(om->colMajor) {
		index = col * om->originalRows + row;
	} else {
		index = row * om->originalCols + col;
	}
	return om->data[index];
}

void omxMarkDirty(omxMatrix *om) { om->isDirty = TRUE; }

unsigned short omxMatrixNeedsUpdate(omxMatrix *om) {
	for(int i = 0; i < om->numPopulateLocations; i++) {
		int index = om->populateFrom[i];
		omxMatrix* sourceMatrix;
		if (index < 0) {
			sourceMatrix = om->currentState->matrixList[~index];
		} else {
			sourceMatrix = om->currentState->algebraList[index];
		}
		if(omxNeedsUpdate(sourceMatrix)) return TRUE;	// Make sure it's up to date
	}
    return FALSE;
};

omxMatrix* omxNewMatrixFromMxMatrix(SEXP mxMatrix, omxState* state) {
/* Creates and populates an omxMatrix with details from an R MxMatrix object. */
	omxMatrix *om = NULL;
	om = omxInitMatrix(NULL, 0, 0, FALSE, state);
	return omxFillMatrixFromMxMatrix(om, mxMatrix, state);
}

omxMatrix* omxNewMatrixFromRPrimitive(SEXP rObject, omxState* state, 
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Creates and populates an omxMatrix with details from an R matrix object. */
	omxMatrix *om = NULL;
	om = omxInitMatrix(NULL, 0, 0, FALSE, state);
	return omxFillMatrixFromRPrimitive(om, rObject, state, hasMatrixNumber, matrixNumber);
}

omxMatrix* omxFillMatrixFromMxMatrix(omxMatrix* om, SEXP mxMatrix, omxState* state) {
/* Populates the fields of a omxMatrix with details from an R Matrix. */
	if(inherits(mxMatrix, "MxMatrix")) {
		if(OMX_DEBUG) { Rprintf("R matrix is Mx Matrix.  Processing.\n"); }
		SEXP matrix;
		PROTECT(matrix = GET_SLOT(mxMatrix,  install("values")));
		om = fillMatrixHelperFunction(om, matrix, state, 0, 0);
		UNPROTECT(1);
	} else {
		error("Recieved unknown matrix type in omxFillMatrixFromMxMatrix.");
	}
	return(om);
}

omxMatrix* omxFillMatrixFromRPrimitive(omxMatrix* om, SEXP rObject, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Populates the fields of a omxMatrix with details from an R object. */
	if(!isMatrix(rObject) && !isVector(rObject)) { // Sanity Check
		error("Recieved unknown matrix type in omxFillMatrixFromRPrimitive.");
	}
	return(fillMatrixHelperFunction(om, rObject, state, hasMatrixNumber, matrixNumber));
}



omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {

	SEXP matrixDims;
	int* dimList;

	if(OMX_DEBUG) { Rprintf("Filling omxMatrix from R matrix.\n"); }

	if(om == NULL) {
		om = omxInitMatrix(NULL, 0, 0, FALSE, state);
	}

	om->data = REAL(AS_NUMERIC(matrix));	// TODO: Class-check first?

	if(isMatrix(matrix)) {
		PROTECT(matrixDims = getAttrib(matrix, R_DimSymbol));
		dimList = INTEGER(matrixDims);
		om->rows = dimList[0];
		om->cols = dimList[1];
		UNPROTECT(1);	// MatrixDims
	} else if (isVector(matrix)) {		// If it's a vector, assume it's a row vector. BLAS doesn't care.
		if(OMX_DEBUG) { Rprintf("Vector discovered.  Assuming rowity.\n"); }
		om->rows = 1;
		om->cols = length(matrix);
	}
	if(OMX_DEBUG) { Rprintf("Matrix connected to (%d, %d) matrix or MxMatrix.\n", om->rows, om->cols); }

	om->localData = FALSE;
	om->colMajor = TRUE;
	om->originalRows = om->rows;
	om->originalCols = om->cols;
	om->originalColMajor = TRUE;
	om->aliasedPtr = NULL;
	om->algebra = NULL;
	om->objective = NULL;
	om->currentState = state;
	om->lastCompute = -1;
	om->lastRow = -1;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;


	if(OMX_DEBUG) { Rprintf("Pre-compute call.\n");}
	omxMatrixCompute(om);
	if(OMX_DEBUG) { Rprintf("Post-compute call.\n");}

	if(OMX_DEBUG) {
		omxPrintMatrix(om, "Finished importing matrix");
	}

	return om;
}

void omxProcessMatrixPopulationList(omxMatrix* matrix, SEXP matStruct) {

	if(OMX_DEBUG) { Rprintf("Processing Population List: %d elements.\n", length(matStruct) - 1); }
	SEXP subList;

	if(length(matStruct) > 1) {
		int numPopLocs = length(matStruct) - 1;
		matrix->numPopulateLocations = numPopLocs;
		matrix->populateFrom = (int*)R_alloc(numPopLocs, sizeof(int));
		matrix->populateFromRow = (int*)R_alloc(numPopLocs, sizeof(int));
		matrix->populateFromCol = (int*)R_alloc(numPopLocs, sizeof(int));
		matrix->populateToRow = (int*)R_alloc(numPopLocs, sizeof(int));
		matrix->populateToCol = (int*)R_alloc(numPopLocs, sizeof(int));
	}

	for(int i = 0; i < length(matStruct)-1; i++) {
		PROTECT(subList = AS_INTEGER(VECTOR_ELT(matStruct, i+1)));

		int* locations = INTEGER(subList);
		if(OMX_DEBUG) { Rprintf("."); } //:::
		matrix->populateFrom[i] = locations[0];
		matrix->populateFromRow[i] = locations[1];
		matrix->populateFromCol[i] = locations[2];
		matrix->populateToRow[i] = locations[3];
		matrix->populateToCol[i] = locations[4];

		UNPROTECT(1); // subList
	}
}

void omxToggleRowColumnMajor(omxMatrix *mat) {

	int i, j;
	int nrows = mat->rows;
	int ncols = mat->cols;
	
	double *newdata = (double*) Calloc(nrows * ncols, double);
	double *olddata = mat->data;

	if (mat->colMajor) {
		for(i = 0; i < ncols; i++) {
			for(j = 0; j < nrows; j++) {
				newdata[i + ncols * j] = olddata[i * nrows + j];
			}
		}
	} else {
		for(i = 0; i < nrows; i++) {
			for(j = 0; j < ncols; j++) {
				newdata[i + nrows * j] = olddata[i * ncols + j];
			}
		}
	}

	if (mat->localData) {
		Free(mat->data);
	}

	mat->localData = TRUE;
	mat->data = newdata;
	mat->colMajor = !mat->colMajor;
}

void omxTransposeMatrix(omxMatrix *mat) {
	mat->colMajor = !mat->colMajor;
	
	if(mat->rows != mat->cols){
        int mid = mat->rows;
        mat->rows = mat->cols;
        mat->cols = mid;
	}
	
	omxMatrixCompute(mat);
}

void omxRemoveRowsAndColumns(omxMatrix *om, int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[])
{
    // TODO: Create short-circuit form of omxRemoveRowsAndCols to remove just rows or just columns.
//	if(OMX_DEBUG_MATRIX) { Rprintf("Removing %d rows and %d columns from 0x%0x.\n", numRowsRemoved, numColsRemoved, om);}

	if(numRowsRemoved < 1 && numColsRemoved < 1) { return; }

	int oldRows, oldCols;

	if(om->aliasedPtr == NULL) {
		if(om->originalRows == 0 || om->originalCols == 0) {
			om->originalRows = om->rows;
			om->originalCols = om->cols;
		}
		oldRows = om->originalRows;
		oldCols = om->originalCols;
	} else {
		oldRows = om->aliasedPtr->rows;
		oldCols = om->aliasedPtr->cols;
	}

	int nextCol = 0;
	int nextRow = 0;

	if(om->rows > om->originalRows || om->cols > om->originalCols) {	// sanity check.
		error("Aliased Matrix is too small for alias.");
	}

	om->rows = oldRows - numRowsRemoved;
	om->cols = oldCols - numColsRemoved;

	// Note:  This really aught to be done using a matrix multiply.  Why isn't it?
	for(int j = 0; j < oldCols; j++) {
		if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("Handling column %d/%d...", j, oldCols);}
		if(colsRemoved[j]) {
			if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("Removed.\n");}
			continue;
		} else {
			nextRow = 0;
			if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("Rows (max %d): ", oldRows); }
			for(int k = 0; k < oldRows; k++) {
				if(rowsRemoved[k]) {
					if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("%d removed....", k);}
					continue;
				} else {
					if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("%d kept....", k);}
					if(om->aliasedPtr == NULL) {
						if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("Self-aliased matrix access.\n");}
						omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j));
					} else {
						if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("Matrix 0x%x re-aliasing to 0x%x.\n", om, om->aliasedPtr);}
						omxSetMatrixElement(om, nextRow, nextCol, omxMatrixElement(om->aliasedPtr, k,  j));
					}
					nextRow++;
				}
			}
			if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { Rprintf("\n");}
			nextCol++;
		}
	}

	omxMatrixCompute(om);
}

/* Function wrappers that switch based on inclusion of algebras */
void omxPrint(omxMatrix *source, char* d) { 					// Pretty-print a (small) matrix
	if(source->algebra != NULL) omxAlgebraPrint(source->algebra, d);
	else if(source->objective != NULL) omxObjectivePrint(source->objective, d);
	else omxPrintMatrix(source, d);
}

unsigned short omxNeedsUpdate(omxMatrix *matrix) {
	unsigned short retval;
	/* Simplest update check: If we're dirty or haven't computed this cycle (iteration or row), we need to. */
	// TODO : Implement a dependency-tree-based dirtiness propagation system
	if(OMX_DEBUG_MATRIX) {Rprintf("Matrix 0x%x NeedsUpdate?", matrix);}

	if(matrix == NULL) {
		if(OMX_DEBUG_MATRIX) {Rprintf("matrix argument is NULL. ");}
		retval = FALSE;		// Not existing means never having to say you need to recompute.
	} else if(matrix->isDirty) {
		if(OMX_DEBUG_MATRIX) {Rprintf("matrix is dirty. ");}
		retval = TRUE;
	} else if(matrix->lastCompute < matrix->currentState->computeCount) {
		if(OMX_DEBUG_MATRIX) {Rprintf("matrix last compute is less than current compute count. ");}
		retval = TRUE;  	// No need to check args if oa's dirty.
	} else if(matrix->lastRow != matrix->currentState->currentRow) {
		if(OMX_DEBUG_MATRIX) {Rprintf("matrix last row is less than current row. ");}
		retval = TRUE;			// Ditto.
	} else if(matrix->algebra != NULL) {
		if(OMX_DEBUG_MATRIX) {Rprintf("checking algebra needs update. ");}
		retval = omxAlgebraNeedsUpdate(matrix->algebra);
	} else if(matrix->objective != NULL) {
		if(OMX_DEBUG_MATRIX) {Rprintf("checking objective function needs update. ");}
		retval = omxObjectiveNeedsUpdate(matrix->objective);
	} else {
		if(OMX_DEBUG_MATRIX) {Rprintf("checking matrix needs update. ");}
		retval = omxMatrixNeedsUpdate(matrix);
	}
	if(OMX_DEBUG_MATRIX && retval) {Rprintf("Yes.\n");}
	if(OMX_DEBUG_MATRIX && !retval) {Rprintf("No.\n");}
	return(retval);
}

void omxRecompute(omxMatrix *matrix) {
	if(!omxNeedsUpdate(matrix)) return;
	if(matrix->algebra != NULL) omxAlgebraCompute(matrix->algebra);
	else if(matrix->objective != NULL) omxObjectiveCompute(matrix->objective);
	else omxMatrixCompute(matrix);
}

void omxCompute(omxMatrix *matrix) {
	if(matrix->algebra != NULL) omxAlgebraCompute(matrix->algebra);
	else if(matrix->objective != NULL) omxObjectiveCompute(matrix->objective);
	else omxMatrixCompute(matrix);
}
