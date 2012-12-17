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
	om->colMajor = (isColMajor ? 1 : 0);

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
	om->fitFunction = NULL;

	om->currentState = os;
	om->isTemporary = FALSE;
	om->name = NULL;
	om->isDirty = FALSE;

	omxMatrixLeadingLagging(om);

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
	dest->isDirty = orig->isDirty;

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

	omxMatrixLeadingLagging(dest);
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
	    Rprintf("Freeing matrix at 0x%0x with data = 0x%x, algebra 0x%x, and fit function 0x%x.\n", 
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

	if(om->fitFunction != NULL) {
		omxFreeFitFunctionArgs(om->fitFunction);
		om->fitFunction = NULL;
	}
	
	if(om->isTemporary) {
		Free(om);
		om = NULL;
	}

}

/**
 * Copies an omxMatrix to a new R matrix object
 *
 * \param om the omxMatrix to copy
 * \return a PROTECT'd SEXP for the R matrix object
 */
SEXP omxExportMatrix(omxMatrix *om) {
	SEXP nextMat;
	PROTECT(nextMat = allocMatrix(REALSXP, om->rows, om->cols));
	for(int row = 0; row < om->rows; row++) {
		for(int col = 0; col < om->cols; col++) {
			REAL(nextMat)[col * om->rows + row] =
				omxMatrixElement(om, row, col);
		}
	}
	return nextMat;
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

	omxMatrixLeadingLagging(om);
}

void omxResetAliasedMatrix(omxMatrix *om) {
	om->rows = om->originalRows;
	om->cols = om->originalCols;
	if(om->aliasedPtr != NULL) {
		memcpy(om->data, om->aliasedPtr->data, om->rows*om->cols*sizeof(double));
		om->colMajor = om->aliasedPtr->colMajor;
	}
	omxMatrixLeadingLagging(om);
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
	free(errstr);  // TODO not reached
}

void setMatrixError(omxMatrix *om, int row, int col, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	static const char *matrixString = "matrix";
	static const char *algebraString = "algebra";
	static const char *fitString = "fit function";
	const char *typeString;
	if (om->algebra != NULL) {
		typeString = algebraString;
	} else if (om->fitFunction != NULL) {
		typeString = fitString;
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
	free(errstr);  // TODO not reached
}

void matrixElementError(int row, int col, int numrow, int numcol) {
	char *errstr = calloc(250, sizeof(char));
	sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.",
		row, col, numrow, numcol);
	error(errstr);
	free(errstr);  // TODO not reached
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
	free(errstr);  // TODO not reached
}

double omxAliasedMatrixElement(omxMatrix *om, int row, int col) {
	int index = 0;
	if(row >= om->originalRows || col >= om->originalCols) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.", 
			row + 1, col + 1, om->originalRows, om->originalCols);
		error(errstr);
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

void omxMarkDirty(omxMatrix *om) { om->isDirty = TRUE; }
void omxMarkClean(omxMatrix *om) { om->isDirty = FALSE; }

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
	om->fitFunction = NULL;
	om->currentState = state;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;
	om->isDirty = FALSE;

	if(OMX_DEBUG) { Rprintf("Pre-compute call.\n");}
	omxMatrixLeadingLagging(om);
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
	
	omxMatrixLeadingLagging(mat);
}

void omxRemoveElements(omxMatrix *om, int numRemoved, int removed[]) {

	if(numRemoved < 1) { return; }

	int oldElements;

	if (om->rows > 1) {
		if(om->aliasedPtr == NULL) {
			if(om->originalRows == 0) {
				om->originalRows = om->rows;
			}
			oldElements = om->originalRows;
		} else {
			oldElements = om->aliasedPtr->rows;
		}
		om->rows = oldElements - numRemoved;
	} else {
		if(om->aliasedPtr == NULL) {
			if(om->originalCols == 0) {
				om->originalCols = om->cols;
			}
			oldElements = om->originalCols;
		} else {
			oldElements = om->aliasedPtr->cols;
		}
		om->cols = oldElements - numRemoved;
	}

	int nextElement = 0;

	for(int j = 0; j < oldElements; j++) {
		if(!removed[j]) {
			if(om->aliasedPtr == NULL) {
				omxUnsafeSetVectorElement(om, nextElement, omxUnsafeVectorElement(om, j));
			} else {
				omxUnsafeSetVectorElement(om, nextElement, omxUnsafeVectorElement(om->aliasedPtr, j));
			}
			nextElement++;
		}
	}

	omxMatrixLeadingLagging(om);
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

	omxMatrixLeadingLagging(om);
}

/* Function wrappers that switch based on inclusion of algebras */
void omxPrint(omxMatrix *source, char* d) { 					// Pretty-print a (small) matrix
    if(source == NULL) Rprintf("%s is NULL.\n", d);
	else if(source->algebra != NULL) omxAlgebraPrint(source->algebra, d);
	else if(source->fitFunction != NULL) omxFitFunctionPrint(source->fitFunction, d);
	else omxPrintMatrix(source, d);
}

void omxPopulateSubstitutions(omxMatrix *om) {
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
		}
	}
}

void omxMatrixLeadingLagging(omxMatrix *om) {
	om->majority = &(omxMatrixMajorityList[(om->colMajor?1:0)]);
	om->minority = &(omxMatrixMajorityList[(om->colMajor?0:1)]);
	om->leading = (om->colMajor?om->rows:om->cols);
	om->lagging = (om->colMajor?om->cols:om->rows);
}

unsigned short omxNeedsUpdate(omxMatrix *matrix) {
	if (OMX_DEBUG_ALGEBRA) 
		Rprintf("Matrix 0x%0x needs update: ", matrix);
	if (matrix->hasMatrixNumber && !matrix->isDirty) {
		if (OMX_DEBUG_ALGEBRA) Rprintf("No\n");
		return(FALSE);
	} else {
		if (OMX_DEBUG_ALGEBRA) Rprintf("Yes\n");
		return(TRUE);
	}
}

void omxRecompute(omxMatrix *matrix) {
	if(matrix->numPopulateLocations > 0) omxPopulateSubstitutions(matrix);
	else if(!omxNeedsUpdate(matrix)) /* do nothing */;
	else if(matrix->algebra != NULL) omxAlgebraRecompute(matrix->algebra);
	else if(matrix->fitFunction != NULL) omxFitFunctionCompute(matrix->fitFunction);
}

void omxCompute(omxMatrix *matrix) {
	if(matrix->numPopulateLocations > 0) omxPopulateSubstitutions(matrix);
	else if (matrix->algebra != NULL) omxAlgebraCompute(matrix->algebra);
	else if(matrix->fitFunction != NULL) omxFitFunctionCompute(matrix->fitFunction);
	else if(matrix->numPopulateLocations > 0) omxPopulateSubstitutions(matrix);
}
