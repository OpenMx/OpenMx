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
#include "omxOpenmpWrap.h"
#include "matrix.h"

// forward declarations
static omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber);

const char omxMatrixMajorityList[3] = "Tn";		// BLAS Column Majority.

void omxPrintMatrix(omxMatrix *source, const char* header)
{
	std::string buf;
	buf += string_snprintf("[%d] %s: (%d x %d) [%s-major] SEXP %p\n%s = matrix(c(",
			       omx_absolute_thread_num(),
			       header, source->rows, source->cols, (source->colMajor?"col":"row"),
			       source->owner, header);

	int first=TRUE;
	if(source->colMajor) {
		for(int j = 0; j < source->rows; j++) {
			buf += "\n";
			for(int k = 0; k < source->cols; k++) {
				if (first) first=FALSE;
				else buf += ",";
				buf += string_snprintf(" %3.6f", source->data[k*source->rows+j]);
			}
		}
	} else {
		for(int j = 0; j < source->cols; j++) {
			buf += "\n";
			for(int k = 0; k < source->rows; k++) {
				if (first) first=FALSE;
				else buf += ",";
				buf += string_snprintf(" %3.6f", source->data[k*source->cols+j]);
			}
		}
	}
	buf += string_snprintf("), byrow=TRUE, nrow=%d, ncol=%d)\n", source->rows, source->cols);
	mxLogBig(buf);
}

omxMatrix* omxInitMatrix(omxMatrix* om, int nrows, int ncols, unsigned short isColMajor, omxState* os) {

	if(om == NULL) om = (omxMatrix*) R_alloc(1, sizeof(omxMatrix));
	if(OMX_DEBUG_MATRIX) { mxLog("Initializing matrix %p to (%d, %d) with state at %p.", om, nrows, ncols, os); }

	om->hasMatrixNumber = 0;
	om->rows = nrows;
	om->cols = ncols;
	om->colMajor = (isColMajor ? 1 : 0);

	om->originalRows = om->rows;
	om->originalCols = om->cols;
	om->originalColMajor=om->colMajor;

	om->owner = NULL;
	if(om->rows == 0 || om->cols == 0) {
		om->data = NULL;
	} else {
		om->data = (double*) Calloc(nrows * ncols, double);
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
	om->version = 1;
	omxMarkClean(om);

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

	if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) {
		const char *oname = "?";
		if (orig->name) oname = orig->name;
		const char *dname = "?";
		if (dest->name) dname = dest->name;
		mxLog("omxCopyMatrix from %s (%p) to %s (%p)", oname, orig, dname, dest);
	}

	int regenerateMemory = TRUE;
	int numPopLocs = orig->numPopulateLocations;

	if(!dest->owner && (dest->originalRows == orig->rows && dest->originalCols == orig->cols)) {
		regenerateMemory = FALSE;				// If it's local data and the right size, we can keep memory.
	}

	dest->rows = orig->rows;
	dest->cols = orig->cols;
	dest->colMajor = orig->colMajor;
	dest->originalRows = dest->rows;
	dest->originalCols = dest->cols;
	dest->originalColMajor = dest->colMajor;

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
	} else {
		if(regenerateMemory) {
			omxFreeMatrixData(dest);											// Free and regenerate memory
			dest->data = (double*) Calloc(dest->rows * dest->cols, double);
		}
		if (dest->data != orig->data) {  // if equal then programmer error? TODO
			memcpy(dest->data, orig->data, dest->rows * dest->cols * sizeof(double));
		}
	}

	dest->aliasedPtr = NULL;

	omxMatrixLeadingLagging(dest);
}

void omxAliasMatrix(omxMatrix *dest, omxMatrix *src) {
	omxCopyMatrix(dest, src);
	dest->aliasedPtr = src;					// Alias now follows back matrix precisely.
}

void omxFreeMatrixData(omxMatrix * om) {

	if(!om->owner && om->data != NULL) {
		if(OMX_DEBUG_MATRIX) { mxLog("Freeing matrix data at %p", om->data); }
		Free(om->data);
	}
	om->owner = NULL;
	om->data = NULL;
}

void omxFreeAllMatrixData(omxMatrix *om) {
    
    if(om == NULL) return;

	if(OMX_DEBUG) { 
	    mxLog("Freeing matrix at %p with data = %p, algebra %p, and fit function %p.", 
		  om, om->data, om->algebra, om->fitFunction);
	}

	omxFreeMatrixData(om);

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
	newMat->name = src->name;
    
    return newMat;    
}

void omxResizeMatrix(omxMatrix *om, int nrows, int ncols, unsigned short keepMemory) {
	// Always Recompute() before you Resize().
	if(OMX_DEBUG_MATRIX) { 
		mxLog("Resizing matrix from (%d, %d) to (%d, %d) (keepMemory: %d)", 
			om->rows, om->cols, 
			nrows, ncols, keepMemory);
	}
	if((keepMemory == FALSE) && (om->rows != nrows || om->cols != ncols)) {
		if(OMX_DEBUG_MATRIX) { mxLog(" and regenerating memory to do it"); }
		omxFreeMatrixData(om);
		om->data = (double*) Calloc(nrows * ncols, double);
	} else if(om->originalRows * om->originalCols < nrows * ncols) {
		warning("Upsizing an existing matrix may cause undefined behavior.\n"); // TODO: Define this behavior?
	}

	if(OMX_DEBUG_MATRIX) { mxLog("."); }
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
	char *errstr = (char*) calloc(250, sizeof(char));
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
	char *errstr = (char*) calloc(250, sizeof(char));
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
	error("Requested improper value (%d, %d) from (%d, %d) matrix", row, col, numrow, numcol);
}

void setVectorError(int index, int numrow, int numcol) {
	char *errstr = (char*) calloc(250, sizeof(char));
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
		char *errstr = (char*) calloc(250, sizeof(char));
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

void omxMarkDirty(omxMatrix *om) { om->version += 1; }
void omxMarkClean(omxMatrix *om) { om->version += 1; om->cleanVersion = om->version; }

omxMatrix* omxNewMatrixFromRPrimitive(SEXP rObject, omxState* state, 
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Creates and populates an omxMatrix with details from an R matrix object. */
	omxMatrix *om = NULL;
	om = omxInitMatrix(NULL, 0, 0, FALSE, state);
	return omxFillMatrixFromRPrimitive(om, rObject, state, hasMatrixNumber, matrixNumber);
}

omxMatrix* omxFillMatrixFromRPrimitive(omxMatrix* om, SEXP rObject, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Populates the fields of a omxMatrix with details from an R object. */
	if(!isMatrix(rObject) && !isVector(rObject)) { // Sanity Check
		error("Recieved unknown matrix type in omxFillMatrixFromRPrimitive.");
	}
	return(fillMatrixHelperFunction(om, rObject, state, hasMatrixNumber, matrixNumber));
}



static omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {

	int* dimList;

	if(OMX_DEBUG) { mxLog("Filling omxMatrix from R matrix."); }

	if(om == NULL) {
		om = omxInitMatrix(NULL, 0, 0, FALSE, state);
	}

	if(isMatrix(matrix)) {
		SEXP matrixDims;
		PROTECT(matrixDims = getAttrib(matrix, R_DimSymbol));
		dimList = INTEGER(matrixDims);
		om->rows = dimList[0];
		om->cols = dimList[1];
		UNPROTECT(1); // matrixDims
	} else if (isVector(matrix)) {		// If it's a vector, assume it's a row vector. BLAS doesn't care.
		if(OMX_DEBUG) { mxLog("Vector discovered.  Assuming rowity."); }
		om->rows = 1;
		om->cols = length(matrix);
	}
	if(OMX_DEBUG) { mxLog("Matrix connected to (%d, %d) matrix or MxMatrix.", om->rows, om->cols); }

	if (TYPEOF(matrix) != REALSXP) {
		// we could avoid a double copy here TODO
		SEXP copy;
		PROTECT(copy = coerceVector(matrix, REALSXP));
		om->data = (double*) Realloc(NULL, om->rows * om->cols, double);
		memcpy(om->data, REAL(copy), om->rows * om->cols * sizeof(double));
	} else {
		om->owner = matrix;
		om->data = REAL(om->owner);
	}

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
	om->version = 1;
	omxMarkClean(om);

	if(OMX_DEBUG) { mxLog("Pre-compute call.");}
	omxMatrixLeadingLagging(om);
	if(OMX_DEBUG) { mxLog("Post-compute call.");}

	if(OMX_DEBUG) {
		omxPrintMatrix(om, "Finished importing matrix");
	}

	return om;
}

void omxProcessMatrixPopulationList(omxMatrix* matrix, SEXP matStruct) {

	if(OMX_DEBUG) { mxLog("Processing Population List: %d elements.", length(matStruct) - 1); }

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
		SEXP subList;
		PROTECT(subList = AS_INTEGER(VECTOR_ELT(matStruct, i+1)));

		int* locations = INTEGER(subList);
		if(OMX_DEBUG) { mxLog("."); } //:::
		matrix->populateFrom[i] = locations[0];
		matrix->populateFromRow[i] = locations[1];
		matrix->populateFromCol[i] = locations[2];
		matrix->populateToRow[i] = locations[3];
		matrix->populateToCol[i] = locations[4];
		UNPROTECT(1); //subList
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

	omxFreeMatrixData(mat);
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
//	if(OMX_DEBUG_MATRIX) { mxLog("Removing %d rows and %d columns from %p.", numRowsRemoved, numColsRemoved, om);}

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
		if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("Handling column %d/%d...", j, oldCols);}
		if(colsRemoved[j]) {
			if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("Removed.");}
			continue;
		} else {
			nextRow = 0;
			if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("Rows (max %d): ", oldRows); }
			for(int k = 0; k < oldRows; k++) {
				if(rowsRemoved[k]) {
					if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("%d removed....", k);}
					continue;
				} else {
					if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("%d kept....", k);}
					if(om->aliasedPtr == NULL) {
						if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("Self-aliased matrix access.");}
						omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j));
					} else {
						if(OMX_DEBUG_MATRIX || OMX_DEBUG_ALGEBRA) { mxLog("Matrix %p re-aliasing to %p.", om, om->aliasedPtr);}
						omxSetMatrixElement(om, nextRow, nextCol, omxMatrixElement(om->aliasedPtr, k,  j));
					}
					nextRow++;
				}
			}
			nextCol++;
		}
	}

	omxMatrixLeadingLagging(om);
}

/* Function wrappers that switch based on inclusion of algebras */
void omxPrint(omxMatrix *source, const char* d) { 					// Pretty-print a (small) matrix
    if(source == NULL) mxLog("%s is NULL.", d);
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
	bool yes;
	if (matrix->hasMatrixNumber && omxMatrixIsClean(matrix)) {
		yes = FALSE;
	} else {
		yes = TRUE;
	}
	const char *name = "?";
	if (matrix->name) name = matrix->name;
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("Matrix %s (%p) %s update", name, matrix, yes? "needs" : "does not need");
	}
	return yes;
}

static void maybeCompute(omxMatrix *matrix, int want)
{
	if(matrix->numPopulateLocations > 0) omxPopulateSubstitutions(matrix);
	else if(!omxNeedsUpdate(matrix)) /* do nothing */;
	else if(matrix->algebra != NULL) omxAlgebraRecompute(matrix->algebra);
	else if(matrix->fitFunction != NULL) {
		omxFitFunctionCompute(matrix->fitFunction, want, NULL);
	}
}

void omxRecompute(omxMatrix *matrix)
{
	maybeCompute(matrix, FF_COMPUTE_FIT);
}

void omxInitialCompute(omxMatrix *matrix)
{
	maybeCompute(matrix, 0);
}

void omxForceCompute(omxMatrix *matrix) {
	if(matrix->numPopulateLocations > 0) omxPopulateSubstitutions(matrix);
	else if (matrix->algebra != NULL) omxAlgebraForceCompute(matrix->algebra);
	else if(matrix->fitFunction != NULL) {
		omxFitFunctionCompute(matrix->fitFunction, FF_COMPUTE_FIT, NULL);
	}
}

/*
 * omxShallowInverse
 * 			Calculates the inverse of (I-A) using an n-step Neumann series
 * Assumes that A reduces to all zeros after numIters steps
 *
 * params:
 * omxMatrix *A				: The A matrix.  I-A will be inverted.  Size MxM.
 * omxMatrix *Z				: On output: Computed (I-A)^-1. MxM.
 * omxMatrix *Ax			: Space for computation. MxM.
 * omxMatrix *I				: Identity matrix. Will not be changed on exit. MxM.
 */

void omxShallowInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I ) {

	omxMatrix* origZ = Z;
    double oned = 1, minusOned = -1.0;

	if(numIters == NA_INTEGER) {
		if(OMX_DEBUG_ALGEBRA) { mxLog("RAM Algebra (I-A) inversion using standard (general) inversion."); }

		/* Z = (I-A)^-1 */
		if(I->colMajor != A->colMajor) {
			omxTransposeMatrix(I);  // transpose I? Hm? TODO
		}
		omxCopyMatrix(Z, A);

		omxDGEMM(FALSE, FALSE, oned, I, I, minusOned, Z);

		Matrix Zmat(Z);
		int info = MatrixInvert1(Zmat);
		if (info) {
			omxRaiseErrorf(A->currentState, "(I-A) is exactly singular (info=%d)", info);
		        return;
		}

		if(OMX_DEBUG_ALGEBRA) {omxPrint(Z, "Z");}

	} else {

		if(OMX_DEBUG_ALGEBRA) { mxLog("RAM Algebra (I-A) inversion using optimized expansion."); }

		/* Taylor Expansion optimized I-A calculation */
		if(I->colMajor != A->colMajor) {
			omxTransposeMatrix(I);
		}

		if(I->colMajor != Ax->colMajor) {
			omxTransposeMatrix(Ax);
		}

		omxCopyMatrix(Z, A);

		/* Optimized I-A inversion: Z = (I-A)^-1 */
		// F77_CALL(omxunsafedgemm)(I->majority, A->majority, &(I->cols), &(I->rows), &(A->rows), &oned, I->data, &(I->cols), I->data, &(I->cols), &oned, Z->data, &(Z->cols));  // Z = I + A = A^0 + A^1
		// omxDGEMM(FALSE, FALSE, 1.0, I, I, 1.0, Z); // Z == A + I

		for(int i = 0; i < A->rows; i++)
			omxSetMatrixElement(Z, i, i, 1);

		for(int i = 1; i <= numIters; i++) { // TODO: Efficiently determine how many times to do this
			// The sequence goes like this: (I + A), I + (I + A) * A, I + (I + (I + A) * A) * A, ...
			// Which means only one DGEMM per iteration.
			if(OMX_DEBUG_ALGEBRA) { mxLog("....RAM: Iteration #%d/%d", i, numIters);}
			omxCopyMatrix(Ax, I);
			// F77_CALL(omxunsafedgemm)(A->majority, A->majority, &(Z->cols), &(Z->rows), &(A->rows), &oned, Z->data, &(Z->cols), A->data, &(A->cols), &oned, Ax->data, &(Ax->cols));  // Ax = Z %*% A + I
			omxDGEMM(FALSE, FALSE, oned, A, Z, oned, Ax);
			omxMatrix* m = Z; Z = Ax; Ax = m;	// Juggle to make Z equal to Ax
		}
		if(origZ != Z) { 	// Juggling has caused Ax and Z to swap
			omxCopyMatrix(Z, Ax);
		}
	}
}

double omxMaxAbsDiff(omxMatrix *m1, omxMatrix *m2)
{
	if (m1->rows != m2->rows || m1->cols != m2->cols) error("Matrices are not the same size");

	double mad = 0;
	int size = m1->rows * m1->cols;
	for (int dx=0; dx < size; ++dx) {
		double mad1 = fabs(m1->data[dx] - m2->data[dx]);
		if (mad < mad1) mad = mad1;
	}
	return mad;
}
