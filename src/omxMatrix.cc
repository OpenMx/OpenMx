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

void omxMatrix::print(char* header) {
	int j, k;
	
	Rprintf("%s: (%d x %d)\n", header, rows, cols);
	
	if(colMajor) {
		for(j = 0; j < rows; j++) {
			for(k = 0; k < cols; k++) {
				Rprintf("\t%3.6f", data[k*rows+j]);
			}
			Rprintf("\n");
		}
	} else {
		for(j = 0; j < cols; j++) {
			for(k = 0; k < rows; k++) {
				Rprintf("\t%3.6f", data[k*cols+j]);
			}
			Rprintf("\n");
		}
	}
}

omxMatrix::omxMatrix() {
	omxMatrix(0,0,FALSE);
}

omxMatrix::omxMatrix(int nrows, int ncols, bool isColMajor) {
	rows = nrows;
	cols = ncols;
	colMajor = (isColMajor?1:0);
	
	originalRows = rows;
	originalCols = cols;
	originalColMajor=colMajor;
	
	localData = TRUE;
	isReused = FALSE;
	majorityList[0] = 'n';
	majorityList[1] = 'T';
	
	data = (double*) Calloc(nrows * ncols, double);
	aliasedPtr = NULL;
	
	recompute();
}

omxMatrix::omxMatrix(const omxMatrix &in) {
	if(OMX_DEBUG) { Rprintf("omxMatrix::Copy Constructor Called, this is %d.\n", this);}
 
	rows = in.rows;
	cols = in.cols;
	colMajor = in.colMajor;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = colMajor;
	
	data = (double*) Calloc(rows * cols, double);
	memcpy(data, in.data, rows * cols * sizeof(double));
	localData = TRUE;
	aliasedPtr = NULL;
	recompute();
}

void omxMatrix::operator=(omxMatrix orig) {
	if(OMX_DEBUG) { Rprintf("Operator =\n"); }
	freeData();

	rows = orig.rows;
	cols = orig.cols;
	colMajor = orig.colMajor;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = colMajor;
	localData = TRUE;
	
	data = (double*) Calloc(rows * cols, double);
	memcpy(data, orig.data, rows * cols * sizeof(double));
	
	aliasedPtr = NULL;
	
	recompute();
}

void omxMatrix::alias(omxMatrix &dM) {
	*this = dM;
	aliasedPtr = dM.data;
}

void omxMatrix::freeData() { 
	if(localData) {
		Free(data);
	}
	localData = FALSE;

}

omxMatrix::~omxMatrix() {
	freeData();
}

void omxMatrix::resize(int nrows, int ncols, bool keepMemory) {
	if(!keepMemory) { 
		freeData();
		data = (double*) Calloc(nrows * ncols, double);
		localData = TRUE;
	} 
	
	rows = nrows;
	cols = ncols;
	recompute();
}

void omxMatrix::reset() {
	rows = originalRows;
	cols = originalCols;
	colMajor = originalColMajor;
	if(aliasedPtr != NULL) {
		if(OMX_DEBUG) { print("I am"); for(int i = 0; i < rows*cols; i++) Rprintf("%3.5f ", aliasedPtr[i]); Rprintf("\n");}
		memcpy(data, aliasedPtr, rows*cols*sizeof(double));
		if(OMX_DEBUG) { print("I am");}
	}
}

void omxMatrix::recompute() {
	if(colMajor) {
		leading = rows;
		lagging = cols;
	} else {
		leading = cols;
		lagging = rows;
	}
	majority = majorityList[colMajor];
}

double* omxMatrix::locationOfElement(int row, int col) {
	int index = 0;
	if(colMajor) {
		index = col * rows + row;
	} else {
		index = row * cols + col;
	}
	return data + index;
}

double omxMatrix::element(int row, int col) {
	int index = 0;
	if(colMajor) {
		index = col * rows + row;
	} else {
		index = row * cols + col;
	}
	return data[index];
}

void omxMatrix::setElement(int row, int col, double value) {
	int index = 0;
	if(colMajor) {
		index = col * rows + row;
	} else {
		index = row * cols + col;
	}
	data[index] = value;
}

void omxMatrix::fillFromS3Matrix(SEXP mxMatrix) {
/* Populates the fields of a omxMatrix with details from an mxMatrix. */ 

	Rprintf("fillDataMatrixFromS3Matrix() Should never be called.\n");

	SEXP objectEnv, valueMatrix, nameString, matrixDims;
	int* dimList;
	double* ans; 

	PROTECT(nameString = getAttrib(mxMatrix, install("class")));
	if(STRING_ELT(nameString,0) == R_NilValue) error("DataMatrix: Not an MxMatrix object.\n");  // TODO: Error System?
	if(strncmp(CHAR(STRING_ELT(nameString,0)), "MxMatrix", 8) != 0) error("DataMatrix: Not an MxMatrix object");
	
    SET_STRING_ELT(nameString, 0, mkChar(".env"));				
	PROTECT(objectEnv = getAttrib(mxMatrix, nameString));							// May duplicate.  TODO: Check.

	SET_STRING_ELT(nameString, 0, mkChar(".values"));
	PROTECT(valueMatrix = findVar(install(CHAR(STRING_ELT(nameString, 0))), objectEnv)); // May duplicate.  TODO: Check.
	data = REAL(valueMatrix);														// Probably duplicates.  TODO: Need to fix.

	PROTECT(matrixDims = getAttrib(valueMatrix, R_DimSymbol));
	dimList = INTEGER(matrixDims);
	rows = dimList[0];
	cols = dimList[1];
	localData = FALSE;
	colMajor = TRUE;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = TRUE;
	isReused = TRUE;

	UNPROTECT(4);
	
	recompute();
	
	return;
}

void omxMatrix::fillFromMatrix(SEXP matrix) {
/* Populates the fields of a omxMatrix with details from an R Matrix. */ 
	
	SEXP className;
	SEXP matrixDims;
	int* dimList;
	
	if(OMX_DEBUG) { Rprintf("Filling omxMatrix from R matrix.\n"); }
	
	/* Sanity Check */
	if(!isMatrix(matrix) && !isVector(matrix)) {
		SEXP values;
		int *rowVec, *colVec;
		double  *dataVec;
		const char *stringName;
		PROTECT(className = getAttrib(matrix, install("class")));
		if(strncmp(CHAR(STRING_ELT(className, 0)), "Symm", 2) != 0) { // Should be "Mx"
			error("omxMatrix::fillFromMatrix--Passed Non-vector, non-matrix SEXP.\n");
		}
		stringName = CHAR(STRING_ELT(className, 0));
		if(strncmp(stringName, "SymmMatrix", 12) == 0) {
			PROTECT(values = GET_SLOT(matrix, install("values")));
			rows = INTEGER(GET_SLOT(values, install("nrow")))[0];
			cols = INTEGER(GET_SLOT(values, install("ncol")))[0];
			
			data = (double*) S_alloc(rows * cols, sizeof(double));		// We can afford to keep through the whole call
			rowVec = INTEGER(GET_SLOT(values, install("rowVector")));
			colVec = INTEGER(GET_SLOT(values, install("colVector")));
			dataVec = REAL(GET_SLOT(values, install("dataVector")));
			for(int j = 0; j < length(GET_SLOT(values, install("dataVector"))); j++) {
				data[(rowVec[j]-1) + (colVec[j]-1) * rows] = dataVec[j];
				data[(rowVec[j]-1) * cols + (colVec[j]-1)] = dataVec[j];  // Symmetric, after all.
			}
			UNPROTECT(1); // values
		}
		UNPROTECT(1); // className
	} else {
		data = REAL(matrix);	// TODO: Class-check first?
		
		if(isMatrix(matrix)) {
			PROTECT(matrixDims = getAttrib(matrix, R_DimSymbol));
			dimList = INTEGER(matrixDims);
			rows = dimList[0];
			cols = dimList[1];
			UNPROTECT(1);	// MatrixDims
		} else if (isVector(matrix)) {			// If it's a vector, assume it's a row vector.
			rows = 1;
			cols = length(matrix);
		}
	}	
	
	localData = FALSE;
	colMajor = TRUE;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = TRUE;
	aliasedPtr = data;
	isReused = TRUE;
	
	recompute();
	
	return;
}

void omxMatrix::removeRowsAndColumns(int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[])
{
	if(aliasedPtr == NULL) {  // This is meant only for aliased matrices.  Maybe Need a subclass?
		error("removeRowsAndColumns intended only for aliased matrices.\n");
	}
	
	
	if(numRowsRemoved < 1 || numColsRemoved < 1) { return; }
		
	int numCols = 0;
	int nextCol = 0;
	int nextRow = 0;
	int oldRows = rows;
	int oldCols = cols;
	int j,k;
	
	rows = rows - numRowsRemoved;
	cols = cols - numColsRemoved;
	
	// Note:  This really aught to be done using a matrix multiply.  Why isn't it?
	if(colMajor) {
		for(int j = 0; j < oldCols; j++) {
			if(colsRemoved[j]) {
				continue;
			} else {
				nextRow = 0;
				for(int k = 0; k <= oldRows; k++) {
					if(rowsRemoved[k]) {
						continue;
					} else {
						setElement(nextRow, nextCol, aliasedPtr[k + j * oldRows]);
						nextRow++;
					}
				}
				nextCol++;
			}
		}
	} else {
		for(int j = 0; j < oldRows; j++) {
			if(rowsRemoved[j]) {
				continue;
			} else {
				nextCol = 0;
				for(int k = 0; k <= oldCols; k++) {
					if(colsRemoved[k]) {
						continue;
					} else {
						setElement(nextRow, nextCol, aliasedPtr[k + j * oldCols]);
						nextCol++;
					}
				}
				nextRow++;
			}
		}
	}

	recompute();
}

//	omxMatrix* censoredMatrix;
//	censoredMatrix = new omxMatrix(rows - numRowsRemoved, cols - numColsRemoved);
//	int numCols;
//	int nextCol;
//	int j,k;
//	
//	// Note:  This really aught to be done using a matrix multiply.  Why isn't it?
//	numCols = 0;
//	nextCol = 0;
//	for(int j = 0; j < cols; j++) {
//		if(rowsRemoved[j]) {
//			continue;
//		} else {
//			nextRow = 0;
//			for(int k = 0; k <= j; k++) {
//				if(colsRemoved[k]) {
//					continue;
//				} else {
//					censoredMatrix->data[nextRow + nextCol*numCols] = data[k + j * expected->cols()];
//					nextRow++;
//				}
//			}
//			nextCol++;
//		}
//	}
//	if(nextRow != numCols || nextCol != numCols) error("Something failed: Matrices are non-working.");
//
//	return censoredDataMatrix;
