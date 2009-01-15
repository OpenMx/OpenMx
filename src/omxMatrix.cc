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

void omxMatrix::init() {
	init(0,0,TRUE);
}

omxMatrix::omxMatrix() {
	init();
}

omxMatrix::omxMatrix(int nrows, int ncols, bool isColMajor) {
	init(nrows, ncols, isColMajor);
}

void omxMatrix::init(int nrows, int ncols, bool isColMajor) {
	if(OMX_DEBUG) { Rprintf("Initializing 0x%0x to (%d, %d).\n", this, nrows, ncols); }
	rows = nrows;
	cols = ncols;
	colMajor = (isColMajor?1:0);
	
	leading = (colMajor?rows:cols);
	
	originalRows = rows;
	originalCols = cols;
	originalColMajor=colMajor;
	
	if(rows == 0 || cols == 0) {
		data = NULL; 
		localData = FALSE; 
	} else {
		data = (double*) Calloc(nrows * ncols, double);
		localData = TRUE;
	}

	aliasedPtr = NULL;
	isAlgebra = false;
	
	compute();
	
}

omxMatrix::omxMatrix(const omxMatrix &in) {
	if(OMX_DEBUG) { Rprintf("omxMatrix::Copy Constructor Called, this is 0x%0x, duplicating 0x%0x.\n", this, &in);}
 
	rows = in.rows;
	cols = in.cols;
	colMajor = in.colMajor;
	leading = in.leading;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = colMajor;
	
	if(rows ==0 || cols == 0) {
		data = NULL;
		localData = FALSE;
	} else {
		data = (double*) Calloc(rows * cols, double);
		memcpy(data, in.data, rows * cols * sizeof(double));
		localData = TRUE;
	}
	
	aliasedPtr = NULL;
	compute();
}

void omxMatrix::operator=(omxMatrix &orig) {
	if(OMX_DEBUG) { Rprintf("Operator = : this is 0x%0x, duplicating 0x%0x.\n", this, &orig); }
	freeData();

	rows = orig.rows;
	cols = orig.cols;
	colMajor = orig.colMajor;
	leading = orig.leading;
	originalRows = rows;
	originalCols = cols;
	originalColMajor = colMajor;

	if(rows == 0 || cols == 0) {
		data = NULL;
		localData=FALSE;
	} else {
		data = (double*) Calloc(rows * cols, double);
		memcpy(data, orig.data, rows * cols * sizeof(double));
		localData = TRUE;
	}

	aliasedPtr = NULL;
	
	compute();
}

void omxMatrix::alias(omxMatrix &dM) {
	*this = dM;
	aliasedPtr = dM.data;
	isAlgebra = false;
}

void omxMatrix::freeData() {
 
	if(localData && data != NULL) {
		if(OMX_DEBUG) { Rprintf("Freeing 0x%0x. Localdata = %d.\n", data, localData); }
		Free(data);
		localData = FALSE;
	}

}

omxMatrix::~omxMatrix() {
	freeData();
}

void omxMatrix::resize(int nrows, int ncols, bool keepMemory) {
	if(OMX_DEBUG) { Rprintf("Resizing matrix from (%d, %d) to (%d, %d) (keepMemory: %d)", rows, cols, nrows, ncols, keepMemory); }
	if(keepMemory == false) { 
		if(OMX_DEBUG) { Rprintf(" and regenerating memory to do it"); }
		freeData();
		data = (double*) Calloc(nrows * ncols, double);
		localData = TRUE;
	} else if(originalRows * originalCols < nrows * ncols) {
		error("Cannot yet keep memory smaller than target while resizing omxMatrix.\n"); // TODO: Allow expansion using memcopy?
	}

	if(OMX_DEBUG) { Rprintf(".\n"); }
	rows = nrows;
	cols = ncols;
	leading = (colMajor?rows:cols);
	recompute();
}

void omxMatrix::reset() {
	rows = originalRows;
	cols = originalCols;
	colMajor = originalColMajor;
	leading = (colMajor?rows:cols);
	if(aliasedPtr != NULL) {
		if(OMX_DEBUG) { print("I was"); for(int i = 0; i < rows*cols; i++) Rprintf("%3.5f ", aliasedPtr[i]); Rprintf("\n");}
		memcpy(data, aliasedPtr, rows*cols*sizeof(double));
		if(OMX_DEBUG) { print("I am");}
	}
}

void omxMatrix::recompute() {
 	if(isDirty) omxMatrix::compute(); 
}

void omxMatrix::compute() {
	if(OMX_DEBUG) { Rprintf("Matrix compute: 0x%0x (needed: %d).\n", 1,1); }
	majority = &(majorityList[colMajor]);
	majority = &(majorityList[(colMajor?0:1)]);
	isDirty = false;
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
	leading = (colMajor?rows:cols);
	originalRows = rows;
	originalCols = cols;
	originalColMajor = TRUE;
	isAlgebra=false;

	UNPROTECT(4);
	
	compute();
	
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
	leading = (colMajor?rows:cols);
	aliasedPtr = data;
	isAlgebra = false;
	
	if(OMX_DEBUG) { Rprintf("Pre-compute call.\n");}
	compute();
	if(OMX_DEBUG) { Rprintf("Post-compute call.\n");}
	
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
	leading = (colMajor?rows:cols);			// Could possibly do this better keeping leading the same. TODO: Improve.
	
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

const char omxMatrix::majorityList[3] = "Tn";
