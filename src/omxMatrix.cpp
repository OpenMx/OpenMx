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
#include "matrix.h"
#include "unsupported/Eigen/MatrixFunctions"
#include "omxState.h"

// forward declarations
static omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber);

static const char *omxMatrixMajorityList[] = {"T", "n"};		// BLAS Column Majority.

// For background, see
// http://epubs.siam.org/doi/abs/10.1137/090768539

void logm_eigen(int n, double *rz, double *out)
{
    Eigen::Map< Eigen::MatrixXd > inMat(rz, n, n);
    Eigen::Map< Eigen::MatrixXd > outMat(out, n, n);
    outMat = inMat.log();
}

void expm_eigen(int n, double *rz, double *out)
{
    Eigen::Map< Eigen::MatrixXd > inMat(rz, n, n);
    Eigen::Map< Eigen::MatrixXd > outMat(out, n, n);
    outMat = inMat.exp();
}

void omxPrintMatrix(omxMatrix *source, const char* header) // make static TODO
{
	EigenMatrixAdaptor Esrc(source);
	if (!header) header = source->name();
	if (!header) header = "?";
	mxPrintMat(header, Esrc);
}

omxMatrix* omxInitMatrix(int nrows, int ncols, unsigned short isColMajor, omxState* os) {

	if (!isColMajor) Rf_error("All matrices are created column major");

	omxMatrix* om = new omxMatrix;

	om->hasMatrixNumber = 0;
	om->rows = nrows;
	om->cols = ncols;
	om->colMajor = (isColMajor ? 1 : 0);

	om->originalRows = om->rows;
	om->originalCols = om->cols;

	om->owner = NULL;
	if(om->rows == 0 || om->cols == 0) {
		om->data = NULL;
	} else {
		om->data = (double*) Calloc(nrows * ncols, double);
	}

	om->algebra = NULL;
	om->fitFunction = NULL;

	om->currentState = os;
	om->nameStr = "?";
	om->version = 1;
	om->cleanVersion = 0;

	omxMatrixLeadingLagging(om);

	return om;

}

static void omxFreeInternalMatrixData(omxMatrix * om)
{
	if(!om->owner && om->data != NULL) {
		Free(om->data);
	}
	om->owner = NULL;
	om->data = NULL;
}

void omxCopyMatrix(omxMatrix *dest, omxMatrix *orig) {
	/* Copy a matrix.  NOTE: Matrix maintains its algebra bindings. */

	int regenerateMemory = TRUE;

	if(!dest->owner && (dest->originalRows == orig->rows && dest->originalCols == orig->cols)) {
		regenerateMemory = FALSE;				// If it's local data and the right size, we can keep memory.
	}

	dest->rows = orig->rows;
	dest->cols = orig->cols;
	dest->colMajor = orig->colMajor;
	dest->originalRows = dest->rows;
	dest->originalCols = dest->cols;
	dest->populate = orig->populate;

	if(dest->rows == 0 || dest->cols == 0) {
		omxFreeInternalMatrixData(dest);
		dest->data = NULL;
	} else {
		if(regenerateMemory) {
			omxFreeInternalMatrixData(dest);											// Free and regenerate memory
			dest->data = (double*) Calloc(dest->rows * dest->cols, double);
		}
		if (dest->data != orig->data) {  // if equal then programmer Rf_error? TODO
			memcpy(dest->data, orig->data, dest->rows * dest->cols * sizeof(double));
		}
	}

	omxMatrixLeadingLagging(dest);
}

void omxFreeMatrix(omxMatrix *om) {
    
    if(om == NULL) return;

	omxFreeInternalMatrixData(om);

	if(om->algebra != NULL) {
		omxFreeAlgebraArgs(om->algebra);
		om->algebra = NULL;
	}

	if(om->fitFunction != NULL) {
		omxFreeFitFunctionArgs(om->fitFunction);
		om->fitFunction = NULL;
	}
	
	if (!om->hasMatrixNumber) delete om;
}

/**
 * Copies an omxMatrix to a new R matrix object
 *
 * \param om the omxMatrix to copy
 * \return a Rf_protect'd SEXP for the R matrix object
 */
SEXP omxExportMatrix(omxMatrix *om) {
	SEXP nextMat;
	Rf_protect(nextMat = Rf_allocMatrix(REALSXP, om->rows, om->cols));
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
		omxResizeMatrix(om, 0, 0);
	}
}

omxMatrix* omxNewIdentityMatrix(int nrows, omxState* state) {
	omxMatrix* newMat = NULL;
	int l,k;

	newMat = omxInitMatrix(nrows, nrows, TRUE, state);
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
	newMat = omxInitMatrix(src->rows, src->cols, TRUE, newState);
	omxCopyMatrix(newMat, src);
	newMat->hasMatrixNumber = src->hasMatrixNumber;
	newMat->matrixNumber    = src->matrixNumber;
	newMat->nameStr = src->nameStr;
    
	newMat->rownames = src->rownames;
	newMat->colnames = src->colnames;

    return newMat;    
}

void omxResizeMatrix(omxMatrix *om, int nrows, int ncols)
{
	// Always Recompute() before you Resize().
	if(OMX_DEBUG_MATRIX) { 
		mxLog("Resizing matrix from (%d, %d) to (%d, %d)",
			om->rows, om->cols, nrows, ncols);
	}

	if( (om->rows != nrows || om->cols != ncols)) {
		omxFreeInternalMatrixData(om);
		om->data = (double*) Calloc(nrows * ncols, double);
	}

	om->rows = nrows;
	om->cols = ncols;
	om->originalRows = om->rows;
	om->originalCols = om->cols;

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
		int Rf_length = (numrow > 1) ? numrow : numcol;
		sprintf(errstr, "Requested improper index (%d) from vector of Rf_length (%d).", 
			index, Rf_length);
	}
	Rf_error(errstr);
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
	sprintf(errstr, "Attempted to set row and column (%d, %d) in %s \"%s\" with dimensions %d x %d.", 
		row, col, typeString, om->name(), numrow, numcol);
	Rf_error(errstr);
	free(errstr);  // TODO not reached
}

void matrixElementError(int row, int col, int numrow, int numcol) {
	Rf_error("Requested improper value (%d, %d) from (%d, %d) matrix", row, col, numrow, numcol);
}

void setVectorError(int index, int numrow, int numcol) {
	char *errstr = (char*) calloc(250, sizeof(char));
	if ((numrow > 1) && (numcol > 1)) {
		sprintf(errstr, "Attempting to set improper index (%d) from a malformed vector of dimensions (%d, %d).", 
			index, numrow, numcol);
	} else {
		int Rf_length = (numrow > 1) ? numrow : numcol;
		sprintf(errstr, "Setting improper index (%d) from vector of Rf_length %d.", 
			index, Rf_length);
	}
	Rf_error(errstr);
	free(errstr);  // TODO not reached
}

static double omxAliasedMatrixElement(omxMatrix *om, int row, int col)
{
	int index = 0;
	if(row >= om->originalRows || col >= om->originalCols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.", 
			row + 1, col + 1, om->originalRows, om->originalCols);
		Rf_error(errstr);
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

void omxMarkDirty(omxMatrix *om) {
	om->version += 1;
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("Marking %s %s dirty", om->getType(), om->name());
	}
}

void omxMarkClean(omxMatrix *om)
{
	om->version += 1;
	om->cleanVersion = om->version;
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("Marking %s %s clean", om->getType(), om->name());
	}
}

omxMatrix* omxNewMatrixFromRPrimitive(SEXP rObject, omxState* state, 
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Creates and populates an omxMatrix with details from an R matrix object. */
	omxMatrix *om = NULL;
	om = omxInitMatrix(0, 0, TRUE, state);
	return omxFillMatrixFromRPrimitive(om, rObject, state, hasMatrixNumber, matrixNumber);
}

omxMatrix* omxFillMatrixFromRPrimitive(omxMatrix* om, SEXP rObject, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Populates the fields of a omxMatrix with details from an R object. */
	if(rObject && !Rf_isMatrix(rObject) && !Rf_isVector(rObject)) { // Sanity Check
		Rf_error("Recieved unknown matrix type in omxFillMatrixFromRPrimitive.");
	}
	return(fillMatrixHelperFunction(om, rObject, state, hasMatrixNumber, matrixNumber));
}



static omxMatrix* fillMatrixHelperFunction(omxMatrix* om, SEXP matrix, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber) {

	int* dimList;

	if(OMX_DEBUG) { mxLog("Filling omxMatrix from R matrix."); }

	if (!om) Rf_error("fillMatrixHelperFunction: matrix must be allocated already");

	if (matrix) {
		if(Rf_isMatrix(matrix)) {
			SEXP matrixDims;
			ScopedProtect p1(matrixDims, Rf_getAttrib(matrix, R_DimSymbol));
			dimList = INTEGER(matrixDims);
			om->rows = dimList[0];
			om->cols = dimList[1];
		} else if (Rf_isVector(matrix)) {		// If it's a vector, assume it's a row vector. BLAS doesn't care.
			if(OMX_DEBUG) { mxLog("Vector discovered.  Assuming rowity."); }
			om->rows = 1;
			om->cols = Rf_length(matrix);
		}
		if(OMX_DEBUG) { mxLog("Matrix connected to (%d, %d) matrix or MxMatrix.", om->rows, om->cols); }

		if (TYPEOF(matrix) != REALSXP) Rf_error("matrix is of type '%s'; only type double is accepted",
							Rf_type2char(TYPEOF(matrix)));

		om->owner = matrix;
		om->data = REAL(om->owner);

		SEXP dimnames;
		ScopedProtect pdn(dimnames, Rf_getAttrib(matrix, R_DimNamesSymbol));
		if (!Rf_isNull(dimnames) && Rf_length(dimnames) == 2) {
			SEXP names;
			{ScopedProtect p1(names, VECTOR_ELT(dimnames, 0));
			int nlen = Rf_length(names);
			om->rownames.resize(nlen);
			for (int nx=0; nx < nlen; ++nx) {
				om->rownames[nx] = CHAR(STRING_ELT(names, nx));
			}
			}

			{ScopedProtect p1(names, VECTOR_ELT(dimnames, 1));
				int nlen = Rf_length(names);
			om->colnames.resize(nlen);
			for (int nx=0; nx < nlen; ++nx) {
				om->colnames[nx] = CHAR(STRING_ELT(names, nx));
			}
			}
		}
	}

	om->colMajor = TRUE;
	om->originalRows = om->rows;
	om->originalCols = om->cols;
	om->algebra = NULL;
	om->fitFunction = NULL;
	om->currentState = state;
	om->hasMatrixNumber = hasMatrixNumber;
	om->matrixNumber = matrixNumber;
	om->version = 1;
	om->cleanVersion = 0;

	omxMatrixLeadingLagging(om);

	return om;
}

void omxMatrix::omxProcessMatrixPopulationList(SEXP matStruct)
{
	if(Rf_length(matStruct) <= 1) return;

	const int numPopLocs = Rf_length(matStruct) - 1;

	if(OMX_DEBUG) { mxLog("Processing Population List: %d elements.", numPopLocs); }

	setNotConstant();

	populate.resize(numPopLocs);

	for(int i = 0; i < numPopLocs; i++) {
		SEXP subList;
		ScopedProtect p1(subList, VECTOR_ELT(matStruct, i+1));

		int* locations = INTEGER(subList);
		populateLocation &pl = populate[i];
		pl.from = locations[0];
		pl.srcRow = locations[1];
		pl.srcCol = locations[2];
		pl.destRow = locations[3];
		pl.destCol = locations[4];
	}
}

void omxMatrix::setNotConstant()
{
	if (!owner) return;

	double *copy = (double*) Realloc(NULL, rows * cols, double);
	memcpy(copy, data, rows * cols * sizeof(double));
	data = copy;
	owner = NULL;
}

void omxEnsureColumnMajor(omxMatrix *mat)
{
	if (mat->colMajor) return;
	omxToggleRowColumnMajor(mat);
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

	omxFreeInternalMatrixData(mat);
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
		if(om->originalRows == 0) Rf_error("Not allocated");
		oldElements = om->originalRows;
		om->rows = oldElements - numRemoved;
	} else {
		if(om->originalCols == 0) Rf_error("Not allocated");
		oldElements = om->originalCols;
		om->cols = oldElements - numRemoved;
	}

	int nextElement = 0;

	for(int j = 0; j < oldElements; j++) {
		if(!removed[j]) {
			omxUnsafeSetVectorElement(om, nextElement, omxUnsafeVectorElement(om, j));
			nextElement++;
		}
	}

	omxMatrixLeadingLagging(om);
}

void omxRemoveRowsAndColumns(omxMatrix *om, int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[])
{
    // TODO: Create short-circuit form of omxRemoveRowsAndCols to remove just rows or just columns.

	if(numRowsRemoved < 1 && numColsRemoved < 1) { return; }

	int oldRows, oldCols;

	if(om->originalRows == 0 || om->originalCols == 0) Rf_error("Not allocated");
	if (om->rows != om->originalRows || om->cols != om->originalCols) {
		// Feasible, but the code is currently not robust to this case
		Rf_error("Can only omxRemoveRowsAndColumns once");
	}

	oldRows = om->originalRows;
	oldCols = om->originalCols;

	int nextCol = 0;
	int nextRow = 0;

	om->rows = oldRows - numRowsRemoved;
	om->cols = oldCols - numColsRemoved;

	if (om->colMajor) {
		for(int j = 0; j < oldCols; j++) {
			if(colsRemoved[j]) continue;
			nextRow = 0;
			for(int k = 0; k < oldRows; k++) {
				if(rowsRemoved[k]) continue;
				omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j));
				nextRow++;
			}
			nextCol++;
		}
	} else {
		for(int k = 0; k < oldRows; k++) {
			if(rowsRemoved[k]) continue;
			nextCol = 0;
			for(int j = 0; j < oldCols; j++) {
				if(colsRemoved[j]) continue;
				omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j));
				nextCol++;
			}
			nextRow++;
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

void omxMatrix::omxPopulateSubstitutions(int want, FitContext *fc)
{
	if (populate.size() == 0) return;
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("omxPopulateSubstitutions %s, %d locations", name(), (int) populate.size());
	}
	for (size_t pi = 0; pi < populate.size(); pi++) {
		populateLocation &pl = populate[pi];
		int index = pl.from;
		omxMatrix* sourceMatrix;
		if (index < 0) {
			sourceMatrix = currentState->matrixList[~index];
		} else {
			sourceMatrix = currentState->algebraList[index];
		}

		omxRecompute(sourceMatrix, fc);

		if (want & FF_COMPUTE_INITIAL_FIT) {
			if (pl.srcRow >= sourceMatrix->rows || pl.srcCol >= sourceMatrix->cols) {
				// may not be properly dimensioned yet
				continue;
			}
		}

		if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_FIT)) {
			double value = omxMatrixElement(sourceMatrix, pl.srcRow, pl.srcCol);
			omxSetMatrixElement(this, pl.destRow, pl.destCol, value);
			if (OMX_DEBUG_ALGEBRA) {
				mxLog("copying %.2f from %s[%d,%d] to %s[%d,%d]",
				      value, sourceMatrix->name(), pl.srcRow, pl.srcCol,
				      name(), pl.destRow, pl.destCol);
			}
		}
	}
}

void omxMatrixLeadingLagging(omxMatrix *om) {
	om->majority = omxMatrixMajorityList[(om->colMajor?1:0)];
	om->minority = omxMatrixMajorityList[(om->colMajor?0:1)];
	om->leading = (om->colMajor?om->rows:om->cols);
	om->lagging = (om->colMajor?om->cols:om->rows);
}

static bool omxNeedsUpdate(omxMatrix *matrix)
{
	bool yes;
	if (matrix->hasMatrixNumber && omxMatrixIsClean(matrix)) {
		yes = FALSE;
	} else {
		yes = TRUE;
	}
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("%s %s is %s", matrix->getType(), matrix->name(), yes? "dirty" : "clean");
	}
	return yes;
}

void omxRecompute(omxMatrix *matrix, FitContext *fc)
{
	int want = matrix->currentState->getWantStage();
	matrix->omxPopulateSubstitutions(want, fc); // could be an algebra!

	if(!omxNeedsUpdate(matrix)) /* do nothing */;
	else if(matrix->algebra != NULL) omxAlgebraRecompute(matrix->algebra, fc);
	else if(matrix->fitFunction != NULL) {
		omxFitFunctionCompute(matrix->fitFunction, want, fc);
	}

	if (want & FF_COMPUTE_FIT) {
		omxMarkClean(matrix);
	}
}

void omxMatrix::transposePopulate()
{
	// This is sub-optimal. If the rows & cols were stored as vectors
	// then we could just swap them.
	for (size_t px=0; px < populate.size(); ++px) populate[px].transpose();
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

void omxShallowInverse(FitContext *fc, int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I )
{
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
			Z->data[0] = nan("singular");
			if (fc) fc->recordIterationError("(I-A) is exactly singular (info=%d)", info);
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
	if (m1->rows != m2->rows || m1->cols != m2->cols) Rf_error("Matrices are not the same size");

	double mad = 0;
	int size = m1->rows * m1->cols;
	for (int dx=0; dx < size; ++dx) {
		double mad1 = fabs(m1->data[dx] - m2->data[dx]);
		if (mad < mad1) mad = mad1;
	}
	return mad;
}

void checkIncreasing(omxMatrix* om, int column, int count, FitContext *fc)
{
	double previous = - INFINITY;
	double current;
	int threshCrossCount = 0;
	if(count > om->rows) {
		count = om->rows;
	}
	for(int j = 0; j < count; j++ ) {
		current = omxMatrixElement(om, j, column);
		if(std::isnan(current) || current == NA_INTEGER) {
			continue;
		}
		if(j == 0) {
			continue;
		}
		previous = omxMatrixElement(om, j-1, column);
		if(std::isnan(previous) || previous == NA_INTEGER) {
			continue;
		}
		if(current <= previous) {
			threshCrossCount++;
		}
	}
	if(threshCrossCount > 0) {
		fc->recordIterationError("Found %d thresholds out of order in column %d.", threshCrossCount, column+1);
	}
}

void omxStandardizeCovMatrix(omxMatrix* cov, double* corList, double* weights) {
	// Maybe coerce this into an algebra or sequence of algebras?

	if(OMX_DEBUG) { mxLog("Standardizing matrix."); }

	int rows = cov->rows;

	for(int i = 0; i < rows; i++) {
		weights[i] = sqrt(omxMatrixElement(cov, i, i));
	}

	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < i; j++) {
			corList[((i*(i-1))/2) + j] = omxMatrixElement(cov, i, j) / (weights[i] * weights[j]);
		}
	}
}

void omxMatrixHorizCat(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	int totalRows = 0, totalCols = 0, currentCol=0;

	if(numArgs == 0) return;

	totalRows = matList[0]->rows;			// Assumed constant.  Assert this below.

	for(int j = 0; j < numArgs; j++) {
		if(totalRows != matList[j]->rows) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Non-conformable matrices in horizontal concatenation (cbind). First argument has %d rows, and argument #%d has %d rows.", totalRows, j + 1, matList[j]->rows);
			omxRaiseError(errstr);
			free(errstr);
			return;
		}
		totalCols += matList[j]->cols;
	}

	if(result->rows != totalRows || result->cols != totalCols) {
		if(OMX_DEBUG_ALGEBRA) { mxLog("ALGEBRA: HorizCat: resizing result.");}
		omxResizeMatrix(result, totalRows, totalCols);
	}

	int allArgumentsColMajor = result->colMajor;
	for(int j = 0; j < numArgs && allArgumentsColMajor; j++) {
		if (!matList[j]->colMajor) allArgumentsColMajor = 0;
	}

	if (allArgumentsColMajor) {
		int offset = 0;
		for(int j = 0; j < numArgs; j++) {	
			omxMatrix* current = matList[j];
			int size = current->rows * current->cols;
			memcpy(result->data + offset, current->data, size * sizeof(double));
			offset += size;
		}
	} else {
		for(int j = 0; j < numArgs; j++) {
			for(int k = 0; k < matList[j]->cols; k++) {
				for(int l = 0; l < totalRows; l++) {		// Gotta be a faster way to do this.
					omxSetMatrixElement(result, l, currentCol, omxMatrixElement(matList[j], l, k));
				}
				currentCol++;
			}
		}
	}
}

void omxMatrixVertCat(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	int totalRows = 0, totalCols = 0, bindRow=0;

	if(numArgs == 0) return;

	totalCols = matList[0]->cols;			// Assumed constant.  Assert this below.

	for(int j = 0; j < numArgs; j++) {
		if(totalCols != matList[j]->cols) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Non-conformable matrices in vertical concatenation (rbind). First argument has %d cols, and argument #%d has %d cols.", totalCols, j + 1, matList[j]->cols);
			omxRaiseError(errstr);
			free(errstr);
			return;
		}
		totalRows += matList[j]->rows;
	}

	if(result->rows != totalRows || result->cols != totalCols) {
		omxResizeMatrix(result, totalRows, totalCols);
	}

	int allArgumentsRowMajor = !result->colMajor;
	for(int j = 0; j < numArgs && allArgumentsRowMajor; j++) {
		if (matList[j]->colMajor) allArgumentsRowMajor = 0;
	}

	if (allArgumentsRowMajor) {
		int offset = 0;
		for(int j = 0; j < numArgs; j++) {	
			omxMatrix* current = matList[j];
			int size = current->rows * current->cols;	
			memcpy(result->data + offset, current->data, size * sizeof(double));
			offset += size;
		}
	} else {
		for(int j = 0; j < numArgs; j++) {
			for(int k = 0; k < matList[j]->rows; k++) {
				for(int l = 0; l < totalCols; l++) {		// Gotta be a faster way to do this.
					omxSetMatrixElement(result, bindRow, l, omxMatrixElement(matList[j], k, l));
				}
				bindRow++;
			}
		}
	}

}

void omxMatrixTrace(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != numArgs && result->cols != numArgs) {
		omxResizeMatrix(result, numArgs, 1);
	}

    for(int i = 0; i < numArgs; i++) {
    	double trace = 0.0;
    	omxMatrix* inMat = matList[i];
        double* values = inMat->data;
        int nrow  = inMat->rows;
        int ncol  = inMat->cols;

    	if(nrow != ncol) {
    		char *errstr = (char*) calloc(250, sizeof(char));
    		sprintf(errstr, "Non-square matrix in Trace().\n");
    		omxRaiseError(errstr);
    		free(errstr);
            return;
    	}

    	/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
        for(int j = 0; j < nrow; j++)
           trace += values[j * nrow + j];

    	omxSetVectorElement(result, i, trace);
    }
}

