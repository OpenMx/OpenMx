/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
#include "omxState.h"
#include <limits>
#include <Eigen/SVD>
#include <Eigen/LU>
#include "Compute.h"
#include "minicsv.h"
#include "EnableWarnings.h"

// forward declarations
static const char *omxMatrixMajorityList[] = {"T", "n"};		// BLAS Column Majority.

void logm_eigen(int n, double *rz, double *out);
void expm_eigen(int n, double *rz, double *out);

std::string stringifyDimnames(omxMatrix *source)
{
	std::string buf;
	auto &rownames = source->rownames;
	auto &colnames = source->colnames;
	if (rownames.size() || colnames.size()) {
		buf += "dimnames=list(";
		if (!rownames.size()) {
			buf += "NULL";
		} else {
			buf += "c(";
			for (auto &nn : rownames) buf += string_snprintf("'%s', ", nn);
			buf += ")";
		}
		buf += ", ";
		if (!colnames.size()) {
			buf += "NULL";
		} else {
			buf += "c(";
			for (auto &nn : colnames) buf += string_snprintf("'%s', ", nn);
			buf += ")";
		}
		buf += ")";
	}
	return buf;
}

void omxPrintMatrix(omxMatrix *source, const char* header) // make static TODO
{
	EigenMatrixAdaptor Esrc(source);
	if (!header) header = source->name();
	if (!header) header = "?";

	std::string buf;
	auto &rownames = source->rownames;
	auto &colnames = source->colnames;
	if (rownames.size() || colnames.size()) {
		buf += ", ";
		buf += stringifyDimnames(source);
	}

	mxPrintMatX(header, Esrc, buf);
}

omxMatrix* omxInitMatrix(int nrows, int ncols, unsigned short isColMajor, omxState* os) {

	if (!isColMajor) mxThrow("All matrices are created column major");

	omxMatrix* om = new omxMatrix;

	om->hasMatrixNumber = 0;
	om->rows = nrows;
	om->cols = ncols;
	om->colMajor = (isColMajor ? 1 : 0);

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

omxMatrix *omxCreateCopyOfMatrix(omxMatrix *orig, omxState *os)
{
	if (!orig) return 0;
	omxMatrix *copy = omxInitMatrix(0, 0, 1, os);
	omxCopyMatrix(copy, orig);
	return copy;
}

void omxCopyMatrix(omxMatrix *dest, omxMatrix *orig) {
	/* Copy a matrix.  NOTE: Matrix maintains its algebra bindings. */

	int regenerateMemory = TRUE;

	if(!dest->owner && (dest->rows == orig->rows && dest->cols == orig->cols)) {
		regenerateMemory = FALSE;				// If it's local data and the right size, we can keep memory.
	}

	dest->rows = orig->rows;
	dest->cols = orig->cols;
	dest->colMajor = orig->colMajor;
	dest->populate = orig->populate;

	if(dest->rows == 0 || dest->cols == 0) {
		omxFreeInternalMatrixData(dest);
		dest->setData(0);
	} else {
		if(regenerateMemory) {
			omxFreeInternalMatrixData(dest);											// Free and regenerate memory
			dest->setData((double*) Calloc(dest->rows * dest->cols, double));
		}
		if (dest->data != orig->data) {  // if equal then programmer stop? TODO
			memcpy(dest->data, orig->data, dest->rows * dest->cols * sizeof(double));
		}
	}

	omxMatrixLeadingLagging(dest);
}

void omxMatrix::setData(double *ptr)
{
	if (allocationLock) mxThrow("Cannot change allocation of matrix '%s'", name());
	data = ptr;
}

void omxMatrix::take(omxMatrix *orig)
{
	omxFreeInternalMatrixData(this);

	this->rows = orig->rows;
	this->cols = orig->cols;
	this->colMajor = orig->colMajor;
	this->populate = orig->populate;

	setData(orig->data);
	this->owner = orig->owner;
	orig->data = 0;
	orig->owner = 0;

	omxMatrixLeadingLagging(this);

	omxMarkDirty(orig);
}

bool omxMatrix::canDiscard()
{
	if (hasMatrixNumber ||             // can be referenced from the front-end
	    populate.size() ||             // has populations from elsewhere
	    (!algebra && !fitFunction) ||  // any simple matrix that we can't recreate
	    (algebra && algebra->oate == &(omxAlgebraSymbolTable[62])))   // broadcast
		return false;
	return true;
}

void omxFreeMatrix(omxMatrix *om) {
    
    if(om == NULL) return;

	omxFreeInternalMatrixData(om);

	if(om->algebra != NULL) {
		omxFreeAlgebraArgs(om->algebra);
		om->algebra = NULL;
	}

	if(om->fitFunction != NULL) {
		delete om->fitFunction;
		om->fitFunction = NULL;
	}

	if (om->freeColnames) for (auto cn : om->colnames) free((void*)cn);
	if (om->freeRownames) for (auto rn : om->rownames) free((void*)rn);
	
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

void omxMatrix::copyAttr(omxMatrix *src)
{
	joinKey = src->joinKey;
	if (src->joinModel) {
		joinModel = omxExpectationFromIndex(src->joinModel->expNum, currentState);
	}
}

// could build index and binary search
int omxMatrix::lookupColumnByName(const char *target)
{
	for (int cx=0; cx < int(colnames.size()); ++cx) {
		if (strEQ(colnames[cx], target)) return cx;
	}
	return -1;
}

void omxResizeMatrix(omxMatrix *om, int nrows, int ncols)
{
	// Always Recompute() before you Resize().
	if(false && OMX_DEBUG) { 
		mxLog("Resizing matrix from (%d, %d) to (%d, %d)",
			om->rows, om->cols, nrows, ncols);
	}

	if( (om->rows != nrows || om->cols != ncols)) {
		omxFreeInternalMatrixData(om);
		om->setData((double*) Calloc(nrows * ncols, double));
	}

	om->rows = nrows;
	om->cols = ncols;

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
	if ((numrow > 1) && (numcol > 1)) {
		mxThrow("Requested improper index (%d) from a malformed vector of dimensions (%d, %d)", 
			index, numrow, numcol);
	} else {
		int Rf_length = (numrow > 1) ? numrow : numcol;
		mxThrow("Requested improper index (%d) from vector of Rf_length (%d)", 
			index, Rf_length);
	}
}

void setMatrixError(omxMatrix *om, int row, int col, int numrow, int numcol) {
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
	mxThrow("Attempted to set row and column (%d, %d) in %s \"%s\" with dimensions %d x %d.", 
		row, col, typeString, om->name(), numrow, numcol);
}

void matrixElementError(int row, int col, omxMatrix *om) {
	mxThrow("Requested improper value (%d, %d) from (%d, %d) matrix '%s'",
		 row, col, om->rows, om->cols, om->name());
}

void setVectorError(int index, int numrow, int numcol) {
	if ((numrow > 1) && (numcol > 1)) {
		mxThrow("Attempting to set improper index (%d) from a malformed vector of dimensions (%d, %d)", 
			index, numrow, numcol);
	} else {
		int Rf_length = (numrow > 1) ? numrow : numcol;
		mxThrow("Setting improper index (%d) from vector of Rf_length %d", 
			index, Rf_length);
	}
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

omxMatrix* omxNewMatrixFromRPrimitive0(SEXP rObject, omxState* state, 
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Creates and populates an omxMatrix with details from an R matrix object. */
	if (Rf_isMatrix(rObject) && Rf_length(rObject) == 1 &&
	    !std::isfinite(REAL(rObject)[0])) return 0;
	omxMatrix *om = NULL;
	om = omxInitMatrix(0, 0, TRUE, state);
	return omxFillMatrixFromRPrimitive(om, rObject, state, hasMatrixNumber, matrixNumber);
}

omxMatrix* omxNewMatrixFromRPrimitive(SEXP rObject, omxState* state, 
	unsigned short hasMatrixNumber, int matrixNumber) {
/* Creates and populates an omxMatrix with details from an R matrix object. */
	omxMatrix *om = NULL;
	om = omxInitMatrix(0, 0, TRUE, state);
	return omxFillMatrixFromRPrimitive(om, rObject, state, hasMatrixNumber, matrixNumber);
}

omxMatrix* omxFillMatrixFromRPrimitive(omxMatrix* om, SEXP rObject, omxState* state,
	unsigned short hasMatrixNumber, int matrixNumber)
{
	if (!om) mxThrow("fillMatrixHelperFunction: matrix must be allocated already");

	if (rObject) {
		if(Rf_isMatrix(rObject)) {
			SEXP matrixDims;
			ScopedProtect p1(matrixDims, Rf_getAttrib(rObject, R_DimSymbol));
			int* dimList = INTEGER(matrixDims);
			om->rows = dimList[0];
			om->cols = dimList[1];
		} else if (Rf_isVector(rObject)) {		// If it's a vector, assume it's a row vector. BLAS doesn't care.
			if(OMX_DEBUG) { mxLog("Vector discovered.  Assuming rowity."); }
			om->rows = 1;
			om->cols = Rf_length(rObject);
		} else {
			mxThrow("Recieved unknown matrix type in omxFillMatrixFromRPrimitive.");
		}
		if(OMX_DEBUG) { mxLog("Matrix connected to (%d, %d) matrix or MxMatrix.", om->rows, om->cols); }

		if (TYPEOF(rObject) != REALSXP) mxThrow("matrix is of type '%s'; only type double is accepted",
							Rf_type2char(TYPEOF(rObject)));

		om->owner = rObject;
		om->data = REAL(om->owner);

		SEXP dimnames;
		ScopedProtect pdn(dimnames, Rf_getAttrib(rObject, R_DimNamesSymbol));
		om->loadDimnames(dimnames);
	}

	om->colMajor = TRUE;
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

void omxMatrix::loadDimnames(SEXP dimnames)
{
	if (!dimnames || Rf_isNull(dimnames)) return;

	if (rownames.size() || colnames.size()) {
		mxThrow("Attempt to load dimnames more than once for %s", name());
	}

	if (Rf_length(dimnames) >= 1) {
		ProtectedSEXP names(VECTOR_ELT(dimnames, 0));
		loadCharVecFromR(name(), names, rownames);
	}
	if (Rf_length(dimnames) >= 2) {
		ProtectedSEXP names(VECTOR_ELT(dimnames, 1));
		loadCharVecFromR(name(), names, colnames);
	}

	if (OMX_DEBUG) {
		std::string buf = stringifyDimnames(this);
		mxLog("matrix %s loaded %s", name(), buf.c_str());
	}
}

void omxMatrix::setJoinInfo(SEXP Rmodel, SEXP Rkey)
{
	int modelIndex = Rf_asInteger(Rmodel);
	if (modelIndex != NA_INTEGER) {
		joinModel = currentState->expectationList[modelIndex - 1];
	}

	int fk = Rf_asInteger(Rkey);
	if (fk != NA_INTEGER) {
		joinKey = fk - 1;
	}
}

void omxMatrix::omxProcessMatrixPopulationList(SEXP matStruct)
{
	setJoinInfo(VECTOR_ELT(matStruct, 1), VECTOR_ELT(matStruct, 2));
	shape = Rf_asInteger(VECTOR_ELT(matStruct, 3));

	const int offsetToPopulationList = 4;
	const int numPopLocs = Rf_length(matStruct) - offsetToPopulationList;

	if(OMX_DEBUG) { mxLog("Processing Population List: %d elements.", numPopLocs); }

	unshareMemoryWithR();

	populate.resize(numPopLocs);

	for(int i = 0; i < numPopLocs; i++) {
		ProtectedSEXP subList(VECTOR_ELT(matStruct, i+offsetToPopulationList));

		int* locations = INTEGER(subList);
		populateLocation &pl = populate[i];
		pl.from = locations[0];
		pl.srcRow = locations[1];
		pl.srcCol = locations[2];
		pl.destRow = locations[3];
		pl.destCol = locations[4];
	}
}

void omxMatrix::addPopulate(omxMatrix *from, int srcRow, int srcCol, int destRow, int destCol)
{
	if (!from->hasMatrixNumber) mxThrow("omxMatrix::addPopulate %s must have matrix number",
					     from->name());
	populate.emplace_back(from->matrixNumber, srcRow, srcCol, destRow, destCol);
}

void omxMatrix::unshareMemoryWithR()
{
	if (!owner) return;

	double *copy = (double*) Realloc(NULL, rows * cols, double);
	memcpy(copy, data, rows * cols * sizeof(double));
	setData(copy);
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
	mat->setData(newdata);  // can probably avoid this, if needed
	mat->colMajor = !mat->colMajor;
}

void omxTransposeMatrix(omxMatrix *mat)
{
	std::swap(mat->colnames, mat->rownames);
	mat->colMajor = !mat->colMajor;
	
	if(mat->rows != mat->cols){
        int mid = mat->rows;
        mat->rows = mat->cols;
        mat->cols = mid;
	}
	
	omxMatrixLeadingLagging(mat);
}

void omxRemoveElements(omxMatrix *om, int removed[])
{
	int oldElements = om->rows * om->cols;
	int nextElement = 0;

	for(int j = 0; j < oldElements; j++) {
		if(removed[j]) continue;
		omxUnsafeSetVectorElement(om, nextElement, omxUnsafeVectorElement(om, j));
		nextElement++;
	}

	if (om->rows > 1) {
		om->rows = nextElement;
	} else {
		om->cols = nextElement;
	}
}

void omxRemoveRowsAndColumns(omxMatrix *om, int rowsRemoved[], int colsRemoved[])
{
	int origRows = om->rows;
	int origCols = om->cols;

	int newRows = origRows;
	for(int j = 0; j < om->rows; j++) {
#if OMX_DEBUG
		if (rowsRemoved[j] != (rowsRemoved[j] & 1)) mxThrow("Removed flag can only be 0 or 1");
#endif
		newRows -= rowsRemoved[j];
	}
	int newCols = origCols;
	for(int j = 0; j < om->cols; j++) {
#if OMX_DEBUG
		if (colsRemoved[j] != (colsRemoved[j] & 1)) mxThrow("Removed flag can only be 0 or 1");
#endif
		newCols -= colsRemoved[j];
	}

	om->rows = newRows;
	om->cols = newCols;
	
	if (om->colMajor) {
		int nextCol = 0;
		for(int j = 0; j < origCols; j++) {
			if(colsRemoved[j]) continue;
			int nextRow = 0;
			for(int k = 0; k < origRows; k++) {
				if(rowsRemoved[k]) continue;
				double val = om->data[j * origRows + k];
				omxSetMatrixElement(om, nextRow, nextCol, val);
				nextRow++;
			}
			nextCol++;
		}
	} else {
		int nextRow = 0;
		for(int k = 0; k < origRows; k++) {
			if(rowsRemoved[k]) continue;
			int nextCol = 0;
			for(int j = 0; j < origCols; j++) {
				if(colsRemoved[j]) continue;
				double val = om->data[k * origCols + j];
				omxSetMatrixElement(om, nextRow, nextCol, val);
				nextCol++;
			}
			nextRow++;
		}
	}

	omxMatrixLeadingLagging(om);  // should be unnecessary, modulo grotesque hacks TODO
}

/* Function wrappers that switch based on inclusion of algebras */
void omxPrint(omxMatrix *source, const char* d) { 					// Pretty-print a (small) matrix
    if(source == NULL) mxLog("%s is NULL.", d);
	else if(source->algebra != NULL) omxAlgebraPrint(source->algebra, d);
	else omxPrintMatrix(source, d);
}

void omxMatrix::markPopulatedEntries()
{
	if (populate.size() == 0) return;

	for (size_t pi = 0; pi < populate.size(); pi++) {
		populateLocation &pl = populate[pi];
		omxSetMatrixElement(this, pl.destRow, pl.destCol, 1.0);
	}
}

void omxMatrix::omxPopulateSubstitutions(int want, FitContext *fc)
{
	if (populate.size() == 0) return;
	if (OMX_DEBUG_ALGEBRA) {
		mxLog("omxPopulateSubstitutions %s, %d locations", name(), (int) populate.size());
	}
	bool modified = false;
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
			if (sourceMatrix->dependsOnParameters()) {
				setDependsOnParameters();
			}
			if (sourceMatrix->dependsOnDefinitionVariables()) {
				setDependsOnDefinitionVariables();
			}
			if (pl.srcRow >= sourceMatrix->rows || pl.srcCol >= sourceMatrix->cols) {
				// may not be properly dimensioned yet
				continue;
			}
		}

		if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_FIT)) {
			double value = omxMatrixElement(sourceMatrix, pl.srcRow, pl.srcCol);
			if (omxMatrixElement(this, pl.destRow, pl.destCol) != value) {
				omxSetMatrixElement(this, pl.destRow, pl.destCol, value);
				modified = true;
				if (OMX_DEBUG_ALGEBRA) {
					mxLog("copying %.2f from %s[%d,%d] to %s[%d,%d]",
					      value, sourceMatrix->name(), pl.srcRow, pl.srcCol,
					      name(), pl.destRow, pl.destCol);
				}
			}
		}
	}
	if (modified) omxMarkClean(this);
}

void omxMatrixLeadingLagging(omxMatrix *om) {
	om->majority = omxMatrixMajorityList[(om->colMajor?1:0)];
	om->minority = omxMatrixMajorityList[(om->colMajor?0:1)];
	om->leading = (om->colMajor?om->rows:om->cols);
	om->lagging = (om->colMajor?om->cols:om->rows);
}

bool omxNeedsUpdate(omxMatrix *matrix)
{
	bool yes;
	if (matrix->hasMatrixNumber && omxMatrixIsClean(matrix)) {
		yes = FALSE;
	} else {
		yes = TRUE;
	}
	if (OMX_DEBUG_ALGEBRA) {
		//mxLog("omxNeedsUpdate %s %s is %s", matrix->getType(), matrix->name(), yes? "dirty" : "clean");
	}
	return yes;
}

void CheckAST(omxMatrix *matrix, FitContext *fc)
{
	if (matrix->algebra) CheckAST(matrix->algebra, fc);
}

void omxRecompute(omxMatrix *matrix, FitContext *fc)
{
	int want = matrix->currentState->getWantStage();
	matrix->omxPopulateSubstitutions(want, fc); // could be an algebra!

	if(!omxNeedsUpdate(matrix)) return;

	if(matrix->algebra) omxAlgebraRecompute(matrix, want, fc);
	else if(matrix->fitFunction) {
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

void MatrixInvert1(omxMatrix *target)
{
	EigenMatrixAdaptor Ein(target);
	// inplace not working as of Aug 2019, maybe try again later
	//Eigen::FullPivLU< Eigen::Ref<Eigen::MatrixXd> > lu(Ein);
	//Eigen::FullPivLU< Eigen::MatrixXd > lu(Ein);
	Eigen::PartialPivLU< Eigen::MatrixXd > lu(Ein);
	// will fail silently if not invertible?
	Ein.derived() = lu.inverse();
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

void omxShallowInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I )
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

		MatrixInvert1(Z);

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
			std::swap(Z, Ax);
			//omxPrint(Z, "Z");
		}
		if(origZ != Z) { 	// Juggling has caused Ax and Z to swap
			omxCopyMatrix(Z, Ax);
		}
	}
}

double omxMaxAbsDiff(omxMatrix *m1, omxMatrix *m2)
{
	if (m1->rows != m2->rows || m1->cols != m2->cols) mxThrow("Matrices are not the same size");

	double mad = 0;
	int size = m1->rows * m1->cols;
	for (int dx=0; dx < size; ++dx) {
		double mad1 = fabs(m1->data[dx] - m2->data[dx]);
		if (mad < mad1) mad = mad1;
	}
	return mad;
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

int omxMatrix::numNonConstElements() const
{
	switch(shape) {
	case 1: //Diag
		return rows;

	case 2: //Full
		return rows * cols;
		
	case 4: //Lower
	case 7: //Symm
		return triangleLoc1(rows);

	case 5: //Sdiag
	case 6: //Stand
		return triangleLoc1(rows - 1);

	case 8: //Unit
	case 9: //Zero
	case 3: //Iden
		return 0;

	default:
		mxThrow("loadFromStream: matrix '%s' with shape %d is unimplemented",
			 name(), shape);
		break;
	}
}

void MoorePenroseInverse(Eigen::Ref<Eigen::MatrixXd> mat)
{
	Eigen::BDCSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeFullV | Eigen::ComputeFullU);
	Eigen::VectorXd sv = svd.singularValues();
	for (int v1=0; v1 < sv.size(); ++v1) {
		if (sv[v1] > 1e-6) sv[v1] = 1.0/sv[v1];
		else sv[v1] = 0;
	}
	mat.derived() = svd.matrixV() * sv.asDiagonal() * svd.matrixU().transpose();
}


SEXP omxMatrix::asR()
{
	int m = rows, n = cols;
	ProtectedSEXP ans(Rcpp::wrap(data, data + m * n));
	ProtectedSEXP dd(Rf_allocVector(INTSXP, 2));
	int *d = INTEGER(dd);
	d[0] = m;
	d[1] = n;
	Rf_setAttrib(ans, R_DimSymbol, dd);
	bool validColnames = int(colnames.size()) == cols;
	bool validRownames = int(rownames.size()) == rows;
	if (validRownames || validColnames) {
		ProtectedSEXP dimnames(Rf_allocVector(VECSXP, 2));
		if (validRownames) {
			ProtectedSEXP names(Rf_allocVector(STRSXP, rows));
			for (int dx=0; dx < rows; ++dx) {
				SET_STRING_ELT(names, dx, Rf_mkChar(rownames[dx]));
			}
			SET_VECTOR_ELT(dimnames, 0, names);
		}
		if (validColnames) {
			ProtectedSEXP names(Rf_allocVector(STRSXP, cols));
			for (int dx=0; dx < cols; ++dx) {
				SET_STRING_ELT(names, dx, Rf_mkChar(colnames[dx]));
			}
			SET_VECTOR_ELT(dimnames, 1, names);
		}
		Rf_setAttrib(ans, R_DimNamesSymbol, dimnames);
	}
	return ans;
}

bool omxMatrix::sameDimnames(omxMatrix *other) const
{
	if (!sameSize(other)) return false;
	if (hasDimnames() != other->hasDimnames()) return false;
	if (!hasDimnames()) return true;
	for (int nx=0; nx < rows; ++nx) if (!strEQ(rownames[nx], other->rownames[nx])) return false;
	for (int nx=0; nx < cols; ++nx) if (!strEQ(colnames[nx], other->colnames[nx])) return false;
	return true;
}
