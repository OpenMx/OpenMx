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
 *  omxData.cc
 *
 *  Created: Timothy R. Brick 	Date: 2009-07-15
 *
 *	Contains code for the omxData class
 *   omxData objects hold data in whatever form it takes
 *
 **********************************************************/
#include "omxData.h"
#include "glue.h"
#include "omxState.h"

omxData* omxInitData(omxState* os) {

	omxData *od = Calloc(1, omxData);

	od->currentState = os;

	if (os != globalState) Rf_error("Too late to create omxData");
	os->dataList.push_back(od);

	return od;

}

omxData* omxDataLookupFromState(SEXP dataObject, omxState* state) {
	int dataIdx = INTEGER(dataObject)[0];

	return state->dataList[dataIdx];
}

omxData* omxNewDataFromMxData(SEXP dataObject, omxState* state) {
	if(OMX_DEBUG) {mxLog("Initializing data Element.");}
	if(dataObject == NULL) {
		Rf_error("Null Data Object detected.  This is an internal Rf_error, and should be reported on the forums.\n");
	}

	omxData* od = omxInitData(state);

	SEXP dataLoc, dataVal;
	int numCols;

	// PARSE MxData Structure
	if(OMX_DEBUG) {mxLog("Processing Data Type.");}
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("type")));
	if(dataLoc == NULL) { Rf_error("Data has no type.  Sorry.\nThis is an internal Rf_error, and should be reported on the forums.\n");}
	Rf_protect(dataVal = STRING_ELT(dataLoc,0));
	od->_type = CHAR(dataVal);
	if(OMX_DEBUG) {mxLog("Element is type %s.", od->_type);}

	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("observed")));
	if(OMX_DEBUG) {mxLog("Processing Data Elements.");}
	if(Rf_isFrame(dataLoc)) {
		if(OMX_DEBUG) {mxLog("Data is a frame.");}
		// Process Data Frame into Columns
		od->cols = Rf_length(dataLoc);
		if(OMX_DEBUG) {mxLog("Data has %d columns.", od->cols);}
		numCols = od->cols;
		od->realData = (double**) R_alloc(numCols, sizeof(double*));
		OMXZERO(od->realData, numCols);
		od->intData = (int**) R_alloc(numCols, sizeof(int*));
		OMXZERO(od->intData, numCols);
		od->location = (int*) R_alloc(numCols, sizeof(int));
		SEXP colnames;
		Rf_protect(colnames = Rf_getAttrib(dataLoc, R_NamesSymbol));
		for(int j = 0; j < numCols; j++) {
			SEXP rcol;
			Rf_protect(rcol = VECTOR_ELT(dataLoc, j));
			if(Rf_isFactor(rcol)) {
				if (Rf_isUnordered(rcol)) {
					Rf_warning("Data[%d] '%s' must be an ordered factor. Please use mxFactor()",
						j+1, CHAR(STRING_ELT(colnames, j)));
				}
				if(OMX_DEBUG) {mxLog("Column %d is a factor.", j);}
				od->intData[j] = INTEGER(rcol);
				od->location[j] = ~j;
				od->numFactor++;
			} else if (Rf_isInteger(rcol)) {
				Rf_error("Internal Rf_error: Column %d is in integer format.", j);
			} else {
				if(OMX_DEBUG) {mxLog("Column %d is a numeric.", j);}
				od->realData[j] = REAL(rcol);
				od->location[j] = j;
				od->numNumeric++;
			}
		}
		od->rows = Rf_length(VECTOR_ELT(dataLoc, 0));
		if(OMX_DEBUG) {mxLog("And %d rows.", od->rows);}
	} else {
		if(OMX_DEBUG) {mxLog("Data contains a matrix.");}
		od->dataMat = omxNewMatrixFromRPrimitive(dataLoc, od->currentState, 0, 0);
		
		if (od->dataMat->colMajor && strncmp(od->_type, "raw", 3) == 0) { 
			omxToggleRowColumnMajor(od->dataMat);
		}
		od->cols = od->dataMat->cols;
		od->rows = od->dataMat->rows;
		od->numNumeric = od->cols;
	}

	if(OMX_DEBUG) {mxLog("Processing Means Matrix.");}
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("means")));
	od->meansMat = omxNewMatrixFromRPrimitive(dataLoc, od->currentState, 0, 0);
	if(od->meansMat->rows == 1 && od->meansMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->meansMat, 0, 0)) ||
	    !isfinite(omxMatrixElement(od->meansMat, 0, 0)))) {
		omxFreeMatrixData(od->meansMat); // Clear just-allocated memory.
		od->meansMat = NULL;  // 1-by-1 matrix of NAs is a null means matrix.
                // FIXME: The above check may cause problems for dynamic data if the means
                //          originally is a 1x1 that has not yet been calculated.  This should be
                //          adjusted.
	}
	
	if(OMX_DEBUG) {
	        if(od->meansMat == NULL) {mxLog("No means found.");}
		else {omxPrint(od->meansMat, "Means Matrix is:");}
	}

	if(OMX_DEBUG) {mxLog("Processing Asymptotic Covariance Matrix.");}
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("acov")));
	od->acovMat = omxNewMatrixFromRPrimitive(dataLoc, od->currentState, 0, 0);
	if(od->acovMat->rows == 1 && od->acovMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->acovMat, 0, 0)) ||
	    !isfinite(omxMatrixElement(od->acovMat, 0, 0)))) {
		omxFreeMatrixData(od->acovMat); // Clear just-allocated memory.
		od->acovMat = NULL;
	}

	if(OMX_DEBUG) {mxLog("Processing Observed Thresholds Matrix.");}
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("thresholds")));
	od->obsThresholdsMat = omxNewMatrixFromRPrimitive(dataLoc, od->currentState, 0, 0);
	if(od->obsThresholdsMat->rows == 1 && od->obsThresholdsMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->obsThresholdsMat, 0, 0)) ||
	    !isfinite(omxMatrixElement(od->obsThresholdsMat, 0, 0)))) {
		omxFreeMatrixData(od->obsThresholdsMat); // Clear just-allocated memory.
		od->obsThresholdsMat = NULL;
	} else {
        int nCol = od->obsThresholdsMat->cols;
		/* Match threshold column names and build ThresholdCols structure */
		od->thresholdCols = (omxThresholdColumn*) R_alloc(nCol, sizeof(omxThresholdColumn));
        Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("thresholdColumns")));
        int *columns = INTEGER(dataLoc);
        Rf_protect(dataVal = R_do_slot(dataObject, Rf_install("thresholdLevels")));
        int *levels = INTEGER(dataVal);
        for(int i = 0; i < od->obsThresholdsMat->cols; i++) {
            od->thresholdCols[i].matrix = od->obsThresholdsMat;
            od->thresholdCols[i].column = columns[i];
            od->thresholdCols[i].numThresholds = levels[i];
			if(OMX_DEBUG) { mxLog("omxData observed thresholds SOURCE column %d.", i); }
			if(OMX_DEBUG) {omxPrintMatrix(od->thresholdCols[i].matrix, "...Thresholds matrix"); }
			if(OMX_DEBUG) { mxLog("omxData observed thresholds DEST column %d.", od->thresholdCols[i].column); }
			if(OMX_DEBUG) { mxLog("omxData observed thresholds numThresholds %d.", od->thresholdCols[i].numThresholds); }
        }
	}

	if(strncmp(od->_type, "raw", 3) != 0) {
		if(OMX_DEBUG) {mxLog("Processing Observation Count.");}
		Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("numObs")));
		od->numObs = REAL(dataLoc)[0];
	} else {
		od->numObs = od->rows;
		if(OMX_DEBUG) {mxLog("Processing presort metadata.");}
		/* For raw data, process sorting metadata. */
		// Process unsorted indices:  // TODO: Generate reverse lookup table
		Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("indexVector")));
		od->indexVector = INTEGER(dataLoc);
		if(Rf_length(dataLoc) == 0 || od->indexVector[0] == R_NaInt) od->indexVector = NULL;
		// Process pre-computed identicality checks
		if(OMX_DEBUG) {mxLog("Processing definition variable identicality.");}
		Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("identicalDefVars")));
		od->identicalDefs = INTEGER(dataLoc);
		if(Rf_length(dataLoc) == 0 || od->identicalDefs[0] == R_NaInt) od->identicalDefs = NULL;
		if(OMX_DEBUG) {mxLog("Processing missingness identicality.");}
		Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("identicalMissingness")));
		od->identicalMissingness = INTEGER(dataLoc);
		if(Rf_length(dataLoc) == 0 || od->identicalMissingness[0] == R_NaInt) od->identicalMissingness = NULL;
		if(OMX_DEBUG) {mxLog("Processing row identicality.");}
		Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("identicalRows")));
		od->identicalRows = INTEGER(dataLoc);
		if(Rf_length(dataLoc) == 0 || od->identicalRows[0] == R_NaInt) od->identicalRows = NULL;
	}
	return od;
}

void resetDefinitionVariables(double *oldDefs, int numDefs) {
	int nextDef;

	for(nextDef = 0; nextDef < numDefs; nextDef++) {
		oldDefs[nextDef] = NA_REAL;					// Def Vars default to NA
	}

}

void omxFreeData(omxData* od) {
	omxFreeAllMatrixData(od->dataMat);
	omxFreeAllMatrixData(od->meansMat);
	omxFreeAllMatrixData(od->acovMat);
	omxFreeAllMatrixData(od->obsThresholdsMat);
	Free(od);
}

double omxDoubleDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		return omxMatrixElement(od->dataMat, row, col);
	}
	int location = od->location[col];
	if(location < 0) {
		return (double)(od->intData[~location][row]);
	} else {
		return od->realData[location][row];
	}
}

int omxIntDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		Rf_error("Use a data frame for factor data");
	}

	int location = od->location[col];
	if(location < 0) {
		return (od->intData[~location][row]);
	} else {
		return (int)(od->realData[location][row]);
	}
}

omxMatrix* omxDataMatrix(omxData *od, omxMatrix* om) {
	double dataElement;

	if(od->dataMat != NULL) {		// Data was entered as a matrix.
		if(om != NULL) {			// It stays as such
			omxCopyMatrix(om, od->dataMat);
			return om;
		}
		return od->dataMat;
	}
	// Otherwise, we must construct the matrix.
	int numRows = od->rows, numCols = od->cols;

	if(om == NULL) {
		om = omxInitMatrix(om, numRows, numCols, TRUE, od->currentState);
	}

	if(om->rows != numRows || om->cols != numCols) {
		omxResizeMatrix(om, numRows, numCols, FALSE);
	}

	for(int j = 0; j < numCols; j++) {
		for(int k = 0; k < numRows; k++) {
			int location = od->location[j];
			if(location < 0) {
				dataElement = (double) od->intData[~location][k];
			} else {
				dataElement = od->realData[location][k];
			}
			omxSetMatrixElement(om, k, j, dataElement);
		}
	}
	return om;
}

omxMatrix* omxDataAcov(omxData *od, omxMatrix* om) {
	if(od->acovMat != NULL) {		// acov was entered as a matrix.
		if(om != NULL) {			// It stays as such
			omxAliasMatrix(om, od->acovMat);
			return om;
		}
		return od->acovMat;
	}
	// Otherwise, we must construct the matrix.
	int numRows = ( (od->rows)*(od->rows + 1) ) / 2;
	
	if(om == NULL) {
		om = omxInitMatrix(om, numRows, numRows, TRUE, od->currentState);
	}
	omxCopyMatrix(om, od->acovMat);//omxAliasMatrix(om, od->acovMat); // Could also be done with omxCopyMatrix.
	return om;
}

unsigned short int omxDataColumnIsFactor(omxData *od, int col) {
	if(od->dataMat != NULL) return FALSE;
	if(col <= od->cols) return (od->location[col] < 0);

	Rf_error("Attempted to access column %d of a %d-column data object", col, od->cols);
	return 0; // not reached
}

omxMatrix* omxDataMeans(omxData *od, omxMatrix* colList, omxMatrix* om) {

	if(od->meansMat == NULL) return NULL;

	if(colList == NULL) {
		if(om == NULL) return od->meansMat;
		omxCopyMatrix(om, od->meansMat);
		return om;
	}

	int cols = colList->cols;

	if(colList == NULL || cols == 0 || cols > od->cols) {
		cols = od->cols;
		if(om == NULL) return od->meansMat;
		omxCopyMatrix(om, od->meansMat);
		return om;
	}

	if(om == NULL) {
		om = omxInitMatrix(om, 1, cols, TRUE, od->currentState);
	}

	for(int i = 0; i < cols; i++) {
		omxSetMatrixElement(om, 1, i, omxVectorElement(od->meansMat, omxVectorElement(colList, i)));
	}

	return om;
}

omxThresholdColumn* omxDataThresholds(omxData *od) {
    return od->thresholdCols;
}

void omxSetContiguousDataColumns(omxContiguousData* contiguous, omxData* data, omxMatrix* colList) {

	contiguous->isContiguous = FALSE;   // Assume not contiguous

	if (data->dataMat == NULL) return; // Data has no matrix elements, so skip.

	omxMatrix* dataMat = data->dataMat;
	if (dataMat->colMajor) return;      // If data matrix is column-major, there's no continuity
	
	int colListLength = colList->cols;              // # of columns in the cov matrix
	double start = omxVectorElement(colList, 0);    // Data col of first column of the covariance
	contiguous->start = (int) start;                // That's our starting point.
	contiguous->length = colListLength;             // And the length is ncol(cov)
	
	for(int i = 1; i < colListLength; i++) {        // Make sure that the col list is 
		double next = omxVectorElement(colList, i); // contiguously increasing in column number
		if (next != (start + i)) return;            // If it isn't, it's not contiguous data
	}
	
	contiguous->isContiguous = TRUE;    // Passed.  This is contiguous.
}

void omxContiguousDataRow(omxData *od, int row, int start, int Rf_length, omxMatrix* om) {
	// TODO: Might be better to combine this with omxDataRow to make a single accessor omxDataRow with a second signature that accepts an omxContiguousData argument.
	if(row >= od->rows) Rf_error("Invalid row");

	if(om == NULL) Rf_error("Must provide an output matrix");
	
	int numcols = od->cols;
	omxMatrix* dataMat = od->dataMat;
	double *dest = om->data;
	double *source = dataMat->data + row * numcols + start;
	memcpy(dest, source, sizeof(double) * Rf_length);
}

void omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om) {

	if(colList == NULL || row >= od->rows) Rf_error("Invalid row or colList");

	if(om == NULL) Rf_error("Must provide an output matrix");

	int numcols = om->cols;
	if(od->dataMat != NULL) { // Matrix Object
		omxMatrix* dataMat = od->dataMat;
		for(int j = 0; j < numcols; j++) {
			omxSetMatrixElement(om, 0, j, omxMatrixElement(dataMat, row, 
								       omxVectorElement(colList, j)));
		}
	} else {		// Data Frame object
		double dataElement;
		int* locations = od->location;
		int** intDataColumns = od->intData;
		double **realDataColumns = od->realData;
		for(int j = 0; j < numcols; j++) {
			int location = locations[(int)omxVectorElement(colList, j)];
			if(location < 0) {
				dataElement = (double) intDataColumns[~location][row];
			} else {
				dataElement = realDataColumns[location][row];
			}
			omxSetMatrixElement(om, 0, j, dataElement);
		}
	}
}

int omxDataIndex(omxData *od, int row) {
	if(od->indexVector != NULL)
		return od->indexVector[row];
	else return row;
}

int omxDataNumIdenticalRows(omxData *od, int row) {
	if(od->identicalRows != NULL)
		return od->identicalRows[row];
	else return 1;
}
int omxDataNumIdenticalMissingness(omxData *od, int row) {
	if(od->identicalMissingness != NULL)
		return od->identicalMissingness[row];
	else return 1;
}

int omxDataNumIdenticalDefs(omxData *od, int row){
	if(od->identicalDefs != NULL)
		return od->identicalDefs[row];
	else return 1;
}

int omxDataNumIdenticalContinuousRows(omxData *od, int row) {
	if(od->numNumeric <= 0) {
		return od->rows;
	}
	return omxDataNumIdenticalRows(od, row);
}

int omxDataNumIdenticalContinuousMissingness(omxData *od, int row) {
	if(od->numNumeric <= 0) {
		return od->rows;
	}
	return omxDataNumIdenticalMissingness(od, row);
}

int omxDataNumIdenticalOrdinalRows(omxData *od, int row) {
	if(od->numFactor <= 0) {
		return od->rows;
	}
	return omxDataNumIdenticalRows(od, row);
}

int omxDataNumIdenticalOrdinalMissingness(omxData *od, int row) {
	if(od->numFactor <= 0) {
		return od->rows;
	}
	return omxDataNumIdenticalMissingness(od, row);
}


double omxDataNumObs(omxData *od) {
	return od->numObs;
}

int omxDataNumFactor(omxData *od) {
	return od->numFactor;
}

int omxDataNumNumeric(omxData *od) {
	return od->numNumeric;
}

const char *omxDataType(omxData *od) {
	return od->_type;
}

int elementEqualsDataframe(SEXP column, int offset1, int offset2) {
	switch (TYPEOF(column)) {
	case REALSXP:
		if(ISNA(REAL(column)[offset1])) return ISNA(REAL(column)[offset2]);
		if(ISNA(REAL(column)[offset2])) return ISNA(REAL(column)[offset1]);
		return(REAL(column)[offset1] == REAL(column)[offset2]);
	case LGLSXP:
	case INTSXP:
		return(INTEGER(column)[offset1] == INTEGER(column)[offset2]);		
	}
	return(0);
}

int testRowDataframe(SEXP data, int numrow, int numcol, int i, int *row, int base) {
	SEXP column;
	int j, equal = TRUE;

	if (i == numrow) {
		equal = FALSE;
	} else {
		for(j = 0; j < numcol && equal; j++) {
			column = VECTOR_ELT(data, j);
			equal = elementEqualsDataframe(column, base, i);
		}
	}

	if (!equal) {
		int gap = i - base;
		for(j = 0; j < gap; j++) {
			row[base + j] = gap - j;
		}
		base = i;
	}
	return(base);
}

int elementEqualsMatrix(SEXP data, int row1, int row2, int numrow, int col) {
	int coloffset = col * numrow;
	switch (TYPEOF(data)) {
	case REALSXP:
		if(ISNA(REAL(data)[row1 + coloffset])) return ISNA(REAL(data)[row2 + coloffset]);
		if(ISNA(REAL(data)[row2 + coloffset])) return ISNA(REAL(data)[row1 + coloffset]);
		return(REAL(data)[row1 + coloffset] == REAL(data)[row2 + coloffset]);
	case LGLSXP:
	case INTSXP:
		return(INTEGER(data)[row1 + coloffset] == INTEGER(data)[row2 + coloffset]);
	}
	return(0);
}

int testRowMatrix(SEXP data, int numrow, int numcol, int i, int *row, int base) {
	int j, equal = TRUE;

	if (i == numrow) {
		equal = FALSE;
	} else {
		for(j = 0; j < numcol && equal; j++) {
			equal = elementEqualsMatrix(data, i, base, numrow, j);
		}
	}

	if (!equal) {
		int gap = i - base;
		for(j = 0; j < gap; j++) {
			row[base + j] = gap - j;
		}
		base = i;
	}
	return(base);
}

SEXP findIdenticalMatrix(SEXP data, SEXP missing, SEXP defvars,
			 SEXP skipMissingExp, SEXP skipDefvarsExp) {

	SEXP retval, identicalRows, identicalMissing, identicalDefvars;
	int i, numrow, numcol, defvarcol;
	int *irows, *imissing, *idefvars;
	int baserows, basemissing, basedefvars;
	int skipMissing, skipDefvars;

	skipMissing = LOGICAL(skipMissingExp)[0];
	skipDefvars = LOGICAL(skipDefvarsExp)[0];
	numrow = Rf_nrows(data);
	numcol = Rf_ncols(data);
	defvarcol = Rf_ncols(defvars);
	Rf_protect(retval = Rf_allocVector(VECSXP, 3));
	Rf_protect(identicalRows = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalMissing = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalDefvars = Rf_allocVector(INTSXP, numrow));
	irows = INTEGER(identicalRows);
	imissing = INTEGER(identicalMissing);
	idefvars = INTEGER(identicalDefvars);
	if (skipMissing) {
		for(i = 0; i < numrow; i++) {
			imissing[i] = numrow - i;
		}
	}
	if (skipDefvars) {
		for(i = 0; i < numrow; i++) {
			idefvars[i] = numrow - i;
		}
	}
	baserows = 0;
	basemissing = 0;
	basedefvars = 0;
	for(i = 1; i <= numrow; i++) {
		baserows = testRowMatrix(data, numrow, numcol, i, irows, baserows); 
		if (!skipMissing) {
			basemissing = testRowMatrix(missing, numrow, numcol, i, imissing, basemissing); 
		}
		if (!skipDefvars) {
			basedefvars = testRowMatrix(defvars, numrow, defvarcol, i, idefvars, basedefvars);
		}
	}
	SET_VECTOR_ELT(retval, 0, identicalRows);
	SET_VECTOR_ELT(retval, 1, identicalMissing);
	SET_VECTOR_ELT(retval, 2, identicalDefvars);
	Rf_unprotect(4); // retval, identicalRows, identicalMissing, identicalDefvars
	return retval;
}

SEXP findIdenticalDataFrame(SEXP data, SEXP missing, SEXP defvars,
			    SEXP skipMissingExp, SEXP skipDefvarsExp) {

	SEXP retval, identicalRows, identicalMissing, identicalDefvars;
	int i, numrow, numcol, defvarcol;
	int *irows, *imissing, *idefvars;
	int baserows, basemissing, basedefvars;
	int skipMissing, skipDefvars;

	skipMissing = LOGICAL(skipMissingExp)[0];
	skipDefvars = LOGICAL(skipDefvarsExp)[0];
	numrow = Rf_length(VECTOR_ELT(data, 0));
	numcol = Rf_length(data);
	defvarcol = Rf_length(defvars);
	Rf_protect(retval = Rf_allocVector(VECSXP, 3));
	Rf_protect(identicalRows = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalMissing = Rf_allocVector(INTSXP, numrow));
	Rf_protect(identicalDefvars = Rf_allocVector(INTSXP, numrow));
	irows = INTEGER(identicalRows);
	imissing = INTEGER(identicalMissing);
	idefvars = INTEGER(identicalDefvars);
	if (skipMissing) {
		for(i = 0; i < numrow; i++) {
			imissing[i] = numrow - i;
		}
	}
	if (skipDefvars) {
		for(i = 0; i < numrow; i++) {
			idefvars[i] = numrow - i;
		}
	}
	baserows = 0;
	basemissing = 0;
	basedefvars = 0;
	for(i = 1; i <= numrow; i++) {
		baserows = testRowDataframe(data, numrow, numcol, i, irows, baserows); 
		if (!skipMissing) {
			basemissing = testRowMatrix(missing, numrow, numcol, i, imissing, basemissing);
		}
		if (!skipDefvars) {
			basedefvars = testRowDataframe(defvars, numrow, defvarcol, i, idefvars, basedefvars);
		}
	}
	SET_VECTOR_ELT(retval, 0, identicalRows);
	SET_VECTOR_ELT(retval, 1, identicalMissing);
	SET_VECTOR_ELT(retval, 2, identicalDefvars);
	Rf_unprotect(4); // retval, identicalRows, identicalMissing, identicalDefvars
	return retval;
}

SEXP findIdenticalRowsData2(SEXP data, SEXP missing, SEXP defvars,
			   SEXP skipMissingness, SEXP skipDefvars) {
	if (Rf_isMatrix(data)) {
		return(findIdenticalMatrix(data, missing, defvars, skipMissingness, skipDefvars));
	} else {
		return(findIdenticalDataFrame(data, missing, defvars, skipMissingness, skipDefvars));
	}
}

SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
			   SEXP skipMissingness, SEXP skipDefvars)
{
	omxManageProtectInsanity protectManager;

	try {
		return findIdenticalRowsData2(data, missing, defvars,
					      skipMissingness, skipDefvars);
	} catch( std::exception& __ex__ ) {
		exception_to_try_Rf_error( __ex__ );
	} catch(...) {
		string_to_try_Rf_error( "c++ exception (unknown reason)" );
	}
}


void omxPrintData(omxData *od, const char *header, int maxRows)
{
	if (!header) header = "Default data";

	if (!od) {
		mxLog("%s: NULL", header);
		return;
	}

	std::string buf;
	buf += string_snprintf("%s(%s): %f observations %d x %d\n", header, od->_type, od->numObs,
			       od->rows, od->cols);
	buf += string_snprintf("Row consists of %d numeric, %d ordered factor:", od->numNumeric, od->numFactor);

        int upto = od->numObs;
        if (maxRows >= 0 && maxRows < upto) upto = maxRows;

	if (od->location) {
		for(int j = 0; j < od->cols; j++) {
			int loc = od->location[j];
			if (loc < 0) {
				buf += " I";
			} else {
				buf += " N";
			}
		}
		buf += "\n";

		for (int vx=0; vx < upto; vx++) {
			for (int j = 0; j < od->cols; j++) {
				int loc = od->location[j];
				if (loc < 0) {
					int *val = od->intData[~loc];
					if (val[vx] == NA_INTEGER) {
						buf += " NA,";
					} else {
						buf += string_snprintf(" %d,", val[vx]);
					}
				} else {
					double *val = od->realData[loc];
					if (val[vx] == NA_REAL) {
						buf += " NA,";
					} else {
						buf += string_snprintf(" %.3g,", val[vx]);
					}
				}
			}
			buf += "\n";
		}
	}

	if (od->identicalRows) {
		buf += "DUPS\trow\tmissing\tdefvars\n";
		for(int j = 0; j < upto; j++) {
			buf += string_snprintf("%d\t%d\t%d\t%d\n", j, od->identicalRows[j],
					       od->identicalMissingness[j], od->identicalDefs[j]);
		}
	}
	mxLogBig(buf);

	if (od->dataMat) omxPrintMatrix(od->dataMat, "dataMat");
	if (od->meansMat) omxPrintMatrix(od->meansMat, "meansMat");
}

void omxPrintData(omxData *od, const char *header)
{
        omxPrintData(od, header, -1);
}

