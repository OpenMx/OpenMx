/*
 *  Copyright 2007-2013 The OpenMx Project
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
#include "npsolWrap.h"

omxData* omxInitData(omxState* os) {

	omxData *od = Calloc(1, omxData);

	od->currentState = os;
	
	if(OMX_DEBUG) {Rprintf("Data's state object is at 0x%x.\n", od->currentState);}

	return od;

}

omxData* omxDataLookupFromState(SEXP dataObject, omxState* state) {
	int dataIdx = INTEGER(dataObject)[0];

	return state->dataList[dataIdx];
}

omxData* omxNewDataFromMxData(SEXP dataObject, omxState* state) {
	if(OMX_DEBUG) {Rprintf("Initializing data Element.\n");}
	if(dataObject == NULL) {
		error("Null Data Object detected.  This is an internal error, and should be reported on the forums.\n");
	}

	omxData* od = omxInitData(state);

	SEXP dataLoc, dataVal;
	int numCols;
	int numInts=0;
	int numReals=0;

	// PARSE MxData Structure
	if(OMX_DEBUG) {Rprintf("Processing Data Type.\n");}
	PROTECT(dataLoc = GET_SLOT(dataObject, install("type")));
	if(dataLoc == NULL) { error("Data has no type.  Sorry.\nThis is an internal error, and should be reported on the forums.\n");}
	PROTECT(dataVal = STRING_ELT(dataLoc,0));
	od->_type = CHAR(dataVal);
	if(OMX_DEBUG) {Rprintf("Element is type %s.\n", od->_type);}

	PROTECT(dataLoc = GET_SLOT(dataObject, install("observed")));
	if(OMX_DEBUG) {Rprintf("Processing Data Elements.\n");}
	if(isFrame(dataLoc)) {
		if(OMX_DEBUG) {Rprintf("Data is a frame.\n");}
		// Process Data Frame into Columns
		od->cols = length(dataLoc);
		if(OMX_DEBUG) {Rprintf("Data has %d columns.\n", od->cols);}
		numCols = od->cols;
		od->realData = (double**) R_alloc(numCols, sizeof(double*));
		od->intData = (int**) R_alloc(numCols, sizeof(int*));
		od->location = (int*) R_alloc(numCols, sizeof(int));
		for(int j = 0; j < numCols; j++) {
			SEXP rcol;
			PROTECT(rcol = VECTOR_ELT(dataLoc, j));
			if(isFactor(rcol)) {
				if(OMX_DEBUG) {Rprintf("Column %d is a factor.\n", j);}
				od->intData[numInts] = INTEGER(rcol);
				od->location[j] = ~(numInts++);
				od->numFactor++;
			} else if (isInteger(rcol)) {
				error("Internal error: Column %d is in integer format.", j);
			} else {
				if(OMX_DEBUG) {Rprintf("Column %d is a numeric.\n", j);}
				od->realData[numReals] = REAL(rcol);
				od->location[j] = (numReals++);
				od->numNumeric++;
			}
		}
		od->rows = length(VECTOR_ELT(dataLoc, 0));
		if(OMX_DEBUG) {Rprintf("And %d rows.\n", od->rows);}
	} else {
		if(OMX_DEBUG) {Rprintf("Data contains a matrix.\n");}
		od->dataMat = omxNewMatrixFromRPrimitive(dataLoc, od->currentState, 0, 0);
		
		if (od->dataMat->colMajor && strncmp(od->_type, "raw", 3) == 0) { 
			omxToggleRowColumnMajor(od->dataMat);
		}
		od->cols = od->dataMat->cols;
		od->rows = od->dataMat->rows;
		od->numNumeric = od->cols;
	}

	if(OMX_DEBUG) {Rprintf("Processing Means Matrix.\n");}
	PROTECT(dataLoc = GET_SLOT(dataObject, install("means")));
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
	        if(od->meansMat == NULL) {Rprintf("No means found.\n");}
		else {omxPrint(od->meansMat, "Means Matrix is:");}
        }

	if(strncmp(od->_type, "raw", 3) != 0) {
		if(OMX_DEBUG) {Rprintf("Processing Observation Count.\n");}
		PROTECT(dataLoc = GET_SLOT(dataObject, install("numObs")));
		od->numObs = REAL(dataLoc)[0];
	} else {
		od->numObs = od->rows;
		if(OMX_DEBUG) {Rprintf("Processing presort metadata.\n");}
		/* For raw data, process sorting metadata. */
		// Process unsorted indices:  // TODO: Generate reverse lookup table
		PROTECT(dataLoc = GET_SLOT(dataObject, install("indexVector")));
		od->indexVector = INTEGER(dataLoc);
		if(od->indexVector[0] == R_NaInt) od->indexVector = NULL;
		// Process pre-computed identicality checks
		if(OMX_DEBUG) {Rprintf("Processing definition variable identicality.\n");}
		PROTECT(dataLoc = GET_SLOT(dataObject, install("identicalDefVars")));
		od->identicalDefs = INTEGER(dataLoc);
		if(od->identicalDefs[0] == R_NaInt) od->identicalDefs = NULL;
		if(OMX_DEBUG) {Rprintf("Processing missingness identicality.\n");}
		PROTECT(dataLoc = GET_SLOT(dataObject, install("identicalMissingness")));
		od->identicalMissingness = INTEGER(dataLoc);
		if(od->identicalMissingness[0] == R_NaInt) od->identicalMissingness = NULL;
		if(OMX_DEBUG) {Rprintf("Processing row identicality.\n");}
		PROTECT(dataLoc = GET_SLOT(dataObject, install("identicalRows")));
		od->identicalRows = INTEGER(dataLoc);
		if(od->identicalRows[0] == R_NaInt) od->identicalRows = NULL;
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
		error("Use a data frame for factor data");
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

unsigned short int omxDataColumnIsFactor(omxData *od, int col) {
	if(od->dataMat != NULL) return FALSE;
	if(col <= od->cols) return (od->location[col] < 0);

	error("Attempted to access column %d of a %d-column data object", col, od->cols);
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

omxMatrix* omxContiguousDataRow(omxData *od, int row, int start, int length, omxMatrix* om) {
	// TODO: Might be better to combine this with omxDataRow to make a single accessor omxDataRow with a second signature that accepts an omxContiguousData argument.
	if(row > od->rows) return NULL;	// Sanity check

	if(om == NULL) {
		om = omxInitMatrix(om, 1, od->cols, TRUE, od->currentState);
	}
	
	int numcols = od->cols;
	omxMatrix* dataMat = od->dataMat;
	double *dest = om->data;
	double *source = dataMat->data + row * numcols + start;
	memcpy(dest, source, sizeof(double) * length);
	return(om);
}

omxMatrix* omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om) {

	if(colList == NULL || row > od->rows) return NULL;	// Sanity check

	if(om == NULL) {
		om = omxInitMatrix(om, 1, od->cols, TRUE, od->currentState);
	}

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
	return om;
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
	numrow = nrows(data);
	numcol = ncols(data);
	defvarcol = ncols(defvars);
	PROTECT(retval = allocVector(VECSXP, 3));
	PROTECT(identicalRows = allocVector(INTSXP, numrow));
	PROTECT(identicalMissing = allocVector(INTSXP, numrow));
	PROTECT(identicalDefvars = allocVector(INTSXP, numrow));
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
	UNPROTECT(4); // retval, identicalRows, identicalMissing, identicalDefvars
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
	numrow = length(VECTOR_ELT(data, 0));
	numcol = length(data);
	defvarcol = length(defvars);
	PROTECT(retval = allocVector(VECSXP, 3));
	PROTECT(identicalRows = allocVector(INTSXP, numrow));
	PROTECT(identicalMissing = allocVector(INTSXP, numrow));
	PROTECT(identicalDefvars = allocVector(INTSXP, numrow));
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
	UNPROTECT(4); // retval, identicalRows, identicalMissing, identicalDefvars
	return retval;
}

SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
			   SEXP skipMissingness, SEXP skipDefvars) {
	if (isMatrix(data)) {
		return(findIdenticalMatrix(data, missing, defvars, skipMissingness, skipDefvars));
	} else {
		return(findIdenticalDataFrame(data, missing, defvars, skipMissingness, skipDefvars));
	}
}

void omxPrintData(omxData *od, const char *header) {
	if (!header) error("omxPrintData: header is NULL");

	if (!od) {
		Rprintf("%s: NULL\n", header);
		return;
	}

	Rprintf("%s(%s): %f observations %d x %d\n", header, od->_type, od->numObs,
		od->rows, od->cols);
	Rprintf("numNumeric %d numFactor %d\n", od->numNumeric, od->numFactor);

	if (od->location) {
		for(int j = 0; j < od->cols; j++) {
			int loc = od->location[j];
			if (loc < 0) {
				Rprintf("Integer[%d]:", j);
				int *val = od->intData[~loc];
				for (int vx=0; vx < od->numObs; vx++) {
					Rprintf(" %d", val[vx]);
				}
			} else {
				Rprintf("Numeric[%d]:", j);
				double *val = od->realData[loc];
				for (int vx=0; vx < od->numObs; vx++) {
					Rprintf(" %.3g", val[vx]);
				}
			}
			Rprintf("\n");
		}
	}

	if (od->identicalRows) {
		Rprintf("\trow\tmissing\tdefvars\n");
		for(int j = 0; j < od->rows; j++) {
			Rprintf("[%d]\t%d\t%d\t%d\n", j, od->identicalRows[j],
				od->identicalMissingness[j], od->identicalDefs[j]);
		}
	}

	if (od->dataMat) omxPrintMatrix(od->dataMat, "dataMat");
	if (od->meansMat) omxPrintMatrix(od->meansMat, "meansMat");
}

