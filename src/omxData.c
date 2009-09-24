/*
 *  Copyright 2007-2009 The OpenMx Project
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

omxData* omxInitData(omxData* od, omxState* os) {

	if(od == NULL) {
		od = (omxData*)R_alloc(1, sizeof(omxData));
	}

	od->dataObject = NULL;
	od->columns = NULL;
	od->dataMat = NULL;
	od->meansMat = NULL;
	od->rows = 0;
	od->cols = 0;
	od->location = NULL;
	od->intData = NULL;
	od->realData = NULL;
	od->currentState = os;
	if(OMX_DEBUG) {Rprintf("Data's state object is at 0x%x.\n", od->currentState);}

	return od;

}

omxData* omxNewDataFromMxDataPtr(SEXP dataObject, omxState* state) {
	int dataIdx = INTEGER(dataObject)[0];
	
	return state->dataList[dataIdx];
}

omxData* omxNewDataFromMxData(omxData* data, SEXP dataObject, omxState* state) {
	if(OMX_DEBUG) {Rprintf("Initializing data Element.\n");}
	if(dataObject == NULL) {
		error("Null Data Object detected.\n");
		// Badness has occurred.  This data object does not exist.
		return NULL;
	}

	omxData* od = data;
	if(od == NULL) {
		od = omxInitData(data, state);
	}
	
	SEXP dataLoc, dataVal;
	int numCols;
	int numInts=0;
	int numReals=0;
		
	// PARSE MxData Structure
	if(OMX_DEBUG) {Rprintf("Processing Data Type.\n");}
	PROTECT(dataLoc = GET_SLOT(dataObject, install("type")));
	if(dataLoc == NULL) {error("Data has no type.  Sorry.\n");}
	PROTECT(dataVal = STRING_ELT(dataLoc,0));
		strncpy(od->type, CHAR(dataVal), 249);
		od->type[249] = '\0';
	UNPROTECT(2); // dataLoc, dataVec
	if(OMX_DEBUG) {Rprintf("Element is type %s.\n", od->type);}
	
	PROTECT(dataLoc = GET_SLOT(dataObject, install("observed")));
	if(OMX_DEBUG) {Rprintf("Processing Data Elements.\n");}
	if(isFrame(dataLoc)) {
		if(OMX_DEBUG) {Rprintf("Data is a frame.\n");}
		// Process Data Frame into Columns
		od->cols = length(dataLoc);
		if(OMX_DEBUG) {Rprintf("Data has %d columns.\n", od->cols);}
		numCols = od->cols;
		od->columns = (SEXP*) R_alloc(numCols, sizeof(SEXP));
		od->realData = (double**) R_alloc(numCols, sizeof(double*));
		od->intData = (int**) R_alloc(numCols, sizeof(int*));
		od->location = (int*) R_alloc(numCols, sizeof(int));
		for(int j = 0; j < numCols; j++) {
			PROTECT(od->columns[j] = VECTOR_ELT(dataLoc, j));
			if(isFactor(od->columns[j])) {
				if(OMX_DEBUG) {Rprintf("Column %d is a factor.\n", j);}
				od->intData[numInts] = INTEGER(od->columns[j]);
				od->location[j] = ~(numInts++);
			} else if (isInteger(od->columns[j])) {
				if(OMX_DEBUG) {Rprintf("Column %d is an integer.\n", j);}
				od->realData[numReals] = REAL(AS_NUMERIC(od->columns[j])); // May need PROTECTion.
				od->location[j] = (numReals++);
			} else {
				if(OMX_DEBUG) {Rprintf("Column %d is a numeric.\n", j);}
				od->realData[numReals] = REAL(od->columns[j]);
				od->location[j] = (numReals++);
			}
			UNPROTECT(1); // columns[j]
		}
		od->rows = length(VECTOR_ELT(dataLoc, 0));
		if(OMX_DEBUG) {Rprintf("And %d rows.\n", od->rows);}
	} else {
		if(OMX_DEBUG) {Rprintf("Data contains a matrix.\n");}
		od->dataMat = omxNewMatrixFromMxMatrix(dataLoc, od->currentState);
		od->cols = od->dataMat->cols;
		od->rows = od->dataMat->rows;
	}	
	UNPROTECT(1); // dataLoc
	
	if(OMX_DEBUG) {Rprintf("Processing Means Matrix.\n");}
	PROTECT(dataLoc = GET_SLOT(dataObject, install("means")));
		od->meansMat = omxNewMatrixFromMxMatrix(dataLoc, od->currentState);
	UNPROTECT(1); // dataLoc
	
	if(strncmp(od->type, "raw", 3) != 0) {
		if(OMX_DEBUG) {Rprintf("Processing Observation Count.\n");}
		PROTECT(dataLoc = GET_SLOT(dataObject, install("numObs")));
		od->numObs = REAL(dataLoc)[0];
		UNPROTECT(1); // dataLoc
	} else {
		od->numObs = od->rows;
	}
	
	return od;
}

void omxFreeData(omxData* od) {
	if(od->dataMat != NULL) omxFreeAllMatrixData(od->dataMat);
	if(od->meansMat != NULL) omxFreeAllMatrixData(od->meansMat);
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
		return (int) omxMatrixElement(od->dataMat, row, col);
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
	char errstr[250];
	sprintf(errstr, "Attempted to access column %d of a %d-column data object.\n", col, od->cols);
	omxRaiseError(od->currentState, -1, errstr);
    return FALSE;
}

omxMatrix* omxDataMeans(omxData *od, omxMatrix* colList, omxMatrix* om) {
	
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

omxMatrix* omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om) {

	if(colList == NULL || row > od->rows) return NULL;	// Sanity check
	
	if(om == NULL) {
		om = omxInitMatrix(om, 1, od->cols, TRUE, od->currentState);
	}

	if(od->dataMat != NULL) { // Matrix Object
		for(int j = 0; j < om->cols; j++) {
			omxSetMatrixElement(om, 0, j, omxDoubleDataElement(od, row, omxVectorElement(colList, j)));
		}
	} else {		// Data Frame object 
		double dataElement;
		for(int j = 0; j < om->cols; j++) {
			int location = od->location[(int)omxVectorElement(colList, j)];
			if(location < 0) {
				dataElement = (double) od->intData[~location][row];
			} else {
				dataElement = od->realData[location][row];
			}
			omxSetMatrixElement(om, 0, j, dataElement);
		}
	}
	return om;
}

double omxDataNumObs(omxData *od) {
	return od->numObs;
}

char* omxDataType(omxData *od) {
	return od->type;
}

void omxPrintData(omxData *source, char* d) {
	Rprintf("NYI: Data Printing not yet implemented.\n");
}

