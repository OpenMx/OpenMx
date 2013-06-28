/*
 *  Copyright 2007-2013 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */

/***********************************************************
 * 
 *  omxData.h
 *
 *  Created: Timothy R. Brick 	Date: 2009-07-15
 *
 *	Contains header information for the omxData class
 *   omxData objects keep data objects in whatever form
 *   they might take.
 *
 **********************************************************/

#ifndef _OMXDATA_H_
#define _OMXDATA_H_

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h> 
#include "omxDefines.h"

typedef struct omxData omxData;
typedef struct omxContiguousData omxContiguousData;
typedef struct omxThresholdColumn omxThresholdColumn;


#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxState.h"

struct omxContiguousData {
	int isContiguous;
	int start;
	int length;
};

struct omxThresholdColumn {		 	// Threshold

	omxMatrix* matrix;		// Which Matrix/Algebra it comes from
	int column;				// Which column has the thresholds
	int numThresholds;		// And how many thresholds

};

struct omxData {						// A matrix
	//TODO: Improve encapsulation
	omxMatrix* dataMat;  // DO NOT ACCESS THIS FIELD DIRECTLY; USE ACCESSOR METHODS
	omxMatrix* meansMat;				// The means, as an omxMatrixObject
	omxMatrix* acovMat;					// The asymptotic covariance, as an omxMatrixObject, added for ordinal WLS
	omxMatrix* obsThresholdsMat;		// The observed thresholds, added for ordinal WLS
	omxThresholdColumn* thresholdCols;  // Wrapper structure for thresholds
	double numObs;						// Number of observations
	const char *_type;
	
	int* location;						// Which of the following contains the data column
	double** realData;					// The actual data objects if numeric
	int** intData;						// The actual data objects if ordinal
	int* indexVector;						// The "official" index into the data set
	int* identicalDefs;					// Number of consecutive rows with identical def. vars
	int* identicalMissingness;			// Number of consecutive rows with identical missingness patterns
	int* identicalRows;					// Number of consecutive rows with identical data
	int numFactor, numNumeric;  // Number of ordinal and continuous columns
	short isDynamic;             // Reserved for when there's actually dynamic data.
	/* Useful Members */
	int rows, cols;						// Matrix size 

	omxState* currentState;				// The Current Optimizer State 	// Might not want this.
};

/* Initialize and Destroy */
omxData* omxInitData(omxState* os);
omxData* omxNewDataFromMxData(SEXP dataObject, omxState* state);
omxData* omxDataLookupFromState(SEXP dataObject, omxState* state);	// Retrieves a data object from the state
void omxFreeData(omxData* od);					// Release any held data.
void omxSetContiguousDataColumns(omxContiguousData* contiguous, omxData* data, omxMatrix* colList);

/* Getters 'n Setters */
double omxDoubleDataElement(omxData *od, int row, int col);					// Returns one data object as a double
int omxIntDataElement(omxData *od, int row, int col);						// Returns one data object as an integer
omxMatrix* omxDataMatrix(omxData *od, omxMatrix* om);						// Populates a matrix with the data (use for covariance matrices)
omxMatrix* omxDataMeans(omxData *od, omxMatrix* colList, omxMatrix* om);	// Populates a matrix with data means
omxMatrix* omxDataAcov(omxData *od, omxMatrix* om); //TODO check this definition
omxThresholdColumn* omxDataThresholds(omxData *od); //TODO define this
// PROBLEM: omxThresholdColumn does not exist because data are prior to expectations.
void omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om);// Populates a matrix with a single data row
void omxContiguousDataRow(omxData *od, int row, int start, int length, omxMatrix* om);// Populates a matrix with a contiguous data row
int omxDataIndex(omxData *od, int row);										// Returns the unsorted (original) index of the current row
int omxDataNumIdenticalRows(omxData *od, int row);							// Returns the number of rows identical to this one in the data set
int omxDataNumIdenticalMissingness(omxData *od, int row);					// Returns the number of rows with definition variables and missingness identical to this one in the data set
int omxDataNumIdenticalContinuousRows(omxData *od, int row);                // Number of rows with continuous variables remaining, or Inf if no continous vars
int omxDataNumIdenticalContinuousMissingness(omxData *od, int row);         // Number of rows with continuous variables remaining, or Inf if no continous vars
int omxDataNumIdenticalOrdinalRows(omxData *od, int row);
int omxDataNumIdenticalOrdinalMissingness(omxData *od, int row);
int omxDataNumIdenticalDefs(omxData *od, int row);							// Returns the number of rows with definition variables identical to this one in the data set

static int OMXINLINE
omxIntDataElementUnsafe(omxData *od, int row, int col)
{
	// All columns must be factors for this to work. You have been warned.
	//	int location = od->location[col];
	return od->intData[col][row];
}
	
double omxDataNumObs(omxData *od);											// Returns number of obs in the dataset
unsigned short int omxDataColumnIsFactor(omxData *od, int col);				// Returns 0 if column is not a factor
const char *omxDataType(omxData *od);			      // TODO: Should this be an ENUM?
	
int omxDataNumNumeric(omxData *od);                   // Number of numeric columns in the data set
int omxDataNumFactor(omxData *od);                    // Number of factor columns in the data set

void resetDefinitionVariables(double *oldDefs, int numDefs);

/* Function wrappers that switch based on inclusion of algebras */

void omxPrintData(omxData *od, const char *header);

#endif /* _OMXDATA_H_ */
