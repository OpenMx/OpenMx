/*
 *  Copyright 2007-2017 The OpenMx Project
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

#include "omxDefines.h"
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h> 

typedef struct omxData omxData;
typedef struct omxContiguousData omxContiguousData;
typedef struct omxThresholdColumn omxThresholdColumn;

#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxState.h"

#include <map>

struct omxDefinitionVar {
	int column;		// in data
	int row, col;           // in matrix
	int matrix;		// matrix number
	int  numDeps;           // number of algebra/matrix dependencies
	int* deps;              // indices of algebra/matrix dependencies

	bool loadData(omxState *state, double val);
};

struct omxContiguousData {
	int isContiguous;
	int start;
	int length;
};

struct omxThresholdColumn {
	int dColumn;			// Which column in the matrix/data.frame
	int column;			// Which column in the thresholds matrix
	int numThresholds;		// And how many thresholds

	// for continuous variables, numThresholds=0
	omxThresholdColumn() : dColumn(-1), column(0), numThresholds(0) {};

	void log() { mxLog("dCol=%d col=%d #thr=%d", dColumn, column, numThresholds); }
};

enum ColumnDataType {
	COLUMNDATA_INVALID,
	COLUMNDATA_ORDERED_FACTOR,
	COLUMNDATA_UNORDERED_FACTOR,
	COLUMNDATA_INTEGER,
	COLUMNDATA_NUMERIC
};

struct ColumnData {
	const char *name;
	ColumnDataType type;
	// exactly one of these is non-null
	double *realData;
	int    *intData;
	SEXP levels;       // factors only
};

class omxData {
 private:
	SEXP rownames;
	void addDynamicDataSource(omxExpectation *ex);
	int primaryKey;   // column of primary key
	int weightCol;
	double *currentWeightColumn;

 public: // move everything to private TODO
	bool hasPrimaryKey() const { return primaryKey >= 0; };
	bool hasWeight() const { return weightCol >= 0; };
	int lookupRowOfKey(int key);
	int primaryKeyOfRow(int row);
	void omxPrintData(const char *header, int maxRows, int *permute);
	void omxPrintData(const char *header, int maxRows);
	void omxPrintData(const char *header);
	void assertColumnIsData(int col);

	const char *name;
	SEXP dataObject;                                // only used for dynamic data
	omxMatrix* dataMat;                             // do not use directly
	omxMatrix* meansMat;				// The means, as an omxMatrixObject
	omxMatrix* acovMat;					// The asymptotic covariance, as an omxMatrixObject, added for ordinal WLS
	omxMatrix* obsThresholdsMat;		// The observed thresholds, added for ordinal WLS
	std::vector< omxThresholdColumn > thresholdCols;
	double numObs;						// Number of observations (sum of rowWeight)
	const char *_type;
	const char *getType() const { return _type; };

	// type=="raw"
	std::vector<ColumnData> rawCols;
	int numFactor, numNumeric;			// Number of ordinal and continuous columns
	bool needSort;

	std::vector<omxDefinitionVar> defVars;
 public:
	int rows, cols;						// Matrix size 
	int verbose;
	std::map<int,int> primaryKeyIndex;

	void loadFakeData(omxState *state, double fake);
	bool hasDefinitionVariables() { return defVars.size() != 0; };
	bool CompareDefVarInMatrix(int lrow, int rrow, omxMatrix *mat, bool &mismatch);
	bool loadDefVars(omxState *state, int row); // prefer omxExpectation member fn

	// Used when the expectation provides the observed data (DataDynamic)
	std::vector<struct omxExpectation *> expectation;   // weak pointers
	int version;

	omxData();
	void newDataStatic(omxState *, SEXP dataObject);
	void connectDynamicData(omxState *currentState);
	void recompute();
	friend void omxDataKeysCompatible(omxData *upper, omxData *lower, int foreignKey);
	double *getWeightColumn() { return currentWeightColumn; };
	void setWeightColumn(double *wc) { currentWeightColumn = wc; }
	double *getOriginalWeightColumn();
	double getRowWeight(int row) {
		if (!hasWeight()) return 1.0;
		return getWeightColumn()[row];
	}
	int numRawRows();
};

omxData* omxNewDataFromMxData(SEXP dataObject, const char *name);

omxData* omxDataLookupFromState(SEXP dataObject, omxState* state);	// Retrieves a data object from the state
void omxFreeData(omxData* od);					// Release any held data.

template <typename T>
void omxSetContiguousDataColumns(omxContiguousData* contiguous, omxData* data,
				 Eigen::MatrixBase<T> &colList)
{
	contiguous->isContiguous = FALSE;   // Assume not contiguous

	if (data->dataMat == NULL) return; // Data has no matrix elements, so skip.

	omxMatrix* dataMat = data->dataMat;
	if (dataMat->colMajor) return;      // If data matrix is column-major, there's no continuity
	
	int colListLength = colList.size();             // # of columns in the cov matrix
	double start = colList[0];                      // Data col of first column of the covariance
	contiguous->start = (int) start;                // That's our starting point.
	contiguous->length = colListLength;             // And the length is ncol(cov)
	
	for(int i = 1; i < colListLength; i++) {        // Make sure that the col list is 
		double next = colList[i];               // contiguously increasing in column number
		if (next != (start + i)) return;            // If it isn't, it's not contiguous data
	}
	
	contiguous->isContiguous = TRUE;    // Passed.  This is contiguous.
}

/* Getters 'n Setters */
int omxDataGetNumFactorLevels(omxData *od, int col);
double omxDoubleDataElement(omxData *od, int row, int col);
double *omxDoubleDataColumn(omxData *od, int col);
int omxIntDataElement(omxData *od, int row, int col);						// Returns one data object as an integer

bool omxDataElementMissing(omxData *od, int row, int col);

inline int omxKeyDataElement(omxData *od, int row, int col)
{
	ColumnData &cd = od->rawCols[col];
	return cd.intData[row];
}

inline bool omxDataHasMatrix(omxData *od) { return od->dataMat != 0; }
omxMatrix* omxDataCovariance(omxData *od);
omxMatrix* omxDataMeans(omxData *od);
omxMatrix* omxDataAcov(omxData *od);

std::vector<omxThresholdColumn> &omxDataThresholds(omxData *od);

void omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om);// Populates a matrix with a single data row

template <typename T>
void omxDataRow(omxData *od, int row, const Eigen::MatrixBase<T> &colList, omxMatrix* om)
{
	if (row >= od->rows) Rf_error("Invalid row");

	if(om == NULL) Rf_error("Must provide an output matrix");

	int numcols = colList.size();
	if(od->dataMat != NULL) {
		omxMatrix* dataMat = od->dataMat;
		for(int j = 0; j < numcols; j++) {
			omxSetMatrixElement(om, 0, j, omxMatrixElement(dataMat, row, colList[j]));
		}
	} else {
		for(int j = 0; j < numcols; j++) {
			omxSetMatrixElement(om, 0, j, omxDoubleDataElement(od, row, colList[j]));
		}
	}
};

void omxContiguousDataRow(omxData *od, int row, int start, int length, omxMatrix* om);// Populates a matrix with a contiguous data row

static OMXINLINE int
omxIntDataElementUnsafe(omxData *od, int row, int col)
{
	return od->rawCols[col].intData[row];
}

static OMXINLINE int *omxIntDataColumnUnsafe(omxData *od, int col)
{
	return od->rawCols[col].intData;
}

double omxDataNumObs(omxData *od);											// Returns number of obs in the dataset
bool omxDataColumnIsFactor(omxData *od, int col);
bool omxDataColumnIsKey(omxData *od, int col);
const char *omxDataColumnName(omxData *od, int col);
const char *omxDataType(omxData *od);			      // TODO: convert to ENUM
	
int omxDataNumNumeric(omxData *od);                   // Number of numeric columns in the data set
int omxDataNumFactor(omxData *od);                    // Number of factor columns in the data set

/* Function wrappers that switch based on inclusion of algebras */

double omxDataDF(omxData *od);

#endif /* _OMXDATA_H_ */
