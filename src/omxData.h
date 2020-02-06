/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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

class omxData;
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

union dataPtr {
	double *realData;
	int *intData;
	dataPtr(double *_p) : realData(_p) {};
	dataPtr(int *_p) : intData(_p) {};
	void clear() { realData=0; intData=0; };
};

struct ColumnData {
	const char *name;
	ColumnDataType type;
	dataPtr ptr;
	std::vector<std::string> levels;       // factors only

	const char *typeName();
};

typedef Eigen::Matrix<int, Eigen::Dynamic, 1> DataColumnIndexVector;

struct WLSVarData {
	Eigen::ArrayXd theta;
	// OLS
	Eigen::ArrayXd resid;
	int contOffset;
	int thrOffset;
};

struct cstrCmp {
	bool operator() (const char *s1, const char *s2) const
	{ return strcmp(s1,s2) < 0; }
};

typedef std::map< const char *, int, cstrCmp > ColMapType;

class obsSummaryStats {
 public:
	std::vector<const char *> dc;
	std::vector<int> exoPred;
	Eigen::ArrayXXi exoFree; // observed by exo matrix
	int totalExoFree;
	const char *wlsType;
	const char *continuousType;
	bool wantFullWeight;
	std::vector<int> index; // rowMult.rows() == index.size()
	Eigen::ArrayXd rowMult;

	bool partial;
	bool output;
	double totalWeight;
	int numOrdinal;  // == thresholdMat->cols or 0 if null
	int numContinuous;
	omxMatrix* covMat;
	omxMatrix *slopeMat; // manifest by exo predictor matrix
	omxMatrix* meansMat;
	omxMatrix* acovMat;			// The asymptotic covariance
	omxMatrix *fullWeight;
	omxMatrix* thresholdMat;
	std::vector< omxThresholdColumn > thresholdCols; // size() == covMat.cols()
	ColMapType colMap;

	// prep
	std::vector<int> contMap;
	int totalThr;
	std::vector<int> thStart;
	std::vector< WLSVarData > perVar;
	Eigen::ArrayXXd SC_VAR;
	Eigen::ArrayXXd SC_SL;
	Eigen::ArrayXXd SC_TH;
	Eigen::ArrayXXd SC_COR;
	Eigen::MatrixXd A21;
	Eigen::ArrayXXd H22;
	Eigen::ArrayXXd H21;

 	obsSummaryStats() :
		wlsType(0), continuousType(0), wantFullWeight(true),
		partial(false), output(false), totalWeight(0), numOrdinal(0), numContinuous(0),
		covMat(0), slopeMat(0), meansMat(0),
		acovMat(0), fullWeight(0), thresholdMat(0), totalThr(0) {};
	~obsSummaryStats();
	void setDimnames(omxData *data);
	void permute(omxData *data);
	void log();
	void loadExoFree(SEXP Ref);
};

class omxData {
 private:
	void addDynamicDataSource(omxExpectation *ex);
	int primaryKey;   // column of primary key
	int weightCol;
	double *currentWeightColumn;
	int freqCol;
	int *currentFreqColumn;
	std::unique_ptr< obsSummaryStats > oss;
	bool parallel;
	bool noExoOptimize;
	bool modified;
	double minVariance;
	bool warnNPDacov;
	std::vector<int> algebra;

	void estimateObservedStats();
	void _prepObsStats(omxState *state, const std::vector<const char *> &dc,
			   std::vector<int> &exoPred, const char *type,
			  const char *continuousType, bool fullWeight);
	bool regenObsStats(const std::vector<const char *> &dc, const char *wlsType);
	void wlsAllContinuousCumulants(omxState *state);
	void convertToDataFrame();

 public:
	bool hasPrimaryKey() const { return primaryKey >= 0; };
	bool hasWeight() const { return weightCol >= 0; };
	bool hasFreq() const { return freqCol >= 0 || currentFreqColumn; };
	int lookupRowOfKey(int key);
	int primaryKeyOfRow(int row);
	void omxPrintData(const char *header, int maxRows, int *permute);
	void omxPrintData(const char *header, int maxRows);
	void omxPrintData(const char *header);
	void assertColumnIsData(int col);
	void setModified() { modified=true; };
	bool isModified() { return modified; };
	double getMinVariance() const { return minVariance; };
	void evalAlgebras(FitContext *fc);

	const char *name;
	SEXP dataObject;                                // only used for dynamic data
	omxMatrix* dataMat;                             // do not use directly
	omxMatrix* meansMat;				// The means, as an omxMatrixObject
	double numObs;						// Number of observations (sum of rowWeight)
	const char *_type;
	const char *getType() const { return _type; };

	// type=="raw"
	struct RawData {
		std::vector<ColumnData> rawCols;
		std::vector<bool> hasNa;
		int rows;
		bool owner;
		RawData() : rows(0), owner(false) {}
		void clear();
		void clearColumn(int col);
		~RawData();
		void refreshHasNa();
		void assertColumnIsData(int col);
	};
	RawData filtered;
	RawData unfiltered;
	enum NaActionType { NA_PASS, NA_FAIL, NA_OMIT, NA_EXCLUDE };
	NaActionType naAction;
	ColMapType rawColMap;
	int numFactor, numNumeric;			// Number of ordinal and continuous columns
	bool needSort;

	SEXP owner;	// The R object owning data or NULL if we own it.

	std::vector<omxDefinitionVar> defVars;
 public:
	void prep();
	void sanityCheck();
	RawData &getUnfilteredRawData() { return unfiltered; }
	bool isRaw() const { return filtered.rawCols.size() != 0; }
	ColumnData &rawCol(int cx) { return filtered.rawCols[cx]; }
	std::vector<ColumnData> &getRawCols() { return filtered.rawCols; }
	int nrows() const { return filtered.rows; }
	int cols;						// dataMat->cols or rawCols.size()
	int verbose;
	std::map<int,int> primaryKeyIndex;

	void loadFakeData(omxState *state, double fake);
	bool hasDefinitionVariables() { return defVars.size() != 0; };
	bool loadDefVars(omxState *state, int row); // prefer omxExpectation member fn

	// Used when the expectation provides the observed data (DataDynamic)
	std::vector<class omxExpectation *> expectation;   // weak pointers
	int version;

	omxData();
	~omxData();
	void reportResults(MxRList &out);
	void newDataStatic(omxState *, SEXP dataObject);
	void connectDynamicData(omxState *currentState);
	void recompute();
	friend void omxDataKeysCompatible(omxData *upper, omxData *lower, int foreignKey);

	double *getWeightColumn();
	double getRowWeight(int row) {
		if (!hasWeight()) return 1.0;
		return getWeightColumn()[row];
	}
	int *getFreqColumn() { return currentFreqColumn; };
	void setFreqColumn(int *wc) { currentFreqColumn = wc; }
	int *getOriginalFreqColumn();
	int getRowFreq(int row) {
		if (!hasFreq()) return 1;
		return getFreqColumn()[row];
	}
	int numRawRows();
	double rowMultiplier(int rx);
	bool containsNAs(int col);
	void prohibitFactor(int col);
	void prohibitNAdefVar(int col);
	void freeInternal();
	bool isDynamic() { return expectation.size() != 0; };
	template <typename T> void visitObsStats(T visitor) {
		visitor(*oss);
	}
	obsSummaryStats &getSingleObsSummaryStats() {
		if (!oss) stop("No observed summary stats");
		return *oss;
	};
	const char *columnName(int col);
	bool columnIsFactor(int col);
	bool hasSummaryStats() { return dataMat != 0 || oss; }
	void prepObsStats(omxState *state, const std::vector<const char *> &dc,
			  std::vector<int> &exoPred, const char *type,
			  const char *continuousType, bool fullWeight);
	template <typename T1>
	void recalcRowWeights(Eigen::ArrayBase<T1> &rowMult, std::vector<int> &index);
	void invalidateCache();
	void invalidateColumnsCache(const std::vector< int > &columns);
	// When FALSE, PolychoricCor uses raw data.
	// When TRUE, PolychoricCor uses summary data (if possible).
	// Both should obtain the same result when no exogenous covariates.
	bool getNoExoOptimize() const { return noExoOptimize; };
};

omxData* omxNewDataFromMxData(SEXP dataObject, const char *name);

omxData* omxDataLookupFromState(SEXP dataObject, omxState* state);	// Retrieves a data object from the state
void omxFreeData(omxData* od);					// Release any held data.

/* Getters 'n Setters */
int omxDataGetNumFactorLevels(omxData *od, int col);
double omxDoubleDataElement(omxData *od, int row, int col);
double *omxDoubleDataColumn(omxData *od, int col);
int omxIntDataElement(omxData *od, int row, int col);						// Returns one data object as an integer

bool omxDataElementMissing(omxData *od, int row, int col);

inline int omxKeyDataElement(omxData *od, int row, int col)
{
	ColumnData &cd = od->rawCol(col);
	return cd.ptr.intData[row];
}

omxMatrix* omxDataCovariance(omxData *od);
omxMatrix* omxDataMeans(omxData *od);

void omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om);// Populates a matrix with a single data row

template <typename T>
void omxDataRow(omxData *od, int row, const Eigen::MatrixBase<T> &colList, omxMatrix* om)
{
	if (row >= od->nrows()) stop("Invalid row");

	if(om == NULL) stop("Must provide an output matrix");

	int numcols = colList.size();
	if(od->isRaw()) {
		for(int j = 0; j < numcols; j++) {
			omxSetMatrixElement(om, 0, j, omxDoubleDataElement(od, row, colList[j]));
		}
	} else {
		omxMatrix* dataMat = od->dataMat;
		for(int j = 0; j < numcols; j++) {
			omxSetMatrixElement(om, 0, j, omxMatrixElement(dataMat, row, colList[j]));
		}
	}
};

static OMXINLINE int
omxIntDataElementUnsafe(omxData *od, int row, int col)
{
	return od->rawCol(col).ptr.intData[row];
}

static OMXINLINE int *omxIntDataColumnUnsafe(omxData *od, int col)
{
	return od->rawCol(col).ptr.intData;
}

double omxDataNumObs(omxData *od);											// Returns number of obs in the dataset
bool omxDataColumnIsKey(omxData *od, int col);
const char *omxDataColumnName(omxData *od, int col);
const char *omxDataType(omxData *od);			      // TODO: convert to ENUM
	
int omxDataNumNumeric(omxData *od);                   // Number of numeric columns in the data set
int omxDataNumFactor(omxData *od);                    // Number of factor columns in the data set

/* Function wrappers that switch based on inclusion of algebras */

double omxDataDF(omxData *od);

inline bool omxDataColumnIsFactor(omxData *od, int col) { return od->columnIsFactor(col); }

// Should not be stored long-term because freq are updated by bootstrap
template <typename T1>
void omxData::recalcRowWeights(Eigen::ArrayBase<T1> &rowMult, std::vector<int> &index)
{
	int rows = nrows();
	index.clear();
	index.reserve(rows);
	rowMult.derived().resize(rows);
	double *rowWeight = getWeightColumn();
	int *rowFreq = getFreqColumn();
	for (int rx=0; rx < rows; ++rx) {
		double ww = 1.0;
		if (rowWeight) ww *= rowWeight[rx];
		if (rowFreq) ww *= rowFreq[rx];
		rowMult[rx] = ww;
		if (ww == 0.0) continue;
		index.push_back(rx);
	}
}

#endif /* _OMXDATA_H_ */
