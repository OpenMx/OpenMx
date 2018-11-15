/*
 *  Copyright 2007-2018 by the individuals mentioned in the source code history
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
#include "omxExpectationBA81.h"  // improve encapsulation TODO
#include "matrix.h"
#include <fstream>
#include <Rmath.h>
#include "nr.h"
#include "omxBVN.h"
#include "omxNLopt.h"
#include "EnableWarnings.h"

omxData::omxData() : primaryKey(NA_INTEGER), weightCol(NA_INTEGER), currentWeightColumn(0),
		     freqCol(NA_INTEGER), currentFreqColumn(0), permuted(false), oss(0),
		     wlsType(0), wlsContinuousType(0),
		     dataObject(0), dataMat(0), meansMat(0), 
		     numObs(0), _type(0), numFactor(0), numNumeric(0),
		     rows(0), cols(0), expectation(0)
{}

omxData::~omxData()
{
	if (oss) delete oss;
}

omxData* omxDataLookupFromState(SEXP dataObject, omxState* state) {
	int dataIdx = INTEGER(dataObject)[0];
	if (dataIdx == NA_INTEGER) return NULL;

	return state->dataList[dataIdx];
}

static void newDataDynamic(SEXP dataObject, omxData *od)
{
	SEXP dataLoc;
	ScopedProtect p1(dataLoc, R_do_slot(dataObject, Rf_install("type")));
	od->_type = CHAR(STRING_ELT(dataLoc,0));
	od->dataObject = dataObject;
	if (!strEQ(od->getType(), "cov")) {
		omxRaiseErrorf("Don't know how to create dynamic data with type '%s'", od->getType());
	}
}

void omxData::addDynamicDataSource(omxExpectation *ex)
{
	expectation.push_back(ex);
	ex->dynamicDataSource = true;
}

void omxData::connectDynamicData(omxState *currentState)
{
	if (!dataObject) return;

	if (expectation.size()) {
		Rf_error("omxData::connectDynamicData called more than once");
	}

	SEXP dataLoc;
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("expectation")));
	if (Rf_length(dataLoc) == 0) {
		omxRaiseError("mxDataDynamic is not connected to a data source");
		return;
	}

	if (Rf_length(dataLoc) == 1) {
		omxExpectation *ex = omxExpectationFromIndex(INTEGER(dataLoc)[0], currentState);
		BA81Expect *other = (BA81Expect *) ex;
		numObs = other->freqSum;
		addDynamicDataSource(ex);
		// nothing special to do
	} else {
		int num = Rf_length(dataLoc);
		expectation.reserve(num);
		int *evec = INTEGER(dataLoc);

		omxExpectation *refE = NULL;
		BA81Expect *refBA81;
		double freqSum = 0;

		for (int sx=0; sx < num; ++sx) {
			omxExpectation *ex = omxExpectationFromIndex(evec[sx], currentState);
			if (!strEQ(ex->expType, "MxExpectationBA81")) {
				omxRaiseErrorf("MxDataDynamic: type='cov' is only valid for MxExpectationBA81, not '%s'",
					       ex->expType);
				continue;
			}
			BA81Expect *other = (BA81Expect *) ex;
			freqSum += other->freqSum;
			if (!refE) {
				refE = ex;
				refBA81 = other;
			} else {
				const char *why = refBA81->getLatentIncompatible(other);
				if (why) {
					omxRaiseErrorf("MxDataDynamic: '%s' is not compatible with '%s' because of %s",
						       ex->name, refE->name, why);
					continue;
				}
			}

			addDynamicDataSource(ex);
		}
		numObs = freqSum;
		if (!refE) return;

		int dims = refBA81->grp.quad.abilities();
		dataMat = omxNewIdentityMatrix(dims, currentState);
		meansMat = omxInitMatrix(dims, 1, TRUE, currentState);
		for (int mx=0; mx < dims; ++mx) omxSetVectorElement(meansMat, mx, 0);
		version = 0;
	}
}

void omxData::recompute()
{
	int num = (int) expectation.size();
	if (num <= 1) return;

	int oldVersion = version;
	ba81AggregateDistributions(expectation, &version, meansMat, dataMat);
	if (oldVersion != version && verbose >= 1) {
		mxLog("MxData: recompute %s", name);
		omxPrint(meansMat, "mean");
		omxPrint(dataMat, "cov");
	}
}

static int carefulMinusOne(int val)
{
	return val != NA_INTEGER? val - 1 : val;
}

static void importDataFrame(SEXP dataLoc, std::vector<ColumnData> &rawCols,
			    int &numNumeric, int &numFactor)
{
	int numCols = Rf_length(dataLoc);
	rawCols.clear();
	rawCols.reserve(numCols);
	numNumeric = 0;
	numFactor = 0;
	ProtectedSEXP colnames(Rf_getAttrib(dataLoc, R_NamesSymbol));

	for(int j = 0; j < numCols; j++) {
		const char *colname = CHAR(STRING_ELT(colnames, j));
		ColumnData cd = { colname, COLUMNDATA_INVALID, NULL, NULL, {} };
		ProtectedSEXP rcol(VECTOR_ELT(dataLoc, j));
		if(Rf_isFactor(rcol)) {
			cd.type = Rf_isUnordered(rcol)? COLUMNDATA_UNORDERED_FACTOR : COLUMNDATA_ORDERED_FACTOR;
			if(OMX_DEBUG) {mxLog("Column[%d] %s is a factor.", j, colname);}
			cd.intData = INTEGER(rcol);
			ProtectedSEXP Rlevels(Rf_getAttrib(rcol, R_LevelsSymbol));
			for (int lx=0; lx < Rf_length(Rlevels); ++lx) {
				cd.levels.push_back(R_CHAR(STRING_ELT(Rlevels, lx)));
			}
			numFactor++;
		} else if (Rf_isInteger(rcol)) {
			cd.intData = INTEGER(rcol);
			cd.type = COLUMNDATA_INTEGER;
		} else if (Rf_isNumeric(rcol)) {
			if(OMX_DEBUG) {mxLog("Column[%d] %s is numeric.", j, colname);}
			cd.realData = REAL(rcol);
			cd.type = COLUMNDATA_NUMERIC;
			numNumeric++;
		} else {
			if(OMX_DEBUG) {mxLog("Column[%d] %s is type %s (ignored)",
					     j, colname, Rf_type2char(TYPEOF(rcol)));}
			cd.type = COLUMNDATA_INVALID;
		}
		rawCols.push_back(cd);
	}
}

void omxData::newDataStatic(omxState *state, SEXP dataObj)
{
	owner = dataObj;
	omxData *od = this;
	SEXP dataLoc;

	// PARSE MxData Structure
	if(OMX_DEBUG) {mxLog("Processing Data '%s'", od->name);}

	{
		ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("type")));
		od->_type = CHAR(STRING_ELT(dataLoc,0));
		if(OMX_DEBUG) {mxLog("Element is type %s.", od->_type);}

		ProtectedSEXP needsort(R_do_slot(dataObj, Rf_install(".needSort")));
		od->needSort = Rf_asLogical(needsort);

		ScopedProtect p2(dataLoc, R_do_slot(dataObj, Rf_install("primaryKey")));
		primaryKey = carefulMinusOne(Rf_asInteger(dataLoc));

		ProtectedSEXP Rweight(R_do_slot(dataObj, Rf_install("weight")));
		weightCol = carefulMinusOne(Rf_asInteger(Rweight));

		ProtectedSEXP Rfrequency(R_do_slot(dataObj, Rf_install("frequency")));
		freqCol = carefulMinusOne(Rf_asInteger(Rfrequency));
	}
	{ScopedProtect pdl(dataLoc, R_do_slot(dataObj, Rf_install("observed")));
	if(OMX_DEBUG) {mxLog("Processing Data Elements.");}
	if (Rf_isFrame(dataLoc)) {
		od->rows = Rf_length(VECTOR_ELT(dataLoc, 0));
		od->cols = Rf_length(dataLoc);
		importDataFrame(dataLoc, od->rawCols, od->numNumeric, od->numFactor);
	} else {
		if(OMX_DEBUG) {mxLog("Data contains a matrix.");}
		od->dataMat = omxNewMatrixFromRPrimitive0(dataLoc, state, 0, 0);
		
		if (od->dataMat->colMajor && strEQ(od->_type, "raw")) {
			omxToggleRowColumnMajor(od->dataMat);
		}
		od->cols = od->dataMat->cols;
		od->rows = od->dataMat->rows;
		od->numNumeric = od->cols;
	}
	}

	ProtectedSEXP RwlsType(R_do_slot(dataObj, Rf_install(".wlsType")));
	wlsType = CHAR(STRING_ELT(RwlsType,0));
	ProtectedSEXP RwlsContType(R_do_slot(dataObj, Rf_install(".wlsContinuousType")));
	wlsContinuousType = CHAR(STRING_ELT(RwlsContType,0));

	if (od->hasPrimaryKey() && !(od->rawCols.size() && od->rawCols[primaryKey].intData)) {
		Rf_error("%s: primary key must be an integer or factor column in raw observed data", od->name);
	}

	if (od->hasWeight() && od->rawCols.size() && od->rawCols[weightCol].type != COLUMNDATA_NUMERIC) {
		Rf_error("%s: weight must be a numeric column in raw observed data", od->name);
	}

	if (od->hasFreq() && od->rawCols.size() && od->rawCols[freqCol].type != COLUMNDATA_INTEGER) {
		Rf_error("%s: frequency must be an integer column in raw observed data", od->name);
	}

	if(OMX_DEBUG) {mxLog("Processing Means Matrix.");}
	{ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("means")));
	od->meansMat = omxNewMatrixFromRPrimitive0(dataLoc, state, 0, 0);
	}

	if(OMX_DEBUG) {
	        if(od->meansMat == NULL) {mxLog("No means found.");}
		else {omxPrint(od->meansMat, "Means Matrix is:");}
	}

	if(OMX_DEBUG) {mxLog("Processing Asymptotic Covariance Matrix.");}
	ProtectedSEXP Racov(R_do_slot(dataObj, Rf_install("acov")));
	omxMatrix *acovMat = omxNewMatrixFromRPrimitive0(Racov, state, 0, 0);
	if (acovMat) {
		oss = new obsSummaryStats;
		auto &o1 = *oss;
		o1.covMat = dataMat;
		dataMat = 0;
		o1.meansMat = meansMat;
		meansMat = 0;
		o1.acovMat = acovMat;

		ProtectedSEXP Rfw(R_do_slot(dataObj, Rf_install("fullWeight")));
		o1.fullWeight = omxNewMatrixFromRPrimitive0(Rfw, state, 0, 0);

		if(OMX_DEBUG) {mxLog("Processing Observed Thresholds Matrix.");}
		ProtectedSEXP Rthr(R_do_slot(dataObj, Rf_install("thresholds")));
		o1.thresholdMat = omxNewMatrixFromRPrimitive0(Rthr, state, 0, 0);

		if(o1.thresholdMat) {
			o1.numOrdinal = o1.thresholdMat->cols;
			o1.thresholdCols.reserve(o1.thresholdMat->cols);
			ProtectedSEXP Rtc(R_do_slot(dataObj, Rf_install("thresholdColumns")));
			int *columns = INTEGER(Rtc);
			ProtectedSEXP Rtl(R_do_slot(dataObj, Rf_install("thresholdLevels")));
			int *levels = INTEGER(Rtl);

			for(int i = 0; i < o1.covMat->cols; i++) {
				omxThresholdColumn tc;
				tc.dColumn = i;
				if (levels[i] != NA_INTEGER) {
					tc.column = columns[i];
					tc.numThresholds = levels[i];
				}
				o1.thresholdCols.push_back(tc);
				if(OMX_DEBUG && levels[i] != NA_INTEGER) {
					mxLog("%s: column %d is ordinal with %d thresholds in threshold column %d.", 
					      name, i, levels[i], columns[i]);
				}
			}
		}
	}

	{
		ProtectedSEXP RdataLoc(R_do_slot(dataObj, Rf_install("numObs")));
		od->numObs = Rf_asReal(RdataLoc);
	}

	if (hasPrimaryKey()) {
		for (int rx=0; rx < rows; ++rx) {
			int key = primaryKeyOfRow(rx);
			std::pair< std::map<int,int>::iterator, bool> ret =
				primaryKeyIndex.insert(std::pair<int,int>(key, rx));
			if (!ret.second) {
				Rf_error("%s: primary keys are not unique (examine rows with key=%d)", od->name, key);
			}
		}
	}

	currentWeightColumn = getWeightColumn();
	currentFreqColumn = getOriginalFreqColumn();
	
	if (currentFreqColumn) {
		for (int rx=0; rx < rows; ++rx) {
			if (currentFreqColumn[rx] >= 0) continue;
			Rf_error("%s: cannot proceed with non-positive weight %d for row %d",
				 name, currentFreqColumn[rx], 1+rx);
		}
	}
}

double *omxData::getWeightColumn()
{
	if (!hasWeight()) return 0;
	if (rawCols.size()) {
		return rawCols[weightCol].realData;
	} else {
		if (dataMat->colMajor) {
			return omxMatrixColumn(dataMat, weightCol);
		} else {
			auto *col = (double*) R_alloc(dataMat->rows, sizeof(double));
			EigenMatrixAdaptor dm(dataMat);
			Eigen::Map< Eigen::VectorXd > Ecol(col, dataMat->rows);
			Ecol.derived() = dm.col(weightCol);
			return col;
		}
	}
}

int *omxData::getOriginalFreqColumn()
{
	if (!hasFreq()) return 0;
	if (rawCols.size()) {
		return rawCols[freqCol].intData;
	} else {
		auto *col = (int*) R_alloc(dataMat->rows, sizeof(int));
		EigenMatrixAdaptor dm(dataMat);
		Eigen::Map< Eigen::VectorXi > Ecol(col, dataMat->rows);
		Ecol.derived() = dm.col(freqCol).cast<int>();
		return col;
	}
}

int omxData::numRawRows()
{
	if (strEQ(getType(), "raw")) return rows;
	return 0;
}

omxData* omxState::omxNewDataFromMxData(SEXP dataObj, const char *name)
{
	if(dataObj == NULL) {
		Rf_error("Null Data Object detected.  This is an internal Rf_error, and should be reported on the forums.\n");
	}

	const char* dclass;
	ProtectedSEXP DataClass(STRING_ELT(Rf_getAttrib(dataObj, R_ClassSymbol), 0));
	dclass = CHAR(DataClass);
	if(OMX_DEBUG) {mxLog("Initializing %s element", dclass);}
	omxData* od = new omxData();
	od->name = name;
	ProtectedSEXP Rverbose(R_do_slot(dataObj, Rf_install("verbose")));
	od->verbose = Rf_asInteger(Rverbose);
	dataList.push_back(od);
	if (strcmp(dclass, "MxDataStatic")==0) od->newDataStatic(this, dataObj);
	else if (strcmp(dclass, "MxDataDynamic")==0) newDataDynamic(dataObj, od);
	else Rf_error("Unknown data class %s", dclass);
	return od;
}

void omxData::freeInternal()
{
	if (owner) {
		owner = 0;
		for (auto &cd : rawCols) {
			cd.realData = 0;
			cd.intData = 0;
		}
	} else {
		for (auto &cd : rawCols) {
			if (cd.realData) delete [] cd.realData;
			if (cd.intData) delete [] cd.intData;
			cd.realData = 0;
			cd.intData = 0;
		}
	}
}

obsSummaryStats::~obsSummaryStats()
{
	omxFreeMatrix(covMat);
	omxFreeMatrix(meansMat);
	omxFreeMatrix(acovMat);
	if (acovMat != fullWeight) omxFreeMatrix(fullWeight);
	omxFreeMatrix(thresholdMat);
}

void omxFreeData(omxData* od)
{
	od->freeInternal();
	omxFreeMatrix(od->dataMat);
	omxFreeMatrix(od->meansMat);
	delete od;
}

bool omxDataElementMissing(omxData *od, int row, int col)
{
	if(od->dataMat != NULL) {
		return std::isnan(omxMatrixElement(od->dataMat, row, col));
	}
	ColumnData &cd = od->rawCols[col];
	if (cd.realData) return std::isnan(cd.realData[row]);
	else return cd.intData[row] == NA_INTEGER;
}

double omxDoubleDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		return omxMatrixElement(od->dataMat, row, col);
	}
	ColumnData &cd = od->rawCols[col];
	if (cd.realData) return cd.realData[row];
	else return cd.intData[row];
}

double *omxDoubleDataColumn(omxData *od, int col)
{
	ColumnData &cd = od->rawCols[col];
	if (!cd.realData) Rf_error("Column '%s' is integer, not real", cd.name);
	else return cd.realData;
}

int omxDataGetNumFactorLevels(omxData *od, int col)
{
	ColumnData &cd = od->rawCols[col];
	if (cd.levels.size() == 0) Rf_error("omxDataGetNumFactorLevels attempt on non-factor");
	return cd.levels.size();
}

int omxIntDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		return (int) omxMatrixElement(od->dataMat, row, col);
	}

	ColumnData &cd = od->rawCols[col];
	if (cd.realData) return cd.realData[row];
	else return cd.intData[row];
}

omxMatrix* omxDataCovariance(omxData *od)
{
	if (od->dataMat) return od->dataMat;

	if (od->expectation.size()) {
		omxExpectation *ex = od->expectation[0];
		return omxGetExpectationComponent(ex, "covariance");
	}

	Rf_error("%s: type='%s' data must be in matrix storage", od->name, od->_type);
}

bool omxData::columnIsFactor(int col)
{
	if(dataMat != NULL) return FALSE;
	ColumnData &cd = rawCols[col];
	return cd.intData && cd.levels.size();
}

bool omxDataColumnIsKey(omxData *od, int col)
{
	if(od->dataMat != NULL) return FALSE;
	ColumnData &cd = od->rawCols[col];
	return cd.intData != 0;
}

void omxData::assertColumnIsData(int col)
{
	if (dataMat) return;
	ColumnData &cd = rawCols[col];
	switch (cd.type) {
	case COLUMNDATA_ORDERED_FACTOR:
	case COLUMNDATA_NUMERIC:
		return;
	case COLUMNDATA_UNORDERED_FACTOR:
		if (++Global->dataTypeWarningCount < 5) {
			Rf_warning("In data '%s', column '%s' must be an ordered factor. "
				   "Please use mxFactor()", name, cd.name);
		}
		return;
	case COLUMNDATA_INTEGER:
		cd.type = COLUMNDATA_NUMERIC;
		cd.realData = (double*) R_alloc(rows, sizeof(double));
		for (int rx=0; rx < rows; ++rx) {
			if (cd.intData[rx] == NA_INTEGER) {
				cd.realData[rx] = NA_REAL;
			} else {
				cd.realData[rx] = cd.intData[rx];
			}
		}
		cd.intData = 0;
		return;
	default:
		Rf_error("In data '%s', column '%s' is an unknown data type", name, cd.name);
	}
}

int omxData::primaryKeyOfRow(int row)
{
	if(dataMat != NULL) Rf_error("%s: only raw data can have a primary key", name);
	ColumnData &cd = rawCols[primaryKey];
	return cd.intData[row];
}

int omxData::lookupRowOfKey(int key)
{
	const std::map<int,int>::iterator it = primaryKeyIndex.find(key);
	if (it == primaryKeyIndex.end()) {
		if (!hasPrimaryKey()) {
			Rf_error("%s: attempt to lookup key=%d but no primary key", name, key);
		}
		ColumnData &cd = rawCols[primaryKey];
		Rf_error("%s: key %d not found in column '%s'", name, key, cd.name);
	}
	return it->second;
}

const char *omxDataColumnName(omxData *od, int col)
{
	return od->columnName(col);
}

const char *omxData::columnName(int col)
{
	if(dataMat) {
		auto &cn = dataMat->colnames;
		if (col < int(cn.size())) return cn[col];
		else return "?";
	}
	if (oss) {
		auto &o1 = *oss;
		auto &cn = o1.covMat->colnames;
		if (col < int(cn.size())) return cn[col];
		else return "?";
	}
	ColumnData &cd = rawCols[col];
	return cd.name;
}

static const char *ColumnDataTypeToString(enum ColumnDataType cdt)
{
	switch (cdt) {
	case COLUMNDATA_INVALID: return "invalid";
	case COLUMNDATA_ORDERED_FACTOR: return "ordered factor";
	case COLUMNDATA_UNORDERED_FACTOR: return "unordered factor";
	case COLUMNDATA_INTEGER: return "integer";
	case COLUMNDATA_NUMERIC: return "real";
	default: Rf_error("type %d unknown", cdt);
	}
}

void omxDataKeysCompatible(omxData *upper, omxData *lower, int foreignKey)
{
	ColumnData &lcd = lower->rawCols[foreignKey];
	if (!upper->hasPrimaryKey()) {
		Rf_error("Attempt to join foreign key '%s' in %s of type '%s' with"
			 " %s which has no primary key declared",
			 lcd.name, lower->name, ColumnDataTypeToString(lcd.type), upper->name);
	}
	ColumnData &ucd = upper->rawCols[upper->primaryKey];
	if (ucd.type != lcd.type) {
		Rf_error("Primary key '%s' in %s of type '%s' cannot be joined with"
			 " foreign key '%s' in %s of type '%s'",
			 ucd.name, upper->name, ColumnDataTypeToString(ucd.type),
			 lcd.name, lower->name, ColumnDataTypeToString(lcd.type));
	}
	if (ucd.type == COLUMNDATA_ORDERED_FACTOR || ucd.type == COLUMNDATA_UNORDERED_FACTOR) {
		if (ucd.levels.size() != lcd.levels.size()) {
			Rf_error("Primary key '%s' in %s has a different number of factor"
				 " levels compared to foreign key '%s' in %s",
				 ucd.name, upper->name, lcd.name, lower->name);
		}
		for (int lx=0; lx < int(ucd.levels.size()); ++lx) {
			auto &ul = ucd.levels[lx];
			auto &ll = lcd.levels[lx];
			if (ul == ll) continue;
			Rf_error("Primary key '%s' in %s has different factor levels ('%s' != '%s')"
				 " compared to foreign key '%s' in %s",
				 ucd.name, upper->name, ul.c_str(), ll.c_str(), lcd.name, lower->name);
		}
	}
}

omxMatrix* omxDataMeans(omxData *od)
{
	if (od->meansMat) return od->meansMat;
	if (od->expectation.size()) {
		omxExpectation *ex = od->expectation[0];
		omxMatrix *mat = omxGetExpectationComponent(ex, "mean");
		if (!mat) return NULL;
		if (mat->rows != 1) omxTransposeMatrix(mat);
		return mat;
	}
	return NULL;
}

void omxContiguousDataRow(omxData *od, int row, int start, int len, omxMatrix* om) {
	// TODO: Might be better to combine this with omxDataRow to make a single accessor omxDataRow with a second signature that accepts an omxContiguousData argument.
	if(row >= od->rows) Rf_error("Invalid row");

	if(om == NULL) Rf_error("Must provide an output matrix");
	
	if (om->cols < len) Rf_error("omxContiguousDataRow: output matrix is too small");
	int numcols = od->cols;
	omxMatrix* dataMat = od->dataMat;
	double *dest = om->data;
	double *source = dataMat->data + row * numcols + start;
	memcpy(dest, source, sizeof(double) * len);
}

void omxDataRow(omxData *od, int row, omxMatrix* colList, omxMatrix* om)
{
	EigenVectorAdaptor ecl(colList);
	omxDataRow(od, row, ecl, om);
}

double omxDataNumObs(omxData *od)
{
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

void omxData::omxPrintData(const char *header, int maxRows, int *permute)
{
	if (!header) header = "Default data";

	omxData *od = this;
	if (!od) {
		mxLog("%s: NULL", header);
		return;
	}

	std::string buf;
	buf += string_snprintf("%s(%s): %f observations %d x %d\n", header, od->_type, od->numObs,
			       od->rows, od->cols);
	if (hasPrimaryKey()) {
		buf += string_snprintf("primaryKey %d\n", od->primaryKey);
	}

	buf += string_snprintf("Row consists of %d numeric, %d ordered factor:", od->numNumeric, od->numFactor);

        int upto = od->numObs;
        if (maxRows >= 0 && maxRows < upto) upto = maxRows;

	if (od->rawCols.size()) {
		for (auto &cd : od->rawCols) {
			buf += " ";
			buf += cd.name;
			if (cd.intData) {
				buf += "[I]";
			} else {
				buf += "[N]";
			}
		}
		buf += "\n";

		for (int vxv=0; vxv < upto; vxv++) {
			int vx = permute? permute[vxv] : vxv;
			for (int j = 0; j < od->cols; j++) {
				ColumnData &cd = od->rawCols[j];
				if (cd.intData) {
					int *val = cd.intData;
					if (val[vx] == NA_INTEGER) {
						buf += " NA,";
					} else {
						buf += string_snprintf(" %d,", val[vx]);
					}
				} else {
					double *val = cd.realData;
					if (!std::isfinite(val[vx])) {
						buf += " NA,";
					} else {
						buf += string_snprintf(" %.3g,", val[vx]);
					}
				}
			}
			buf += string_snprintf("\t# %d \n", vxv);
		}
	}

	mxLogBig(buf);

	if (od->dataMat) omxPrintMatrix(od->dataMat, "dataMat");
	if (od->meansMat) omxPrintMatrix(od->meansMat, "meansMat");
}

void omxData::omxPrintData(const char *header)
{
        omxPrintData(header, -1, 0);
}

void omxData::omxPrintData(const char *header, int maxRows)
{
        omxPrintData(header, maxRows, 0);
}

double omxDataDF(omxData *od)
{
	const char *type = od->_type;
	if (strEQ(type, "cov") || strEQ(type, "sscp")) {
		omxMatrix *cov = omxDataCovariance(od);
		int df = triangleLoc1(cov->rows);
		omxMatrix *mm = omxDataMeans(od);
		if (mm) df += mm->rows * mm->cols;
		return df;
	} else if (strEQ(type, "cor")) {
		omxMatrix *cov = omxDataCovariance(od);
		int df = triangleLoc1(cov->rows - 1);
		omxMatrix *mm = omxDataMeans(od);
		if (mm) df += mm->rows * mm->cols;
		return df;
	}
	return NA_REAL;
}

static void markDefVarDependencies(omxState* os, omxDefinitionVar* defVar)
{
	int numDeps = defVar->numDeps;
	int *deps = defVar->deps;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(os->matrixList[~value]);
		} else {
			omxMarkDirty(os->algebraList[value]);
		}
	}
}

bool omxData::loadDefVars(omxState *state, int row)
{
	bool changed = false;
	for (int k=0; k < int(defVars.size()); ++k) {
		double newDefVar;
		if (defVars[k].column == weightCol) {
			newDefVar = getRowWeight(row);
		} else if (defVars[k].column == freqCol) {
			newDefVar = getRowFreq(row);
		} else {
			newDefVar = omxDoubleDataElement(this, row, defVars[k].column);
		}
		changed |= defVars[k].loadData(state, newDefVar);
	}
	if (changed && OMX_DEBUG_ROWS(row)) { mxLog("Processing Definition Vars for row %d", row); }
	return changed;
}

void omxData::prohibitNAs(int col)
{
	if(dataMat != NULL) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(omxMatrixElement(dataMat, rx, col))) continue;
			Rf_error("%s: NA in definition variable '%s' row %d",
				 name, omxDataColumnName(this, col), 1+rx);
		}
		return;
	}
	if (col == weightCol) {
		double *wc = getWeightColumn();
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(wc[rx])) continue;
			Rf_error("%s: NA in row weight %d", name, 1+rx);
		}
		return;
	}
	if (col == freqCol) {
		int *wc = getFreqColumn();
		for (int rx=0; rx < rows; ++rx) {
			if (wc[rx] != NA_INTEGER) continue;
			Rf_error("%s: NA in row frequency %d", name, 1+rx);
		}
		return;
	}
	ColumnData &cd = rawCols[col];
	if (cd.realData) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(cd.realData[rx])) continue;
			Rf_error("%s: NA in definition variable '%s' row %d",
				 name, omxDataColumnName(this, col), 1+rx);
		}
	} else {
		if (cd.type != COLUMNDATA_INTEGER) {
			Rf_warning("%s: definition variable '%s' is a factor;"
				   " note that it will be treated as integer (as is done by ?unclass)."
				   " Is this really what you want to do? Really?",
				   name, omxDataColumnName(this, col));
		}
		for (int rx=0; rx < rows; ++rx) {
			if (cd.intData[rx] != NA_INTEGER) continue;
			Rf_error("%s: NA in definition variable '%s' row %d",
				 name, omxDataColumnName(this, col), 1+rx);
		}
	}
}

void omxData::loadFakeData(omxState *state, double fake)
{
	for (int dx=0; dx < int(defVars.size()); ++dx) {
		defVars[dx].loadData(state, fake);
	}
}

bool omxDefinitionVar::loadData(omxState *state, double val)
{
	omxMatrix *mat = state->matrixList[matrix];
	if (val == omxMatrixElement(mat, row, col)) return false;
	omxSetMatrixElement(mat, row, col, val);
	if (OMX_DEBUG) {
		mxLog("Load data %f into %s[%d,%d], state[%d]",
		      val, mat->name(), row, col, state->getId());
	}
	omxMarkClean(mat);
	markDefVarDependencies(state, this);
	return true;
}

struct cstrCmp {
	bool operator() (const char *s1, const char *s2) const
	{ return strcmp(s1,s2) < 0; }
};

typedef std::map< const char *, int, cstrCmp > OrigColMapType;

static int plookup(OrigColMapType &map, const char *str)
{
	auto it = map.find(str);
	if (it == map.end()) Rf_error("Can't find '%s'", str);
	return it->second;
}

void obsSummaryStats::permute(omxData *data, const Eigen::Ref<const DataColumnIndexVector> &dc)
{
	Eigen::VectorXi invDataColumns(dc.size()); // data -> expectation order
	bool needPermute = false;
	for (int cx=0; cx < int(dc.size()); ++cx) {
		//mxLog("%d %s", cx, omxDataColumnName(data, dc[cx]));
		invDataColumns[dc[cx]] = cx;
		if (dc[cx] != cx) needPermute = true;
	}
	if (!needPermute) return;

	covMat->unshareMemoryWithR();
	if (meansMat) meansMat->unshareMemoryWithR();
	acovMat->unshareMemoryWithR();
	if (fullWeight) fullWeight->unshareMemoryWithR();

	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> p1(dc);

	OrigColMapType colMap;
	for (int cx=0; cx < int(acovMat->colnames.size()); ++cx) {
		OrigColMapType::value_type v(acovMat->colnames[cx], cx);
		//mxLog("%s -> %d", acovMat->colnames[cx], cx);
		colMap.insert(v);
	}

	auto &thresh = thresholdCols;

	Eigen::PermutationMatrix<Eigen::Dynamic> p2(acovMat->cols);
	int px = 0;

	if (thresh.size() == 0 && meansMat) {
		for (int cx=0; cx < int(dc.size()); ++cx) {
			p2.indices()[px++] = plookup(colMap, omxDataColumnName(data, dc[cx]));
		}
	} else if (thresh.size()) {
		for (int cx=0; cx < int(dc.size()); ++cx) {
			auto &t1 = thresh[ dc[cx] ];
			auto cn = omxDataColumnName(data, dc[cx]);
			if (t1.numThresholds) {
				for (int tx=1; tx <= t1.numThresholds; ++tx) {
					auto s1 = string_snprintf("%st%d", cn, tx);
					p2.indices()[px++] = plookup(colMap, s1.c_str());
				}
			} else {
				p2.indices()[px++] = plookup(colMap, cn);
			}
		}
	}

	for (int cx=0; cx < int(dc.size()); ++cx) {
		auto &t1 = thresh[ dc[cx] ];
		if (t1.numThresholds) continue;
		auto cn = omxDataColumnName(data, dc[cx]);
		auto s1 = string_snprintf("var_%s", cn);
		p2.indices()[px++] = plookup(colMap, s1.c_str());
	}
	
	for (int cx=0; cx < int(dc.size())-1; ++cx) {
		auto cn = omxDataColumnName(data, dc[cx]);
		for (int rx=cx+1; rx < int(dc.size()); ++rx) {
			auto rn = omxDataColumnName(data, dc[rx]);
			auto s1 = string_snprintf("poly_%s_%s", rn, cn);
			auto s2 = string_snprintf("poly_%s_%s", cn, rn);
			auto it = colMap.find(s1.c_str());
			if (it == colMap.end()) {
				it = colMap.find(s2.c_str());
			}
			if (it == colMap.end()) Rf_error("Can't find '%s' or '%s'",
							 s1.c_str(), s2.c_str());
			p2.indices()[px++] = it->second;
		}
	}

	EigenVectorAdaptor Emean(meansMat);
	EigenMatrixAdaptor Ecov(covMat);
	EigenMatrixAdaptor Eacov(acovMat);
	EigenMatrixAdaptor Efw(fullWeight);

	//mxPrintMat("p2", p2.indices());

	Emean.derived() = (p1.transpose() * Emean).eval();
	Ecov.derived() = (p1.transpose() * Ecov * p1).eval();
	Eacov.derived() = (p2.transpose() * Eacov * p2).eval();
	Efw.derived() = (p2.transpose() * Efw * p2).eval();

	for (auto &th : thresh) th.dColumn = invDataColumns[th.dColumn];
	std::sort(thresh.begin(), thresh.end(),
		  [](const omxThresholdColumn &a, const omxThresholdColumn &b) -> bool
		  { return a.dColumn < b.dColumn; });
}

void obsSummaryStats::log()
{
	mxLog("numObs %d numOrdinal %d", numObs, numOrdinal);
	if (covMat) omxPrint(covMat, "cov");
	if (meansMat) omxPrint(meansMat, "mean");
	if (acovMat) omxPrint(acovMat, "acov");
	if (fullWeight && acovMat != fullWeight) omxPrint(fullWeight, "full");
	for (auto &th : thresholdCols) { th.log(); }
	if (thresholdMat) omxPrint(thresholdMat, "thr");
}

void omxData::permute(const Eigen::Ref<const DataColumnIndexVector> &dc)
{
	if (!dc.size()) return;
	if (permuted) Rf_error("Cannot permute '%s' two different ways", name);
	permuted = true;

	getSingleObsSummaryStats().permute(this, dc);
}

template <typename T>
void getContRow(std::vector<ColumnData> &df,
		int row,
		const Eigen::Ref<const DataColumnIndexVector> &dc,
		Eigen::MatrixBase<T> &out)
{
	for (int cx=0; cx < dc.size(); ++cx) {
		auto &cd = df[ dc[cx] ];
		out[cx] = cd.realData[row];
	}
}

void omxData::wlsAllContinuousCumulants(omxState *state,
					const Eigen::Ref<const DataColumnIndexVector> &dc)
{
	// permit mxData(type='raw') TODO

	int numCols = dc.size();
	int numColsStar = numCols*(numCols+1)/2;
	auto &o1 = *oss;

	o1.covMat = omxInitMatrix(numCols, numCols, state);
	//o1.meansMat = omxInitMatrix(1, numCols, state);
	o1.fullWeight = omxInitMatrix(numColsStar, numColsStar, state);
	EigenMatrixAdaptor Ecov(o1.covMat);
	Eigen::VectorXd Emean(numCols);
	Ecov.setZero();
	Emean.setZero();

	Eigen::ArrayXXd data(int(numObs), numCols);
	Eigen::VectorXd r1(numCols);
	for (int rx=0; rx < numObs; ++rx) {
		getContRow(rawCols, rx, dc, r1);
		Emean += r1;
		Ecov += r1 * r1.transpose();
		data.row(rx) = r1;
	}
	Emean /= numObs;
	Ecov -= numObs * Emean * Emean.transpose();
	Ecov /= numObs - 1;

	data.rowwise() -= Emean.array().transpose();
	Eigen::MatrixXd Vmat = Ecov * (numObs-1.) / numObs;
	EigenMatrixAdaptor Umat(o1.fullWeight);

	Eigen::ArrayXXd M(numColsStar, 2);
	for (int x1=0, mx=0; x1 < numCols; ++x1) {
		for (int x2=x1; x2 < numCols; ++x2) {
			M(mx, 0) = x2;
			M(mx, 1) = x1;
			mx += 1;
		}
	}

	Eigen::Array<double, 4, 1> ind;
	for (int jx=0; jx < numColsStar; ++jx) {
		for (int ix=jx; ix < numColsStar; ++ix) {
			ind.segment(0,2) = M.row(ix);
			ind.segment(2,2) = M.row(jx);
			Umat(ix,jx) = (data.col(ind[0]) * data.col(ind[1]) *
				       data.col(ind[2]) * data.col(ind[3])).sum() / numObs -
				Vmat(ind[0],ind[1]) * Vmat(ind[2],ind[3]);
		}
	}

	int singular = InvertSymmetricPosDef(Umat, 'L');
	if (singular) {
		omxRaiseErrorf("%s: cannot invert full weight matrix (%d)", name, singular);
		return;
	}
	Umat.triangularView<Eigen::Upper>() = Umat.transpose().triangularView<Eigen::Upper>();
	Umat *= numObs;

	if (strEQ(wlsType, "WLS")) {
		o1.acovMat = o1.fullWeight;
	} else {
		o1.acovMat = omxInitMatrix(numColsStar, numColsStar, state);
		EigenMatrixAdaptor acov(o1.acovMat);
		if (strEQ(wlsType, "ULS")) {
			acov.setIdentity();
		} else { // DWLS
			acov.setZero();
			for (int ix=0; ix < numColsStar; ++ix) {
				acov(ix,ix) = 1./((data.col(M(ix, 0)) * data.col(M(ix, 0)) *
						  data.col(M(ix, 1)) * data.col(M(ix, 1))).sum() / numObs -
						  Vmat(M(ix, 0), M(ix, 1)) * Vmat(M(ix, 0), M(ix, 1)));
			}
		}
		acov *= numObs;
	}
}

template <typename T1, typename T2>
void tabulate(Eigen::MatrixBase<T1> &data, Eigen::MatrixBase<T2> &out)
{
	out.setZero();
	for (int rx=0; rx < data.rows(); ++rx) {
		if (data[rx] == NA_INTEGER) continue;
		out[ data[rx] - 1 ] += 1;
	}
}

template <typename T1>
void cumsum(Eigen::MatrixBase<T1> &data)
{
	for (int rx=data.rows()-2; rx >= 0; --rx) {
		data.segment(rx+1, data.size() - rx - 1).array() += data[rx];
	}
}

struct OLSRegression {
	omxData &data;
	int rows;
	ColumnData *response;
	Eigen::MatrixXd pred;
	Eigen::MatrixXd predCov;
	Eigen::ArrayXd resid;
	Eigen::VectorXd beta;
	Eigen::MatrixXd scores;
	double var;
	OLSRegression(omxData *_d, const Eigen::Ref<const Eigen::MatrixXd> _pred);
	void setResponse(ColumnData &response);
	void calcScores();
};

OLSRegression::OLSRegression(omxData *_d, const Eigen::Ref<const Eigen::MatrixXd> _pred)
	: data(*_d), rows(int(data.numObs))
{
	resid.resize(rows);
	pred.resize(rows, 1 + _pred.cols());
	pred.col(0).setConstant(1.0);
	pred.block(0,1,rows,_pred.cols()) = _pred;
	predCov = pred.transpose() * pred;
	int singular = InvertSymmetricPosDef(predCov, 'L');
	if (singular) {
		omxRaiseErrorf("%s: cannot invert exogenous predictors (%d)",
			       data.name, singular);
		return;
	}
}

void OLSRegression::setResponse(ColumnData &cd)
{
	response = &cd;
	Eigen::Map< Eigen::VectorXd > ycol(cd.realData, rows);
	if (pred.cols()) {
		beta = predCov.selfadjointView<Eigen::Lower>() * pred.transpose() * ycol;
		resid = ycol - pred * beta;
	} else {
		beta.resize(1);
		beta[0] = ycol.mean();
		resid = ycol.array() - beta[0];
	}
	var = resid.square().sum() / rows;
}

void OLSRegression::calcScores()
{
	Eigen::Map< Eigen::VectorXd > ycol(response->realData, rows);
	scores.resize(rows, 1 + pred.cols());  // mean pred var
	scores.block(0,0,rows,pred.cols()) = (pred.array().colwise() * resid) / var;
	scores.col(pred.cols()) = -1./(2*var) + 1./(2*var*var) * resid * resid;
}

struct ProbitRegression : NewtonRaphsonObjective {
	omxData &data;
	int rows;
	int numThr;
	ColumnData *response;
	std::vector<int> exoPred;
	const Eigen::Ref<const Eigen::MatrixXd> pred;
	int verbose;
	Eigen::VectorXd param;
	std::vector<std::string> pnames;
	double fit;
	Eigen::ArrayXd pr;
	bool stale;
	Eigen::ArrayXXd zi;
	Eigen::ArrayXXd dzi;
	Eigen::ArrayXXd scores;
	Eigen::VectorXd grad;
	Eigen::ArrayXXd dxa;
	Eigen::ArrayXXd Y1, Y2;
	Eigen::MatrixXd hess;

	ProbitRegression(omxData *_d, std::vector<int> &_exoPred,
			 const Eigen::Ref<const Eigen::MatrixXd> _pred);
	void setResponse(ColumnData &response);
	virtual double getFit() { return fit; }
	virtual const char *paramIndexToName(int px)
	{ return pnames[px].c_str(); }
	void evaluate0();
	void calcScores();
	virtual void evaluateFit();
	virtual void evaluateDerivs(int want);
	virtual double *getParamVec() { return param.data(); }
	virtual double *getGrad() { return grad.data(); }
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir);
};

ProbitRegression::ProbitRegression(omxData *_d, std::vector<int> &_exoPred,
				   const Eigen::Ref<const Eigen::MatrixXd> _pred) :
	data(*_d), rows(int(data.numObs)), numThr(0),
	response(0), exoPred(_exoPred), pred(_pred), verbose(data.verbose), stale(true)
{
	zi.resize(rows, 2);
	dzi.resize(rows, 2);
	pr.resize(rows);
}

void ProbitRegression::setResponse(ColumnData &_r)
{
	response = &_r;
	numThr = response->levels.size()-1;

	Eigen::Map< Eigen::VectorXi > ycol(response->intData, rows);
	Eigen::VectorXi tab(response->levels.size());
	tabulate(ycol, tab);
	Eigen::VectorXd prop = (tab.cast<double>() / double(tab.sum())).
		segment(0, numThr);
	cumsum(prop);
	param.resize(prop.size() + pred.cols());
	pnames.clear();
	for (int px=0; px < prop.size(); ++px) {
		param[px] = Rf_qnorm5(prop[px], 0., 1., 1, 0);
		if (verbose >= 1) pnames.push_back(string_snprintf("th%d", 1+px));
	}
	if (verbose >= 1) {
		for (int cx=0; cx < pred.cols(); ++cx) {
			auto &c1 = data.rawCols[ exoPred[cx] ];
			pnames.push_back(c1.name);
		}
	}
	param.segment(numThr, pred.cols()).array() = 0;

	dxa.resize(rows, numThr);

	Y1.resize(rows, numThr);
	Y2.resize(rows, numThr);
	Y1.setZero();
	Y2.setZero();
	for (int rx=0; rx < rows; ++rx) {
		if (ycol[rx]-2 >= 0)     Y2(rx, ycol[rx]-2) = 1;
		if (ycol[rx]-1 < numThr) Y1(rx, ycol[rx]-1) = 1;
	}

	lbound.resize(param.size());
	lbound.setConstant(NEG_INF);
	ubound.resize(param.size());
	ubound.setConstant(INF);
	scores.resize(rows, param.size());
	hess.resize(param.size(), param.size());
	stale = true;
}

void ProbitRegression::evaluate0()
{
	Eigen::Map< Eigen::VectorXi > ycol(response->intData, rows);
	Eigen::VectorXd th(1 + response->levels.size());
	th.segment(1, numThr) = param.segment(0, numThr);
	th[0] = -std::numeric_limits<double>::infinity();
	th[response->levels.size()] = std::numeric_limits<double>::infinity();

	for (int rx=0; rx < rows; ++rx) {
		double eta = 0;
		if (pred.cols()) eta = pred.row(rx).matrix() * param.segment(numThr, pred.cols());
		zi(rx,0) = std::min(INF, th[ycol[rx]] - eta);
		zi(rx,1) = std::max(NEG_INF, th[ycol[rx]-1] - eta);
		pr[rx] = Rf_pnorm5(zi(rx,0), 0., 1., 1, 0) - Rf_pnorm5(zi(rx,1), 0., 1., 1, 0);
	}
	stale = false;
}

void ProbitRegression::evaluateFit()
{
	evaluate0();
	fit = -pr.array().log().sum();
}

void ProbitRegression::calcScores()
{
	if (stale) evaluate0();
	Eigen::Map< Eigen::VectorXi > ycol(response->intData, rows);
	dxa.setZero();
	for (int rx=0; rx < rows; ++rx) {
		dzi(rx,0) = Rf_dnorm4(zi(rx,0), 0., 1., 0);
		dzi(rx,1) = Rf_dnorm4(zi(rx,1), 0., 1., 0);
		if (ycol[rx]-2 >= 0) {
			dxa(rx,ycol[rx]-2) -= dzi(rx,1);
		}
		if (ycol[rx]-1 < numThr) {
			dxa(rx, ycol[rx]-1) += dzi(rx,0);
		}
	}
	scores.block(0,0,rows,numThr) = dxa.colwise() / pr;
	scores.block(0,numThr,rows,pred.cols()) =
		pred.array().colwise() * ((dzi.col(1) - dzi.col(0)) / pr);
}

void ProbitRegression::evaluateDerivs(int want)
{
	if (want & FF_COMPUTE_FIT) {
		evaluateFit();
	} else {
		evaluate0();
	}

	calcScores();
	grad = -scores.colwise().sum();
	
	Eigen::ArrayXXd gdzi = (dzi * zi).colwise() / -pr;
	Eigen::ArrayXd pr2 = pr * pr;
	pr2 = 1.0/pr2;

	hess.block(0,0,numThr,numThr) =
		dxa.transpose().matrix() * (dxa.colwise() * pr2).matrix() -
		((Y1.colwise() * gdzi.col(0)).transpose().matrix() * Y1.matrix() -
		 (Y2.colwise() * gdzi.col(1)).transpose().matrix() * Y2.matrix());

	Eigen::ArrayXXd dxb = pred.array().colwise() * (dzi.col(0) - dzi.col(1));

	hess.block(numThr,numThr,pred.cols(),pred.cols()) =
		dxb.transpose().matrix() * (dxb.colwise() * pr2).matrix() -
		((pred.array().colwise() * gdzi.col(0)).transpose().matrix() * pred.matrix() -
		 (pred.array().colwise() * gdzi.col(1)).transpose().matrix() * pred.matrix());

	hess.block(0,numThr,numThr,pred.cols()) =
		-(dxa.transpose().matrix() * (dxb.colwise() * pr2).matrix() -
		  ((Y1.colwise() * gdzi.col(0)).transpose().matrix() * pred.matrix() -
		   (Y2.colwise() * gdzi.col(1)).transpose().matrix() * pred.matrix()));
	
	hess.block(numThr,0,pred.cols(),numThr) =
		hess.block(0,numThr,numThr,pred.cols()).transpose();

	if (want & FF_COMPUTE_HESSIAN) ; // TODO
}

void ProbitRegression::setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir)
{
	Eigen::MatrixXd ihess = hess;
	int singular = InvertSymmetricPosDef(ihess, 'U');
	if (singular) {
		singular = InvertSymmetricIndef(ihess, 'U');
		if (singular) ihess = Eigen::DiagonalMatrix<double, Eigen::Dynamic>(param.size());
	}
	searchDir = ihess.selfadjointView<Eigen::Upper>() * grad;
}

template <typename T1, typename T2>
void regressOrdinalThresholds(const Eigen::MatrixBase<T1> &pred,
			     ColumnData &oc, WLSVarData &ov,
			     Eigen::ArrayBase<T2> &zi)
{
	int rows = pred.rows();
	zi.derived().resize(rows, 2);

	int numThr = oc.levels.size() - 1;
	Eigen::VectorXd th(2 + numThr);
	th.segment(1, numThr) = ov.theta.segment(0, numThr);
	th[0] = -std::numeric_limits<double>::infinity();
	th[numThr + 1] = std::numeric_limits<double>::infinity();

	Eigen::Map< Eigen::VectorXi > ycol(oc.intData, rows);
	for (int rx=0; rx < rows; ++rx) {
		double eta = 0;
		if (pred.cols()) pred.row(rx) * ov.theta.matrix().segment(numThr, pred.cols());
		zi(rx,0) = std::min(INF, th[ycol[rx]] - eta);
		zi(rx,1) = std::max(NEG_INF, th[ycol[rx]-1] - eta);
	}
}

struct PolyserialCor : UnconstrainedObjective {
	// continuous
	double var;
	Eigen::ArrayXd &resid;
	Eigen::ArrayXd zee;
	// ordinal
	Eigen::ArrayXXd zi;
	Eigen::ArrayXXd dzi;
	// other stuff
	omxData &data;
	int rows;
	int numThr;
	ColumnData &oc;
	WLSVarData &ov;
	double rho;
	double param;
	double R;
	const Eigen::Ref<const Eigen::MatrixXd> pred;
	Eigen::ArrayXXd tau;
	Eigen::ArrayXXd tauj;
	Eigen::ArrayXd pr;
	Eigen::ArrayXXd scores;

	PolyserialCor(omxData *_d, WLSVarData &cv, ColumnData &_oc, WLSVarData &_ov,
		      const Eigen::Ref<const Eigen::MatrixXd> _pred) :
		resid(cv.resid), data(*_d), rows(int(data.numObs)), oc(_oc), ov(_ov), pred(_pred)
	{
		lbound.resize(1);
		lbound.setConstant(NEG_INF);
		ubound.resize(1);
		ubound.setConstant(INF);

		var = cv.theta[cv.theta.size()-1];
		zee = resid / sqrt(var);

		pr.resize(rows);
		dzi.resize(rows, 2);

		regressOrdinalThresholds(pred, oc, ov, zi);

		Eigen::Map< Eigen::VectorXi > ycol(oc.intData, rows);
		numThr = oc.levels.size() - 1;
		double den = 0;
		for (int tx=0; tx < numThr; ++tx) den += Rf_dnorm4(ov.theta[tx], 0., 1., 0);
		rho = (zee * ycol.cast<double>().array()).sum() / (rows * sqrt(var) * den);

		if (fabs(rho) >= 1.0) rho = 0;
		param = atanh(rho);
	}
	virtual double *getParamVec() { return &param; };
	virtual double getFit(const double *_x)
	{
		rho = tanh(_x[0]);
		R = sqrt(1 - rho * rho);
		tau = (zi.colwise() - rho * zee) / R;

		for (int rx=0; rx < rows; ++rx) {
			pr[rx] = Rf_pnorm5(tau(rx,0), 0., 1., 1, 0) - Rf_pnorm5(tau(rx,1), 0., 1., 1, 0);
		}
		pr.max(std::numeric_limits<double>::epsilon());
		double fit = -pr.log().sum();
		// py1 <- dnorm(Y1, mean=y1.ETA, sd=y1.SD) <-- is constant
		return fit;
	}
	virtual void getGrad(const double *_x, double *grad)
	{
		// can assume getFit just called
		for (int rx=0; rx < rows; ++rx) {
			dzi(rx,0) = Rf_dnorm4(tau(rx,0), 0., 1., 0);
			dzi(rx,1) = Rf_dnorm4(tau(rx,1), 0., 1., 0);
		}

		tauj = dzi * ((zi * rho).colwise() - zee);
		double dx_rho = (1./(R*R*R*pr) * (tauj.col(0) - tauj.col(1))).sum();

		double cosh_x = cosh(_x[0]);
		grad[0] = -dx_rho * 1./(cosh_x * cosh_x);
	}
	void calcScores()
	{
		// mu1 var1 th2 beta1 beta2 rho
		scores.resize(rows, 2 + numThr + pred.cols() * 2 + 1);
		scores.block(0,2,rows,numThr).setZero();

		double R3 = R*R*R;
		Eigen::Map< Eigen::VectorXi > ycol(oc.intData, rows);
		for (int rx=0; rx < rows; ++rx) {
			double irpr = 1.0 / (R * pr[rx]);
			scores(rx,0) = 1.0/sqrt(var) *
				(zee[rx] + rho * irpr * (dzi(rx,0)-dzi(rx,1)));
			scores(rx,1) =
				1.0/(2*var) * ((zee[rx]*zee[rx] - 1.0) +
					       rho*zee[rx] * irpr * (dzi(rx,0)-dzi(rx,1)));
			if (ycol(rx)-1 < numThr)
				scores(rx, 2 + ycol(rx)-1) = dzi(rx,0) * irpr;
			if (ycol(rx)-2 >= 0)
				scores(rx, 2 + ycol(rx)-2) = -dzi(rx,1) * irpr;
			scores.row(rx).segment(2+numThr, pred.cols()) =
				scores(rx,0) * pred.row(rx);
			scores.row(rx).segment(2+numThr+pred.cols(), pred.cols()) =
				-pred.row(rx) * (dzi(rx,0)-dzi(rx,1)) * irpr;
			scores(rx, 2 + numThr + pred.cols() * 2) =
				1./(R3*pr[rx]) * (tauj(rx,0) - tauj(rx,1));
		}
	}
};

struct PolychoricCor : UnconstrainedObjective {
	typedef UnconstrainedObjective super;
	omxData &data;
	int rows;
	ColumnData &c1;
	WLSVarData &v1;
	ColumnData &c2;
	WLSVarData &v2;
	const Eigen::Ref<const Eigen::MatrixXd> pred;
	int numThr1;
	int numThr2;
	Eigen::ArrayXXd z1;
	Eigen::ArrayXXd z2;
	Eigen::ArrayXd pr;
	double rho;
	double param;
	Eigen::ArrayXXd scores;

	PolychoricCor(omxData *_d, ColumnData &_c1, WLSVarData &_v1,
		      ColumnData &_c2, WLSVarData &_v2,
		      const Eigen::Ref<const Eigen::MatrixXd> _pred)
		: data(*_d), rows(int(data.numObs)), c1(_c1), v1(_v1), c2(_c2), v2(_v2), pred(_pred)
	{
		lbound.resize(1);
		lbound.setConstant(NEG_INF);
		ubound.resize(1);
		ubound.setConstant(INF);
		numThr1 = c1.levels.size()-1;
		numThr2 = c2.levels.size()-1;

		// when exoPred is empty, massive speedups are possible TODO

		regressOrdinalThresholds(pred.derived(), c1, v1, z1);
		regressOrdinalThresholds(pred.derived(), c2, v2, z2);

		pr.resize(rows);

		Eigen::Map< Eigen::ArrayXi > y1(c1.intData, rows);
		Eigen::Map< Eigen::ArrayXi > y2(c2.intData, rows);
		Eigen::ArrayXd y1c = y1.cast<double>() - y1.cast<double>().mean();
		Eigen::ArrayXd y2c = y2.cast<double>() - y2.cast<double>().mean();
		rho = (y1c * y2c).sum() / (sqrt((y1c*y1c).sum() * (y2c*y2c).sum()));

		if (fabs(rho) >= 1.0) rho = 0;
		param = atanh(rho);
	}

	virtual double *getParamVec() { return &param; };
	virtual double getFit(const double *_x)
	{
		rho = tanh(_x[0]);

		for (int rx=0; rx < rows; ++rx) {
			pr[rx] = pbivnorm(z1(rx,1), z2(rx,1), z1(rx,0), z2(rx,0), rho);
		}

		return -pr.log().sum();
	}
	virtual void getGrad(const double *_x, double *grad)
	{
		double dx = 0;
		for (int rx=0; rx < rows; ++rx) {
			dx += dbivnorm(z1(rx,1), z2(rx,1), z1(rx,0), z2(rx,0), rho) / pr[rx];
		}
		double cosh_x = cosh(_x[0]);
		grad[0] = -dx / (cosh_x * cosh_x);
	}
	void calcScores()
	{
		// th1 th2 beta1 beta2 rho
		Eigen::Map< Eigen::VectorXi > y1(c1.intData, rows);
		Eigen::Map< Eigen::VectorXi > y2(c2.intData, rows);
		Eigen::ArrayXXd Z1(rows, 2);
		Eigen::ArrayXXd Z2(rows, 2);

		scores.resize(rows, numThr1 + numThr2 + pred.cols() * 2 + 1);
		scores.setZero();

		double R = sqrt(1 - rho*rho);
		for (int rx=0; rx < rows; ++rx) {
			double ipr = 1.0 / pr[rx];
			Z1(rx,0) = Rf_dnorm4(z1(rx,0), 0., 1., 0) *
				(Rf_pnorm5((z2(rx,0) - rho*z1(rx,0))/R, 0., 1., 1, 0) -
				 Rf_pnorm5((z2(rx,1) - rho*z1(rx,0))/R, 0., 1., 1, 0)) * ipr;
			Z1(rx,1) = Rf_dnorm4(z1(rx,1), 0., 1., 0) *
				(Rf_pnorm5((z2(rx,0) - rho*z1(rx,1))/R, 0., 1., 1, 0) -
				 Rf_pnorm5((z2(rx,1) - rho*z1(rx,1))/R, 0., 1., 1, 0)) * ipr;
			Z2(rx,0) = Rf_dnorm4(z2(rx,0), 0., 1., 0) *
				(Rf_pnorm5((z1(rx,0) - rho*z2(rx,0))/R, 0., 1., 1, 0) -
				 Rf_pnorm5((z1(rx,1) - rho*z2(rx,0))/R, 0., 1., 1, 0)) * ipr;
			Z2(rx,1) = Rf_dnorm4(z2(rx,1), 0., 1., 0) *
				(Rf_pnorm5((z1(rx,0) - rho*z2(rx,1))/R, 0., 1., 1, 0) -
				 Rf_pnorm5((z1(rx,1) - rho*z2(rx,1))/R, 0., 1., 1, 0)) * ipr;

			if (y1(rx)-1 < numThr1) scores(rx, y1(rx)-1) = Z1(rx,0);
			if (y1(rx)-2 >= 0)      scores(rx, y1(rx)-2) = -Z1(rx,1);
			if (y2(rx)-1 < numThr2) scores(rx, numThr1 + y2(rx)-1) = Z2(rx,0);
			if (y2(rx)-2 >= 0)      scores(rx, numThr1 + y2(rx)-2) = -Z2(rx,1);
			scores.row(rx).segment(numThr1+numThr2, pred.cols()) =
				(Z1(rx,1)-Z1(rx,0)) * pred.row(rx).array();
			scores.row(rx).segment(numThr1+numThr2+pred.cols(), pred.cols()) =
				(Z2(rx,1)-Z2(rx,0)) * pred.row(rx).array();
			scores(rx,numThr1+numThr2+2*pred.cols()) =
				dbivnorm(z1(rx,1), z2(rx,1), z1(rx,0), z2(rx,0), rho) * ipr;
		}
	}
};

struct PearsonCor {
	double rho;
	Eigen::ArrayXXd scores;

	PearsonCor(WLSVarData &pv1, WLSVarData &pv2,
		   const Eigen::Ref<const Eigen::MatrixXd> pred)
	{
		int rows = pv1.resid.size();
		rho = 2.*pv1.resid.matrix().dot(pv2.resid.matrix()) /
			(pv1.resid.square().sum()+pv2.resid.square().sum());
		double R = (1 - rho*rho);
		double i2r = 1./(2.*R);
		double var_y1 = pv1.theta[pv1.theta.size()-1];
		double sd_y1 = sqrt(var_y1);
		double var_y2 = pv2.theta[pv2.theta.size()-1];
		double sd_y2 = sqrt(var_y2);

		scores.resize(rows, 4 + pred.cols()*2 + 1);

		scores.col(0) = (2*pv1.resid/var_y1 - 2*rho*pv2.resid/(sd_y1*sd_y2)) * i2r;
		scores.col(1) = (2*pv2.resid/var_y2 - 2*rho*pv1.resid/(sd_y1*sd_y2)) * i2r;
		scores.col(2) = -(.5/var_y1 - ((pv1.resid*pv1.resid)/(var_y1*var_y1) -
					       rho*pv1.resid*pv2.resid/(var_y1*sd_y1*sd_y2)) * i2r);
		scores.col(3) = -(.5/var_y2 - ((pv2.resid*pv2.resid)/(var_y2*var_y2) -
					       rho*pv1.resid*pv2.resid/(var_y2*sd_y1*sd_y2)) * i2r);
		scores.block(0,4,rows,pred.cols()) = pred.array().colwise() * scores.col(0);
		scores.block(0,4+pred.cols(),rows,pred.cols()) = pred.array().colwise() * scores.col(1);

		Eigen::ArrayXd zee = (pv1.resid*pv1.resid/var_y1
				      -2*rho*pv1.resid*pv2.resid/(sd_y1*sd_y2) +
				      pv2.resid*pv2.resid/var_y2);
		scores.col(4+2*pred.cols()) = rho/R + pv1.resid*pv2.resid/(sd_y1*sd_y2*R) - zee*rho/(R*R);
	}
};

template <typename T1, typename T2, typename T3>
void copyBlockwise(const Eigen::MatrixBase<T1> &in, Eigen::MatrixBase<T3> &out, T2 includeTest)
{
	for (int cx=0; cx < out.cols(); ++cx) {
		if (!includeTest(cx)) continue;
		for (int rx=0; rx < out.rows(); ++rx) {
			if (!includeTest(rx)) continue;
			out(rx,cx) = in(rx,cx);
		}
	}
}

void omxData::recalcWLSStats(omxState *state, const Eigen::Ref<const DataColumnIndexVector> &dc,
			     std::vector<int> &exoPred)
{
	if (strEQ(getType(), "acov")) {
		permute(dc);
		return;
	}

	// ensure no missing data TODO

	int numCols = dc.size();
	int numColsStar = numCols*(numCols+1)/2;
	if (numObs-1 < numColsStar) {
		Rf_error("%s: too few observations (%d) for the number of columns (%d).\n"
			 "For WLS, you need at least n*(n+1)/2 + 1 = %d observations.\n"
			 "Better start rubbing two pennies together.",
			 name, numObs, numCols, numColsStar+1);
	}

	permuted = true;

	if (oss) delete oss;
	oss = new obsSummaryStats;
	auto &o1 = *oss;
	o1.numObs = int(numObs);

	// Maybe still applicable if exoPred is empty? TODO
	//wlsAllContinuousCumulants(state, dc);

	Eigen::MatrixXd pred(o1.numObs, exoPred.size());
	for (int cx=0; cx < int(exoPred.size()); ++cx) {
		auto &e1 = rawCols[ exoPred[cx] ];
		Eigen::Map< Eigen::VectorXd > vec(e1.realData, o1.numObs);
		pred.col(cx) = vec;
	}

	o1.covMat = omxInitMatrix(numCols, numCols, state);
	o1.meansMat = omxInitMatrix(1, numCols, state);
	EigenMatrixAdaptor Ecov(o1.covMat);
	EigenVectorAdaptor Emean(o1.meansMat);

	std::vector<int> contMap;
	int numContinuous = 0;
	int totalThr = 0;
	int maxNumThr = 0;
	std::vector<int> thStart(numCols);
	for (int yy=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		thStart[yy] = totalThr;
		if (cd.type == COLUMNDATA_NUMERIC) {
			contMap.push_back(numContinuous);
			totalThr += 1;  // mean
			numContinuous += 1;
		} else {
			contMap.push_back(-1);
			int numThr = cd.levels.size() - 1;
			totalThr += numThr;
			maxNumThr = std::max(maxNumThr, numThr);
		}
	}

	o1.perVar.resize(numCols);
	o1.SC_VAR.resize(rows, numContinuous);
	o1.SC_SL.resize(rows, numCols * pred.cols());
	o1.SC_TH.resize(rows, totalThr);
	o1.numOrdinal = 0;
	double eps = sqrt(std::numeric_limits<double>::epsilon());
	OLSRegression olsr(this, pred);
	ProbitRegression pr(this, exoPred, pred);

	// based on lav_samplestats_step1.R, lavaan 0.6-2
	for (int yy=0, contOffset=0, thrOffset=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		WLSVarData &pv = o1.perVar[yy];
		if (cd.type == COLUMNDATA_NUMERIC) {
			omxThresholdColumn tc;
			tc.dColumn = dc[yy];
			tc.column = -1;
			tc.numThresholds = 0;
			o1.thresholdCols.push_back(tc);

			olsr.setResponse(cd);
			olsr.calcScores();
			pv.resid = olsr.resid;
			pv.theta.resize(olsr.beta.size() + 1);
			pv.theta.segment(0, olsr.beta.size()) = olsr.beta;
			pv.theta[olsr.beta.size()] = olsr.var;
			Ecov(yy,yy) = olsr.var;
			Emean[yy] = pv.theta[0];
			o1.SC_TH.block(0,thrOffset,rows,1) = 
				olsr.scores.block(0,0,rows,1);
			for (int px=0; px < pred.cols(); ++px)
				o1.SC_SL.col(yy+numCols*px) = olsr.scores.col(1+px);
			o1.SC_VAR.block(0,contOffset,rows,1) = 
				olsr.scores.block(0,1+pred.cols(),rows,1);
			contOffset += 1;
			thrOffset += 1;
		} else {
			pr.setResponse(cd);
			if (exoPred.size()) {
				NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
				nro(pr);
			} else {
				pr.calcScores();
			}
			pv.theta = pr.param;
			omxThresholdColumn tc;
			tc.dColumn = dc[yy];
			tc.column = o1.numOrdinal++;
			tc.numThresholds = pr.numThr;
			o1.thresholdCols.push_back(tc);
			Ecov(yy,yy) = 1.;
			Emean[yy] = 0.;
			o1.SC_TH.block(0,thrOffset,rows,pr.numThr) = 
				pr.scores.block(0,0,rows,pr.numThr);
			for (int px=0; px < pred.cols(); ++px)
				o1.SC_SL.col(yy+numCols*px) = pr.scores.col(pr.numThr+px);
			thrOffset += pr.numThr;
		}
	}

	o1.thresholdMat = omxInitMatrix(maxNumThr, o1.numOrdinal, state);
	EigenMatrixAdaptor Ethr(o1.thresholdMat);
	Ethr.setConstant(nan("uninit"));
	for (int yy=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		if (cd.type == COLUMNDATA_NUMERIC) continue;
		WLSVarData &pv = o1.perVar[yy];
		auto &tc = o1.thresholdCols[yy];
		Ethr.block(0,tc.column,tc.numThresholds,1) = pv.theta.segment(0,tc.numThresholds);
	}

	int pstar = triangleLoc1(numCols-1);
	o1.SC_COR.resize(rows, pstar);
	int A11_size = o1.SC_TH.cols() + o1.SC_SL.cols() + o1.SC_VAR.cols();
	Eigen::MatrixXd A21(pstar, A11_size);
	A21.setZero();
	Eigen::ArrayXXd H22(pstar, pstar);
	H22.setZero();
	Eigen::ArrayXXd H21(pstar, A11_size);
	H21.setZero();

	// based on lav_samplestats_step2.R, lavaan 0.6-2
	for (int jj=0; jj < numCols-1; ++jj) {
		for (int ii=jj+1; ii < numCols; ++ii) {
			int pstar_idx = ii-(jj+1) + pstar - triangleLoc1(numCols - jj - 1);
			ColumnData &cd1 = rawCols[ dc[jj] ];
			ColumnData &cd2 = rawCols[ dc[ii] ];
			WLSVarData &pv1 = o1.perVar[jj];
			WLSVarData &pv2 = o1.perVar[ii];
			//mxLog("consider %s %s [%d]", cd1.name, cd2.name, pstar_idx);
			double rho;
			if (cd1.type == COLUMNDATA_NUMERIC && cd2.type == COLUMNDATA_NUMERIC) {
				PearsonCor pc(pv2, pv1, pred);
				o1.SC_COR.col(pstar_idx) = pc.scores.col(4+2*pred.cols());
				A21(pstar_idx,thStart[ii]) = (o1.SC_COR.col(pstar_idx) * pc.scores.col(0)).sum();
				A21(pstar_idx,thStart[jj]) = (o1.SC_COR.col(pstar_idx) * pc.scores.col(1)).sum();
				for (int px=0; px < pred.cols(); ++px) {
					A21(pstar_idx, totalThr + ii+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(4+px)).sum();
					A21(pstar_idx, totalThr + jj+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(4+pred.cols()+px)).sum();
				}
				A21(pstar_idx, totalThr + pred.cols()*numCols + contMap[ii]) =
					(o1.SC_COR.col(pstar_idx) * pc.scores.col(2)).sum();
				A21(pstar_idx, totalThr + pred.cols()*numCols + contMap[jj]) =
					(o1.SC_COR.col(pstar_idx) * pc.scores.col(3)).sum();
				double sd1 = sqrt(pv1.theta[pv1.theta.size()-1]);
				double sd2 = sqrt(pv2.theta[pv2.theta.size()-1]);
				H21(pstar_idx, totalThr + pred.cols()*numCols + contMap[ii]) =
					sd1 * pc.rho / (2. * sd2);
				H21(pstar_idx, totalThr + pred.cols()*numCols + contMap[jj]) =
					sd2 * pc.rho / (2. * sd1);
				H22(pstar_idx,pstar_idx) = sd1 * sd2;
				rho = pc.rho * H22(pstar_idx,pstar_idx);
			} else if (cd1.type == COLUMNDATA_NUMERIC) {
				PolyserialCor ps(this, pv1, cd2, pv2, pred);
				UnconstrainedSLSQPOptimizer uo(name, 100, eps, verbose);
				uo(ps);
				ps.calcScores();
				o1.SC_COR.col(pstar_idx) = ps.scores.col(2 + ps.numThr + 2*pred.cols());
				A21(pstar_idx, thStart[jj]) = (o1.SC_COR.col(pstar_idx) * ps.scores.col(0)).sum();
				for (int tx=0; tx < ps.numThr; ++tx)
					A21(pstar_idx, thStart[ii]+tx) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+tx)).sum();
				for (int px=0; px < pred.cols(); ++px) {
					A21(pstar_idx, totalThr + jj+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+ps.numThr+px)).sum();
					A21(pstar_idx, totalThr + ii+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+ps.numThr+pred.cols()+px)).sum();
				}
				A21(pstar_idx, totalThr + pred.cols()*numCols + contMap[jj]) =
					(o1.SC_COR.col(pstar_idx) * ps.scores.col(1)).sum();
				double sd1 = sqrt(pv1.theta[pv1.theta.size()-1]);
				H21(pstar_idx, totalThr + pred.cols()*numCols + contMap[jj]) =
					ps.rho / (2. * sd1);
				H22(pstar_idx,pstar_idx) = sd1;
				rho = ps.rho * sd1;
			} else if (cd2.type == COLUMNDATA_NUMERIC) {
				PolyserialCor ps(this, pv2, cd1, pv1, pred);
				UnconstrainedSLSQPOptimizer uo(name, 100, eps, verbose);
				uo(ps);
				ps.calcScores();
				o1.SC_COR.col(pstar_idx) = ps.scores.col(2 + ps.numThr + 2*pred.cols());
				A21(pstar_idx, thStart[ii]) = (o1.SC_COR.col(pstar_idx) * ps.scores.col(0)).sum();
				for (int tx=0; tx < ps.numThr; ++tx)
					A21(pstar_idx, thStart[jj]+tx) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+tx)).sum();
				for (int px=0; px < pred.cols(); ++px) {
					A21(pstar_idx, totalThr + ii+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+ps.numThr+px)).sum();
					A21(pstar_idx, totalThr + jj+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * ps.scores.col(2+ps.numThr+pred.cols()+px)).sum();
				}
				A21(pstar_idx, totalThr + pred.cols()*numCols + contMap[ii]) =
					(o1.SC_COR.col(pstar_idx) * ps.scores.col(1)).sum();
				double sd1 = sqrt(pv2.theta[pv2.theta.size()-1]);
				H21(pstar_idx, totalThr + pred.cols()*numCols + contMap[ii]) =
					ps.rho / (2. * sd1);
				H22(pstar_idx,pstar_idx) = sd1;
				rho = ps.rho * sd1;
			} else {
				PolychoricCor pc(this, cd2, pv2, cd1, pv1, pred);
				UnconstrainedSLSQPOptimizer uo(name, 100, eps, verbose);
				uo(pc);
				H22(pstar_idx,pstar_idx) = 1.0;
				rho = pc.rho;
				pc.calcScores();
				o1.SC_COR.col(pstar_idx) = pc.scores.col(pc.numThr1 + pc.numThr2 + 2*pred.cols());
				for (int tx=0; tx < pc.numThr1; ++tx)
					A21(pstar_idx, thStart[ii]+tx) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(tx)).sum();
				for (int tx=0; tx < pc.numThr2; ++tx)
					A21(pstar_idx, thStart[jj]+tx) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(pc.numThr1 + tx)).sum();
				int numThr = pc.numThr1 + pc.numThr2;
				for (int px=0; px < pred.cols(); ++px) {
					A21(pstar_idx, totalThr + ii+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(numThr+px)).sum();
					A21(pstar_idx, totalThr + jj+px*numCols) =
						(o1.SC_COR.col(pstar_idx) * pc.scores.col(numThr+pred.cols()+px)).sum();
				}
			}
			Ecov(ii,jj) = rho;
			Ecov(jj,ii) = rho;
		}
	}

	// mxPrintMat("SC_TH", o1.SC_TH.block(0,0,4,o1.SC_TH.cols())); // good
	// mxPrintMat("SC_SL", o1.SC_SL.block(0,0,4,o1.SC_SL.cols())); // good
	// mxPrintMat("SC_VAR", o1.SC_VAR.block(0,0,4,o1.SC_VAR.cols())); // good
	// mxPrintMat("SC_COR", o1.SC_COR.block(0,0,4,o1.SC_COR.cols())); // good
	// mxPrintMat("A21", A21); // good
	// mxPrintMat("H22", H22); // good
	// mxPrintMat("H21", H21); // good

	// based on lav_muthen1984, lavaan 0.6-2
	int acov_size = A11_size + o1.SC_COR.cols();
	Eigen::MatrixXd SC(rows, acov_size);
	SC.block(0,0,rows,o1.SC_TH.cols()) = o1.SC_TH;
	SC.block(0,o1.SC_TH.cols(),rows,o1.SC_SL.cols()) = o1.SC_SL;
	SC.block(0,o1.SC_TH.cols()+o1.SC_SL.cols(),rows,o1.SC_VAR.cols()) = o1.SC_VAR;
	SC.block(0,o1.SC_TH.cols()+o1.SC_SL.cols()+o1.SC_VAR.cols(),rows,o1.SC_COR.cols()) = o1.SC_COR;
	Eigen::MatrixXd INNER = SC.transpose() * SC;
	// mxPrintMat("INNER", INNER); // good

	Eigen::MatrixXd A11(A11_size,A11_size);
	A11.setZero();
	for (int yy=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		std::vector<bool> mask(A11_size, false);
		int numThr = cd.type == COLUMNDATA_NUMERIC? 1 : cd.levels.size() - 1;
		for (int tx=0; tx < numThr; ++tx) mask[thStart[yy] + tx] = true;
		for (int px=0; px < pred.cols(); ++px) mask[totalThr + yy+numCols*px] = true;
		if (cd.type == COLUMNDATA_NUMERIC)
			mask[totalThr + numCols * pred.cols() + contMap[yy]] = true;
		copyBlockwise(INNER, A11, [&mask](int xx){ return mask[xx]; });
	}
	// mxPrintMat("A11", A11); // good

	if (InvertSymmetricPosDef(A11, 'L')) {
		MoorePenroseInverse(A11);
	}

	Eigen::MatrixXd A22(pstar, pstar);
	A22.setZero();
	for (int ii=0; ii < pstar; ++ii) {
		double val = o1.SC_COR.col(ii).square().sum();
		if (val != 0) A22(ii,ii) = 1.0/val;
	}

	Eigen::MatrixXd A21i = -(A22 * A21 * A11.selfadjointView<Eigen::Lower>());
	Eigen::MatrixXd Bi(acov_size,acov_size);
	Bi.setZero();
	Bi.block(0,0,A11.rows(),A11.cols()) = A11.selfadjointView<Eigen::Lower>();
	Bi.block(A11.rows(),0,A21i.rows(),A21i.cols()) = A21i;
	Bi.block(A11.rows(),A11.cols(), A22.rows(), A22.cols()) = A22;
	// mxPrintMat("Bi", Bi); // good

	o1.fullWeight = omxInitMatrix(acov_size, acov_size, state);
	EigenMatrixAdaptor Efw(o1.fullWeight);
	Efw.derived() = Bi * INNER * Bi.transpose();

	//mxPrintMat("Efw", Efw); //good

	if (numContinuous) {
		Eigen::MatrixXd H(acov_size,acov_size);
		H.setZero();
		H.block(0,0,A11_size,A11_size) = Eigen::MatrixXd::Identity(A11_size,A11_size);
		H.block(A11_size,0,H21.rows(),H21.cols()) = H21;
		H.block(A11_size,A11_size,H22.rows(),H22.cols()) = H22;
		//mxPrintMat("H", H); //good

		Efw.derived() = (H * Efw * H.transpose()).eval();
	}

	//mxPrintMat("Efw", Efw);

	if (strEQ(wlsType, "WLS")) {
		o1.acovMat = o1.fullWeight;
	} else {
		o1.acovMat = omxInitMatrix(acov_size, acov_size, state);
		EigenMatrixAdaptor acov(o1.acovMat);
		if (strEQ(wlsType, "ULS")) {
			// can leave acov unallocated? TODO
			acov.setIdentity();
		} else { // DWLS
			acov.setZero();
			for (int ix=0; ix < acov_size; ++ix) {
				acov(ix,ix) = 1.0 / Efw(ix,ix);
			}
		}
	}
	if (InvertSymmetricPosDef(Efw, 'L')) Rf_error("Attempt to invert acov failed");

	// lavaan divides Efw by rows, we don't
	Efw.derived() = Efw.selfadjointView<Eigen::Lower>();
	
	//mxPrintMat("Efw", Efw);

	//o1.log();
	//Rf_error("Not implemented yet");
}
