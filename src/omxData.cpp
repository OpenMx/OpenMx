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
#include "EnableWarnings.h"

omxData::omxData() : primaryKey(NA_INTEGER), weightCol(NA_INTEGER), currentWeightColumn(0),
		     freqCol(NA_INTEGER), currentFreqColumn(0), permuted(false),
		     wlsType(0), wlsContinuousType(0),
		     dataObject(0), dataMat(0), meansMat(0), 
		     numObs(0), _type(0), numFactor(0), numNumeric(0),
		     rows(0), cols(0), expectation(0)
{}

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

	// temporary TODO
	ProtectedSEXP RwlsData(R_do_slot(dataObj, Rf_install(".rawData")));
	if (Rf_isFrame(RwlsData)) {
		importDataFrame(RwlsData, wlsCols, numNumeric, numFactor);
		
		ProtectedSEXP RwlsType(R_do_slot(dataObj, Rf_install(".wlsType")));
		wlsType = CHAR(STRING_ELT(RwlsType,0));
		ProtectedSEXP RwlsContType(R_do_slot(dataObj, Rf_install(".wlsContinuousType")));
		wlsContinuousType = CHAR(STRING_ELT(RwlsContType,0));
	}

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
		obsStatsVec.resize(1);
		auto &o1 = obsStatsVec[0];
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
				if (levels[i] == NA_INTEGER) continue;
				omxThresholdColumn tc;
				tc.dColumn = i;
				tc.column = columns[i];
				tc.numThresholds = levels[i];
				o1.thresholdCols.push_back(tc);
				if(OMX_DEBUG) {
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
	if (obsStatsVec.size() == 1) {
		auto &o1 = obsStatsVec[0];
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

void obsSummaryStats::permute(const Eigen::Ref<const DataColumnIndexVector> &dc)
{
	std::vector< omxThresholdColumn > &origThresh = thresholdCols;
	std::vector< omxThresholdColumn > oThresh = origThresh;

	covMat->unshareMemoryWithR();
	if (meansMat) meansMat->unshareMemoryWithR();
	acovMat->unshareMemoryWithR();
	if (fullWeight) fullWeight->unshareMemoryWithR();

	Eigen::VectorXi invDataColumns(dc.size()); // data -> expectation order
	for (int cx=0; cx < int(dc.size()); ++cx) {
		invDataColumns[dc[cx]] = cx;
	}
	//mxPrintMat("invDataColumns", invDataColumns);
	//Eigen::VectorXi invDataColumns = dc;
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pm(invDataColumns);
	EigenMatrixAdaptor Ecov(covMat);
	Ecov.derived() = (pm * Ecov * pm.transpose()).eval();
	if (meansMat) {
		EigenVectorAdaptor Emean(meansMat);
		Emean.derived() = (pm * Emean).eval();
	}

	Eigen::MatrixXi mm(dc.size(), dc.size());
	for (int cx=0, en=0; cx < dc.size(); ++cx) {
		for (int rx=cx; rx < dc.size(); ++rx) {
			mm(rx,cx) = en;
			en += 1;
		}
	}
	mm = mm.selfadjointView<Eigen::Lower>();
	mm = (pm * mm * pm.transpose()).eval();
	//mxPrintMat("mm", mm);

	Eigen::VectorXi tstart(origThresh.size() + 1);
	tstart[0] = 0;
	int totalThresholds = 0;
	for (int tx=0; tx < int(origThresh.size()); ++tx) {
		totalThresholds += origThresh[tx].numThresholds;
		tstart[tx+1] = totalThresholds;
	}

	int wpermSize = triangleLoc1(dc.size()) + totalThresholds;
	if (meansMat) wpermSize += dc.size();
	Eigen::VectorXi wperm(wpermSize);

	for (int cx=0, en=0; cx < dc.size(); ++cx) {
		for (int rx=cx; rx < dc.size(); ++rx) {
			wperm[en] = mm(rx,cx);
			en += 1;
		}
	}

	if (meansMat) {
		wperm.segment(triangleLoc1(dc.size()), dc.size()) = dc.array() + triangleLoc1(dc.size());
	}

	std::vector<int> newOrder;
	newOrder.reserve(origThresh.size());
	for (int xx=0; xx < int(origThresh.size()); ++xx) newOrder.push_back(xx);

	std::sort(newOrder.begin(), newOrder.end(),
		  [&](const int &a, const int &b) -> bool
		  { return invDataColumns[origThresh[a].dColumn] < invDataColumns[origThresh[b].dColumn]; });

	//for (auto &order : newOrder) mxLog("new order %d lev %d", order, origThresh[order].numThresholds);

	int thStart = triangleLoc1(dc.size());
	if (meansMat) thStart += dc.size();
	for (int t1=0, dest=0; t1 < int(newOrder.size()); ++t1) {
		int oldIndex = newOrder[t1];
		auto &th = oThresh[oldIndex];
		for (int t2=0; t2 < th.numThresholds; ++t2) {
			wperm[thStart + dest] = thStart + tstart[oldIndex] + t2;
			dest += 1;
		}
	}

	for (auto &th : oThresh) th.dColumn = invDataColumns[th.dColumn];
	std::sort(oThresh.begin(), oThresh.end(),
		  [](const omxThresholdColumn &a, const omxThresholdColumn &b) -> bool
		  { return a.dColumn < b.dColumn; });

	origThresh = oThresh;

	//mxPrintMat("wperm", wperm);
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> wpm(wperm);
	EigenMatrixAdaptor Eweights(acovMat);
	Eweights.derived() = (wpm.transpose() * Eweights * wpm).eval();

	if (fullWeight) {
		EigenMatrixAdaptor Efw(fullWeight);
		Efw.derived() = (wpm.transpose() * Efw * wpm).eval();
	}
	//mxPrintMat("ew", Eweights);
}

void obsSummaryStats::log()
{
	mxLog("numObs %d numOrdinal %d", numObs, numOrdinal);
	if (covMat) omxPrint(covMat, "cov");
	if (meansMat) omxPrint(meansMat, "mean");
	if (acovMat) omxPrint(acovMat, "acov");
	if (fullWeight && acovMat != fullWeight) omxPrint(fullWeight, "full");
	if (thresholdMat) omxPrint(thresholdMat, "thr");
}

void omxData::permute(const Eigen::Ref<const DataColumnIndexVector> &dc)
{
	if (!dc.size()) return;
	if (permuted) Rf_error("Cannot permute '%s' two different ways", name);
	permuted = true;

	if (obsStatsVec.size() != 1) Rf_error("obsStatsVec.size() != 1");

	obsStatsVec[0].permute(dc);
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
	auto &o1 = obsStatsVec[0];

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
		getContRow(wlsCols, rx, dc, r1);
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
	std::vector<int> exoPred;
	ColumnData *response;
	Eigen::MatrixXd pred;
	Eigen::MatrixXd predCov;
	Eigen::ArrayXd resid;
	Eigen::VectorXd beta;
	Eigen::MatrixXd scores;
	double var;
	OLSRegression(omxData *_d, std::vector<int> &_exoPred);
	void setResponse(ColumnData &response);
	void calcScores();
};

OLSRegression::OLSRegression(omxData *_d, std::vector<int> &_exoPred)
	: data(*_d), rows(int(data.numObs)), exoPred(_exoPred)
{
	resid.resize(rows);
	pred.resize(rows, 1 + exoPred.size());
	pred.col(0).setConstant(1.0);
	for (int cx=0; cx < int(exoPred.size()); ++cx) {
		auto &c1 = data.rawCols[ exoPred[cx] ];
		Eigen::Map< Eigen::VectorXd > vec(c1.realData, rows);
		pred.col(1+cx) = vec;
	}
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
	if (exoPred.size()) {
		beta = predCov.selfadjointView<Eigen::Lower>() * pred.transpose() * ycol;
		resid = ycol - pred * beta;
	} else {
		beta.resize(1);
		beta[0] = ycol.mean();
		resid = ycol.array() - beta.col(0).array();
	}
	var = resid.square().sum() / rows;
}

void OLSRegression::calcScores()
{
	Eigen::Map< Eigen::VectorXd > ycol(response->realData, rows);
	scores.resize(rows, 1 + pred.cols());
	scores.block(0,0,rows,pred.cols()) = (pred.array().colwise() * resid) / var;
	scores.col(pred.cols()) = -1./(2*var) + 1./(2*var*var) * resid * resid;
}

struct ProbitRegression : NewtonRaphsonObjective {
	omxData &data;
	int rows;
	int numThr;
	ColumnData *response;
	std::vector<int> exoPred;
	int verbose;
	Eigen::VectorXd param;
	std::vector<std::string> pnames;
	Eigen::ArrayXXd pred;
	double fit;
	Eigen::ArrayXd pr;
	Eigen::ArrayXXd zi;
	Eigen::ArrayXXd dzi;
	Eigen::ArrayXXd scores;
	Eigen::VectorXd grad;
	Eigen::ArrayXXd dxa;
	Eigen::ArrayXXd Y1, Y2;
	Eigen::MatrixXd hess;

	ProbitRegression(omxData *_d, std::vector<int> &_exoPred);
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

ProbitRegression::ProbitRegression(omxData *_d, std::vector<int> &_exoPred) :
	data(*_d), rows(int(data.numObs)), numThr(0),
	response(0), exoPred(_exoPred), verbose(data.verbose)
{
	pred.resize(rows, exoPred.size());
	for (int cx=0; cx < int(exoPred.size()); ++cx) {
		auto &c1 = data.rawCols[ exoPred[cx] ];
		Eigen::Map< Eigen::VectorXd > vec(c1.realData, rows);
		pred.col(cx) = vec;
	}

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
	param.resize(prop.size() + exoPred.size());
	pnames.clear();
	for (int px=0; px < prop.size(); ++px) {
		param[px] = Rf_qnorm5(prop[px], 0., 1., 1, 0);
		if (verbose >= 1) pnames.push_back(string_snprintf("th%d", 1+px));
	}
	if (verbose >= 1) {
		for (int cx=0; cx < int(exoPred.size()); ++cx) {
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
}

void ProbitRegression::evaluate0()
{
	Eigen::Map< Eigen::VectorXi > ycol(response->intData, rows);
	Eigen::VectorXd th(1 + response->levels.size());
	th.segment(1, numThr) = param.segment(0, numThr);
	th[0] = -std::numeric_limits<double>::infinity();
	th[response->levels.size()] = std::numeric_limits<double>::infinity();

	for (int rx=0; rx < rows; ++rx) {
		double eta = pred.row(rx).matrix() * param.segment(numThr, pred.cols());
		zi(rx,0) = std::min(100., th[ycol[rx]] - eta);
		zi(rx,1) = std::max(-100., th[ycol[rx]-1] - eta);
		pr[rx] = Rf_pnorm5(zi(rx,0), 0., 1., 1, 0) - Rf_pnorm5(zi(rx,1), 0., 1., 1, 0);
	}
}

void ProbitRegression::evaluateFit()
{
	evaluate0();
	fit = -pr.array().log().sum();
}

void ProbitRegression::calcScores()
{
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
		pred.colwise() * ((dzi.col(1) - dzi.col(0)) / pr);
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

	Eigen::ArrayXXd dxb = pred.colwise() * (dzi.col(0) - dzi.col(1));

	hess.block(numThr,numThr,pred.cols(),pred.cols()) =
		dxb.transpose().matrix() * (dxb.colwise() * pr2).matrix() -
		((pred.colwise() * gdzi.col(0)).transpose().matrix() * pred.matrix() -
		 (pred.colwise() * gdzi.col(1)).transpose().matrix() * pred.matrix());

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
	obsStatsVec.clear();
	obsStatsVec.resize(1);
	auto &o1 = obsStatsVec[0];

	// Maybe still applicable if exoPred is empty? TODO
	//wlsAllContinuousCumulants(state, dc);

	o1.covMat = omxInitMatrix(numCols, numCols, state);
	o1.meansMat = omxInitMatrix(1, numCols, state);
	EigenMatrixAdaptor Ecov(o1.covMat);
	EigenVectorAdaptor Emean(o1.meansMat);

	o1.perVar.resize(numCols);
	OLSRegression olsr(this, exoPred);
	ProbitRegression pr(this, exoPred);

	// based on lav_samplestats_step1.R, lavaan 0.6-2
	int maxNumThr = 0;
	for (int yy=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		WLSVarData &pv = o1.perVar[yy];
		if (cd.type == COLUMNDATA_NUMERIC) {
			olsr.setResponse(cd);
			olsr.calcScores();
			pv.resid = olsr.resid;
			pv.theta.resize(olsr.beta.size() + 1);
			pv.theta.segment(0, olsr.beta.size()) = olsr.beta;
			pv.theta[olsr.beta.size()] = olsr.var;
			pv.scores = olsr.scores;
			Ecov(yy,yy) = olsr.var;
			Emean[yy] = pv.theta[0];
		} else {
			pr.setResponse(cd);
			if (exoPred.size()) {
				double eps = sqrt(std::numeric_limits<double>::epsilon());
				NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
				nro(pr);
			} else {
				pr.calcScores();
			}
			pv.theta = pr.param;
			pv.scores = pr.scores;
			omxThresholdColumn tc;
			tc.dColumn = dc[yy];
			tc.column = o1.thresholdCols.size();
			tc.numThresholds = pr.numThr;
			o1.thresholdCols.push_back(tc);
			maxNumThr = std::max(maxNumThr, pr.numThr);
			Ecov(yy,yy) = 1.;
			Emean[yy] = 0.;
		}
	}

	o1.numObs = int(numObs);
	o1.numOrdinal = o1.thresholdCols.size();
	o1.thresholdMat = omxInitMatrix(maxNumThr, o1.numOrdinal, state);
	EigenMatrixAdaptor Ethr(o1.thresholdMat);
	Ethr.setConstant(nan("uninit"));
	for (int yy=0, tx=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCols[ dc[yy] ];
		if (cd.type == COLUMNDATA_NUMERIC) continue;
		WLSVarData &pv = o1.perVar[yy];
		auto &tc = o1.thresholdCols[tx];
		Ethr.block(0,tx,tc.numThresholds,1) = pv.theta.segment(0,tc.numThresholds);
		tx += 1;
	}

	// based on lav_samplestats_step2.R, lavaan 0.6-2
	for (int jj=0; jj < numCols-1; ++jj) {
		ColumnData &cd1 = rawCols[ dc[jj] ];
		for (int ii=jj+1; ii < numCols; ++ii) {
			ColumnData &cd2 = rawCols[ dc[ii] ];
			if (cd1.type == COLUMNDATA_NUMERIC && cd2.type == COLUMNDATA_NUMERIC) {
				WLSVarData &pv1 = o1.perVar[jj];
				WLSVarData &pv2 = o1.perVar[ii];
				double tmp = pv1.resid.matrix().dot(pv2.resid.matrix());
				double cor = tmp / (o1.numObs *
						    sqrt(pv1.theta[pv1.theta.size()-1]) *
						    sqrt(pv2.theta[pv2.theta.size()-1]));
				Ecov(ii,jj) = cor;
				Ecov(jj,ii) = cor;
			} else if (cd1.type == COLUMNDATA_NUMERIC) {
				
			} else if (cd2.type == COLUMNDATA_NUMERIC) {
			} else {
			}
		}
	}

	o1.log();

	Rf_error("Not implemented yet");
}
