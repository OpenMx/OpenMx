/*
 *  Copyright 2007-2017 The OpenMx Project
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
#include "EnableWarnings.h"

omxData::omxData() : rownames(0), primaryKey(-1), weightCol(-1), currentWeightColumn(0),
		     dataObject(0), dataMat(0), meansMat(0), acovMat(0), obsThresholdsMat(0),
		     thresholdCols(0), numObs(0), _type(0), numFactor(0), numNumeric(0),
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

	{ScopedProtect p2(dataLoc, R_do_slot(dataObject, Rf_install("verbose")));
	od->verbose = Rf_asInteger(dataLoc);
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
		numObs = other->weightSum;
		addDynamicDataSource(ex);
		// nothing special to do
	} else {
		int num = Rf_length(dataLoc);
		expectation.reserve(num);
		int *evec = INTEGER(dataLoc);

		omxExpectation *refE = NULL;
		BA81Expect *refBA81;
		double weightSum = 0;

		for (int sx=0; sx < num; ++sx) {
			omxExpectation *ex = omxExpectationFromIndex(evec[sx], currentState);
			if (!strEQ(ex->expType, "MxExpectationBA81")) {
				omxRaiseErrorf("MxDataDynamic: type='cov' is only valid for MxExpectationBA81, not '%s'",
					       ex->expType);
				continue;
			}
			BA81Expect *other = (BA81Expect *) ex;
			weightSum += other->weightSum;
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
		numObs = weightSum;
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

void omxData::newDataStatic(omxState *state, SEXP dataObj)
{
	omxData *od = this;
	SEXP dataLoc, dataVal;
	int numCols;

	// PARSE MxData Structure
	if(OMX_DEBUG) {mxLog("Processing Data '%s'", od->name);}

	{
		ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("type")));
		od->_type = CHAR(STRING_ELT(dataLoc,0));
		if(OMX_DEBUG) {mxLog("Element is type %s.", od->_type);}

		ProtectedSEXP needsort(R_do_slot(dataObj, Rf_install(".needSort")));
		od->needSort = Rf_asLogical(needsort);

		ScopedProtect p2(dataLoc, R_do_slot(dataObj, Rf_install("primaryKey")));
		int pk = Rf_asInteger(dataLoc);
		if (pk != NA_INTEGER) {
			primaryKey = pk - 1;
		}

		ProtectedSEXP Rweight(R_do_slot(dataObj, Rf_install("weight")));
		weightCol = Rf_asInteger(Rweight) - 1;
	}
	{ScopedProtect pdl(dataLoc, R_do_slot(dataObj, Rf_install("observed")));
	if(OMX_DEBUG) {mxLog("Processing Data Elements.");}
	if (Rf_isFrame(dataLoc)) {
		if(OMX_DEBUG) {mxLog("Data is a frame.");}
		// Process Data Frame into Columns
		od->rownames = Rf_getAttrib(dataLoc, R_RowNamesSymbol);
		od->cols = Rf_length(dataLoc);
		if(OMX_DEBUG) {mxLog("Data has %d columns.", od->cols);}
		numCols = od->cols;
		SEXP colnames;
		ScopedProtect p2(colnames, Rf_getAttrib(dataLoc, R_NamesSymbol));
		od->rawCols.reserve(numCols);
		for(int j = 0; j < numCols; j++) {
			const char *colname = CHAR(STRING_ELT(colnames, j));
			ColumnData cd = { colname, COLUMNDATA_INVALID, NULL, NULL, NULL };
			SEXP rcol;
			ScopedProtect p1(rcol, VECTOR_ELT(dataLoc, j));
			if(Rf_isFactor(rcol)) {
				cd.type = Rf_isUnordered(rcol)? COLUMNDATA_UNORDERED_FACTOR : COLUMNDATA_ORDERED_FACTOR;
				if(OMX_DEBUG) {mxLog("Column[%d] %s is a factor.", j, colname);}
				cd.intData = INTEGER(rcol);
				cd.levels = Rf_getAttrib(rcol, R_LevelsSymbol);
				od->numFactor++;
			} else if (Rf_isInteger(rcol)) {
				if (j == primaryKey) {
					if(OMX_DEBUG) {mxLog("Column[%d] %s is the primary key", j, colname);}
				} else {
					if(OMX_DEBUG) {mxLog("Column[%d] %s is a foreign key", j, colname);}
				}
				cd.intData = INTEGER(rcol);
				cd.type = COLUMNDATA_INTEGER;
			} else if (Rf_isNumeric(rcol)) {
				if(OMX_DEBUG) {mxLog("Column[%d] %s is numeric.", j, colname);}
				cd.realData = REAL(rcol);
				cd.type = COLUMNDATA_NUMERIC;
				od->numNumeric++;
			} else {
				if(OMX_DEBUG) {mxLog("Column[%d] %s is type %s (ignored)",
						     j, colname, Rf_type2char(TYPEOF(rcol)));}
				cd.type = COLUMNDATA_INVALID;
			}
			od->rawCols.push_back(cd);
		}
		od->rows = Rf_length(VECTOR_ELT(dataLoc, 0));
		if(OMX_DEBUG) {mxLog("And %d rows.", od->rows);}
	} else {
		if(OMX_DEBUG) {mxLog("Data contains a matrix.");}
		od->dataMat = omxNewMatrixFromRPrimitive(dataLoc, state, 0, 0);
		
		if (od->dataMat->colMajor && strEQ(od->_type, "raw")) {
			omxToggleRowColumnMajor(od->dataMat);
		}
		od->cols = od->dataMat->cols;
		od->rows = od->dataMat->rows;
		od->numNumeric = od->cols;
	}
	}

	if (od->hasPrimaryKey() && !(od->rawCols.size() && od->rawCols[primaryKey].intData)) {
		Rf_error("%s: primary key must be an integer or factor column in raw observed data", od->name);
	}

	if (od->hasWeight() && od->rawCols.size() && od->rawCols[weightCol].type != COLUMNDATA_NUMERIC) {
		Rf_error("%s: weight must be a numeric column in raw observed data", od->name);
	}

	if(OMX_DEBUG) {mxLog("Processing Means Matrix.");}
	{ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("means")));
	od->meansMat = omxNewMatrixFromRPrimitive(dataLoc, state, 0, 0);
	}

	if(od->meansMat->rows == 1 && od->meansMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->meansMat, 0, 0)) ||
	    !std::isfinite(omxMatrixElement(od->meansMat, 0, 0)))) {
		omxFreeMatrix(od->meansMat); // Clear just-allocated memory.
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
	{ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("acov")));
	od->acovMat = omxNewMatrixFromRPrimitive(dataLoc, state, 0, 0);
	}

	if(od->acovMat->rows == 1 && od->acovMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->acovMat, 0, 0)) ||
	    !std::isfinite(omxMatrixElement(od->acovMat, 0, 0)))) {
		omxFreeMatrix(od->acovMat); // Clear just-allocated memory.
		od->acovMat = NULL;
	}

	if(OMX_DEBUG) {mxLog("Processing Observed Thresholds Matrix.");}
	{ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("thresholds")));
	od->obsThresholdsMat = omxNewMatrixFromRPrimitive(dataLoc, state, 0, 0);
	if (!strEQ(od->obsThresholdsMat->getType(), "matrix")) {
		Rf_error("Observed thresholds must be constant");
	}
	}

	if(od->obsThresholdsMat->rows == 1 && od->obsThresholdsMat->cols == 1 && 
	   (!R_finite(omxMatrixElement(od->obsThresholdsMat, 0, 0)) ||
	    !std::isfinite(omxMatrixElement(od->obsThresholdsMat, 0, 0)))) {
		omxFreeMatrix(od->obsThresholdsMat); // Clear just-allocated memory.
		od->obsThresholdsMat = NULL;
	} else {
		od->thresholdCols.reserve(od->obsThresholdsMat->cols);
		int *columns;
		{
			ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("thresholdColumns")));
			columns = INTEGER(dataLoc);
		}

		int *levels;
		{
			ScopedProtect p1(dataVal, R_do_slot(dataObj, Rf_install("thresholdLevels")));
			levels = INTEGER(dataVal);
		}

		//for(int i = 0; i < od->obsThresholdsMat->cols; i++) {
		for(int i = 0; i < od->cols; i++) {
			if (levels[i] == NA_INTEGER) continue;
			omxThresholdColumn tc;
			tc.dColumn = i;
			tc.column = columns[i];
			tc.numThresholds = levels[i];
			od->thresholdCols.push_back(tc);
			od->numFactor++; //N.B. must increment numFactor when data@type=='raw' (above) AND when data@type=='acov' (here)
			if(OMX_DEBUG) {
				mxLog("%s: column %d is ordinal with %d thresholds in threshold column %d.", 
				      name, i, levels[i], columns[i]);
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

	currentWeightColumn = getOriginalWeightColumn();
}

double *omxData::getOriginalWeightColumn()
{
	if (!hasWeight()) return 0;
	if (rawCols.size()) {
		return rawCols[weightCol].realData;
	} else {
		return omxMatrixColumn(dataMat, weightCol);
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

	SEXP DataClass;
	const char* dclass;
	{ScopedProtect p1(DataClass, STRING_ELT(Rf_getAttrib(dataObj, R_ClassSymbol), 0));
		dclass = CHAR(DataClass);
	}
	if(OMX_DEBUG) {mxLog("Initializing %s element", dclass);}
	omxData* od = new omxData();
	od->name = name;
	dataList.push_back(od);
	if (strcmp(dclass, "MxDataStatic")==0) od->newDataStatic(this, dataObj);
	else if (strcmp(dclass, "MxDataDynamic")==0) newDataDynamic(dataObj, od);
	else Rf_error("Unknown data class %s", dclass);
	return od;
}

void omxFreeData(omxData* od) {
	omxFreeMatrix(od->dataMat);
	omxFreeMatrix(od->meansMat);
	omxFreeMatrix(od->acovMat);
	omxFreeMatrix(od->obsThresholdsMat);
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
	if (!cd.levels) Rf_error("omxDataGetNumFactorLevels attempt on non-factor");
	return Rf_length(cd.levels);
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

omxMatrix* omxDataAcov(omxData *od)
{
	return od->acovMat;  // can be NULL
}

bool omxDataColumnIsFactor(omxData *od, int col)
{
	if(od->dataMat != NULL) return FALSE;
	ColumnData &cd = od->rawCols[col];
	return cd.intData && cd.levels;
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
	if(od->dataMat) return od->dataMat->colnames[col];
	ColumnData &cd = od->rawCols[col];
	return cd.name;
}

void omxDataKeysCompatible(omxData *upper, omxData *lower, int foreignKey)
{
	ColumnData &ucd = upper->rawCols[upper->primaryKey];
	ColumnData &lcd = lower->rawCols[foreignKey];
	if (ucd.type != lcd.type) {
		Rf_error("Primary key '%s' in %s is not the same type"
			 " as foreign key '%s' in %s",
			 ucd.name, upper->name, lcd.name, lower->name);
	}
	if (ucd.type == COLUMNDATA_ORDERED_FACTOR || ucd.type == COLUMNDATA_UNORDERED_FACTOR) {
		if (Rf_length(ucd.levels) != Rf_length(lcd.levels)) {
			Rf_error("Primary key '%s' in %s has a different number of factor"
				 " levels compared to foreign key '%s' in %s",
				 ucd.name, upper->name, lcd.name, lower->name);
		}
		for (int lx=0; lx < Rf_length(ucd.levels); ++lx) {
			const char *ul = CHAR(STRING_ELT(ucd.levels, lx));
			const char *ll = CHAR(STRING_ELT(lcd.levels, lx));
			if (strEQ(ul, ll)) continue;
			Rf_error("Primary key '%s' in %s has different factor levels ('%s' != '%s')"
				 " compared to foreign key '%s' in %s",
				 ucd.name, upper->name, ul, ll, lcd.name, lower->name);
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

std::vector<omxThresholdColumn> &omxDataThresholds(omxData *od)
{
    return od->thresholdCols;
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
		double newDefVar = omxDoubleDataElement(this, row, defVars[k].column);
		if(ISNA(newDefVar)) {
			Rf_error("Error: NA value for a definition variable is Not Yet Implemented.");
		}
		changed |= defVars[k].loadData(state, newDefVar);
	}
	if (changed && OMX_DEBUG_ROWS(row)) { mxLog("Processing Definition Vars for row %d", row); }
	return changed;
}

void omxData::loadFakeData(omxState *state, double fake)
{
	for (int dx=0; dx < int(defVars.size()); ++dx) {
		defVars[dx].loadData(state, fake);
	}
}

bool omxData::CompareDefVarInMatrix(int lrow, int rrow, omxMatrix *mat, bool &mismatch)
{
	int mnum = ~mat->matrixNumber;
	mismatch = true;
	for (int dx=0; dx < int(defVars.size()); ++dx) {
		omxDefinitionVar &dv = defVars[dx];
		if (dv.matrix != mnum) continue;
		double lval = omxDoubleDataElement(this, lrow, dv.column);
		double rval = omxDoubleDataElement(this, rrow, dv.column);
		if (lval != rval) return lval < rval;
	}
	mismatch = false;
	return false;
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

