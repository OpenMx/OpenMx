/*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
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
#include <RcppEigenCholmod.h>
#include <RcppEigenWrap.h>
#include "CovEntrywisePar.h"
#include "Compute.h"
#include "EnableWarnings.h"

static double clamp(double v, double lo, double hi) { // replace with std::clamp C++-17
  return v < lo ? lo : hi < v ? hi : v;
}

omxData::omxData() : primaryKey(NA_INTEGER), weightCol(NA_INTEGER), currentWeightColumn(0),
                     freqCol(NA_INTEGER), currentFreqColumn(0), numEstimatedEntries(0),
                     parallel(true),
		     noExoOptimize(true), modified(false), minVariance(0), warnNPDuseWeight(true),
		     dataObject(0), dataMat(0), meansMat(0),
										 numObs(0), u_type(0), naAction(NA_PASS), numFactor(0), numNumeric(0),
		     cols(0), expectation(0)
{}

omxData::~omxData()
{
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
	od->u_type = CHAR(STRING_ELT(dataLoc,0));
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
		mxThrow("omxData::connectDynamicData called more than once");
	}

	SEXP dataLoc;
	Rf_protect(dataLoc = R_do_slot(dataObject, Rf_install("expectation")));
	if (Rf_length(dataLoc) == 0) {
		omxRaiseErrorf("mxDataDynamic is not connected to a data source");
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
			if (!strEQ(ex->name, "MxExpectationBA81")) {
				omxRaiseErrorf("MxDataDynamic: type='cov' is only valid for MxExpectationBA81, not '%s'",
					       ex->name);
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
	const int debug = 0;
	int numCols = Rf_length(dataLoc);
	rawCols.clear();
	rawCols.reserve(numCols);
	numNumeric = 0;
	numFactor = 0;
	ProtectedSEXP colnames(Rf_getAttrib(dataLoc, R_NamesSymbol));

	for(int j = 0; j < numCols; j++) {
		const char *colname = CHAR(STRING_ELT(colnames, j));
		ColumnData cd(colname);
		ProtectedSEXP rcol(VECTOR_ELT(dataLoc, j));
		if(Rf_isFactor(rcol)) {
			cd.type = Rf_isUnordered(rcol)? COLUMNDATA_UNORDERED_FACTOR : COLUMNDATA_ORDERED_FACTOR;
			if(debug+OMX_DEBUG) {mxLog("Column[%d] %s is a factor.", j, colname);}
			cd.setBorrow(INTEGER(rcol));
			ProtectedSEXP Rlevels(Rf_getAttrib(rcol, R_LevelsSymbol));
			for (int lx=0; lx < Rf_length(Rlevels); ++lx) {
				cd.levelNames.push_back(R_CHAR(STRING_ELT(Rlevels, lx)));
			}
      cd.setMaxValueFromLevels();
			numFactor++;
		} else if (Rf_isInteger(rcol) || Rf_isLogical(rcol)) {
      // Otherwise logical will be handled by isNumeric
			if(debug+OMX_DEBUG) {mxLog("Column[%d] %s is integer.", j, colname);}
			cd.setBorrow(INTEGER(rcol));
			cd.type = COLUMNDATA_INTEGER;
      cd.setMinValue(0);
		} else if (Rf_isNumeric(rcol)) {
			if(debug+OMX_DEBUG) {mxLog("Column[%d] %s is numeric.", j, colname);}
			cd.setBorrow(REAL(rcol));
			cd.type = COLUMNDATA_NUMERIC;
			numNumeric++;
		} else {
			if(debug+OMX_DEBUG) {mxLog("Column[%d] %s is type %s (ignored)",
					     j, colname, Rf_type2char(TYPEOF(rcol)));}
			cd.type = COLUMNDATA_INVALID;
		}
		rawCols.emplace_back(cd);
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
		od->u_type = CHAR(STRING_ELT(dataLoc,0));
		if(OMX_DEBUG) {mxLog("Element is type %s.", od->u_type);}
    if (strEQ(od->u_type, "none")) od->u_type = "raw";

		if (R_has_slot(dataObj, Rf_install("naAction"))) {
			ProtectedSEXP RactStr(R_do_slot(dataObj, Rf_install("naAction")));
			if (Rf_length(RactStr)) {
				const char *naActStr = as<const char *>(RactStr);
				if (strEQ(naActStr, "pass")) ; // is default
				else if (strEQ(naActStr, "fail")) naAction = NA_FAIL;
				else if (strEQ(naActStr, "omit")) naAction = NA_OMIT;
				else if (strEQ(naActStr, "exclude")) naAction = NA_EXCLUDE;
				else mxThrow("%s: unknown naAction '%s'", name, naActStr);
			}
		}
    fitTolerance = sqrt(std::numeric_limits<double>::epsilon());
		if (R_has_slot(dataObj, Rf_install("fitTolerance"))) {
			ProtectedSEXP Rft(R_do_slot(dataObj, Rf_install("fitTolerance")));
      fitTolerance = Rf_asReal(Rft);
    }
    gradientTolerance = 0.01;
		if (R_has_slot(dataObj, Rf_install("gradientTolerance"))) {
			ProtectedSEXP Rgt(R_do_slot(dataObj, Rf_install("gradientTolerance")));
      gradientTolerance = Rf_asReal(Rgt);
    }

		ProtectedSEXP needsort(R_do_slot(dataObj, Rf_install(".needSort")));
		od->needSort = Rf_asLogical(needsort);

		ScopedProtect p2(dataLoc, R_do_slot(dataObj, Rf_install("primaryKey")));
		primaryKey = carefulMinusOne(Rf_asInteger(dataLoc));

		ProtectedSEXP Rweight(R_do_slot(dataObj, Rf_install("weight")));
		weightCol = carefulMinusOne(Rf_asInteger(Rweight));

		ProtectedSEXP Rfrequency(R_do_slot(dataObj, Rf_install("frequency")));
		freqCol = carefulMinusOne(Rf_asInteger(Rfrequency));
	}
	if (R_has_slot(dataObj, Rf_install(".parallel"))) {
		ProtectedSEXP Rpar(R_do_slot(dataObj, Rf_install(".parallel")));
		parallel = Rf_asLogical(Rpar);
	}
	if (R_has_slot(dataObj, Rf_install(".noExoOptimize"))) {
		ProtectedSEXP Rneo(R_do_slot(dataObj, Rf_install(".noExoOptimize")));
		noExoOptimize = Rf_asLogical(Rneo);
	}
	if (R_has_slot(dataObj, Rf_install("minVariance"))) {
		ProtectedSEXP Rmv(R_do_slot(dataObj, Rf_install("minVariance")));
		minVariance = Rf_asReal(Rmv);
	}
	if (R_has_slot(dataObj, Rf_install("warnNPDuseWeight"))) {
		ProtectedSEXP Rnpd(R_do_slot(dataObj, Rf_install("warnNPDuseWeight")));
		warnNPDuseWeight = Rf_asLogical(Rnpd);
	}

	{ScopedProtect pdl(dataLoc, R_do_slot(dataObj, Rf_install("observed")));

    if (dataLoc == R_NilValue) {
      // OK, probably observedStats only
    } else if (Rf_isFrame(dataLoc)) {
      auto &rd = unfiltered;
      rd.rows = Rf_length(VECTOR_ELT(dataLoc, 0));
      importDataFrame(dataLoc, rd.rawCols, od->numNumeric, od->numFactor);
      cols = int(rd.rawCols.size());
      for (int cx=0; cx < cols; ++cx) {
        rawColMap.emplace(rd.rawCols[cx].name, cx);
      }

      if (hasPrimaryKey()) {
        int kt = rd.rawCols[primaryKey].type;
        if (kt != COLUMNDATA_INTEGER &&
            kt != COLUMNDATA_ORDERED_FACTOR && kt != COLUMNDATA_UNORDERED_FACTOR) {
          mxThrow("%s: primary key must be an integer or factor column in raw observed data", name);
        }
      }

      if (hasWeight() && rd.rawCols[weightCol].type != COLUMNDATA_NUMERIC) {
        mxThrow("%s: weight must be a numeric column in raw observed data", name);
      }

      if (hasFreq() && rd.rawCols[freqCol].type != COLUMNDATA_INTEGER) {
        mxThrow("%s: frequency must be an integer column in raw observed data", name);
      }
    } else {
      if(OMX_DEBUG) {mxLog("Data contains a matrix.");}
      od->dataMat = omxNewMatrixFromRPrimitive0(dataLoc, state, 0, 0);

      if (od->dataMat->colMajor && strEQ(od->u_type, "raw")) {
        omxToggleRowColumnMajor(od->dataMat);
      }
      od->cols = od->dataMat->cols;
      unfiltered.rows = od->dataMat->rows;
      od->numNumeric = od->cols;

      // Raw data is already stored like a data frame
      if (strEQ(u_type, "raw")) convertToDataFrame();
    }
	}

	if (!strEQ(u_type, "raw") && (hasPrimaryKey() || hasWeight() || hasFreq())) {
		mxThrow("%s: weight, freq, or primary key require raw data", name);
	}

	if(OMX_DEBUG) {mxLog("Processing Means Matrix.");}
	{ScopedProtect p1(dataLoc, R_do_slot(dataObj, Rf_install("means")));
	od->meansMat = omxNewMatrixFromRPrimitive0(dataLoc, state, 0, 0);
	}

	if(OMX_DEBUG) {
	        if(od->meansMat == NULL) {mxLog("No means found.");}
		else {omxPrint(od->meansMat, "Means Matrix is:");}
	}
	{
		ProtectedSEXP RdataLoc(R_do_slot(dataObj, Rf_install("numObs")));
		od->numObs = Rf_asReal(RdataLoc);
	}

	if(OMX_DEBUG) {mxLog("Processing Asymptotic Covariance Matrix.");}
	if (strEQ(od->u_type, "acov")) {  // old style for backward compatibility
		oss = std::unique_ptr< obsSummaryStats >(new obsSummaryStats);
		auto &o1 = *oss;
		o1.covMat = omxCreateCopyOfMatrix(dataMat, state);
		o1.covMat->colnames = dataMat->colnames;
		o1.covMat->rownames = dataMat->rownames;
		o1.meansMat = meansMat;
		meansMat = 0;
		// slopeMat unimplemented for legacy version
		ProtectedSEXP Racov(R_do_slot(dataObj, Rf_install("acov")));
		o1.useWeight = omxNewMatrixFromRPrimitive(Racov, state, 0, 0);
		ProtectedSEXP Rfw(R_do_slot(dataObj, Rf_install("fullWeight")));
		o1.asymCov = omxNewMatrixFromRPrimitive0(Rfw, state, 0, 0);
		ProtectedSEXP Rthr(R_do_slot(dataObj, Rf_install("thresholds")));
		o1.thresholdMat = omxNewMatrixFromRPrimitive0(Rthr, state, 0, 0);
	}
	if (R_has_slot(dataObj, Rf_install("observedStats"))) {
    if (!std::isfinite(numObs))
      mxThrow("%s: numObs is required when using observedStats", name);
		ProtectedSEXP RobsStats(R_do_slot(dataObj, Rf_install("observedStats")));
		ProtectedSEXP RobsStatsName(Rf_getAttrib(RobsStats, R_NamesSymbol));
		if (Rf_length(RobsStats)) oss = std::unique_ptr< obsSummaryStats >(new obsSummaryStats);
    int apiVersion = NA_INTEGER;
		for (int ax=0; ax < Rf_length(RobsStats); ++ax) {
			const char *key = R_CHAR(STRING_ELT(RobsStatsName, ax));
			auto &o1 = *oss;
			if (strEQ(key, "cov")) {
				o1.covMat = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
			} else if (strEQ(key, "slope")) {
				o1.slopeMat = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
				if (int(o1.slopeMat->colnames.size()) != o1.slopeMat->cols)
					mxThrow("%s: observedStats$slope must have colnames", name);
			} else if (strEQ(key, "means")) {
				o1.meansMat = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
			} else if (strEQ(key, "asymCov")) {
        if (apiVersion == NA_INTEGER || apiVersion == 1) apiVersion = 1;
        else mxThrow("%s: observedStats contains 'asymCov'; i'm confused", name);
				o1.asymCov = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
        o1.asymCov->unshareMemoryWithR();
        EigenMatrixAdaptor ac(o1.asymCov);
        ac *= od->numObs;
			} else if (strEQ(key, "useWeight")) {
        if (apiVersion == NA_INTEGER || apiVersion == 1) apiVersion = 1;
        else mxThrow("%s: observedStats contains 'useWeight'; i'm confused", name);
				o1.useWeight = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
			} else if (strEQ(key, "acov")) {
        if (apiVersion == NA_INTEGER || apiVersion == 0) apiVersion = 0;
        else mxThrow("%s: observedStats contains 'acov'; i'm confused", name);
				o1.useWeight = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
        Rf_warning("%s: observedStats item 'acov' is deprecated", name);
			} else if (strEQ(key, "fullWeight")) {
        if (apiVersion == NA_INTEGER || apiVersion == 0) apiVersion = 0;
        else mxThrow("%s: observedStats contains 'fullWeight'; i'm confused", name);
				o1.asymCov = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
        Rf_warning("%s: observedStats item 'fullWeight' is deprecated", name);
			} else if (strEQ(key, "thresholds")) {
				o1.thresholdMat = omxNewMatrixFromRPrimitive(VECTOR_ELT(RobsStats, ax), state, 0, 0);
			} else if (strEQ(key, "numEstimatedEntries")) {
        // ignore, just diagnostic output
      } else if (strstr(key, ".vcov")) {
        // Ignore for now. Is the vcov regenerated?
			} else {
				Rf_warning("%s: observedStats key '%s' ignored", name, key);
			}
		}
    if (oss && !oss->covMat) oss = 0;
	}
	if (oss) {
		auto &o1 = *oss;
    o1.totalWeight = numObs; // may not have data available to recalc
		if (o1.thresholdMat) o1.numOrdinal = o1.thresholdMat->cols;
		if (!o1.covMat) mxThrow("%s: observedStats must include a covariance matrix", name);
		if (int(o1.covMat->colnames.size()) != o1.covMat->cols)
			mxThrow("%s: observedStats$cov must have colnames", name);
		if (o1.numOrdinal) {
			if (int(o1.thresholdMat->colnames.size()) != o1.thresholdMat->cols)
				mxThrow("%s: observedStats$thresholds must have colnames", name);
			EigenMatrixAdaptor Ethr(o1.thresholdMat);
			ColMapType thrMap;
			for (int cx=0; cx < int(o1.thresholdMat->colnames.size()); ++cx) {
				thrMap.emplace(o1.thresholdMat->colnames[cx], cx);
			}
			int foundOrd = 0;
			for (int cx=0; cx < int(o1.covMat->colnames.size()); ++cx) {
				const char *cn = o1.covMat->colnames[cx];
				auto it = thrMap.find(cn);
				omxThresholdColumn tc;
				tc.isDiscrete = false;
				tc.dataColumn = cx;
				if (it != thrMap.end()) {
					tc.column = it->second;
					int numThr = 0;
					while (numThr < o1.thresholdMat->rows &&
					       std::isfinite(Ethr(numThr, it->second))) ++numThr;
					tc.numThresholds = numThr;
					foundOrd += 1;
				}
				o1.thresholdCols.push_back(tc);
			}
			if (foundOrd != o1.numOrdinal) mxThrow("%s: cannot match all threshold columns", name);
		}
	}

	if (R_has_slot(dataObj, Rf_install("algebra"))) {
		ProtectedSEXP Ralg(R_do_slot(dataObj, Rf_install("algebra")));
		int len = Rf_length(Ralg);
		if (len) {
			algebra.reserve(len);
			int *aptr = INTEGER(Ralg);
			for (int ax=0; ax < len; ++ax) algebra.push_back(aptr[ax]);
		}
	}
}

void omxData::prep()
{
  //mxLog("omxData::prep naAction=%d", naAction);
	switch (naAction) {
	case NA_PASS:
		filtered = unfiltered;
		break;
	case NA_FAIL:
		unfiltered.refreshHasNa();
		if (std::any_of(unfiltered.hasNa.begin(),
										unfiltered.hasNa.end(), [](bool na){ return na; })) {
			mxThrow("%d: contains at least one NA and naAction='fail'", name);
		}
		filtered = unfiltered;
		break;
	case NA_OMIT:{ // remove rows
		if (filtered.rawCols.size() != unfiltered.rawCols.size()) {
			filtered = unfiltered;
		}
		unfiltered.refreshHasNa();
		int rows = std::count(unfiltered.hasNa.begin(),
													unfiltered.hasNa.end(), false);
		if (filtered.rows != rows || filtered.hasNa != unfiltered.hasNa) {
			if (verbose >= 1) mxLog("omit: NA pattern changed, clearing cache");
			filtered.rows = rows;
			filtered.hasNa = unfiltered.hasNa;
			for (auto &cd : filtered.rawCols) {
				cd.clear();
			}
			oss.reset();
		}
		int filterCount = 0;
		for (int cx=0; cx < int(unfiltered.rawCols.size()); ++cx) {
			auto &src = unfiltered.rawCols[cx];
			auto &dest = filtered.rawCols[cx];
			if (dest.d()) continue;
			filterCount += 1;
      dest.setMinValue(src.getMinValue());
      dest.setMaxValue(src.getMaxValueUnsafe());
			switch (src.type) {
			case COLUMNDATA_ORDERED_FACTOR:
			case COLUMNDATA_UNORDERED_FACTOR:
			case COLUMNDATA_INTEGER:
				dest.setOwn(new int[rows]);
				for (int sx=0, dx=0; sx < unfiltered.rows; ++sx) {
					if (unfiltered.hasNa[sx]) continue;
					dest.i()[dx++] = src.i()[sx];
				}
				break;
			case COLUMNDATA_NUMERIC:
				dest.setOwn(new double[rows]);
				for (int sx=0, dx=0; sx < unfiltered.rows; ++sx) {
					if (unfiltered.hasNa[sx]) continue;
					dest.d()[dx++] = src.d()[sx];
				}
				break;
			default: ;
        // continue is fine here
			}
		}
		if (verbose >= 1) mxLog("omit: filtered %d columns", filterCount);
		break;}
	case NA_EXCLUDE:{ // set to zero frequency
		unfiltered.refreshHasNa();
    filtered = unfiltered; // resets minValue; can usually avoid? TODO
		if (!hasFreq() || int(filtered.rawCols.size()) <= freqCol) {
			if (verbose >= 1) mxLog("exclude: adding freq column");
			freqCol = filtered.rawCols.size();
			auto *newFreq = new int[filtered.rows];
			filtered.rawCols.emplace_back("freq", COLUMNDATA_INTEGER, newFreq);
			for (int rx=0; rx < filtered.rows; ++rx) {
				if (unfiltered.hasNa[rx]) newFreq[rx] = 0;
				else newFreq[rx] = 1;
			}
		} else {
			int *oldFreq = unfiltered.rawCols[freqCol].i();
      auto &dc = filtered.rawCols[freqCol];
			auto *newFreq = new int[filtered.rows];
			for (int rx=0; rx < filtered.rows; ++rx) {
				if (unfiltered.hasNa[rx]) newFreq[rx] = 0;
        else if (oldFreq) newFreq[rx] = oldFreq[rx];
				else newFreq[rx] = 1;
			}
      dc.setOwn(newFreq);
		}
    oss.reset();
		break;}
	default: mxThrow("unknown naAction %d", naAction);
	}

	currentWeightColumn = getWeightColumn();
	currentFreqColumn = getOriginalFreqColumn();

  if (OMX_DEBUG) {
    for (auto &col : filtered.rawCols) {
      // This check is too strict. It is legitimate to for the
      // data.frame to contain non-data columns.
      if (!col.i()) OOPS;
    }
  }

	if (verbose >= 2) omxPrintData("prep", 10);
}

void omxData::sanityCheck()
{
	if (hasPrimaryKey()) {
		for (int rx=0; rx < nrows(); ++rx) {
			int key = primaryKeyOfRow(rx);
			std::pair< std::map<int,int>::iterator, bool> ret =
				primaryKeyIndex.insert(std::pair<int,int>(key, rx));
			if (ret.second) continue;
			mxThrow("%s: primary keys are not unique (examine rows with key=%d)", name, key);
		}
	}

	if (currentFreqColumn) {
		for (int rx=0; rx < nrows(); ++rx) {
			if (currentFreqColumn[rx] >= 0) continue;
			mxThrow("%s: cannot proceed with non-positive frequency %d for row %d",
				 name, currentFreqColumn[rx], 1+rx);
		}
	}
}

double *omxData::getWeightColumn()
{
	if (!hasWeight()) return 0;
	if (isRaw()) {
		return rawCol(weightCol).d();
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
	if (freqCol < 0) return 0;
	return rawCol(freqCol).i();
}

int omxData::numRawRows()
{
	if (strEQ(getType(), "raw")) return nrows();
	return 0;
}

omxData* omxState::omxNewDataFromMxData(SEXP dataObj, const char *name)
{
	if(dataObj == NULL) {
		mxThrow("Null Data Object detected.  This is an internal error, and should be reported on the forums.\n");
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
	if (strEQ(dclass, "MxDataStatic") || strEQ(dclass, "MxDataLegacyWLS")) od->newDataStatic(this, dataObj);
	else if (strcmp(dclass, "MxDataDynamic")==0) newDataDynamic(dataObj, od);
	else mxThrow("Unknown data class %s", dclass);
	od->prep();
	od->sanityCheck();
	return od;
}

void omxData::RawData::refreshHasNa()
{
	hasNa.resize(rows);
	for (int rx=0; rx < rows; ++rx) {
		bool na = false;
		for (auto &cd : rawCols) {
			switch (cd.type) {
			case COLUMNDATA_INVALID: continue;
			case COLUMNDATA_ORDERED_FACTOR:
			case COLUMNDATA_UNORDERED_FACTOR:
			case COLUMNDATA_INTEGER:
				na |= cd.i()[rx] == NA_INTEGER;
				break;
			case COLUMNDATA_NUMERIC:
				na |= !std::isfinite(cd.d()[rx]);
				break;
			}
			hasNa[rx] = na;
		}
	}
}

void omxData::RawData::clearColumn(int col)
{
	rawCols[col].clear();
}

void omxData::RawData::clear()
{
	rawCols.clear();
	rows = 0;
}

omxData::RawData::~RawData()
{ clear(); }

obsSummaryStats::~obsSummaryStats()
{
	omxFreeMatrix(covMat);
	omxFreeMatrix(slopeMat);
	omxFreeMatrix(meansMat);
	omxFreeMatrix(asymCov);
	omxFreeMatrix(useWeight);
	omxFreeMatrix(thresholdMat);
}

void omxFreeData(omxData* od)
{
	omxFreeMatrix(od->dataMat);
	omxFreeMatrix(od->meansMat);
	delete od;
}

bool omxDataElementMissing(omxData *od, int row, int col)
{
	if(od->dataMat != NULL) {
		return std::isnan(omxMatrixElement(od->dataMat, row, col));
	}
	ColumnData &cd = od->rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) return std::isnan(cd.d()[row]);
	else return cd.i()[row] == NA_INTEGER;
}

double omxDoubleDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		return omxMatrixElement(od->dataMat, row, col);
	}
	ColumnData &cd = od->rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) return cd.d()[row];
	else return cd.i()[row];
}

double *omxDoubleDataColumn(omxData *od, int col)
{
	ColumnData &cd = od->rawCol(col);
	if (cd.type != COLUMNDATA_NUMERIC) mxThrow("Column '%s' is integer, not real", cd.name);
	else return cd.d();
}

int omxIntDataElement(omxData *od, int row, int col) {
	if(od->dataMat != NULL) {
		return (int) omxMatrixElement(od->dataMat, row, col);
	}

	ColumnData &cd = od->rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) return cd.d()[row];
	else return cd.i()[row];
}

omxMatrix* omxDataCovariance(omxData *od)
{
	if (od->dataMat) return od->dataMat;

	if (od->expectation.size()) {
		omxExpectation *ex = od->expectation[0];
		return omxGetExpectationComponent(ex, "covariance");
	}

	mxThrow("%s: type='%s' data must be in matrix storage", od->name, od->u_type);
}

bool omxData::columnIsFactor(int col)
{
	if(dataMat != NULL) return FALSE;
	ColumnData &cd = rawCol(col);
	return cd.type == COLUMNDATA_ORDERED_FACTOR || cd.type == COLUMNDATA_INTEGER;
}

bool omxDataColumnIsKey(omxData *od, int col)
{
	if(od->dataMat != NULL) return FALSE;
	ColumnData &cd = od->rawCol(col);
	return cd.type != COLUMNDATA_NUMERIC;
}

void omxData::RawData::assertColumnIsData(int col, OmxDataType dt, bool isUnfiltered)
{
  //mxLog("omxData::RawData::assertColumnIsData(%d, %d, %d)", col, dt, isUnfiltered);
	if (col < 0 || col >= int(rawCols.size())) {
		mxThrow("Column %d requested but only %d columns of data",
						col, int(rawCols.size()));
	}
	ColumnData &cd = rawCols[col];
	switch (cd.type) {
	case COLUMNDATA_ORDERED_FACTOR:
		if (dt == OMXDATA_ORDINAL || dt == OMXDATA_COUNT) {
      if (!isUnfiltered) cd.setZeroMinValue(rows);
      return;
    }
		mxThrow("Don't know how to interpret factor column '%s' as numeric.\n"
						"You may want to specify thresholds for your model like this: "
						"mxThreshold(vars='%s', nThresh=%d)",
						cd.name, cd.name, cd.getNumThresholds());
	case COLUMNDATA_NUMERIC:{
		if (dt == OMXDATA_REAL) return;
    if (dt == OMXDATA_ORDINAL) {
      mxThrow("Don't know how to interpret numeric column '%s' as ordinal", cd.name);
    }
		// convert for dt == OMXDATA_COUNT
		cd.type = COLUMNDATA_INTEGER;
		double *realData = cd.d();
    int *intData = new int[rows];
		for (int rx=0; rx < rows; ++rx) {
      intData[rx] = cast_with_NA(realData[rx]);
		}
    cd.setOwn(intData);
    cd.setMinValue(0);
    cd.verifyMinValue(rows);
    if (!isUnfiltered) cd.setMaxValueFromData(rows);
    break;}
	case COLUMNDATA_UNORDERED_FACTOR:
		if (dt == OMXDATA_ORDINAL) {
			if (isUnfiltered && ++Global->dataTypeWarningCount < 5) {
				Rf_warning("Column '%s' must be an ordered factor. "
									 "Please use mxFactor()", cd.name);
			}
      if (!isUnfiltered) cd.setZeroMinValue(rows);
		} else if (dt == OMXDATA_COUNT) {
			mxThrow("Don't know how to interpret unordered factor '%s' as a count", cd.name);
		} else {
			mxThrow("Don't know how to interpret unordered factor '%s' as numeric", cd.name);
		}
		break;
	case COLUMNDATA_INTEGER:{
		if (dt == OMXDATA_COUNT) {
      cd.verifyMinValue(rows);
      if (!isUnfiltered) cd.setMaxValueFromData(rows);
      return;
    }
		if (dt == OMXDATA_ORDINAL) {
			mxThrow("Don't know how to interpret integer column '%s' as ordinal. "
							"Please use mxFactor()", cd.name);
		}
		// convert for dt == OMXDATA_REAL
		cd.type = COLUMNDATA_NUMERIC;
		int *intData = cd.i();
    double *realData = new double[rows];
		for (int rx=0; rx < rows; ++rx) {
      realData[rx] = cast_with_NA(intData[rx]);
		}
    cd.setOwn(realData);
		return;}
	default:
		mxThrow("Column '%s' is an unknown data type", cd.name);
	}
}

void omxData::assertColumnIsData(int col, OmxDataType dt)
{
	if (dataMat) return;
  if (verbose) {
    ColumnData &cd = unfiltered.rawCols[col];
    mxLog("%s: assertColumnIsData(%s(%d), %d)", name, cd.name, col, dt);
  }
	unfiltered.assertColumnIsData(col, dt, true);
	filtered.assertColumnIsData(col, dt, false);
}

int omxData::primaryKeyOfRow(int row)
{
	if(dataMat != NULL) mxThrow("%s: only raw data can have a primary key", name);
	ColumnData &cd = rawCol(primaryKey);
	return cd.i()[row];
}

int omxData::lookupRowOfKey(int key)
{
	const std::map<int,int>::iterator it = primaryKeyIndex.find(key);
	if (it == primaryKeyIndex.end()) {
		if (!hasPrimaryKey()) {
			mxThrow("%s: attempt to lookup key=%d but no primary key", name, key);
		}
		ColumnData &cd = rawCol(primaryKey);
		mxThrow("%s: key %d not found in column '%s'", name, key, cd.name);
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
	ColumnData &cd = rawCol(col);
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
	default: mxThrow("type %d unknown", cdt);
	}
}

const char *ColumnData::typeName()
{
	return ColumnDataTypeToString(type);
}

ColumnData ColumnData::clone() const
{
  ColumnData ret(name);
  ret.type = type;
  ret.setBorrow(ptr);
  ret.levelNames = levelNames;
  ret.minValue = minValue;
  ret.maxValue = maxValue;
  return ret;
}

void ColumnData::verifyMinValue(int nrows)
{
  int m0 = std::numeric_limits<int>::max();
  for (int *it=i(); it < i()+nrows; ++it) {
    if (*it == NA_INTEGER) continue;
    if (*it < m0) m0 = *it;
  }
  if (m0 != minValue) {
    mxThrow("column %s: minimum value is %d not %d", name, m0, minValue);
  }
}

void ColumnData::setMaxValueFromData(int nrows)
{
  // NA_INTEGER is negative so can ignore
  auto it = std::max_element(i(), i() + nrows);
  maxValue = *it;
}

void ColumnData::setZeroMinValue(int rows)
{
  if (minValue == 0) return;
  if (type == COLUMNDATA_NUMERIC)
    mxThrow("ColumnData::setZeroMinValue not implemented for numeric data");
  if (minValue != 1) OOPS;
  int *oldIntData = i();
  int *id = new int[rows];
  for (int xx=0; xx < rows; ++xx) {
    if (oldIntData[xx] == NA_INTEGER) {
      id[xx] = NA_INTEGER;
    } else {
      id[xx] = oldIntData[xx] - 1;
    }
  }
  setOwn(id);
  minValue = 0;
  if (maxValue == NA_INTEGER) OOPS;
  maxValue -= 1;
}

void omxData::RawData::operator=(const RawData &other)
{
  rawCols.clear();
  for (auto &rc : other.rawCols) {
    rawCols.push_back(rc.clone());
  }
  hasNa = other.hasNa;
  rows = other.rows;
}

void omxDataKeysCompatible(omxData *upper, omxData *lower, int foreignKey)
{
	ColumnData &lcd = lower->rawCol(foreignKey);
	if (!upper->hasPrimaryKey()) {
		mxThrow("Attempt to join foreign key '%s' in %s of type '%s' with"
			 " %s which has no primary key declared",
			 lcd.name, lower->name, ColumnDataTypeToString(lcd.type), upper->name);
	}
	ColumnData &ucd = upper->rawCol(upper->primaryKey);
	if (ucd.type != lcd.type) {
		mxThrow("Primary key '%s' in %s of type '%s' cannot be joined with"
			 " foreign key '%s' in %s of type '%s'",
			 ucd.name, upper->name, ColumnDataTypeToString(ucd.type),
			 lcd.name, lower->name, ColumnDataTypeToString(lcd.type));
	}
	if (ucd.type == COLUMNDATA_ORDERED_FACTOR || ucd.type == COLUMNDATA_UNORDERED_FACTOR) {
		if (ucd.getMaxValue() != lcd.getMaxValue()) {
			mxThrow("Primary key '%s' in %s has a different number of factor"
				 " levels compared to foreign key '%s' in %s",
				 ucd.name, upper->name, lcd.name, lower->name);
		}
		for (int lx=0; lx < int(ucd.levelNames.size()); ++lx) {
			auto &ul = ucd.levelNames[lx];
			auto &ll = lcd.levelNames[lx];
			if (ul == ll) continue;
			mxThrow("Primary key '%s' in %s has different factor levels ('%s' != '%s')"
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
	if(row >= od->nrows()) mxThrow("Invalid row");

	if(om == NULL) mxThrow("Must provide an output matrix");

	if (om->cols < len) mxThrow("omxContiguousDataRow: output matrix is too small");
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
	return od->u_type;
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
	buf += string_snprintf("%s(%s): %f observations %d x %d\n", header, od->u_type, od->numObs,
												 nrows(), cols);
	if (hasPrimaryKey()) {
		buf += string_snprintf("primaryKey %d\n", od->primaryKey);
	}

	buf += string_snprintf("Row consists of %d numeric, %d ordered factor:", od->numNumeric, od->numFactor);

	int upto = nrows();
	if (maxRows >= 0) upto = std::min(upto, maxRows);

	if (isRaw()) {
		auto &rd = filtered;
		for (auto &cd : rd.rawCols) {
			buf += " ";
			buf += cd.name;
			if (cd.type != COLUMNDATA_NUMERIC) {
				buf += "[I]";
			} else {
				buf += "[N]";
			}
		}
		buf += "\n";

		for (int vxv=0; upto > 0; vxv++) {
			int vx = permute? permute[vxv] : vxv;
			if (hasFreq() && getRowFreq(vx) == 0) continue;
			upto -= 1;
			for (int j = 0; j < int(rd.rawCols.size()); j++) {
				ColumnData &cd = rd.rawCols[j];
				if (cd.type == COLUMNDATA_INVALID) continue;
				if (cd.type != COLUMNDATA_NUMERIC) {
					int *val = cd.i();
					if (!val || val[vx] == NA_INTEGER) {
						buf += " NA,";
					} else {
						buf += string_snprintf(" %d,", val[vx]);
					}
				} else {
					double *val = cd.d();
					if (!val || !std::isfinite(val[vx])) {
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
		auto &dv = defVars[k];
		double newDefVar;
		if (dv.column == weightCol) {
			newDefVar = getRowWeight(row);
		} else if (dv.column == freqCol) {
			newDefVar = getRowFreq(row);
		} else {
			newDefVar = omxDoubleDataElement(this, row, dv.column);
		}
		changed |= dv.loadData(state, newDefVar);
	}
	if (changed && OMX_DEBUG_ROWS(row)) { mxLog("%s: load def vars for row %d", name, row); }
	return changed;
}

double omxData::rowMultiplier(int rx)
{
	double ww = 1.0;
	double *rowWeight = getWeightColumn();
	int *rowFreq = getFreqColumn();
	if (rowWeight) ww *= rowWeight[rx];
	if (rowFreq) ww *= rowFreq[rx];
	return ww;
}

bool omxData::containsNAs(int col)
{
	int rows = nrows();
	if(dataMat != NULL) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(omxMatrixElement(dataMat, rx, col))) continue;
			return true;
		}
		return false;
	}

	if (col == weightCol || col == freqCol) return false;

	ColumnData &cd = rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(cd.d()[rx]) ||
					rowMultiplier(rx) == 0) continue;
			return true;
		}
	} else {
		for (int rx=0; rx < rows; ++rx) {
			if (cd.i()[rx] != NA_INTEGER ||
					rowMultiplier(rx) == 0) continue;
			return true;
		}
	}
	return false;
}

double omxData::countObs(int col)
{
  double obs = 0;
	int rows = nrows();
	if (dataMat) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(omxMatrixElement(dataMat, rx, col))) obs += 1;
		}
		return obs;
	}

	if (col == weightCol || col == freqCol) return false;

	ColumnData &cd = rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) {
		for (int rx=0; rx < rows; ++rx) {
			if (std::isfinite(cd.d()[rx])) obs += rowMultiplier(rx);
		}
	} else {
		for (int rx=0; rx < rows; ++rx) {
			if (cd.i()[rx] != NA_INTEGER) obs += rowMultiplier(rx);
		}
	}
	return obs;
}

void omxData::prohibitFactor(int col)
{
	if (!isRaw()) return;
	if (col == weightCol || col == freqCol) return;
	auto &cd = rawCol(col);
	if (cd.type == COLUMNDATA_NUMERIC) return;
	if (cd.type != COLUMNDATA_INTEGER) {
		Rf_warning("%s: definition variable '%s' is of type '%s';"
			   " note that it will be treated as integer (as is done by ?unclass)."
			   " Is this really what you want to do? Really?",
			   name, columnName(col), cd.typeName());
	}
}

void omxData::prohibitNAdefVar(int col)
{
	if (!containsNAs(col)) return;
	if (!dataMat) {
		if (col == weightCol) {
			mxThrow("%s: NA in row weights", name);
		}
		if (col == freqCol) {
			mxThrow("%s: NA in row frequencies", name);
		}
	}
	mxThrow("%s: NA in definition variable '%s'",
		 name, omxDataColumnName(this, col));
}

void omxData::loadFakeData(omxState *state, double fake)
{
	for (int dx=0; dx < int(defVars.size()); ++dx) {
		auto &dv = defVars[dx];
		dv.loadData(state, fake);
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

static int plookup(ColMapType &map, const char *str)
{
	auto it = map.find(str);
	if (it == map.end()) mxThrow("Can't find '%s'", str);
	return it->second;
}

void obsSummaryStats::setDimnames(omxData *data)
{
	colMap.clear();
	for (int cx=0; cx < int(dc.size()); ++cx) colMap.emplace(dc[cx], cx);

	if (int(dc.size()) != covMat->cols)
		mxThrow("%s: internal error; dc.size() %d != covMat->cols %d",
			 data->name, int(dc.size()), covMat->cols);
	covMat->colnames.resize(covMat->cols);
	covMat->rownames.resize(covMat->cols);
	for (int cx=0; cx < covMat->cols; ++cx) {
		covMat->colnames[cx] = dc[cx];
		covMat->rownames[cx] = dc[cx];
	}

	if (slopeMat) {
		slopeMat->colnames.resize(exoPred.size());
		for (int cx=0; cx < int(exoPred.size()); ++cx) {
			slopeMat->colnames[cx] = data->columnName(exoPred[cx]);
		}
		slopeMat->rownames.resize(covMat->cols);
		for (int cx=0; cx < covMat->cols; ++cx) {
			slopeMat->rownames[cx] = dc[cx];
		}
	}

	if (thresholdMat) {
		thresholdMat->colnames.resize(thresholdMat->cols);
		for (auto &th : thresholdCols) {
			if (!th.numThresholds) continue;
			thresholdMat->colnames[th.column] = dc[th.dataColumn];
		}
	}

	if (asymCov) {
		asymCov->clearDimnames();
		asymCov->colnames.reserve(asymCov->cols);
		if (thresholdMat || meansMat) {
			for (auto &tc : thresholdCols) {
				if (tc.numThresholds == 0) {
					asymCov->colnames.push_back(strdup(dc[tc.dataColumn]));
				} else {
					for (int th=1; th <= tc.numThresholds; ++th) {
						auto str = string_snprintf("%st%d", dc[tc.dataColumn], th);
						asymCov->colnames.push_back(strdup(str.c_str()));
					}
				}
			}
		}
		// slopes TODO
		for (int cx=0; cx < covMat->cols; ++cx) {
			if (thresholdMat && thresholdCols[cx].numThresholds) continue;
			auto str = string_snprintf("var_%s", dc[cx]);
			asymCov->colnames.push_back(strdup(str.c_str()));
		}
		for (int cx=0; cx < covMat->cols-1; ++cx) {
			for (int rx=cx+1; rx < covMat->cols; ++rx) {
				auto str = string_snprintf("poly_%s_%s", dc[rx], dc[cx]);
				asymCov->colnames.push_back(strdup(str.c_str()));
			}
		}
		asymCov->freeColnames = true;
		asymCov->rownames = asymCov->colnames;
	}
}

void obsSummaryStats::permute(omxData *data)
{
	covMat->unshareMemoryWithR();
	if (meansMat) meansMat->unshareMemoryWithR();
	if (asymCov) asymCov->unshareMemoryWithR();
	if (useWeight) useWeight->unshareMemoryWithR();

	ColMapType dataMap;
	for (int cx=0; cx < int(dc.size()); ++cx) dataMap.emplace(dc[cx], cx);

	Eigen::VectorXi invDataColumns(dc.size()); // data -> expectation order
	if (int(covMat->colnames.size()) != covMat->cols) mxThrow("%s: cannot permute without cov dimnames", data->name);
	for (int cx=0; cx < int(covMat->colnames.size()); ++cx) {
		auto it = dataMap.find(covMat->colnames[cx]);
		if (it == dataMap.end()) OOPS;
		invDataColumns[cx] = it->second;
		//mxLog("%d %s", cx, omxDataColumnName(data, dc[cx]));
	}

	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> p1(invDataColumns);

	ColMapType acovMap;
	if (int(asymCov->colnames.size()) != asymCov->cols) mxThrow("%s: cannot permute without acov dimnames", data->name);
	for (int cx=0; cx < int(asymCov->colnames.size()); ++cx) {
		//mxLog("%s -> %d", asymCov->colnames[cx], cx);
		acovMap.emplace(asymCov->colnames[cx], cx);
	}

	auto &thresh = thresholdCols;

	Eigen::PermutationMatrix<Eigen::Dynamic> p2(asymCov->cols);
	int px = 0;

	ColMapType thrMap;
	if (thresh.size() == 0 && meansMat) {
		for (int cx=0; cx < int(dc.size()); ++cx) {
			p2.indices()[px++] = plookup(acovMap, dc[cx]);
		}
	} else if (thresh.size()) {
		for (int cx=0; cx < int(thresh.size()); ++cx) {
			thrMap.emplace(covMat->colnames[cx], cx);
			//mxLog("%s -> %d", covMat->colnames[cx], cx);
		}
		for (int cx=0; cx < int(dc.size()); ++cx) {
			auto cn = dc[cx];
			auto &t1 = thresh[ thrMap[cn] ];
			if (t1.numThresholds) {
				for (int tx=1; tx <= t1.numThresholds; ++tx) {
					auto s1 = string_snprintf("%st%d", cn, tx);
					p2.indices()[px++] = plookup(acovMap, s1.c_str());
				}
			} else {
				p2.indices()[px++] = plookup(acovMap, cn);
			}
		}
	}

	for (int cx=0; cx < int(dc.size()); ++cx) {
		auto &t1 = thresh[ thrMap[dc[cx]] ];
		if (t1.numThresholds) continue;
		auto cn = dc[cx];
		auto s1 = string_snprintf("var_%s", cn);
		p2.indices()[px++] = plookup(acovMap, s1.c_str());
	}

	for (int cx=0; cx < int(dc.size())-1; ++cx) {
		auto cn = dc[cx];
		for (int rx=cx+1; rx < int(dc.size()); ++rx) {
			auto rn = dc[rx];
			auto s1 = string_snprintf("poly_%s_%s", rn, cn);
			auto s2 = string_snprintf("poly_%s_%s", cn, rn);
			auto it = acovMap.find(s1.c_str());
			if (it == acovMap.end()) {
				it = acovMap.find(s2.c_str());
			}
			if (it == acovMap.end()) mxThrow("Can't find '%s' or '%s'",
							 s1.c_str(), s2.c_str());
			p2.indices()[px++] = it->second;
		}
	}

	EigenVectorAdaptor Emean(meansMat);
	EigenMatrixAdaptor Ecov(covMat);
	EigenMatrixAdaptor Eacov(asymCov);

	//Eigen::MatrixXi mp1 = p1;
	//mxPrintMat("p1", mp1);
	//mxPrintMat("p2", p2.indices());

	Emean.derived() = (p1 * Emean).eval();
	Ecov.derived() = (p1 * Ecov * p1.transpose()).eval();
	Eacov.derived() = (p2.transpose() * Eacov * p2).eval();
  if (useWeight) {
    EigenMatrixAdaptor Euw(useWeight);
    Euw.derived() = (p2.transpose() * Euw * p2).eval();
  }

	for (auto &th : thresh) th.dataColumn = invDataColumns[th.dataColumn];
	std::sort(thresh.begin(), thresh.end(),
		  [](const omxThresholdColumn &a, const omxThresholdColumn &b) -> bool
		  { return a.dataColumn < b.dataColumn; });
}

void obsSummaryStats::log()
{
	mxLog("totalWeight %f numOrdinal %d", totalWeight, numOrdinal);
	if (covMat) omxPrint(covMat, "cov");
	if (slopeMat) omxPrint(slopeMat, "slope");
	if (meansMat) omxPrint(meansMat, "mean");
	if (asymCov) {
    EigenMatrixAdaptor Eacov(asymCov);
    if (Eacov.cols() < 30) {
      mxPrintMat("asymCov", Eacov);
    } else {
      mxPrintMat("asymCov (topleft)", Eacov.block(0,0,30,30));
    }
  }
	if (useWeight) {
    EigenMatrixAdaptor Efw(useWeight);
    if (Efw.cols() < 30) {
      mxPrintMat("useW", Efw);
    } else {
      mxPrintMat("useW (topleft)", Efw.block(0,0,30,30));
    }
	}
	for (auto &th : thresholdCols) { th.log(); }
	if (thresholdMat) omxPrint(thresholdMat, "thr");
}

int obsSummaryStats::numPredictors(int vx)
{
    int numTh = thresholdCols[vx].numThresholds;
    if (numTh == 0) numTh = 1;
    int ef = exoFree.row(vx).sum();
    return numTh + ef;
}

std::string omxData::getExoPredictorName(int vx, int nx)
{
  auto &o1 = *oss;

  if (!(0 <= nx && nx < o1.numPredictors(vx))) mxThrow("nx %d out of range for vx %d", nx, vx);

  int numTh = o1.thresholdCols[vx].numThresholds;
  if (numTh == 0) { // continuous
    if (nx == 0) return "(intercept)";
    nx -= 1;
  } else {  // ordinal
    if (nx < numTh) return string_snprintf("th%d", 1+nx);
    nx -= numTh;
  }

  for (int xx=0, px=0; xx < int(o1.exoPred.size()); ++xx) {
    if (!o1.exoFree(vx, xx)) continue;
    if (nx == px) return columnName(o1.exoPred[xx]);
    px += 1;
  }
  return "unknown";
}

void omxData::reportResults(MxRList &out)
{
	out.add("numObs", Rf_ScalarReal(omxDataNumObs(this)));

	auto &rc = filtered.rawCols;
	int numValid = std::count_if(rc.begin(), rc.end(),
				     [](ColumnData &c1)->bool{ return c1.type != COLUMNDATA_INVALID; });
	if (isModified() && numValid) {
		int rows = filtered.rows;
		Rcpp::CharacterVector colNames(numValid);
		Rcpp::List columns(numValid);
		for (int cx=0, dx=0; cx < int(rc.size()); ++cx) {
			auto &c1 = rc[cx];
			if (c1.type == COLUMNDATA_INVALID) continue;
			colNames[dx] = c1.name;
			if (c1.type == COLUMNDATA_NUMERIC) {
				Eigen::Map< Eigen::VectorXd > vec(c1.d(), rows);
				columns[dx] = Rcpp::wrap(vec);
			} else {
				Eigen::Map< Eigen::VectorXi > vec(c1.i(), rows);
				columns[dx] = Rcpp::wrap(vec);
			}
			dx += 1;
		}
		columns.attr("names") = colNames;
		markAsDataFrame(columns, rows);
		out.add("rawData", columns);
	}

	if (!oss) return;
	auto &o1 = *oss;
	if (!o1.output) return;

	if (o1.covMat) out.add("cov", o1.covMat->asR());
	if (o1.meansMat) out.add("means", o1.meansMat->asR());
	if (o1.asymCov) {
    double Ninv = 1 / omxDataNumObs(this);
    EigenMatrixAdaptor ac(o1.asymCov);
    ac *= Ninv;
    out.add("asymCov", o1.asymCov->asR());
    ac *= omxDataNumObs(this);
  }
	if (o1.slopeMat) out.add("slope", o1.slopeMat->asR());
	if (o1.useWeight) out.add("useWeight", o1.useWeight->asR());
	if (o1.thresholdMat) out.add("thresholds", o1.thresholdMat->asR());
  out.add("numEstimatedEntries", Rcpp::wrap(numEstimatedEntries));

  for (int vx=0; vx < int(o1.perVar.size()); ++vx) {
    auto &pv = o1.perVar[vx];
    if (pv.vcov.rows() == 0) continue;
    std::string key = o1.dc[vx];
    key += ".vcov";
    int ef = 0;
    if (o1.exoPred.size()) ef = o1.exoFree.row(vx).sum();
    int numTh = o1.thresholdCols[vx].numThresholds;
    if (numTh == 0) numTh = 1;
    if (pv.vcov.rows() != numTh+ef)
      mxThrow("%s pv.vcov.rows() %d != ef %d", key.c_str(), pv.vcov.rows(), ef);
    Rcpp::NumericMatrix vcov = Rcpp::wrap(pv.vcov);
    Rcpp::CharacterVector ch(pv.vcov.rows());
    int nx=0;
    if (o1.thresholdCols[vx].numThresholds == 0) {
      ch(nx++) = "(intercept)";
    } else {
      for (int xx=0; xx < numTh; ++xx) {
        std::string tname = string_snprintf("th%d", 1+xx);
        ch(nx++) = tname.c_str();
      }
    }
    for (int xx=0; xx < int(o1.exoPred.size()); ++xx) {
      if (!o1.exoFree(vx, xx)) continue;
      ch(nx++) = columnName(o1.exoPred[xx]);
    }
    colnames(vcov) = ch;
    rownames(vcov) = ch;
    out.add(key.c_str(), vcov);
  }
}

template <typename T>
void getContRow(std::vector<ColumnData> &df,
		int row,
		const Eigen::Ref<const DataColumnIndexVector> &dc,
		Eigen::MatrixBase<T> &out)
{
	for (int cx=0; cx < dc.size(); ++cx) {
		auto &cd = df[ dc[cx] ];
		out[cx] = cd.d()[row];
	}
}

void omxData::wlsAllContinuousCumulants(omxState *state)
{
	auto &o1 = *oss;
	const char *wlsType = o1.wlsType;
	const std::vector<const char *> &dc = o1.dc;
	const Eigen::ArrayXd &rowMult = o1.rowMult;
	std::vector<int> &index = o1.index;

	if (verbose >= 1) mxLog("%s: using wlsAllContinuousCumulants type=%s", name, wlsType);

	DataColumnIndexVector dci(dc.size());
	for (int cx=0; cx < int(dc.size()); ++cx) {
		int dx = rawColMap[dc[cx]];
		dci[cx] = dx;
		if (!containsNAs(dx)) continue;
		mxThrow("%s: all continuous data with missingness (column '%s') cannot "
			 "be handled using the cumulants method. Use na.omit(yourDataFrame) "
			 "to remove rows with missing values or use allContinuousMethod='marginals' "
			 "or use maximum likelihood", name, columnName(dx));
	}

	int numCols = dc.size();
	int numColsStar = numCols*(numCols+1)/2;

	double totalWeight = o1.totalWeight;
	o1.covMat = omxInitMatrix(numCols, numCols, state);
	o1.asymCov = omxInitMatrix(numColsStar, numColsStar, state);
	EigenMatrixAdaptor Ecov(o1.covMat);
	Eigen::VectorXd Emean(numCols);
	Ecov.setZero();
	Emean.setZero();

	Eigen::ArrayXXd data(int(index.size()), numCols);
	Eigen::VectorXd r1(numCols);
	for (int rx=0; rx < int(index.size()); ++rx) {
		int ix = index[rx];
		getContRow(filtered.rawCols, ix, dci, r1);
		Emean += r1 * rowMult[rx];
		Ecov += r1 * r1.transpose() * rowMult[rx];
		data.row(rx) = r1;
	}
	Emean /= totalWeight;
	Ecov -= totalWeight * Emean * Emean.transpose();
	Ecov /= totalWeight - 1;
	for (int cx=0; cx < int(dc.size()); ++cx) {
		if (Ecov(cx,cx) > minVariance) continue;
		mxThrow("%s: '%s' has observed variance less than %g",
			name, dc[cx], minVariance);
	}

	data.rowwise() -= Emean.array().transpose();
	Eigen::MatrixXd Vmat = Ecov * (totalWeight-1.) / totalWeight;
	EigenMatrixAdaptor Umat(o1.asymCov);

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
		for (int ix=jx; ix < numColsStar; ++ix) { // lower triangle!
			ind.segment(0,2) = M.row(ix);
			ind.segment(2,2) = M.row(jx);
			Umat(ix,jx) = (data.col(ind[0]) * data.col(ind[1]) *
				       data.col(ind[2]) * data.col(ind[3]) * rowMult).sum() / totalWeight -
			  Vmat(static_cast<int>(ind[0]),static_cast<int>(ind[1])) * Vmat(static_cast<int>(ind[2]),static_cast<int>(ind[3]));
		}
	}
  Umat.triangularView<Eigen::Upper>() = Umat.transpose().triangularView<Eigen::Upper>();

	Eigen::PermutationMatrix<Eigen::Dynamic> p1(o1.asymCov->cols);
	for (int vx=0; vx < Ecov.cols(); ++vx) {
		p1.indices()[vx] = o1.asymCov->cols - triangleLoc0(Ecov.cols() - vx - 1) - 1;
	}
	for (int v1=0, vx=Ecov.cols(); vx < o1.asymCov->cols; ++vx) {
		if (p1.indices()[v1] == vx - Ecov.cols() + v1) v1 += 1;
		p1.indices()[vx] = vx - Ecov.cols() + v1;
	}
	Umat.derived() = (p1.transpose() * Umat * p1).eval();

	if (strEQ(wlsType, "WLS")) {
    o1.useWeight = omxInitMatrix(numColsStar, numColsStar, state);
    EigenMatrixAdaptor uw(o1.useWeight);
    uw = Umat;
    int singular = InvertSymmetricPosDef(uw, 'L');
    if (singular) {
      mxThrow("%s: weight matrix is rank deficient; cannot invert", name);
      return;
    }
    uw.triangularView<Eigen::Upper>() = uw.transpose().triangularView<Eigen::Upper>();
    uw /= totalWeight;
	} else {
		if (strEQ(wlsType, "ULS")) {
			// OK
		} else { // DWLS
			o1.useWeight = omxInitMatrix(numColsStar, numColsStar, state);
			EigenMatrixAdaptor uw(o1.useWeight);
			uw.setZero();
			for (int ix=0; ix < numColsStar; ++ix) {
				uw(ix,ix) = 1/((data.col(M(ix, 0)) * data.col(M(ix, 0)) *
						  data.col(M(ix, 1)) * data.col(M(ix, 1))).sum() / totalWeight -
							 Vmat(static_cast<int>(M(ix, 0)), static_cast<int>(M(ix, 1))) * Vmat(static_cast<int>(M(ix, 0)), static_cast<int>(M(ix, 1))));
			}
			uw.derived() = (p1.transpose() * uw * p1).eval();
		}
	}

  Umat /= totalWeight;
}

template <typename T1, typename T2, typename T3>
void tabulate(Eigen::MatrixBase<T1> &data, const Eigen::ArrayBase<T3> &weight, Eigen::MatrixBase<T2> &out)
{
	out.setZero();
	for (int rx=0; rx < data.rows(); ++rx) {
		if (data[rx] == NA_INTEGER) continue;
		out[ data[rx] ] += weight[rx];
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
	int naCount;
	double totalWeight;
	const Eigen::Ref<const Eigen::ArrayXd> rowMult;
	std::vector<int> &index;
	ColumnData *response;
	Eigen::MatrixXd pred;
	Eigen::ArrayXd resid;
	Eigen::VectorXd beta;
	Eigen::MatrixXd scores;
	double var;
	Eigen::VectorXd ycol;
  Eigen::MatrixXd vcov;
	OLSRegression(omxData *u_d, double u_totalWeight,
		      const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
		      std::vector<int> &u_index);
	void setResponse(ColumnData &cd, WLSVarData &pv,
									 std::vector< Eigen::Ref<Eigen::VectorXd> > &predCols);
	void calcScores();
};

OLSRegression::OLSRegression(omxData *u_d,
			     double u_totalWeight,
			     const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
			     std::vector<int> &u_index)
	: data(*u_d), totalWeight(u_totalWeight), rowMult(u_rowMult), index(u_index)
{
}

void OLSRegression::setResponse(ColumnData &cd, WLSVarData &pv,
																std::vector< Eigen::Ref<Eigen::VectorXd> > &predCols)
{
	pred.resize(rowMult.size(), 1 + predCols.size());
	pred.col(0).setConstant(1.0);
	for (int cx=0; cx < int(predCols.size()); ++cx)
		pred.col(1+cx) = predCols[cx];

	response = &cd;
	Eigen::Map< Eigen::VectorXd > ycolFull(cd.d(), data.nrows());
	ycol.resize(pred.rows());
	subsetVector(ycolFull, index, ycol);
	auto notMissingF = [&](int rx){ return std::isfinite(ycol[rx]); };
  double indicatorWeight = totalWeight;
	naCount = 0;
	for (int rx=0; rx < int(ycol.size()); ++rx) {
		if (!notMissingF(rx)) {
      naCount += 1;
      indicatorWeight -= rowMult[rx];
    }
	}
	Eigen::VectorXd ycolF(ycol.size() - naCount);
	subsetVector(ycol, notMissingF, ycolF);
	Eigen::ArrayXd rowMultF(rowMult.size() - naCount);
	subsetVector(rowMult, notMissingF, rowMultF);
  Eigen::MatrixXd predCov;
	if (pred.cols() > 1) {
		Eigen::DiagonalMatrix<double, Eigen::Dynamic> weightMat(rowMultF.matrix());
		Eigen::MatrixXd predF(pred.rows() - naCount, pred.cols());
		subsetRows(pred, notMissingF, predF);
		predCov = predF.transpose() * weightMat * predF;
		int singular = InvertSymmetricPosDef(predCov, 'L');
		if (singular) {
			omxRaiseErrorf("%s: cannot invert exogenous predictors (%d)",
				       data.name, singular);
			return;
		}
		beta = predCov.selfadjointView<Eigen::Lower>() * predF.transpose() * weightMat * ycolF;
		resid = ycol - pred * beta;
	} else {
    predCov.resize(1, 1);
    predCov(0,0) = 1 / indicatorWeight;
		beta.resize(1);
		beta[0] = (ycolF.array() * rowMultF).sum() / indicatorWeight;
		resid = ycol.array() - beta[0];
	}
	subsetVectorStore(resid, [&](int rx){ return !std::isfinite(ycol[rx]); }, 0.);
  double residSqSum = (resid.square() * rowMult).sum();
  double sigma2 = residSqSum / (indicatorWeight - beta.size());
  vcov = sigma2 * predCov.selfadjointView<Eigen::Lower>();
	var = residSqSum / indicatorWeight;
}

void OLSRegression::calcScores()
{
	scores.resize(index.size(), 1 + pred.cols());  // mean pred var
	scores.block(0,0,index.size(),pred.cols()) = (pred.array().colwise() * resid) / var;
	scores.col(pred.cols()) = -1./(2*var) + 1./(2*var*var) * resid * resid;
	scores.array().colwise() *= rowMult;
}

struct ProbitRegression : NewtonRaphsonObjective {
	omxData &data;
	obsSummaryStats &o1;
	double totalWeight;
	const Eigen::Ref<const Eigen::ArrayXd> rowMult;
	std::vector<int> &index;
	std::vector< Eigen::Ref<Eigen::VectorXd> > &predCols;
	Eigen::ArrayXXd pred;
	int numThr;
	ColumnData *response;
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
	Eigen::VectorXi ycol;

	// Could merge setResponse with constructor
	ProbitRegression(omxData *u_d,
									 obsSummaryStats &u_o1,
			 const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
									 std::vector<int> &u_index,
									 std::vector< Eigen::Ref<Eigen::VectorXd> > &predCols);
	void setResponse(ColumnData &u_r, int yy);
	virtual double getFit() override { return fit; }
	virtual const char *paramIndexToName(int px) override
	{ return pnames[px].c_str(); }
	void evaluate0();
	void calcScores();
	virtual void evaluateFit() override;
	virtual void evaluateDerivs(int want) override;
	virtual void getParamVec(Eigen::Ref<Eigen::VectorXd> out) override
  { out = param;  }
  virtual void setParamVec(const Eigen::Ref<const Eigen::VectorXd> in) override
  { param = in; }
	virtual double *getGrad() override { return grad.data(); }
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir) override;
};

ProbitRegression::
ProbitRegression(omxData *u_d,
								 obsSummaryStats &u_o1,
								 const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
								 std::vector<int> &u_index,
								 std::vector< Eigen::Ref<Eigen::VectorXd> > &u_predCols) :
	data(*u_d), o1(u_o1), totalWeight(u_o1.totalWeight), rowMult(u_rowMult),
	index(u_index), predCols(u_predCols), numThr(0),
	response(0), verbose(data.verbose), stale(true)
{
	zi.resize(index.size(), 2);
	dzi.resize(index.size(), 2);
	pr.resize(index.size());
	pred.resize(rowMult.rows(), predCols.size());
	for (int ex=0; ex < int(predCols.size()); ++ex) {
		pred.col(ex) = predCols[ex].array();
	}
}

void ProbitRegression::setResponse(ColumnData &u_r, int yy)
{
	response = &u_r;
	numThr = response->getNumThresholds();

	Eigen::Map< Eigen::VectorXi > ycolFull(response->i(), data.nrows());
	ycol.resize(pred.rows());
	subsetVector(ycolFull, index, ycol);
	auto notMissingF = [&](int rx){ return ycol[rx] != NA_INTEGER; };
	int naCount=0;
	for (int rx=0; rx < int(index.size()); ++rx) {
		if (!notMissingF(rx)) naCount += 1;
	}
	Eigen::VectorXi ycolF(ycol.size() - naCount);
	subsetVector(ycol, notMissingF, ycolF);
	Eigen::VectorXd tab(response->getNumOutcomes());
	tabulate(ycolF, rowMult.derived(), tab);
	if ((tab.array()==0).any()) {
		int x,y;
		tab.minCoeff(&x,&y);
    std::string outcome = string_snprintf("%d", 1+x);
    if (response->levelNames.size()) outcome = response->levelNames[x];
		mxThrow("%s: variable '%s' has a zero frequency outcome '%s'.\n"
            "Eliminate this level in your mxFactor() or combine categories in some other way.\n"
            "Do not pass go. Do not collect $200.",
            data.name, response->name, outcome.c_str());
	}
	Eigen::VectorXd prop = (tab.cast<double>() / double(tab.sum())).
		segment(0, numThr);
	if (data.verbose >= 3) mxPrintMat("starting prop", prop);
	cumsum(prop);
	param.resize(prop.size() + pred.cols());
	pnames.clear();
	for (int px=0; px < prop.size(); ++px) {
		param[px] = Rf_qnorm5(prop[px], 0., 1., 1, 0);
		if (verbose >= 1) pnames.push_back(string_snprintf("th%d", 1+px));
	}
	if (verbose >= 1) {
		for (int ex=0; ex < int(o1.exoPred.size()); ++ex) {
			if (!o1.exoFree(yy,ex)) continue;
			auto &c1 = data.rawCol( o1.exoPred[ex] );
			pnames.push_back(c1.name);
		}
	}
	param.segment(numThr, pred.cols()).array() = 0;

	dxa.resize(index.size(), numThr);

	Y1.resize(index.size(), numThr);
	Y2.resize(index.size(), numThr);
	Y1.setZero();
	Y2.setZero();
	for (int rx=0; rx < ycol.size(); ++rx) {
		if (ycol[rx] == NA_INTEGER) continue;
		if (ycol[rx]-1 >= 0)   Y2(rx, ycol[rx]-1) = 1;
		if (ycol[rx] < numThr) Y1(rx, ycol[rx]) = 1;
	}

	lbound.resize(param.size());
	lbound.setConstant(NEG_INF);
	ubound.resize(param.size());
	ubound.setConstant(INF);
	scores.resize(index.size(), param.size());
	hess.resize(param.size(), param.size());
	stale = true;
}

template <typename T1, typename T2, typename T3>
void regressOrdinalThresholds(const Eigen::MatrixBase<T3> &ycol,
															const std::vector< Eigen::Ref<T1> > &pred,
			     ColumnData &oc, WLSVarData &ov,
			      Eigen::ArrayBase<T2> &zi)
{
	zi.derived().resize(ycol.size(), 2);

	int numThr = oc.getNumThresholds();
	Eigen::VectorXd th(2 + numThr);  // pass in as argument TODO
	th.segment(1, numThr) = ov.theta.segment(0, numThr);
	th[0] = NEG_INF;
	th[numThr + 1] = INF;

	if (pred.size() == 0) zi.col(0).setZero();
	for (int px=0; px < int(pred.size()); ++px) {
		if (px == 0) {
			zi.col(0) = -pred[px].array() * ov.theta[numThr + px];
		} else {
			zi.col(0) -= pred[px].array() * ov.theta[numThr + px];
		}
	}
	zi.col(1) = zi.col(0);

	for (int rx=0; rx < ycol.size(); ++rx) {
		if (ycol[rx] == NA_INTEGER) {
			zi(rx,0) = INF;
			zi(rx,1) = NEG_INF;
			continue;
		}
		zi(rx,0) += th[ycol[rx]+1];
		zi(rx,1) += th[ycol[rx]];
	}
}

void ProbitRegression::evaluate0()
{
	Eigen::VectorXd th(1+response->getNumOutcomes());
	th.segment(1, numThr) = param.segment(0, numThr);
	th[0] = -std::numeric_limits<double>::infinity();
	th[response->getNumOutcomes()] = std::numeric_limits<double>::infinity();

	for (int rx=0; rx < ycol.size(); ++rx) {
		if (ycol[rx] == NA_INTEGER) {
			zi(rx,0) = INF;
			zi(rx,1) = NEG_INF;
			pr[rx] = 1.;
			continue;
		}
		double eta = 0;
		if (pred.cols()) eta = pred.row(rx).matrix() * param.segment(numThr, pred.cols());
		zi(rx,0) = std::min(INF, th[ycol[rx]+1] - eta);
		zi(rx,1) = std::max(NEG_INF, th[ycol[rx]] - eta);
		pr[rx] = Rf_pnorm5(zi(rx,0), 0., 1., 1, 0) - Rf_pnorm5(zi(rx,1), 0., 1., 1, 0);
	}
	stale = false;
}

void ProbitRegression::evaluateFit()
{
	evaluate0();
	fit = -(pr.array().log() * rowMult).sum();
}

void ProbitRegression::calcScores()
{
	if (stale) evaluate0();
	dxa.setZero();
	for (int rx=0; rx < ycol.size(); ++rx) {
		dzi(rx,0) = Rf_dnorm4(zi(rx,0), 0., 1., 0);
		dzi(rx,1) = Rf_dnorm4(zi(rx,1), 0., 1., 0);
		if (ycol[rx] == NA_INTEGER) continue;
		if (ycol[rx]-1 >= 0) {
			dxa(rx,ycol[rx]-1) -= dzi(rx,1);
		}
		if (ycol[rx] < numThr) {
			dxa(rx, ycol[rx]) += dzi(rx,0);
		}
	}
	scores.block(0,0,index.size(),numThr) = dxa.colwise() / pr;
	scores.block(0,numThr,index.size(),pred.cols()) =
		pred.array().colwise() * ((dzi.col(1) - dzi.col(0)) / pr);
	scores.colwise() *= rowMult;
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
	if (data.verbose >= 3) mxPrintMat("grad", grad);

	Eigen::ArrayXXd gdzi = (dzi * zi).colwise() * (rowMult / -pr);
	Eigen::ArrayXd pr2 = pr * pr;
	pr2 = 1.0/pr2;

	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weightMat(rowMult.matrix());

	hess.block(0,0,numThr,numThr) =
		dxa.transpose().matrix() * weightMat * (dxa.colwise() * pr2).matrix() -
		((Y1.colwise() * gdzi.col(0)).transpose().matrix() * Y1.matrix() -
		 (Y2.colwise() * gdzi.col(1)).transpose().matrix() * Y2.matrix());

	Eigen::ArrayXXd dxb = pred.array().colwise() * (dzi.col(0) - dzi.col(1));

	hess.block(numThr,numThr,pred.cols(),pred.cols()) =
		dxb.transpose().matrix() * weightMat * (dxb.colwise() * pr2).matrix() -
		((pred.array().colwise() * gdzi.col(0)).transpose().matrix() * pred.matrix() -
		 (pred.array().colwise() * gdzi.col(1)).transpose().matrix() * pred.matrix());

	hess.block(0,numThr,numThr,pred.cols()) =
		-(dxa.transpose().matrix() * weightMat * (dxb.colwise() * pr2).matrix() -
		  ((Y1.colwise() * gdzi.col(0)).transpose().matrix() * pred.matrix() -
		   (Y2.colwise() * gdzi.col(1)).transpose().matrix() * pred.matrix()));

	hess.block(numThr,0,pred.cols(),numThr) =
		hess.block(0,numThr,numThr,pred.cols()).transpose().eval();

	if (data.verbose >= 3) mxPrintMat("hess", hess);

	if (want & FF_COMPUTE_HESSIAN) {
		// TODO for accurate logging
	}
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

struct PolyserialCor : NewtonRaphsonObjective {
	double totalWeight;
	const Eigen::Ref<const Eigen::ArrayXd> rowMult;
	std::vector<int> &index;
	// continuous
	double var;
	Eigen::ArrayXd &resid;
	Eigen::ArrayXd zee;
	// ordinal
	Eigen::ArrayXXd zi;
	Eigen::ArrayXXd dzi;
	// other stuff
	omxData &data;
	int numThr;
	ColumnData &oc;
	WLSVarData &ov;
	double param, grad;
	double fit;
	std::vector< Eigen::Ref<Eigen::VectorXd> > &pred1;
	std::vector< Eigen::Ref<Eigen::VectorXd> > &pred2;
	Eigen::ArrayXXd tau;
	Eigen::ArrayXXd tauj;
	Eigen::ArrayXd pr;
	Eigen::ArrayXXd scores;
	Eigen::VectorXd ycol;

	PolyserialCor(omxData *u_d, WLSVarData &cv, ColumnData &u_oc, WLSVarData &u_ov,
								std::vector< Eigen::Ref<Eigen::VectorXd> > &u_pred1,
								std::vector< Eigen::Ref<Eigen::VectorXd> > &u_pred2,
		      double u_totalWeight,
		      const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
		      std::vector<int> &u_index) :
		totalWeight(u_totalWeight), rowMult(u_rowMult), index(u_index),
		resid(cv.resid), data(*u_d), oc(u_oc), ov(u_ov), pred1(u_pred1),
		pred2(u_pred2)
	{
		lbound.resize(1);
		lbound.setConstant(NEG_INF);
		ubound.resize(1);
		ubound.setConstant(INF);

		var = cv.theta[cv.theta.size()-1];
		zee = resid / sqrt(var);

		pr.resize(index.size());
		dzi.resize(index.size(), 2);

		Eigen::Map< Eigen::VectorXi > ycolFull(oc.i(), data.nrows());
		ycol.resize(rowMult.rows());
		subsetVector(ycolFull, index, ycol);

		regressOrdinalThresholds(ycol, pred2, oc, ov, zi);

		numThr = oc.getNumThresholds();

		int naCount=0;
		auto notMissingF = [&](int rx){ return ycol[rx] != NA_INTEGER; };
		for (int rx=0; rx < int(index.size()); ++rx) {
			if (!notMissingF(rx)) naCount += 1;
		}
		Eigen::VectorXi ycolF(ycol.rows() - naCount);
		subsetVector(ycol, notMissingF, ycolF);
		Eigen::ArrayXd zeeF(ycolF.size());
		subsetVector(zee, notMissingF, zeeF);

		Eigen::ArrayXd rowMultF(ycolF.size());
		subsetVector(rowMult, notMissingF, rowMultF);
		double den = 0;
		for (int tx=0; tx < numThr; ++tx) den += Rf_dnorm4(ov.theta[tx], 0., 1., 0);
		double rho = (zeeF * ycolF.cast<double>().array() * rowMultF).sum() /
			((totalWeight - naCount) * sqrt(var) * den);
		if (!std::isfinite(rho)) mxThrow("PolyserialCor starting value not finite");
		if (fabs(rho) >= 1.0) rho = 0;
		if (data.verbose >= 3) mxLog("starting ps rho = %f", rho);
		param = atanh(rho);
	}
	virtual double getFit() override { return fit; };
	virtual const char *paramIndexToName(int px) override { return "rho"; }
	virtual void evaluateFit() override
	{
		double rho = tanh(clamp(param,-100,100));
		double R = sqrt(1 - rho * rho);
		tau = (zi.colwise() - rho * zee) / R;

		for (int rx=0; rx < ycol.rows(); ++rx) {
			pr[rx] = std::max(Rf_pnorm5(tau(rx,0), 0., 1., 1, 0) -
					  Rf_pnorm5(tau(rx,1), 0., 1., 1, 0),
					  std::numeric_limits<double>::epsilon());
		}
		fit = -(pr.log() * rowMult).sum();
	}
	virtual void getParamVec(Eigen::Ref<Eigen::VectorXd> out) override
  { out[0] = param; }
  virtual void setParamVec(const Eigen::Ref<const Eigen::VectorXd> in) override
  { param = in[0]; }
	virtual double *getGrad() override { return &grad; };
	virtual void evaluateDerivs(int want) override
	{
		if (want & FF_COMPUTE_FIT) evaluateFit();

		for (int rx=0; rx < ycol.rows(); ++rx) {
			dzi(rx,0) = Rf_dnorm4(tau(rx,0), 0., 1., 0);
			dzi(rx,1) = Rf_dnorm4(tau(rx,1), 0., 1., 0);
		}

		double rho = tanh(clamp(param,-100,100));
		double R = sqrt(1 - rho * rho);
		tauj = dzi * ((zi * rho).colwise() - zee);
		double dx_rho = (1./(R*R*R*pr) * (tauj.col(0) - tauj.col(1)) * rowMult).sum();

		double cosh_x = cosh(param);
		grad = -dx_rho/(cosh_x * cosh_x);
	}
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir) override
	{
		// Can fix Hessian at 1.0 because only 1 parameter.
		// Line search takes care of scaling.
		searchDir[0] = grad;
	}
	virtual void adjustSpeed(double &speed) override
  { speed = fabs(0.1 / grad); }
	void calcScores()
	{
		// mu1 var1 th2 beta1 beta2 rho
		scores.resize(index.size(), 2 + numThr + pred1.size() + pred2.size() + 1);
		scores.setZero();

		evaluateDerivs(FF_COMPUTE_FIT);
		double rho = tanh(param);
		double R = sqrt(1 - rho * rho);
		double R3 = R*R*R;
		for (int rx=0; rx < ycol.rows(); ++rx) {
			if (ycol[rx] == NA_INTEGER) continue;
			double irpr = 1.0 / (R * pr[rx]);
			scores(rx,0) = 1.0/sqrt(var) *
				(zee[rx] + rho * irpr * (dzi(rx,0)-dzi(rx,1)));
			scores(rx,1) =
				1.0/(2*var) * ((zee[rx]*zee[rx] - 1.0) +
					       rho*zee[rx] * irpr * (dzi(rx,0)-dzi(rx,1)));
			if (ycol(rx) < numThr)
  			        scores(rx, 2 + static_cast<int>(ycol(rx))) = dzi(rx,0) * irpr;
			if (ycol(rx)-1 >= 0)
				scores(rx, 2 + static_cast<int>(ycol(rx))-1) = -dzi(rx,1) * irpr;
			for (int px=0; px < int(pred1.size()); ++px) {
				scores(rx, 2+numThr+px) = scores(rx,0) * pred1[px][rx];
			}
			for (int px=0; px < int(pred2.size()); ++px) {
				scores(rx, 2+numThr+pred1.size()+px) =
					-pred2[px][rx] * (dzi(rx,0)-dzi(rx,1)) * irpr;
			}
			scores(rx, 2 + numThr + pred1.size() + pred2.size()) =
				1./(R3*pr[rx]) * (tauj(rx,0) - tauj(rx,1));
		}
		scores.colwise() *= rowMult;
	}
};

struct PolychoricCor : NewtonRaphsonObjective {
	double totalWeight;
	const Eigen::Ref<const Eigen::ArrayXd> rowMult;
	std::vector<int> &index;
	omxData &data;
	ColumnData &c1;
	WLSVarData &v1;
	ColumnData &c2;
	WLSVarData &v2;
	std::vector< Eigen::Ref<Eigen::VectorXd> > &pred1;
	std::vector< Eigen::Ref<Eigen::VectorXd> > &pred2;
	int numThr1;
	int numThr2;
	Eigen::ArrayXXd z1;
	Eigen::ArrayXXd z2;
	Eigen::ArrayXd pr;
	Eigen::ArrayXd den;
	double param;
	double fit, grad;
	Eigen::ArrayXXd scores;
	Eigen::ArrayXi y1;
	Eigen::ArrayXi y2;
	Eigen::ArrayXd th1;
	Eigen::ArrayXd th2;
	Eigen::ArrayXXd obsTable;

	PolychoricCor(omxData *u_d, ColumnData &u_c1, WLSVarData &u_v1,
		      ColumnData &u_c2, WLSVarData &u_v2,
								std::vector< Eigen::Ref<Eigen::VectorXd> > &u_pred1,
								std::vector< Eigen::Ref<Eigen::VectorXd> > &u_pred2,
		      double u_totalWeight,
		      const Eigen::Ref<const Eigen::ArrayXd> u_rowMult,
		      std::vector<int> &u_index)
		: totalWeight(u_totalWeight), rowMult(u_rowMult), index(u_index),
		  data(*u_d), c1(u_c1), v1(u_v1), c2(u_c2), v2(u_v2), pred1(u_pred1), pred2(u_pred2)
	{
		lbound.resize(1);
		lbound.setConstant(NEG_INF);
		ubound.resize(1);
		ubound.setConstant(INF);
		numThr1 = c1.getNumThresholds();
		numThr2 = c2.getNumThresholds();
		th1.resize(2 + numThr1);
		th1.segment(1, numThr1) = v1.theta.segment(0, numThr1);
		th1[0] = NEG_INF;
		th1[numThr1 + 1] = INF;
		th2.resize(2 + numThr2);
		th2.segment(1, numThr2) = v2.theta.segment(0, numThr2);
		th2[0] = NEG_INF;
		th2[numThr2 + 1] = INF;

		Eigen::Map< Eigen::ArrayXi > y1Full(c1.i(), data.nrows());
		y1.resize(index.size());
		subsetVector(y1Full, index, y1);
		Eigen::Map< Eigen::ArrayXi > y2Full(c2.i(), data.nrows());
		y2.resize(index.size());
		subsetVector(y2Full, index, y2);

    double pairwiseWeight = totalWeight;
		int naCount = 0;
		for (int rx=0; rx < rowMult.rows(); ++rx) {
			if (y1[rx] != NA_INTEGER && y2[rx] != NA_INTEGER) continue;
			naCount += 1;
      pairwiseWeight -= rowMult[rx];
		}
		auto notMissingF = [&](int rx){ return y1[rx] != NA_INTEGER && y2[rx] != NA_INTEGER; };
		Eigen::ArrayXi y1F(rowMult.rows() - naCount);
		Eigen::ArrayXi y2F(rowMult.rows() - naCount);
		Eigen::ArrayXd rowMultF(rowMult.rows() - naCount);
		subsetVector(y1, notMissingF, y1F);
		subsetVector(y2, notMissingF, y2F);
		subsetVector(rowMult, notMissingF, rowMultF);
		Eigen::ArrayXd y1c = y1F.cast<double>() - (y1F.cast<double>() * rowMultF).sum() / pairwiseWeight;
		Eigen::ArrayXd y2c = y2F.cast<double>() - (y2F.cast<double>() * rowMultF).sum() / pairwiseWeight;
		double rho = (y1c * y2c * rowMultF).sum() /
			(sqrt((y1c*y1c*rowMultF).sum() * (y2c*y2c*rowMultF).sum()));
		if (fabs(rho) >= 1.0) rho = 0;
		if (data.verbose >= 3) mxLog("starting rho = %f", rho);
		param = atanh(rho);

		if (pred1.size() || pred2.size() || !data.getNoExoOptimize()) {
			pr.resize(index.size());
			den.resize(index.size());
			regressOrdinalThresholds(y1.matrix(), pred1, c1, v1, z1);
			regressOrdinalThresholds(y2.matrix(), pred2, c2, v2, z2);
		} else {
			obsTable.resize(c1.getNumOutcomes(), c2.getNumOutcomes());
			obsTable.setZero();
			pr.resize(obsTable.size());
			den.resize(obsTable.size());
			for (int rx=0; rx < y1F.rows(); ++rx) {
				obsTable(y1F[rx], y2F[rx]) += rowMultF[rx];
			}
		}
	}

	virtual void getParamVec(Eigen::Ref<Eigen::VectorXd> out) override
  { out[0] = param; }
  virtual void setParamVec(const Eigen::Ref<const Eigen::VectorXd> in) override
  { param = in[0]; }
	virtual double *getGrad() override { return &grad; };
	virtual const char *paramIndexToName(int px) override { return "rho"; };
	virtual double getFit() override { return fit; };
	virtual void evaluateFit() override
	{
		double rho = tanh(clamp(param,-100,100));

		const double eps = std::numeric_limits<double>::epsilon();

		if (pred1.size() || pred2.size() || !data.getNoExoOptimize()) {
			for (int rx=0; rx < int(index.size()); ++rx) {
				pr[rx] = std::max(pbivnorm(z1(rx,1), z2(rx,1), z1(rx,0), z2(rx,0), rho), eps);
			}
			fit = -(pr.log() * rowMult).sum();
		} else {
			fit = 0;
			for (int cx=0; cx < obsTable.cols(); ++cx) {
				for (int rx=0; rx < obsTable.rows(); ++rx) {
					int px = cx * obsTable.rows() + rx;
					pr[px] = std::max(pbivnorm(th1[rx], th2[cx],
								   th1[rx+1], th2[cx+1], rho), eps);
					fit -= log(pr[px]) * obsTable(rx,cx);
				}
			}
		}
	}
	virtual void evaluateDerivs(int want) override
	{
		if (want & FF_COMPUTE_FIT) evaluateFit();

		double rho = tanh(clamp(param,-100,100));
		double dx = 0;

		if (pred1.size() || pred2.size() || !data.getNoExoOptimize()) {
			for (int rx=0; rx < int(index.size()); ++rx) {
				den[rx] = dbivnorm(z1(rx,1), z2(rx,1), z1(rx,0), z2(rx,0), rho);
				dx += rowMult[rx] * den[rx] / pr[rx];
			}
		} else {
			for (int cx=0; cx < obsTable.cols(); ++cx) {
				for (int rx=0; rx < obsTable.rows(); ++rx) {
					int px = cx * obsTable.rows() + rx;
					den[px] = dbivnorm(th1[rx], th2[cx],
							   th1[rx+1], th2[cx+1], rho);
					dx += obsTable(rx,cx) * den[px] / pr[px];
				}
			}
		}
		double cosh_x = cosh(param);
		grad = -dx / (cosh_x * cosh_x);
	}
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir) override
	{
		// Can fix Hessian at 1.0 because only 1 parameter.
		// Line search takes care of scaling.
		searchDir[0] = grad;
	}
	virtual void adjustSpeed(double &speed) override
  { speed = fabs(0.1 / grad); }
	void calcScores()
	{
		// th1 th2 beta1 beta2 rho
		scores.resize(index.size(), numThr1 + numThr2 + pred1.size() + pred2.size() + 1);
		scores.setZero();

		evaluateDerivs(FF_COMPUTE_FIT);
		double rho = tanh(param);
		double R = sqrt(1 - rho*rho);

		Eigen::ArrayXXd Z1;
		Eigen::ArrayXXd Z2;
		bool slow = pred1.size() || pred2.size() || !data.getNoExoOptimize();
		if (slow) {
			Z1.resize(index.size(), 2);
			Z2.resize(index.size(), 2);

			for (int rx=0; rx < rowMult.rows(); ++rx) {
				if (y1[rx] == NA_INTEGER || y2[rx] == NA_INTEGER) continue;
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
			}
		} else {
			Z1.resize(pr.size(), 2);
			Z2.resize(pr.size(), 2);

			for (int cx=0; cx < obsTable.cols(); ++cx) {
				for (int rx=0; rx < obsTable.rows(); ++rx) {
					int px = cx * obsTable.rows() + rx;
					double ipr = 1.0 / pr[px];
					Z1(px,0) = Rf_dnorm4(th1[rx+1], 0., 1., 0) *
						(Rf_pnorm5((th2[cx+1] - rho*th1[rx+1])/R, 0., 1., 1, 0) -
						 Rf_pnorm5((th2[cx] - rho*th1[rx+1])/R, 0., 1., 1, 0)) * ipr;
					Z1(px,1) = Rf_dnorm4(th1[rx], 0., 1., 0) *
						(Rf_pnorm5((th2[cx+1] - rho*th1[rx])/R, 0., 1., 1, 0) -
						 Rf_pnorm5((th2[cx] - rho*th1[rx])/R, 0., 1., 1, 0)) * ipr;
					Z2(px,0) = Rf_dnorm4(th2[cx+1], 0., 1., 0) *
						(Rf_pnorm5((th1[rx+1] - rho*th2[cx+1])/R, 0., 1., 1, 0) -
						 Rf_pnorm5((th1[rx] - rho*th2[cx+1])/R, 0., 1., 1, 0)) * ipr;
					Z2(px,1) = Rf_dnorm4(th2[cx], 0., 1., 0) *
						(Rf_pnorm5((th1[rx+1] - rho*th2[cx])/R, 0., 1., 1, 0) -
						 Rf_pnorm5((th1[rx] - rho*th2[cx])/R, 0., 1., 1, 0)) * ipr;
				}
			}
		}

		for (int rx=0; rx < rowMult.rows(); ++rx) {
			if (y1[rx] == NA_INTEGER || y2[rx] == NA_INTEGER) continue;
			int px = slow? rx : y2[rx] * obsTable.rows() + y1[rx];
			if (y1(rx) < numThr1) scores(rx, y1(rx)) = Z1(px,0);
			if (y1(rx)-1 >= 0)      scores(rx, y1(rx)-1) = -Z1(px,1);
			if (y2(rx) < numThr2) scores(rx, numThr1 + y2(rx)) = Z2(px,0);
			if (y2(rx)-1 >= 0)      scores(rx, numThr1 + y2(rx)-1) = -Z2(px,1);
			for (int ex=0; ex < int(pred1.size()); ++ex) {
				scores(rx, numThr1+numThr2+ex) =
					(Z1(px,1)-Z1(px,0)) * pred1[ex][rx];
			}
			for (int ex=0; ex < int(pred2.size()); ++ex) {
				scores(rx, numThr1+numThr2+pred1.size()+ex) =
					(Z2(px,1)-Z2(px,0)) * pred2[ex][rx];
			}
			scores(rx,numThr1+numThr2+pred1.size()+pred2.size()) = den[px] / pr[px];
		}
		scores.colwise() *= rowMult;
	}
};

struct PearsonCor {
	double rho;
	Eigen::ArrayXXd scores;

	PearsonCor(WLSVarData &pv1, WLSVarData &pv2,
						 std::vector< Eigen::Ref<Eigen::VectorXd> > &pred1,
						 std::vector< Eigen::Ref<Eigen::VectorXd> > &pred2,
		   const Eigen::Ref<const Eigen::ArrayXd> rowMult,
		   std::vector<int> &index)
	{
		int rows = pv1.resid.size();

		double var_y1 = pv1.theta[pv1.theta.size()-1];
		double sd_y1 = sqrt(var_y1);
		double var_y2 = pv2.theta[pv2.theta.size()-1];
		double sd_y2 = sqrt(var_y2);
		double joint_N = ((pv1.resid != 0.).cast<double>() * (pv2.resid != 0.).cast<double>() * rowMult).sum();
		// rho is the correlation, computed as the standardized covariance
		rho = ((pv1.resid * pv2.resid * rowMult).sum() / joint_N) / (sd_y1 * sd_y2);
		double R = (1 - rho*rho);
		double i2r = 1./(2.*R);

		int numPred = pred1.size() + pred2.size();
		scores.resize(rows, 4 + numPred + 1);

		scores.col(0) = (2*pv1.resid/var_y1 - 2*rho*pv2.resid/(sd_y1*sd_y2)) * i2r;
		scores.col(1) = (2*pv2.resid/var_y2 - 2*rho*pv1.resid/(sd_y1*sd_y2)) * i2r;
		scores.col(2) = -(.5/var_y1 - ((pv1.resid*pv1.resid)/(var_y1*var_y1) -
					       rho*pv1.resid*pv2.resid/(var_y1*sd_y1*sd_y2)) * i2r);
		scores.col(3) = -(.5/var_y2 - ((pv2.resid*pv2.resid)/(var_y2*var_y2) -
					       rho*pv1.resid*pv2.resid/(var_y2*sd_y1*sd_y2)) * i2r);
		for (int ex=0; ex < int(pred1.size()); ++ex) {
			scores.col(4+ex) = pred1[ex].array() * scores.col(0);
		}
		for (int ex=0; ex < int(pred2.size()); ++ex) {
			scores.col(4+pred1.size()+ex) = pred2[ex].array() * scores.col(1);
		}

		Eigen::ArrayXd zee = (pv1.resid*pv1.resid/var_y1
				      -2*rho*pv1.resid*pv2.resid/(sd_y1*sd_y2) +
				      pv2.resid*pv2.resid/var_y2);
		scores.col(4+numPred) =
			rho/R + pv1.resid*pv2.resid/(sd_y1*sd_y2*R) - zee*rho/(R*R);
		scores.colwise() *= rowMult;
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

std::string omxData::regenObsStats(const std::vector<const char *> &dc, const char *wlsType)
{
	if (!oss) return string_snprintf("%s: no observed data summary available", name);
	auto &o1 = *oss;
	// implement checks for exoPred (slope matrix) TODO

  if (o1.totalWeight <= 0) mxThrow("%s: o1.totalWeight <= 0", name);

	if (int(dc.size()) != o1.covMat->cols) {
    return string_snprintf("%s: cov is dimension %d but model is dimension %d",
                           name, o1.covMat->cols, int(dc.size()));
	}

	o1.dc = dc;
	o1.wlsType = wlsType;

	ColMapType dataMap;
	for (int cx=0; cx < int(dc.size()); ++cx) {
		//mxLog("%s", dc[cx]);
		dataMap.emplace(dc[cx], cx);
	}

	ColMapType thrMap;
	if (o1.thresholdMat) {
		if (int(o1.thresholdMat->colnames.size()) != o1.thresholdMat->cols) {
      return string_snprintf("%s: thresholdMat has no colnames", name);
		}
		for (int cx=0; cx < int(o1.thresholdMat->colnames.size()); ++cx) {
			thrMap.emplace(o1.thresholdMat->colnames[cx], cx);
		}
	}

	bool permute = false;
	if (int(o1.covMat->colnames.size()) != o1.covMat->cols) {
    return string_snprintf("%s: cov has no colnames", name);
	}
	for (int cx=0; cx < int(o1.covMat->colnames.size()); ++cx) {
		const char *cn = o1.covMat->colnames[cx];
		auto it = dataMap.find(cn);
		if (it == dataMap.end()) {
			return string_snprintf("%s: observedStats don't include column '%s'", name, cn);
		}
		//mxLog("%d %d %s", cx, it->second, cn);
		if (cx != it->second) permute = true;
		auto rci = rawColMap.find(cn);
		if (rci == rawColMap.end()) continue; //hope for the best
		auto &rc = rawCol(rci->second);
		auto it2 = thrMap.find(cn);
		if (it2 == thrMap.end()) {
			if (rc.type != COLUMNDATA_NUMERIC) {
				return string_snprintf("%s: column '%s' is continuous but found %s",
                               name, cn, ColumnDataTypeToString(rc.type));
			}
		} else {
			if (rc.type != COLUMNDATA_ORDERED_FACTOR) {
				return string_snprintf("%s: column '%s' is ordinal data but found %s",
                               name, cn, ColumnDataTypeToString(rc.type));
			}
			o1.thresholdCols[cx].column = it2->second;
		}
	}

	if (o1.thresholdMat) {
		for (int cx=0; cx < int(o1.thresholdMat->colnames.size()); ++cx) {
			auto it = rawColMap.find(o1.thresholdMat->colnames[cx]);
			if (it == rawColMap.end()) continue;
			auto &rc = rawCol(it->second);
			EigenMatrixAdaptor Eth(o1.thresholdMat);
			for (int tx=0; tx <= o1.thresholdMat->rows; ++tx) {
				if (tx < o1.thresholdMat->rows && std::isfinite(Eth(tx,cx))) continue;
				int nthr = rc.getNumThresholds();
				if (tx != nthr) {
          return string_snprintf("%s: threshold '%s' implies %d levels but data has %d levels",
                                 name, o1.thresholdMat->colnames[cx], 1+tx, 1+nthr);
				}
				break;
			}
		}
	}

	if (strEQ(wlsType, "ULS")) {
		omxFreeMatrix(o1.useWeight);
		o1.useWeight = 0;
	} else {
		EigenMatrixAdaptor Eweight(o1.useWeight);
		Eigen::MatrixXd offDiagW = Eweight.triangularView<Eigen::StrictlyUpper>();
		double offDiag = offDiagW.array().abs().sum();
		if (strEQ(wlsType, "DWLS")) {
			if (offDiag > 0) {
        return string_snprintf("%s: DWLS requested but full weight matrix found", name);
			}
		} else {
			if (offDiag == 0.0) {
        return string_snprintf("%s: WLS requested but diagonal weight matrix found", name);
      }
		}
	}

	//omxPrint(o1.covMat, "cov");
	if (permute) {
		if (o1.slopeMat) {
      return string_snprintf("%s: observedStats could be permuted but "
                             "has slopes (not implemented)", name);
		}
		if (int(o1.covMat->colnames.size()) != o1.covMat->cols ||
		    !o1.asymCov || int(o1.asymCov->colnames.size()) != o1.asymCov->cols) {
			return string_snprintf("%s: observedStats could be permuted but dimnames are "
                             "unavailable (cov=%d, acov=%d)", name,
                             int(o1.covMat->colnames.size()),
                             o1.asymCov? int(o1.asymCov->colnames.size()) : -1);
		}
		if (verbose >= 1) mxLog("%s: observedStats needs permutation", name);
		o1.permute(this);
	}
	//omxPrint(o1.covMat, "cov");

	if (verbose >= 1) {
    mxLog("%s: pre-existing observedStats looks good", name);
    if (verbose >= 2) o1.log();
  }

	return "";
}

void omxData::prepObsStats(omxState *state, const std::vector<const char *> &dc,
			   std::vector<int> &exoPred, const char *type,
			  const char *continuousType, bool fullWeight)
{
	if (state->isClone()) mxThrow("omxData::prepObsStats called in a thread context");

	if (strEQ(u_type, "acov")) {
		// ignore request from fit function (legacy, deprecated)
		auto &o1 = *oss;
		if (!o1.thresholdMat && !o1.meansMat) {
			continuousType = "cumulants";
		} else {
			continuousType = "marginals";
		}
		if (!o1.useWeight) {
			type = "ULS";
		} else {
			EigenMatrixAdaptor Eweight(o1.useWeight);
			Eigen::MatrixXd offDiagW = Eweight.triangularView<Eigen::StrictlyUpper>();
			double offDiag = offDiagW.array().abs().sum();
			if (offDiag > 0) type = "WLS";
			else type = "DWLS";
		}
	}

	u_prepObsStats(state, dc, exoPred, type, continuousType, fullWeight);
	if (oss) oss->setDimnames(this);
}

struct sampleStats {
	// based on lav_samplestats_step[12].R, lavaan 0.6-2

	struct FilterPred {  // to finesse constructor initialization order
		omxData &data;
		std::vector<ColumnData> &rawCols;
		std::vector<int> &exoPred;
		std::vector<Eigen::VectorXd> allPred;

		FilterPred(omxData *u_d, obsSummaryStats &u_o1,
							 int rows, std::vector<int> &index) :
			data(*u_d), rawCols(data.filtered.rawCols),
			exoPred(u_o1.exoPred),
			allPred(exoPred.size())
		{
			for (auto &v : allPred) v.resize(rows);
			for (int cx=0; cx < int(exoPred.size()); ++cx) {
				auto &e1 = rawCols[ exoPred[cx] ];
				Eigen::Map< Eigen::VectorXd > vec(e1.d(), data.nrows());
				auto &v1 = allPred[cx];
				for (int ix=0; ix < rows; ++ix) {
					v1[ix] = vec[ index[ix] ];
				}
			}
		}
	};

	omxData &data;
	const std::vector<const char *> &dc;
	Eigen::Ref<Eigen::ArrayXXi> exoOrder;
	obsSummaryStats &o1;
	const Eigen::ArrayXd &rowMult;
	std::vector<int> &index;
	EigenVectorAdaptor Emean;
	EigenMatrixAdaptor Ecov;
	EigenMatrixAdaptor0 Ethr;
	FilterPred fPred;
	int totalExoFree;
	double eps;
	int numCols;
	int pstar;
	int verbose;
	const int rows;
	ColMapType &rawColMap;
	std::vector<ColumnData> &rawCols;
	std::vector<int> &contMap;
	std::vector<int> &thStart;
	const int totalThr;
	Eigen::MatrixXd &A21;
	Eigen::ArrayXXd &H22;
	Eigen::ArrayXXd &H21;
	Eigen::Map< Eigen::ArrayXi > freq;

	sampleStats(omxData *u_d, const std::vector<const char *> &u_dc,
							Eigen::Ref<Eigen::ArrayXXi> u_exoOrder, obsSummaryStats &u_o1) :
		data(*u_d), dc(u_dc), exoOrder(u_exoOrder),
		o1(u_o1), rowMult(o1.rowMult), index(o1.index),
		Emean(o1.meansMat), Ecov(o1.covMat), Ethr(o1.thresholdMat),
		fPred(u_d, u_o1, rowMult.rows(), index),
		totalExoFree(u_o1.totalExoFree),
		rows(data.nrows()),
		rawColMap(data.rawColMap),
		rawCols(data.filtered.rawCols),
		contMap(o1.contMap),
		thStart(o1.thStart),
		totalThr(o1.totalThr),
		A21(o1.A21),
		H22(o1.H22),
		H21(o1.H21),
		freq(data.getFreqColumn(), rows)
	{
		eps = u_d->getFitTolerance();
		numCols = dc.size();
		pstar = triangleLoc1(numCols-1);
		verbose = data.verbose;
	}

	template <typename T1, typename T2>
	void copyScores(Eigen::ArrayBase<T1> &dest, int destCol,
			const Eigen::ArrayBase<T2> &src, int srcCol, int nc=1)
	{
		for (int cx=0; cx < nc; ++cx) {
			if (data.hasFreq()) {
				for (int dx=0, sx=0, ix=0; ix < freq.size(); ++ix) {
					if (freq[ix] == 0) continue;
					double val = src(sx++, srcCol + cx) / double(freq[ix]);
					for (int fx=0; fx < freq[ix]; ++fx) {
						dest(dx++, destCol + cx) = val;
					}
				}
			} else {
				dest.col(destCol + cx) = src.col(srcCol + cx);
			}
		}
	}

	template <typename T1, typename T2>
	double scoreDotProd(const Eigen::ArrayBase<T1> &a1,
			    const Eigen::ArrayBase<T2> &a2)
	{
		if (data.hasFreq()) {
			double result = 0;
			for (int dx=0, sx=0, ix=0; ix < freq.size(); ++ix) {
				if (freq[ix] == 0) continue;
				result += a1(dx) * a2(sx++);
				dx += freq[ix];
			}
			return result;
		} else {
			return (a1 * a2).sum();
		}
	}

	int getNumCols() { return numCols; };
	bool isDone(int rx, int cx) { return std::isfinite(Ecov(rx,cx)); };
	void reportProgress(int numDone) {
		std::string detail = std::to_string(numDone) + "/" + std::to_string(triangleLoc1(numCols));
		Global->reportProgress1(data.name, detail);
	}

	void onDiag(int yy)
	{
		std::vector< Eigen::Ref<Eigen::VectorXd> > pred;
		std::vector<int> predX;
		for (int ex=0; ex < int(o1.exoPred.size()); ++ex) {
			if (!o1.exoFree(yy,ex)) continue;
			predX.push_back(ex);
			pred.emplace_back(fPred.allPred[ex]);
		}
		ColumnData &cd = rawCols[ rawColMap[dc[yy]] ];
		WLSVarData &pv = o1.perVar[yy];
		if (verbose >= 3) mxLog("consider %s", cd.name);
		if (cd.type == COLUMNDATA_NUMERIC) {
			OLSRegression olsr(&data, o1.totalWeight, rowMult, index);
			olsr.setResponse(cd, pv, pred);
			if (isErrorRaised()) return;
			olsr.calcScores();
			pv.resid = olsr.resid;
			pv.theta.resize(olsr.beta.size() + 1);
			pv.theta.segment(0, olsr.beta.size()) = olsr.beta;
			pv.theta[olsr.beta.size()] = olsr.var;
      pv.vcov = olsr.vcov;
			Ecov(yy,yy) = olsr.var;
			if (olsr.var < data.getMinVariance()) {
				omxRaiseErrorf("%s: '%s' has observed variance less than %g",
					       data.name, dc[yy], data.getMinVariance());
				return;
			}
			Emean[yy] = pv.theta[0];
			copyScores(o1.SC_TH, pv.thrOffset, olsr.scores.array(), 0);
			if (pred.size()) {
				EigenMatrixAdaptor Eslope(o1.slopeMat);
				for (int ex=0, sx=0; ex < int(o1.exoPred.size()); ++ex) {
					if (!o1.exoFree(yy,ex)) continue;
					Eslope(yy,ex) = olsr.beta[++sx];
				}
				for (int px=0; px < int(pred.size()); ++px)
					copyScores(o1.SC_SL, exoOrder(yy,predX[px]),
										 olsr.scores.array(), 1+px);
			}
			copyScores(o1.SC_VAR, pv.contOffset, olsr.scores.array(), 1+pred.size());
		} else {
			ProbitRegression pr(&data, o1, rowMult, index, pred);
			pr.setResponse(cd, yy);
			if (pred.size()) {
				NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
        nro.setGradTolerance(data.getGradientTolerance());
				nro(pr);
			} else {
				pr.calcScores();
        pr.evaluateDerivs(FF_COMPUTE_HESSIAN);
			}
      pv.vcov = pr.hess;
      InvertSymmetricIndef(pv.vcov, 'U');
      pv.vcov.triangularView<Eigen::Lower>() = pv.vcov.transpose().triangularView<Eigen::Lower>();
			pv.theta = pr.param;
			auto &tc = o1.thresholdCols[yy];
			Ethr.block(0,tc.column,tc.numThresholds,1) = pv.theta.segment(0,tc.numThresholds);
			Ecov(yy,yy) = 1.;
			Emean[yy] = 0.;
			copyScores(o1.SC_TH, pv.thrOffset, pr.scores, 0, pr.numThr);
			if (pred.size()) {
				EigenMatrixAdaptor Eslope(o1.slopeMat);
				for (int ex=0, sx=0; ex < int(o1.exoPred.size()); ++ex) {
					if (!o1.exoFree(yy,ex)) continue;
					Eslope(yy,ex) = pr.param[pr.numThr + sx++];
				}
				for (int px=0; px < int(pred.size()); ++px)
					copyScores(o1.SC_SL, exoOrder(yy,predX[px]),
										 pr.scores, pr.numThr+px);
			}
		}
	}

	void offDiag(int jj, int ii)
	{
		// assume jj < ii, upper triangle
		std::vector< Eigen::Ref<Eigen::VectorXd> > predj, predi;
		std::vector<int> predXj, predXi;
		for (int ex=0; ex < int(o1.exoPred.size()); ++ex) {
			if (o1.exoFree(jj,ex)) {
				predXj.push_back(ex);
				predj.emplace_back(fPred.allPred[ex]);
			}
			if (o1.exoFree(ii,ex)) {
				predXi.push_back(ex);
				predi.emplace_back(fPred.allPred[ex]);
			}
		}
		int pstar_idx = ii-(jj+1) + pstar - triangleLoc1(numCols - jj - 1);
		ColumnData &cd1 = rawCols[ rawColMap[dc[jj]] ];
		ColumnData &cd2 = rawCols[ rawColMap[dc[ii]] ];
		WLSVarData &pv1 = o1.perVar[jj];
		WLSVarData &pv2 = o1.perVar[ii];
		if (verbose >= 3) mxLog("consider %s %s [%d]", cd1.name, cd2.name, pstar_idx);
		double rho;
		if (cd1.type == COLUMNDATA_NUMERIC && cd2.type == COLUMNDATA_NUMERIC) {
			PearsonCor pc(pv2, pv1, predi, predj, rowMult, index);
			//mxPrintMat("PearsonScores", pc.scores.block(0,0,4,pc.scores.cols()));
			copyScores(o1.SC_COR, pstar_idx, pc.scores, 4+predi.size()+predj.size());
			A21(pstar_idx,thStart[ii]) = scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(0));
			A21(pstar_idx,thStart[jj]) = scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(1));
			for (int px=0; px < int(predi.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(ii,predXi[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(4+px));
			}
			for (int px=0; px < int(predj.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(jj,predXj[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(4+predi.size()+px));
			}
			A21(pstar_idx, totalThr + totalExoFree + contMap[ii]) =
				scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(2));
			A21(pstar_idx, totalThr + totalExoFree + contMap[jj]) =
				scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(3));
			double sd1 = sqrt(pv1.theta[pv1.theta.size()-1]);
			double sd2 = sqrt(pv2.theta[pv2.theta.size()-1]);
			H21(pstar_idx, totalThr + totalExoFree + contMap[ii]) =
				sd1 * pc.rho / (2. * sd2);
			H21(pstar_idx, totalThr + totalExoFree + contMap[jj]) =
				sd2 * pc.rho / (2. * sd1);
			H22(pstar_idx,pstar_idx) = sd1 * sd2;
			rho = pc.rho * H22(pstar_idx,pstar_idx);
		} else if (cd1.type == COLUMNDATA_NUMERIC) {
			PolyserialCor ps(&data, pv1, cd2, pv2, predj, predi,
											 o1.totalWeight, rowMult, index);
			NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
      nro.setGradTolerance(data.getGradientTolerance());
			nro(ps);
			ps.calcScores();
			//mxPrintMat("PolyserialScores", ps.scores.block(0,0,4,ps.scores.cols()));
			copyScores(o1.SC_COR, pstar_idx,
								 ps.scores, 2 + ps.numThr + predj.size() + predi.size());
			A21(pstar_idx, thStart[jj]) = scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(0));
			for (int tx=0; tx < ps.numThr; ++tx)
				A21(pstar_idx, thStart[ii]+tx) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+tx));
			for (int px=0; px < int(predj.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(jj,predXj[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+ps.numThr+px));
			}
			for (int px=0; px < int(predi.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(ii,predXi[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+ps.numThr+predj.size()+px));
			}
			A21(pstar_idx, totalThr + totalExoFree + contMap[jj]) =
				scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(1));
			double std_rho = tanh(ps.param);
			double sd1 = sqrt(pv1.theta[pv1.theta.size()-1]);
			H21(pstar_idx, totalThr + totalExoFree + contMap[jj]) =
				std_rho / (2. * sd1);
			H22(pstar_idx,pstar_idx) = sd1;
			rho = std_rho * sd1;
		} else if (cd2.type == COLUMNDATA_NUMERIC) {
			PolyserialCor ps(&data, pv2, cd1, pv1, predi, predj,
											 o1.totalWeight, rowMult, index);
			NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
      nro.setGradTolerance(data.getGradientTolerance());
			nro(ps);
			ps.calcScores();
			//mxPrintMat("PolyserialScores", ps.scores.block(0,0,4,ps.scores.cols()));
			copyScores(o1.SC_COR, pstar_idx,
								 ps.scores, 2 + ps.numThr + predj.size() + predi.size());
			A21(pstar_idx, thStart[ii]) = scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(0));
			for (int tx=0; tx < ps.numThr; ++tx)
				A21(pstar_idx, thStart[jj]+tx) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+tx));
			for (int px=0; px < int(predi.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(ii,predXi[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+ps.numThr+px));
			}
			for (int px=0; px < int(predj.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(jj,predXj[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(2+ps.numThr+predi.size()+px));
			}
			A21(pstar_idx, totalThr + totalExoFree + contMap[ii]) =
				scoreDotProd(o1.SC_COR.col(pstar_idx), ps.scores.col(1));
			double sd1 = sqrt(pv2.theta[pv2.theta.size()-1]);
			double std_rho = tanh(ps.param);
			H21(pstar_idx, totalThr + totalExoFree + contMap[ii]) =
				std_rho / (2. * sd1);
			H22(pstar_idx,pstar_idx) = sd1;
			rho = std_rho * sd1;
		} else {
			PolychoricCor pc(&data, cd2, pv2, cd1, pv1, predi, predj,
											 o1.totalWeight, rowMult, index);
			NewtonRaphsonOptimizer nro("nr", 100, eps, verbose);
      nro.setGradTolerance(data.getGradientTolerance());
			nro(pc);
			H22(pstar_idx,pstar_idx) = 1.0;
			rho = tanh(pc.param);
			pc.calcScores();
			//mxPrintMat("PolychoricScores", pc.scores.block(0,0,4,pc.scores.cols()));
			copyScores(o1.SC_COR, pstar_idx,
								 pc.scores, pc.numThr1 + pc.numThr2 + predi.size() + predj.size());
			for (int tx=0; tx < pc.numThr1; ++tx)
				A21(pstar_idx, thStart[ii]+tx) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(tx));
			for (int tx=0; tx < pc.numThr2; ++tx)
				A21(pstar_idx, thStart[jj]+tx) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(pc.numThr1 + tx));
			int numThr = pc.numThr1 + pc.numThr2;
			for (int px=0; px < int(predi.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(ii,predXi[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(numThr+px));
			}
			for (int px=0; px < int(predj.size()); ++px) {
				A21(pstar_idx, totalThr + exoOrder(jj,predXj[px])) =
					scoreDotProd(o1.SC_COR.col(pstar_idx), pc.scores.col(numThr+predi.size()+px));
			}
		}
		if (verbose >= 3) mxLog("cov %s %s [%d] -> %f", cd1.name, cd2.name, pstar_idx, rho);
		Ecov(ii,jj) = rho;
		Ecov(jj,ii) = rho;
	}
};

void omxData::convertToDataFrame()
{
	auto &rc = unfiltered.rawCols;
	rc.reserve(cols);
	numNumeric = cols;

	if (!dataMat->colMajor) omxToggleRowColumnMajor(dataMat);

	for(int j = 0; j < cols; j++) {
		const char *colname = dataMat->colnames[j];
		if (j == freqCol || j == primaryKey) {
			int rows = unfiltered.rows;
			auto *col = new int[rows];
			Eigen::Map< Eigen::VectorXi > Ecol(col, rows);
			Eigen::Map< Eigen::VectorXd > dm(omxMatrixColumn(dataMat, j), rows);
			Ecol.derived() = dm.cast<int>();
			rc.emplace_back(colname, COLUMNDATA_INTEGER, col);
		} else {
			ColumnData cd(colname);
      cd.type = COLUMNDATA_NUMERIC;
      cd.setBorrow(omxMatrixColumn(dataMat, j));
			rc.emplace_back(cd);
		}
	}

	rawColMap.clear();
	for (int cx=0; cx < int(rc.size()); ++cx) {
		rawColMap.emplace(rc[cx].name, cx);
	}
}

void obsSummaryStats::loadExoFree(SEXP Ref)
{
	int numRows, numCols;
	getMatrixDims(Ref, &numRows, &numCols);
	exoFree.resize(numRows, numCols);
	int *ef = LOGICAL(Ref);
	for (int cx=0; cx < numCols; ++cx) {
		for (int rx=0; rx < numRows; ++rx) {
			int v = ef[cx * numRows + rx];
			if (v != 0 && v != 1) mxThrow("exoFree matrix cell [%d,%d] "
																		"is not TRUE/FALSE", 1+rx, 1+cx);
			exoFree(rx,cx) = v;
		}
	}
	// check dimnames? TODO
}

void omxData::u_prepObsStats(omxState *state, const std::vector<const char *> &dc,
			   std::vector<int> &exoPred, const char *wlsType,
			  const char *continuousType, bool wantAsymCov)
{
	if (!dc.size()) return;

  std::string reason = regenObsStats(dc, wlsType);
  if (reason.empty()) {
		if (verbose >= 1) mxLog("%s: reusing pre-existing observedStats (partial=%d)",
					name, int(oss->partial));
		if (oss->partial) estimateObservedStats();
		return;
	}

	if (!isRaw()) {
		mxThrow("%s and raw data are also not available", reason);
	} else {
    if (verbose >= 1) mxLog("%s", reason.c_str());
  }

	int numCols = dc.size();
	int numColsStar = triangleLoc1(numCols);

	oss = std::unique_ptr< obsSummaryStats >(new obsSummaryStats);
	auto &o1 = *oss;
	o1.dc = dc;
	o1.exoPred = exoPred;
	if (R_has_slot(owner, Rf_install("exoFree"))) {
		ProtectedSEXP Ref(R_do_slot(owner, Rf_install("exoFree")));
		if (Rf_length(Ref)) o1.loadExoFree(Ref);
	}
	if (o1.exoFree.size() == 0) {
		o1.exoFree.resize(dc.size(), exoPred.size());
		o1.exoFree.setConstant(1);
	} else {
		if (o1.exoFree.rows() != int(dc.size()) || o1.exoFree.cols() != int(exoPred.size())) {
			mxThrow("%s: exoFree must be %dx%d (observed by exogeneous predictor)",
							name, (int)dc.size(), (int)exoPred.size());
		}
	}
	o1.totalExoFree = o1.exoFree.sum();
	o1.wlsType = wlsType;
	o1.continuousType = continuousType;
	o1.wantAsymCov = wantAsymCov;
	o1.output = true;
	Eigen::ArrayXd rowMultFull;
	std::vector<int> &index = o1.index;
	recalcRowWeights(rowMultFull, index);
	Eigen::ArrayXd &rowMult = o1.rowMult;
	rowMult.resize(index.size()); // rowMult.rows() == index.size()
	subsetVector(rowMultFull, index, rowMult);
	o1.totalWeight = rowMult.sum();

	int maxSizeCov = floor(sqrt(2. * std::numeric_limits<int>::max() / double(index.size())));
	if (numCols > maxSizeCov) {
		mxThrow("%s: for %d rows, WLS cannot handle more than %d columns",
			 name, int(index.size()), maxSizeCov);
	}
	if (rowMult.size() < numColsStar) {
		mxThrow("%s: too few observations (%d) for the number of columns (%d).\n"
			 "For WLS, you need at least n*(n+1)/2 + 1 = %d observations.\n"
			 "Better start rubbing two pennies together.",
			 name, rowMult.size(), numCols, numColsStar+1);
	}
	if (numFactor == 0 && strEQ(continuousType, "cumulants")) {
		if (exoPred.size() != 0) {
			mxThrow("%s: allContinuousMethod cumulants does not work "
				 "with exogenous predictors. Use 'marginals' instead", name);
		}
		wlsAllContinuousCumulants(state);
		return;
	}

	if (verbose >= 1) mxLog("%s: computing marginals stats", name);

	if (exoPred.size()) {
		o1.slopeMat = omxInitMatrix(numCols, exoPred.size(), state);
		EigenMatrixAdaptor Eslope(o1.slopeMat);
		Eslope.setZero();
	}

	o1.covMat = omxInitMatrix(numCols, numCols, state);
	{
		EigenMatrixAdaptor Ecov(o1.covMat);
		Ecov.setConstant(nan("unset"));
	}
	o1.meansMat = omxInitMatrix(1, numCols, state);

	std::vector<int> &contMap = o1.contMap;
	int &numContinuous = o1.numContinuous;
	numContinuous = 0;
	int &totalThr = o1.totalThr;
	totalThr = 0;
	int maxNumThr = 0;
	std::vector<int> &thStart = o1.thStart;
	thStart.resize(numCols);
	o1.numOrdinal = 0;
	o1.perVar.resize(numCols);
	for (int yy=0, thrOffset=0; yy < numCols; ++yy) {
		auto &pv = o1.perVar[yy];
    int dc1 = rawColMap[dc[yy]];
		ColumnData &cd = rawCol(dc1);
		thStart[yy] = totalThr;
		if (cd.type == COLUMNDATA_NUMERIC) {
			omxThresholdColumn tc;
			tc.dataColumn = yy;
			tc.column = -1;
			tc.numThresholds = 0;
			tc.isDiscrete = false;
			o1.thresholdCols.push_back(tc);

			pv.contOffset = numContinuous;
			pv.thrOffset = thrOffset++;
			contMap.push_back(numContinuous);
			totalThr += 1;  // mean
			numContinuous += 1;
		} else {
			int numThr = cd.getNumThresholds();

			omxThresholdColumn tc;
			tc.dataColumn = yy;
			tc.column = o1.numOrdinal++;
			tc.numThresholds = numThr;
			tc.isDiscrete = false;
			o1.thresholdCols.push_back(tc);

			pv.contOffset = -1;
			pv.thrOffset = thrOffset;
			contMap.push_back(-1);
			thrOffset += numThr;
			totalThr += numThr;
			maxNumThr = std::max(maxNumThr, numThr);
		}
	}

	if (o1.numOrdinal) {
		o1.thresholdMat = omxInitMatrix(maxNumThr, o1.numOrdinal, state);
		EigenMatrixAdaptor Ethr(o1.thresholdMat);
		Ethr.setConstant(NA_REAL);
	}

	int scoreRows = index.size();
	if (hasFreq()) {
		Eigen::Map< Eigen::ArrayXi > freq(getFreqColumn(), nrows());
		scoreRows = freq.sum();
	}
	if (verbose >= 1) mxLog("%s: orig %d rows; index.size() = %d; scoreRows = %d; totalWeight = %f",
													name, nrows(), int(index.size()), scoreRows, o1.totalWeight);

	int pstar = triangleLoc1(numCols-1);
	o1.SC_VAR.resize(scoreRows, numContinuous);
	o1.SC_SL.resize(scoreRows, o1.totalExoFree);
	o1.SC_TH.resize(scoreRows, totalThr);
	o1.SC_COR.resize(scoreRows, pstar);

	int A11_size = o1.SC_TH.cols() + o1.SC_SL.cols() + o1.SC_VAR.cols();
	int full_acov_size = o1.SC_TH.cols() + o1.exoFree.size() + o1.SC_VAR.cols() + pstar;
	if (!strEQ(wlsType, "ULS")) {
		o1.useWeight = omxInitMatrix(full_acov_size, full_acov_size, state);
	}
	if (wantAsymCov) {
		o1.asymCov = omxInitMatrix(full_acov_size, full_acov_size, state);
	}

	auto &A21 = o1.A21;
	A21.resize(pstar, A11_size);
	A21.setZero();
	auto &H22 = o1.H22;
	H22.resize(pstar, pstar);
	H22.setZero();
	auto &H21 = o1.H21;
	H21.resize(pstar, A11_size);
	H21.setZero();

	estimateObservedStats();
}

void omxData::estimateObservedStats()
{
	auto &o1 = *oss;
	std::vector<const char *> &dc = o1.dc;
	int numCols = dc.size();
	int totalExoFree = o1.totalExoFree;
	Eigen::ArrayXXi exoOrder;
	const char *wlsType = o1.wlsType;
	std::vector<int> &thStart = o1.thStart;
	auto &A21 = o1.A21;
	auto &H22 = o1.H22;
	auto &H21 = o1.H21;

	exoOrder.resizeLike(o1.exoFree);
	exoOrder.setConstant(NA_INTEGER);
	for (int cx=0, dx=0; cx < o1.exoFree.cols(); ++cx) {
		for (int rx=0; rx < o1.exoFree.rows(); ++rx) {
			if (!o1.exoFree(rx,cx)) continue;
			exoOrder(rx,cx) = dx++;
		}
	}

	{
		EigenMatrixAdaptor Ecov(o1.covMat);
    numEstimatedEntries += triangleLoc1(Ecov.rows()) - Ecov.array().isFinite().count();
		sampleStats ss(this, dc, exoOrder, o1);
		int numThr = parallel? Global->numThreads : 1;
		//int numThr = 1;
		CovEntrywiseParallel(numThr, ss);
		if (isErrorRaised()) return;
	}

	if (1) {
		EigenMatrixAdaptor Ecov(o1.covMat);
		SimpCholesky< Eigen::MatrixXd > sc;
		sc.compute(Ecov);
		if (sc.info() != Eigen::Success || !sc.isPositive()) {
			Rf_warning("%s: marginal covariance matrix "
				   "is non-positive definite", name);
		}
	}

	// Small optimization opportunity:
	// We could avoid above score computations if !fullWeight
	if (!o1.wantAsymCov) return;

	//mxPrintMat("SC_TH", o1.SC_TH.block(0,0,4,o1.SC_TH.cols())); // good
	//mxPrintMat("SC_SL", o1.SC_SL.block(0,0,4,o1.SC_SL.cols())); // good
	//mxPrintMat("SC_VAR", o1.SC_VAR.block(0,0,4,o1.SC_VAR.cols())); // good
	//mxPrintMat("SC_COR", o1.SC_COR.block(0,0,4,o1.SC_COR.cols())); // good
	//mxPrintMat("A21", A21); // good
	//mxPrintMat("H22", H22); // good
	//mxPrintMat("H21", H21); // good

	// based on lav_muthen1984, lavaan 0.6-2
	int A11_size = o1.SC_TH.cols() + o1.SC_SL.cols() + o1.SC_VAR.cols();
	int acov_size = A11_size + o1.SC_COR.cols();
	int scoreRows = o1.SC_COR.rows();
	Eigen::MatrixXd SC(scoreRows, acov_size);
	SC.block(0,0,scoreRows,o1.SC_TH.cols()) = o1.SC_TH;
	SC.block(0,o1.SC_TH.cols(),scoreRows,o1.SC_SL.cols()) = o1.SC_SL;
	SC.block(0,o1.SC_TH.cols()+o1.SC_SL.cols(),scoreRows,o1.SC_VAR.cols()) = o1.SC_VAR;
	SC.block(0,o1.SC_TH.cols()+o1.SC_SL.cols()+o1.SC_VAR.cols(),scoreRows,o1.SC_COR.cols()) = o1.SC_COR;
	// mxPrintMat("SC", SC.block(0,0,10,SC.cols())); // good
	Eigen::MatrixXd INNER = SC.transpose() * SC;
	// mxPrintMat("INNER", INNER); // good

	Eigen::MatrixXd A11(A11_size,A11_size);
	A11.setZero();
	int numContinuous = o1.numContinuous;
	int totalThr = o1.totalThr;
	std::vector<int> &contMap = o1.contMap;
	for (int yy=0; yy < numCols; ++yy) {
		ColumnData &cd = rawCol( rawColMap[dc[yy]] );
		std::vector<bool> mask(A11_size, false);
		int numThr = cd.type == COLUMNDATA_NUMERIC? 1 : cd.getNumThresholds();
		for (int tx=0; tx < numThr; ++tx) mask[thStart[yy] + tx] = true;
		for (int px=0; px < o1.exoFree.cols(); ++px) {
			if (o1.exoFree(yy,px)) mask[totalThr + exoOrder(yy,px)] = true;
		}
		if (cd.type == COLUMNDATA_NUMERIC)
			mask[totalThr + totalExoFree + contMap[yy]] = true;
		copyBlockwise(INNER, A11, [&mask](int xx){ return mask[xx]; });
	}
	//mxPrintMat("A11", A11); // good

	if (InvertSymmetricPosDef(A11, 'L')) {
		MoorePenroseInverseSq(A11);
	}

	int pstar = H22.rows();
	Eigen::MatrixXd A22(pstar, pstar);
	A22.setZero();
	for (int ii=0; ii < pstar; ++ii) {
		double val = o1.SC_COR.col(ii).square().sum();
		if (val != 0) A22(ii,ii) = 1.0/val;
	}

	Eigen::MatrixXd A21i;
	if (pstar) A21i = -(A22 * A21 * A11.selfadjointView<Eigen::Lower>());
	Eigen::MatrixXd Bi(acov_size,acov_size);
	Bi.setZero();
	Bi.block(0,0,A11.rows(),A11.cols()) = A11.selfadjointView<Eigen::Lower>();
	if (pstar) {
		Bi.block(A11.rows(),0,A21i.rows(),A21i.cols()) = A21i;
		Bi.block(A11.rows(),A11.cols(), A22.rows(), A22.cols()) = A22;
	}
	//mxPrintMat("Bi", Bi); // good

	Eigen::MatrixXd fw = Bi * INNER * Bi.transpose();

	if (numContinuous) {
		Eigen::MatrixXd H(acov_size,acov_size);
		H.setZero();
		H.block(0,0,A11_size,A11_size) = Eigen::MatrixXd::Identity(A11_size,A11_size);
		H.block(A11_size,0,H21.rows(),H21.cols()) = H21;
		H.block(A11_size,A11_size,H22.rows(),H22.cols()) = H22;
		//mxPrintMat("H", H); //good

		fw = (H * fw * H.transpose()).eval();
	}

	EigenMatrixAdaptor Eac(o1.asymCov);
	if (totalExoFree == o1.exoFree.size()) {
		Eac.derived() = fw;
	} else {
		Eac.setIdentity();
		std::vector<bool> mask(Eac.cols(), true);
		for (int cx=0,dx=0; cx < o1.exoFree.cols(); ++cx) {
			for (int rx=0; rx < o1.exoFree.rows(); ++rx) {
				mask[totalThr + dx++] = o1.exoFree(rx,cx);
			}
		}
		subsetCovarianceStore(Eac, [&mask](int x){ return mask[x]; }, fw);
	}
	//mxPrintMat("Eac", Eac);

	if (strEQ(wlsType, "WLS")) {
			EigenMatrixAdaptor uw(o1.useWeight);
      uw.derived() = Eac;

      if (InvertSymmetricPosDef(uw, 'L')) {
        if (InvertSymmetricIndef(uw, 'L')) {
          std::ostringstream temp;
          Eigen::DiagonalMatrix<double, Eigen::Dynamic>
            isd((1.0 / INNER.diagonal().array().sqrt()).matrix());
          Eigen::MatrixXd innerCor = isd * INNER * isd;
          for (int cx=0; cx < innerCor.cols()-1; ++cx) {
            for (int rx=cx+1; rx < innerCor.rows(); ++rx) {
              if (fabs(fabs(innerCor(rx,cx)) - 1.0) > 1e-6) continue;
              // Better if we reported the names of the summary stats
              temp << " [" << (1+rx) << "," << (1+cx) << "]";
            }
          }
          std::string diag = temp.str();
          mxThrow("%s: the uw matrix is rank deficient as "
                  "determined by LU factorization; perfectly correlated gradients:%s",
                  name, diag.c_str());
        }
        else if (warnNPDuseWeight) Rf_warning("%s: useWeight matrix is not positive definite", name);
      }

      uw.triangularView<Eigen::Upper>() = uw.transpose().triangularView<Eigen::Upper>();
      uw /= o1.totalWeight;
	} else {
		if (strEQ(wlsType, "ULS")) {
			// OK
		} else {
			EigenMatrixAdaptor uw(o1.useWeight);
			uw.setZero();
			for (int ix=0; ix < uw.cols(); ++ix) {
				uw(ix,ix) = 1. / (o1.totalWeight * Eac(ix,ix)); // DWLS
			}
		}
	}
	//mxPrintMat("Eac", Eac);

	o1.partial = false;
	if (verbose >= 2) o1.log();
}

void omxData::invalidateCache()
{
	prep();
	oss.reset();
}

void omxData::invalidateColumnsCache(const std::vector< int > &columns)
{
	if (naAction == NA_OMIT) {
		for (auto col : columns) {
			filtered.clearColumn(col);
		}
	}
	prep();

	if (!oss) return;

	auto &o1 = *oss;
	// If no meansMat then we used the cumulants method
	if (!o1.meansMat || !o1.covMat) { invalidateCache(); return; }

	bool fail = false;
	EigenMatrixAdaptor Ecov(o1.covMat);
	for (auto col : columns) {
		auto it = o1.colMap.find(rawCol(col).name);
		if (it == o1.colMap.end()) {
			if (verbose >= 1) mxLog("%s: column '%s' is not an observed indicator; "
															"must re-estimate all observed stats",
															name, rawCol(col).name);
			fail = true; break;
		}
		Ecov.row(it->second).setConstant(nan("uninit"));
		Ecov.col(it->second).setConstant(nan("uninit"));
		o1.partial = true;
	}
	if (fail) invalidateCache();
  //else mxThrow("partial!");
}

void omxData::evalAlgebras(FitContext *fc)
{
	if (algebra.size()) setModified();
	for (auto ax : algebra) {
		omxMatrix *a1 = fc->state->algebraList[ax];
		if (verbose >= 2) mxLog("%s::evalAlgebras %s(%d)",
														name, a1->name(), ax);
		if (!a1->colnames.size()) {
			mxThrow("%s: algebra '%s' must have colnames", name, a1->name());
		}
		std::vector<int> colMap;
		int numCols = a1->colnames.size();
		for (int cx=0; cx < numCols; ++cx) {
			auto *cn = a1->colnames[cx];
			auto it = rawColMap.find(cn);
			if (it == rawColMap.end()) {
				mxThrow("%s: cannot find column '%s'", name, cn);
			}
			int dc = it->second;
			if (rawCol(dc).type != COLUMNDATA_NUMERIC) {
				mxThrow("%s: column '%s' must be type of numeric not %s",
						 name, cn, ColumnDataTypeToString(rawCol(dc).type));
			}
			//mxLog("%s -> %d", a1->colnames[cx], dc);
			colMap.push_back(dc);
		}
		for (int rx=0; rx < nrows(); ++rx) {
			loadDefVars(fc->state, rx);
			omxRecompute(a1, fc); // should not depend on free parameters TODO
			if (a1->rows != 1) mxThrow("%s: algebra '%s' must evaluate to a row vector "
						   "instead of %dx%d", name, a1->name(), a1->rows, a1->cols);
			if (a1->cols < numCols) mxThrow("%s: algebra '%s' must have at least "
							"%d columns (found %d)",
							name, a1->name(), numCols, a1->cols);
			EigenVectorAdaptor result(a1);
			for (int cx=0; cx < numCols; ++cx) {
				if (verbose >= 3) mxLog("%s::evalAlgebras [%d,%d] <- %f",
							name, 1+rx, 1+cx, result[cx]);
				rawCol( colMap[cx] ).d()[rx] = result[cx];
			}
		}
	}
}
