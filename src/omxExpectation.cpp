/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
*  omxExpectation.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Expectation objects carry distributional expectations
* 		for the model.  Because they have no requirement
*		to produce a single matrix of output, they are
*		not a subclass of mxMatrix, but rather their own
*		strange beast.
*	// TODO:  Create a multi-matrix Algebra type, and make
*	//	MxExpectation a subtype of that.
*
**********************************************************/

#include "omxExpectation.h"
#include "glue.h"
#include "Compute.h"
#include "EnableWarnings.h"

typedef struct omxExpectationTableEntry omxExpectationTableEntry;

struct omxExpectationTableEntry {
	char name[32];
	omxExpectation *(*initFun)(omxState *os, int num);
};

static const omxExpectationTableEntry omxExpectationSymbolTable[] = {
	{"MxExpectationLISREL",			&omxInitLISRELExpectation},
	{"MxExpectationStateSpace",			&omxInitStateSpaceExpectation},
	{"MxExpectationNormal", 		&omxInitNormalExpectation},
	{"MxExpectationRAM",			&omxInitRAMExpectation},
	{"MxExpectationBA81", &omxInitExpectationBA81},
	{"MxExpectationGREML", &omxInitGREMLExpectation},
	{"MxExpectationHiddenMarkov", &InitHiddenMarkovExpectation},
	{"MxExpectationMixture", &InitMixtureExpectation}
};

void omxFreeExpectationArgs(omxExpectation *ox) {
	if(ox==NULL) return;
	delete ox;
}

void omxExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	if (data) data->recompute(); // for dynamic data
	if (thresholdsMat) {
		omxRecompute(thresholdsMat, fc);

		for (auto &th : thresholds) {
			if (!th.numThresholds) continue;
			int column = th.column;
			int count = th.numThresholds;
			int threshCrossCount = 0;
			omxMatrix *om = thresholdsMat;
			if(count > om->rows) {
				mxThrow("Too many thresholds (%d) requested from %dx%d thresholds matrix (in column %d)",
								count, om->rows, om->cols, column);
			}
			for(int j = 1; j < count; j++) {
				double lower = omxMatrixElement(om, j-1, column);
				double upper = omxMatrixElement(om, j, column);
				if (upper - lower < sqrt(std::numeric_limits<double>::epsilon()) * (fabs(lower) + fabs(upper))) {
					threshCrossCount++;
				}
			}
			if(threshCrossCount > 0) {
				fc->recordIterationError("Found %d thresholds too close together in column %d",
																 threshCrossCount, column+1);
			}
		}
	}
	if (discreteMat) {
		omxRecompute(discreteMat, fc);
		EigenMatrixAdaptor dm(discreteMat);
		discreteCache.resize(discreteMat->cols);
		for (int cx=0; cx < discreteMat->cols; ++cx) {
			auto &vec = discreteCache[cx];
			vec.resize(dm(0,cx));
			switch(int(dm(2,cx))) {
			case 1:
				for (int ox = 0; ox < vec.size(); ++ox) vec[ox] = Rf_ppois(ox, dm(3,cx), 1,0);
				break;
			case 2:
				for (int ox = 0; ox < vec.size(); ++ox)
					vec[ox] = Rf_pnbinom(ox, dm(3,cx), dm(4,cx),1,0);
				break;
			case 3:
				for (int ox = 0; ox < vec.size(); ++ox)
					vec[ox] = Rf_pnbinom_mu(ox, dm(3,cx), dm(4,cx),1,0);
				break;
			default:
				mxThrow("%s: unknown discrete distribution %d in '%s' column %d",
								name, dm(2,cx), discreteMat->name(), cx);
			}
			double zif = Rf_plogis(dm(1,cx), 0,1,1,0);
			for(int j = 0; j < vec.size(); j++) {
				vec[j] = Rf_qnorm5(zif + (1-zif) * vec[j], 0, 1, 1, 0);  // p2z
			}
			//mxLog("[%d]", cx);
			//mxPrintMat("DC", vec);
		}
	}
}

omxMatrix* omxGetExpectationComponent(omxExpectation* ox, const char* component)
{
	if(component == NULL) return NULL;

	return(ox->getComponent(component));
}

void omxSetExpectationComponent(omxExpectation* ox, const char* component, omxMatrix* om)
{
	ox->mutate(component, om);
}

omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState) {

	return omxNewIncompleteExpectation(src->rObj, src->expNum, newState);
}

omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os)
{
	omxExpectation* ox = os->expectationList.at(expIndex);
	return ox;
}

void omxExpectation::loadThresholds()
{
	numOrdinal = 0;
	if (!thresholdsMat && !discreteMat) return;

	bool debug = false;
	if (thresholdsMat) CheckAST(thresholdsMat, 0);
	if (discreteMat) CheckAST(discreteMat, 0);

	//omxPrint(thresholdsMat, "loadThr");

	auto dc = base::getDataColumns();
	thresholds.resize(dc.size());
	for(int dx = 0; dx < int(dc.size()); dx++) {
		thresholds[dx].dColumn = dc[dx];
	}

	std::vector<bool> tfound(thresholdsMat? thresholdsMat->cols : 0, false);
	std::vector<bool> dfound(discreteMat? discreteMat->cols : 0, false);

	for(int dx = 0; dx < int(dc.size()); dx++) {
		omxThresholdColumn &col = thresholds[dx];
		int index = col.dColumn;
		const char *colname = data->columnName(index);

		if (thresholdsMat) {
			int tc = thresholdsMat->lookupColumnByName(colname);
			if (tc >= 0) {
				tfound[tc] = true;
				col.column = tc;
				col.isDiscrete = false;
				if (data->isRaw()) {
					col.numThresholds = omxDataGetNumFactorLevels(data, index) - 1;
					data->assertColumnIsData(col.dColumn, OMXDATA_ORDINAL);
				} else {
					// See omxData
				}
				if(debug || OMX_DEBUG) {
					mxLog("%s: column[%d] '%s' is ordinal with %d thresholds in column %d", 
								name, index, colname, col.numThresholds, tc);
				}
				numOrdinal++;
			}
		}
		if (discreteMat) {
			int tc = discreteMat->lookupColumnByName(colname);
			if (tc >= 0) {
				dfound[tc] = true;
				col.column = tc;
				col.isDiscrete = true;
				col.numThresholds = omxMatrixElement(discreteMat, 0, tc);
				if (data->isRaw()) {
					data->assertColumnIsData(col.dColumn, OMXDATA_COUNT);
					ColumnData &cd = data->rawCol(col.dColumn);
					// const auto range =
					// 	std::minmax_element(cd.ptr.intData, cd.ptr.intData + data->nrows());
					// mxLog("%s range %d-%d", colname, *range.first, *range.second);
				} else {
					// See omxData
				}
				if(debug || OMX_DEBUG) {
					mxLog("%s: column[%d] '%s' is discrete with %d thresholds in column %d", 
								name, index, colname, col.numThresholds, tc);
				}
				numOrdinal++;
			}
		}
	}

	if (thresholdsMat) {
		std::string buf;
		for (int cx=0; cx < thresholdsMat->cols; ++cx) {
			if (tfound[cx]) continue;
			buf += string_snprintf(" %s(%d)", thresholdsMat->colnames[cx], 1+cx);
		}
		if (buf.size()) {
			mxThrow("%s: cannot find data for threshold columns:%s\n(Do appropriate threshold column names match data column names?)", name, buf.c_str());
		}
	}
	if (discreteMat) {
		std::string buf;
		for (int cx=0; cx < discreteMat->cols; ++cx) {
			if (dfound[cx]) continue;
			buf += string_snprintf(" %s(%d)", discreteMat->colnames[cx], 1+cx);
		}
		if (buf.size()) {
			mxThrow("%s: cannot find data for discrete columns:%s\n(Do appropriate discrete column names match data column names?)", name, buf.c_str());
		}
	}

	for (auto &th : thresholds) {
		if (th.column >= 0) continue;
		data->assertColumnIsData(th.dColumn, OMXDATA_REAL);
		if(debug || OMX_DEBUG) {
			mxLog("%s: column[%d] '%s' is continuous",
						name, th.dColumn, data->columnName(th.dColumn));
		}
	}
}

void omxExpectation::loadDataColFromR()
{
	if (!rObj || !data) return;

	auto ox = this;

	int numCols=0;
	{
		ProtectedSEXP Rdc(R_do_slot(rObj, Rf_install("dataColumns")));
		numCols = Rf_length(Rdc);
		ox->saveDataColumnsInfo(Rdc);
		if(OMX_DEBUG) mxPrintMat("Variable mapping", base::getDataColumns());
		if (R_has_slot(rObj, Rf_install("dataColumnNames"))) {
			ProtectedSEXP Rdcn(R_do_slot(rObj, Rf_install("dataColumnNames")));
			loadCharVecFromR(name, Rdcn, dataColumnNames);
		}
		if (numCols && !dataColumnNames.size()) {
			// eventually deprecate slot 'dataColumns'
			if (usesDataColumnNames()) {
				Rf_warning("Slot MxData@dataColumnNames is not set up; OpenMx bug? Improvising...");
			}
			auto dc = base::getDataColumns();
			for (int cx=0; cx < int(dc.size()); ++cx) {
				dataColumnNames.push_back(data->columnName(dc[cx]));
			}
		}
	}
}

void omxExpectation::loadThresholdFromR()
{
	if (R_has_slot(rObj, Rf_install("thresholds"))) {
		ProtectedSEXP threshMatrix(R_do_slot(rObj, Rf_install("thresholds")));
		if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
			thresholdsMat = omxMatrixLookupFromState1(threshMatrix, currentState);
		}
	}
	if (R_has_slot(rObj, Rf_install("discrete"))) {
		ProtectedSEXP mat(R_do_slot(rObj, Rf_install("discrete")));
		if(INTEGER(mat)[0] != NA_INTEGER) {
			discreteMat = omxMatrixLookupFromState1(mat, currentState);
		}
	}
	loadThresholds();
	bool isRaw = strEQ(omxDataType(data), "raw");
	if (isRaw) {
		if (thresholds.size() == 0) {
			auto dc = base::getDataColumns();
			for (int cx=0; cx < int(dc.size()); ++cx) {
				int var = dc[cx];
				data->assertColumnIsData(var, OMXDATA_REAL);
			}
		}
	}
}

void omxExpectation::generateData(FitContext *, MxRList &out)
{
	mxThrow("%s: generateData not implemented", name);
}

omxExpectation *
omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os)
{
	const char *name;
	{ProtectedSEXP ExpectationClass(STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
		name = CHAR(ExpectationClass);
	}

	omxExpectation* expect = 0;

	/* Switch based on Expectation type. */ 
	for (size_t ex=0; ex < OMX_STATIC_ARRAY_SIZE(omxExpectationSymbolTable); ex++) {
		const omxExpectationTableEntry *entry = omxExpectationSymbolTable + ex;
		if(strEQ(name, entry->name)) {
			expect = entry->initFun(os, expNum);
			expect->name = entry->name;
			break;
		}
	}

	if (!expect) mxThrow("expectation '%s' not recognized", name);

	expect->canDuplicate = true;
	expect->dynamicDataSource = false;
	expect->rObj = rObj;
	
	ProtectedSEXP Rdata(R_do_slot(rObj, Rf_install("data")));
	if (TYPEOF(Rdata) == INTSXP) {
		expect->data = omxDataLookupFromState(Rdata, os);
	}

	return expect;
}

void omxCompleteExpectation(omxExpectation *ox) {
	
	if(ox->isComplete) return;
	ox->isComplete = TRUE;

	ox->init();

	if (OMX_DEBUG) {
		omxData *od = ox->data;
		omxState *state = ox->currentState;
		std::string msg = string_snprintf("Expectation '%s' has"
						  " %d definition variables:\n", ox->name,
						  int(od->defVars.size()));
		for (int dx=0; dx < int(od->defVars.size()); ++dx) {
			omxDefinitionVar &dv = od->defVars[dx];
			msg += string_snprintf("[%d] column '%s' ->", dx, omxDataColumnName(od, dv.column));
			msg += string_snprintf(" %s[%d,%d]", state->matrixToName(~dv.matrix),
					       dv.row, dv.col);
			msg += "\n  dirty:";
			for (int mx=0; mx < dv.numDeps; ++mx) {
				msg += string_snprintf(" %s", state->matrixToName(dv.deps[mx]));
			}
			msg += "\n";
		}
		mxLogBig(msg);
	}
}

const std::vector<const char *> &omxExpectation::getDataColumnNames() const
{ return dataColumnNames; }

const Eigen::Map<omxExpectation::DataColumnIndexVector> omxExpectation::getDataColumns()
{
	return Eigen::Map<DataColumnIndexVector>(dataColumnsPtr, numDataColumns);
}

std::vector< omxThresholdColumn > &omxExpectation::getThresholdInfo()
{
	return thresholds;
}

double omxExpectation::getThreshold(int r, int c)
{
	auto &allTh = getThresholdInfo();
	auto &th = allTh[c];
	if (th.isDiscrete) {
		auto &vec = discreteCache[c];
		return vec[r];
	} else {
		EigenMatrixAdaptor Eth(thresholdsMat);
		return Eth(r,c);
	}
}

void omxExpectation::print()
{
	mxLog("(Expectation, type %s) ", (name==NULL?"Untyped":name));
}

void omxExpectationPrint(omxExpectation* ox, char* d)
{
	ox->print();
}

void complainAboutMissingMeans(omxExpectation *off)
{
	omxRaiseErrorf("%s: raw data observed but no expected means "
		       "vector was provided. Add something like mxPath(from = 'one',"
		       " to = manifests) to your model.", off->name);
}

bool omxExpectation::loadDefVars(int row)
{
	if (!data) return false;
	return data->loadDefVars(currentState, row);
}

void omxExpectation::loadFakeDefVars()
{
	if (!data) return;
	data->loadFakeData(currentState, 1.0);
}

int omxExpectation::numSummaryStats()
{
	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		mxThrow("%s::numSummaryStats is not implemented", name);
	}

	omxMatrix *mean = getComponent("means");

	int count = 0;

	omxMatrix *slope = getComponent("slope");
	if (slope) count += slope->rows * slope->cols;

	auto &ti = getThresholdInfo();
	if (ti.size() == 0) {
		// all continuous
		count += triangleLoc1(cov->rows);
		if (mean) count += cov->rows;
		return count;
	}

	count += triangleLoc1(cov->rows - 1);  // covariances
	for (auto &th : ti) {
		// mean + variance
		count += th.numThresholds? th.numThresholds : 2;
	}

	return count;
}

void omxExpectation::asVector1(FitContext *fc, int row, Eigen::Ref<Eigen::VectorXd> out)
{
	loadDefVars(row);
	omxExpectationCompute(fc, this, 0);

	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		mxThrow("%s::asVector is not implemented", name);
	}

	normalToStdVector(cov, getComponent("means"), getComponent("slope"),
										[this](int r, int c)->double{ return this->getThreshold(r,c); },
										getThresholdInfo(), out);
}

bool omxExpectation::isClone() const
{ return currentState->isClone(); }
