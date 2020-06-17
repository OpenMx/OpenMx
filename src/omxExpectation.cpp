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
	bool debug = false;
	CheckAST(thresholdsMat, 0);

	//omxPrint(thresholdsMat, "loadThr");

	auto dc = base::getDataColumns();
	thresholds.reserve(dc.size());
	std::vector<bool> found(thresholdsMat->cols, false);
	for(int dx = 0; dx < int(dc.size()); dx++) {
		int index = dc[dx];
		omxThresholdColumn col;
		col.dColumn = index;

		const char *colname = data->columnName(index);
		int tc = thresholdsMat->lookupColumnByName(colname);

		if (tc < 0 || (data->isRaw() && !omxDataColumnIsFactor(data, index))) {	// Continuous variable
			if(debug || OMX_DEBUG) {
				mxLog("%s: column[%d] '%s' is continuous (tc=%d isFactor=%d)",
				      name, index, colname, tc, omxDataColumnIsFactor(data, index));
			}
			thresholds.push_back(col);
		} else {
			found[tc] = true;
			col.column = tc;
			if (data->isRaw()) {
				col.numThresholds = omxDataGetNumFactorLevels(data, index) - 1;
			} else {
				// See omxData::permute
			}

			thresholds.push_back(col);
			if(debug || OMX_DEBUG) {
				mxLog("%s: column[%d] '%s' is ordinal with %d thresholds in threshold column %d.", 
				      name, index, colname, col.numThresholds, tc);
			}
			numOrdinal++;
		}
	}

	if (numOrdinal != thresholdsMat->cols) {
		std::string buf;
		for (int cx=0; cx < thresholdsMat->cols; ++cx) {
			if (found[cx]) continue;
			buf += string_snprintf(" %d", 1+cx);
		}
		omxRaiseErrorf("%s: cannot find data for threshold columns:%s\n(Do appropriate threshold column names match data column names?)", name, buf.c_str());
	}

	if(debug || OMX_DEBUG) {
		mxLog("%d threshold columns processed.", numOrdinal);
	}
}

void omxExpectation::loadFromR()
{
	if (!rObj || !data) return;

	auto ox = this;

	int numCols=0;
	bool isRaw = strEQ(omxDataType(data), "raw");
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

	numOrdinal = 0;

	if (R_has_slot(rObj, Rf_install("thresholds"))) {
		if(OMX_DEBUG) {
			mxLog("Accessing Threshold matrix.");
		}
		ProtectedSEXP threshMatrix(R_do_slot(rObj, Rf_install("thresholds")));

		if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
			ox->thresholdsMat = omxMatrixLookupFromState1(threshMatrix, ox->currentState);
			ox->loadThresholds();
		}
	}
	if (R_has_slot(rObj, Rf_install("discrete"))) {
		// TODO
	}
	if (isRaw) {
		auto dc = base::getDataColumns();
		for (int cx=0; cx < numCols; ++cx) {
			int var = dc[cx];
			data->assertColumnIsData(var);
			// (thresholds[cx].numThresholds>0) == data->columnIsFactor(var)  HOLDS
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

	ox->loadFromR();
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
	EigenMatrixAdaptor Eth(thresholdsMat);
	return Eth(r,c);
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
