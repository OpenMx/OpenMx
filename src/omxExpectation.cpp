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
#include "EnableWarnings.h"

typedef struct omxExpectationTableEntry omxExpectationTableEntry;

struct omxExpectationTableEntry {
	char name[32];
	omxExpectation *(*initFun)();
};

static const omxExpectationTableEntry omxExpectationSymbolTable[] = {
	{"MxExpectationLISREL",			&omxInitLISRELExpectation},
	{"MxExpectationStateSpace",			&omxInitStateSpaceExpectation},
	{"MxExpectationNormal", 		&omxInitNormalExpectation},
	{"MxExpectationRAM",			&omxInitRAMExpectation},
	{"MxExpectationBA81", &omxInitExpectationBA81},
	{"MxExpectationGREML", &omxInitGREMLExpectation},
	{"MxExpectationHiddenMarkov", &InitHiddenMarkovExpectation},
	{"MxExpectationMixture", &InitMixtureExpectation},
};

void omxFreeExpectationArgs(omxExpectation *ox) {
	if(ox==NULL) return;
	delete ox;
}

void omxExpectationRecompute(FitContext *fc, omxExpectation *ox)
{
	omxExpectationCompute(fc, ox, NULL);
}

void omxExpectationCompute(FitContext *fc, omxExpectation *ox, const char *what, const char *how)
{
	if (!ox) return;

	if (ox->data) ox->data->recompute(); // for dynamic data
	if (ox->thresholdsMat) omxRecompute(ox->thresholdsMat, fc);
	ox->compute(fc, what, how);
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
	numOrdinal = 0;

	//omxPrint(thresholdsMat, "loadThr");

	auto dc = base::getDataColumns();
	thresholds.reserve(dc.size());
	std::vector<bool> found(thresholdsMat->cols, false);
	for(int dx = 0; dx < int(dc.size()); dx++) {
		int index = dc[dx];
		omxThresholdColumn col;
		col.dColumn = index;

		const char *colname = omxDataColumnName(data, index);
		int tc = thresholdsMat->lookupColumnByName(colname);

		if (tc < 0 || (data->rawCols.size() && !omxDataColumnIsFactor(data, index))) {	// Continuous variable
			if(debug || OMX_DEBUG) {
				mxLog("%s: column[%d] '%s' is continuous (tc=%d isFactor=%d)",
				      name, index, colname, tc, omxDataColumnIsFactor(data, index));
			}
			thresholds.push_back(col);
		} else {
			found[tc] = true;
			col.column = tc;
			if (data->rawCols.size()) {
				col.numThresholds = omxDataGetNumFactorLevels(data, index) - 1;
			} else {
				// See omxWLSFitFunction. There are lots of complications
				// to permit the order of variables to differ between
				// the observed data and expected model.
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
	if (isRaw || omxDataHasMatrix(data)) {
		ProtectedSEXP Rdc(R_do_slot(rObj, Rf_install("dataColumns")));
		numCols = Rf_length(Rdc);
		ox->saveDataColumnsInfo(Rdc);
		if(OMX_DEBUG) mxPrintMat("Variable mapping", base::getDataColumns());
		if (isRaw) {
			auto dc = base::getDataColumns();
			for (int cx=0; cx < numCols; ++cx) {
				int var = dc[cx];
				data->assertColumnIsData(var);
			}
		}
	}

	if (R_has_slot(rObj, Rf_install("thresholds"))) {
		if(OMX_DEBUG) {
			mxLog("Accessing Threshold matrix.");
		}
		ProtectedSEXP threshMatrix(R_do_slot(rObj, Rf_install("thresholds")));

		if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
			ox->thresholdsMat = omxMatrixLookupFromState1(threshMatrix, ox->currentState);
			ox->loadThresholds();
		} else {
			if (OMX_DEBUG) {
				mxLog("No thresholds matrix; not processing thresholds.");
			}
			ox->numOrdinal = 0;
		}
	}
}

void omxExpectation::generateData(FitContext *, MxRList &out)
{
	Rf_error("%s: generateData not implemented for '%s'", name, expType);
}

omxExpectation* omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os) {

	SEXP ExpectationClass;
	const char *expType;
	{ScopedProtect p1(ExpectationClass, STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
		expType = CHAR(ExpectationClass);
	}

	omxExpectation* expect = omxNewInternalExpectation(expType, os);

	expect->rObj = rObj;
	expect->expNum = expNum;
	
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
		std::string msg = string_snprintf("Expectation '%s' of type '%s' has"
						  " %d definition variables:\n", ox->name, ox->expType,
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

const Eigen::Map<omxExpectation::DataColumnType> omxExpectation::getDataColumns()
{
	return Eigen::Map<DataColumnType>(dataColumnsPtr, numDataColumns);
}

std::vector< omxThresholdColumn > &omxExpectation::getThresholdInfo()
{
	return thresholds;
}

omxExpectation *
omxNewInternalExpectation(const char *expType, omxState* os)
{
	omxExpectation* expect = 0;

	/* Switch based on Expectation type. */ 
	for (size_t ex=0; ex < OMX_STATIC_ARRAY_SIZE(omxExpectationSymbolTable); ex++) {
		const omxExpectationTableEntry *entry = omxExpectationSymbolTable + ex;
		if(strEQ(expType, entry->name)) {
			expect = entry->initFun();
		        expect->expType = entry->name;
			break;
		}
	}

	if (!expect) Rf_error("expectation '%s' not recognized", expType);

	expect->currentState = os;
	expect->canDuplicate = true;
	expect->dynamicDataSource = false;

	return expect;
}

void omxExpectation::print()
{
	mxLog("(Expectation, type %s) ", (expType==NULL?"Untyped":expType));
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
	bool changed = false;
	for (int k=0; k < int(data->defVars.size()); ++k) {
		omxDefinitionVar &dv = data->defVars[k];
		double newDefVar = omxDoubleDataElement(data, row, dv.column);
		if(ISNA(newDefVar)) {
			Rf_error("Error: NA value for a definition variable is Not Yet Implemented.");
		}
		changed |= dv.loadData(currentState, newDefVar);
	}
	if (changed && OMX_DEBUG_ROWS(row)) { mxLog("%s: loading definition vars for row %d", name, row); }
	return changed;
}

int omxExpectation::numSummaryStats()
{
	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		Rf_error("%s::numSummaryStats is not implemented (for object '%s')", expType, name);
	}
	int count = triangleLoc1(cov->rows);
	omxMatrix *mean = getComponent("means");
	if (mean) count += cov->rows;
	
	for (auto &th : getThresholdInfo()) {
		count += th.numThresholds;
	}
	return count;
}

void normalToStdVector(omxMatrix *cov, omxMatrix *mean, omxMatrix *thr,
		       int numOrdinal, std::vector< omxThresholdColumn > &ti,
		       Eigen::Ref<Eigen::VectorXd> out)
{
	EigenMatrixAdaptor Ecov(cov);
	if (numOrdinal == 0) {
		int dx = 0;
		for (int cx=0; cx < cov->cols; ++cx) {
			for (int rx=cx; rx < cov->rows; ++rx) {
				out[dx++] = Ecov(rx,cx);
			}
		}
		if (mean) {
			EigenVectorAdaptor Emean(mean);
			for (int rx=0; rx < cov->cols; ++rx) {
				out[dx++] = Emean(rx);
			}
		}
		return;
	}
	if (!mean) Rf_error("ordinal indicators and no mean vector");

	EigenVectorAdaptor Emean(mean);
	Eigen::VectorXd stdMean = Emean;
	EigenMatrixAdaptor Eth(thr);
	Eigen::VectorXd sdTmp(1.0/Ecov.diagonal().array().sqrt());
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> sd(Emean.size());
	sd.setIdentity();
	
	{
		int tx = triangleLoc1(cov->rows) + cov->rows;
		for (auto &th : ti) {
			for (int t1=0; t1 < th.numThresholds; ++t1) {
				double sd1 = sdTmp[th.dColumn];
				out[tx + t1] = (Eth(t1, th.column) - Emean[th.dColumn]) * sd1;
				sd.diagonal()[th.dColumn] = sd1;
			}
			if (th.numThresholds) stdMean(th.dColumn) = 0.0;
			tx += th.numThresholds;
		}
	}
	
	Eigen::MatrixXd stdCov(sd * Ecov * sd);

	int dx = 0;
	for (int cx=0; cx < cov->cols; ++cx) {
		for (int rx=cx; rx < cov->rows; ++rx) {
			out[dx++] = stdCov(rx,cx);
		}
	}
	out.segment(dx, cov->cols) = stdMean;
}

void omxExpectation::asVector1(FitContext *fc, int row, Eigen::Ref<Eigen::VectorXd> out)
{
	loadDefVars(row);
	omxExpectationCompute(fc, this, 0);

	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		Rf_error("%s::asVector is not implemented (for object '%s')", expType, name);
	}

	normalToStdVector(cov, getComponent("means"), thresholdsMat,
			  numOrdinal, getThresholdInfo(), out);
}
