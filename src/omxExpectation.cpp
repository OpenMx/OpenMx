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
#include "EnableWarnings.h"

typedef struct omxExpectationTableEntry omxExpectationTableEntry;

struct omxExpectationTableEntry {
	char name[32];
	omxExpectation *(*initFun)(omxState *os);
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
	{"MxExpectationPOVRAM",			povRAMExpectationInit}
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

		const char *colname = data->columnName(index);
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
		if (isRaw) {
			auto dc = base::getDataColumns();
			for (int cx=0; cx < numCols; ++cx) {
				int var = dc[cx];
				data->assertColumnIsData(var);
			}
		}
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
	mxThrow("%s: generateData not implemented for '%s'", name, expType);
}

omxExpectation *
omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os)
{
	SEXP ExpectationClass;
	const char *expType;
	{ScopedProtect p1(ExpectationClass, STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
		expType = CHAR(ExpectationClass);
	}

	omxExpectation* expect = 0;

	/* Switch based on Expectation type. */ 
	for (size_t ex=0; ex < OMX_STATIC_ARRAY_SIZE(omxExpectationSymbolTable); ex++) {
		const omxExpectationTableEntry *entry = omxExpectationSymbolTable + ex;
		if(strEQ(expType, entry->name)) {
			expect = entry->initFun(os);
			expect->expType = entry->name;
			break;
		}
	}

	if (!expect) mxThrow("expectation '%s' not recognized", expType);

	expect->canDuplicate = true;
	expect->dynamicDataSource = false;
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
	if (!data) return false;
	return data->loadDefVars(currentState, row);
}

int omxExpectation::numSummaryStats()
{
	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		mxThrow("%s::numSummaryStats is not implemented (for object '%s')", expType, name);
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

void normalToStdVector(omxMatrix *cov, omxMatrix *mean, omxMatrix *slope, omxMatrix *thr,
		       int numOrdinal, std::vector< omxThresholdColumn > &ti,
		       Eigen::Ref<Eigen::VectorXd> out)
{
	// order of elements: (c.f. lav_model_wls, lavaan 0.6-2)
	// 1. thresholds + means (interleaved)
	// 2. slopes (if any, columnwise per exo)
	// 3. variances (continuous indicators only)
	// 4. covariances; not correlations (lower triangle)

	EigenMatrixAdaptor Ecov(cov);
	if (numOrdinal == 0) {
		int dx = 0;
		if (mean) {
			EigenVectorAdaptor Emean(mean);
			for (int rx=0; rx < cov->cols; ++rx) {
				out[dx++] = Emean(rx);
			}
		}
		if (slope) {
			EigenMatrixAdaptor Eslope(slope);
			for (int cx=0; cx < Eslope.cols(); ++cx) {
				for (int rx=0; rx < Eslope.rows(); ++rx) {
					out[dx++] = Eslope(rx,cx);
				}
			}
		}
		for (int cx=0; cx < cov->cols; ++cx) {
			out[dx++] = Ecov(cx,cx);
		}
		for (int cx=0; cx < cov->cols-1; ++cx) {
			for (int rx=cx+1; rx < cov->rows; ++rx) {
				out[dx++] = Ecov(rx,cx);
			}
		}
		return;
	}
	if (!mean) mxThrow("ordinal indicators and no mean vector");

	EigenVectorAdaptor Emean(mean);
	EigenMatrixAdaptor Eth(thr);
	Eigen::VectorXd sdTmp(1.0/Ecov.diagonal().array().sqrt());
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> sd(Emean.size());
	sd.setIdentity();
	
	int dx = 0;
	for (auto &th : ti) {
		for (int t1=0; t1 < th.numThresholds; ++t1) {
			double sd1 = sdTmp[th.dColumn];
			out[dx++] = (Eth(t1, th.column) - Emean[th.dColumn]) * sd1;
			sd.diagonal()[th.dColumn] = sd1;
		}
		if (!th.numThresholds) {
			out[dx++] = Emean[th.dColumn];
		}
	}
	
	if (slope) {
		EigenMatrixAdaptor Eslope(slope);
		for (int cx=0; cx < Eslope.cols(); ++cx) {
			for (int rx=0; rx < Eslope.rows(); ++rx) {
				out[dx++] = Eslope(rx,cx);
			}
		}
	}

	Eigen::MatrixXd stdCov(sd * Ecov * sd);

	for (int cx=0; cx < cov->cols; ++cx) {
		if (ti[cx].numThresholds) continue;
		out[dx++] = stdCov(cx,cx);
	}

	for (int cx=0; cx < cov->cols-1; ++cx) {
		for (int rx=cx+1; rx < cov->rows; ++rx) {
			out[dx++] = stdCov(rx,cx);
		}
	}
}

void omxExpectation::asVector1(FitContext *fc, int row, Eigen::Ref<Eigen::VectorXd> out)
{
	loadDefVars(row);
	omxExpectationCompute(fc, this, 0);

	omxMatrix *cov = getComponent("cov");
	if (!cov) {
		mxThrow("%s::asVector is not implemented (for object '%s')", expType, name);
	}

	normalToStdVector(cov, getComponent("means"), getComponent("slope"), thresholdsMat,
			  numOrdinal, getThresholdInfo(), out);
}
