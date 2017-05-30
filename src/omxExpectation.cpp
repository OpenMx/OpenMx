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

void omxExpectation::loadThresholds(int numCols, int *thresholdColumn, int *thresholdNumber)
{
	numOrdinal = 0;
	thresholds.reserve(numCols);
	for(int index = 0; index < numCols; index++) {
		omxThresholdColumn col;
		col.dColumn = index;
		if(thresholdColumn[index] == NA_INTEGER) {	// Continuous variable
			if(OMX_DEBUG) {
				mxLog("Column %d is continuous.", index);
			}
			thresholds.push_back(col);
		} else {
			col.column = thresholdColumn[index];
			col.numThresholds = thresholdNumber[index];
			thresholds.push_back(col);
			if(OMX_DEBUG) {
				mxLog("%s: column %d is ordinal with %d thresholds in threshold column %d.", 
				      name, index, thresholdNumber[index], thresholdColumn[index]);
			}
			numOrdinal++;
		}
	}
	if(OMX_DEBUG) {
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
			if(OMX_DEBUG) {
				mxLog("Accessing Threshold Mappings.");
			}
        
			ox->thresholdsMat = omxMatrixLookupFromState1(threshMatrix, ox->currentState);

			/* Process the data and threshold mapping structures */
			/* if (threshMatrix == NA_INTEGER), then we could ignore the slot "thresholdColumns"
			 * and fill all the thresholds with {NULL, 0, 0}.
			 * However the current path does not have a lot of overhead. */
			int* thresholdColumn, *thresholdNumber;
			ProtectedSEXP nextMatrix (R_do_slot(rObj, Rf_install("thresholdColumns")));
			thresholdColumn = INTEGER(nextMatrix);
			ProtectedSEXP itemList(R_do_slot(rObj, Rf_install("thresholdLevels")));
			thresholdNumber = INTEGER(itemList);

			ox->loadThresholds(numCols, thresholdColumn, thresholdNumber);
		} else {
			if (OMX_DEBUG) {
				mxLog("No thresholds matrix; not processing thresholds.");
			}
			ox->numOrdinal = 0;
		}
	}
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
