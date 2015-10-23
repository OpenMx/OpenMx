/*
 *  Copyright 2007-2015 The OpenMx Project
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

typedef struct omxExpectationTableEntry omxExpectationTableEntry;

struct omxExpectationTableEntry {
	char name[32];
	void (*initFun)(omxExpectation*);
};

static const omxExpectationTableEntry omxExpectationSymbolTable[] = {
	{"MxExpectationLISREL",			&omxInitLISRELExpectation},
	{"MxExpectationStateSpace",			&omxInitStateSpaceExpectation},
	{"MxExpectationNormal", 		&omxInitNormalExpectation},
	{"MxExpectationRAM",			&omxInitRAMExpectation},
	{"MxExpectationBA81", &omxInitExpectationBA81},
  {"MxExpectationGREML", &omxInitGREMLExpectation}
};

void omxFreeExpectationArgs(omxExpectation *ox) {
	if(ox==NULL) return;
    
	if (ox->destructFun) ox->destructFun(ox);
	omxFreeMatrix(ox->dataColumns);
	Free(ox);
}

void omxExpectationRecompute(omxExpectation *ox) {
	for(int i = 0; i < int(ox->thresholds.size()); i++) {
		if (!ox->thresholds[i].matrix) continue;
		omxRecompute(ox->thresholds[i].matrix, NULL);
	}

	omxExpectationCompute(ox, NULL);
}

void omxExpectationCompute(omxExpectation *ox, const char *what, const char *how)
{
	if (!ox) return;

	ox->data->recompute(); // for dynamic data
	ox->computeFun(ox, what, how);
}

omxMatrix* omxGetExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component)
{
	if(component == NULL) return NULL;

	/* Hard-wired expectation components */
	if(strEQ("dataColumns", component)) {
		return ox->dataColumns;
	}

	if(ox->componentFun == NULL) return NULL;

	return(ox->componentFun(ox, off, component));
	
}

void omxSetExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component, omxMatrix* om) {
	if(!strcmp(ox->expType, "MxExpectationStateSpace")) {
		ox->mutateFun(ox, off, component, om);
	}
}

omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState) {

	return omxNewIncompleteExpectation(src->rObj, src->expNum, newState);
}

omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os)
{
	omxExpectation* ox = os->expectationList.at(expIndex);
	return ox;
}

static void omxExpectationProcessDataStructures(omxExpectation* ox, SEXP rObj)
{
	int index, numCols, numOrdinal=0;
	SEXP nextMatrix, itemList, threshMatrix; 
	
	if(rObj == NULL) return;

	if(OMX_DEBUG) {
		mxLog("Accessing variable mapping structure.");
	}

	if (R_has_slot(rObj, Rf_install("dataColumns"))) {
		{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("dataColumns")));
		ox->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, ox->currentState, 0, 0);
		}

		if(OMX_DEBUG) {
			omxPrint(ox->dataColumns, "Variable mapping");
		}
	
		numCols = ox->dataColumns->cols;

		if (R_has_slot(rObj, Rf_install("thresholds"))) {
			if(OMX_DEBUG) {
				mxLog("Accessing Threshold matrix.");
			}
			ScopedProtect p1(threshMatrix, R_do_slot(rObj, Rf_install("thresholds")));

			if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
				if(OMX_DEBUG) {
					mxLog("Accessing Threshold Mappings.");
				}
        
				/* Process the data and threshold mapping structures */
				/* if (threshMatrix == NA_INTEGER), then we could ignore the slot "thresholdColumns"
				 * and fill all the thresholds with {NULL, 0, 0}.
				 * However the current path does not have a lot of overhead. */
				int* thresholdColumn, *thresholdNumber;
				{ScopedProtect pc(nextMatrix, R_do_slot(rObj, Rf_install("thresholdColumns")));
				thresholdColumn = INTEGER(nextMatrix);
				}
				{ScopedProtect pi(itemList, R_do_slot(rObj, Rf_install("thresholdLevels")));
				thresholdNumber = INTEGER(itemList);
				}
				ox->thresholds.reserve(numCols);
				for(index = 0; index < numCols; index++) {
					if(thresholdColumn[index] == NA_INTEGER) {	// Continuous variable
						if(OMX_DEBUG) {
							mxLog("Column %d is continuous.", index);
						}
						omxThresholdColumn col;
						ox->thresholds.push_back(col);
					} else {
						omxThresholdColumn col;
						col.matrix = omxMatrixLookupFromState1(threshMatrix, ox->currentState);
						col.column = thresholdColumn[index];
						col.numThresholds = thresholdNumber[index];
						ox->thresholds.push_back(col);
						if(OMX_DEBUG) {
							mxLog("Column %d is ordinal with %d thresholds in threshold column %d.", 
								index, thresholdNumber[index], thresholdColumn[index]);
						}
						numOrdinal++;
					}
				}
				if(OMX_DEBUG) {
					mxLog("%d threshold columns processed.", numOrdinal);
				}
				ox->numOrdinal = numOrdinal;
			} else {
				if (OMX_DEBUG) {
					mxLog("No thresholds matrix; not processing thresholds.");
				}
				ox->numOrdinal = 0;
			}
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
	
	SEXP nextMatrix;
	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("data")));
	expect->data = omxDataLookupFromState(nextMatrix, os);
	}

	return expect;
}

void omxCompleteExpectation(omxExpectation *ox) {
	
	if(ox->isComplete) return;

	if (ox->rObj) {
		omxState *os = ox->currentState;
		SEXP rObj = ox->rObj;
		SEXP slot;
		{ScopedProtect(slot, R_do_slot(rObj, Rf_install("container")));
		if (Rf_length(slot) == 1) {
			int ex = INTEGER(slot)[0];
			ox->container = os->expectationList.at(ex);
		}
		}

		{ScopedProtect(slot, R_do_slot(rObj, Rf_install("submodels")));
		if (Rf_length(slot)) {
			int numSubmodels = Rf_length(slot);
			int *submodel = INTEGER(slot);
			for (int ex=0; ex < numSubmodels; ex++) {
				int sx = submodel[ex];
				ox->submodels.push_back(omxExpectationFromIndex(sx, os));
			}
		}
		}
	}

	omxExpectationProcessDataStructures(ox, ox->rObj);

	int numSubmodels = (int) ox->submodels.size();
	for (int ex=0; ex < numSubmodels; ex++) {
		omxCompleteExpectation(ox->submodels[ex]);
	}

	ox->initFun(ox);

	if(ox->computeFun == NULL) {
		if (isErrorRaised()) {
			Rf_error("Failed to initialize '%s' of type %s: %s", ox->name, ox->expType, Global->getBads());
		} else {
			Rf_error("Failed to initialize '%s' of type %s", ox->name, ox->expType);
		}
	}

	if (OMX_DEBUG) {
		omxData *od = ox->data;
		omxState *state = ox->currentState;
		std::string msg = string_snprintf("Expectation '%s' of type '%s' has"
						  " %d definition variables:\n", ox->name, ox->expType,
						  int(od->defVars.size()));
		for (int dx=0; dx < int(od->defVars.size()); ++dx) {
			omxDefinitionVar &dv = od->defVars[dx];
			msg += string_snprintf("[%d] column '%s' ->", dx, omxDataColumnName(od, dv.column));
			for (int lx=0; lx < dv.numLocations; ++lx) {
				msg += string_snprintf(" %s[%d,%d]", state->matrixToName(~dv.matrices[lx]),
						       dv.rows[lx], dv.cols[lx]);
			}
			msg += "\n  dirty:";
			for (int mx=0; mx < dv.numDeps; ++mx) {
				msg += string_snprintf(" %s", state->matrixToName(dv.deps[mx]));
			}
			msg += "\n";
		}
		mxLogBig(msg);
	}

	ox->isComplete = TRUE;
}

static void defaultSetVarGroup(omxExpectation *ox, FreeVarGroup *fvg)
{
	if (OMX_DEBUG && ox->freeVarGroup && ox->freeVarGroup != fvg) {
		Rf_warning("setFreeVarGroup called with different group (%d vs %d) on %s",
			ox->name, ox->freeVarGroup->id[0], fvg->id[0]);
	}
	ox->freeVarGroup = fvg;
}

void setFreeVarGroup(omxExpectation *ox, FreeVarGroup *fvg)
{
	(*ox->setVarGroup)(ox, fvg);
}

omxExpectation *
omxNewInternalExpectation(const char *expType, omxState* os)
{
	omxExpectation* expect = Calloc(1, omxExpectation);
	expect->setVarGroup = defaultSetVarGroup;

	/* Switch based on Expectation type. */ 
	for (size_t ex=0; ex < OMX_STATIC_ARRAY_SIZE(omxExpectationSymbolTable); ex++) {
		const omxExpectationTableEntry *entry = omxExpectationSymbolTable + ex;
		if(strEQ(expType, entry->name)) {
		        expect->expType = entry->name;
			expect->initFun = entry->initFun;
			break;
		}
	}

	if(!expect->initFun) {
		Free(expect);
		Rf_error("Expectation %s not implemented", expType);
	}

	expect->currentState = os;
	expect->canDuplicate = true;
	expect->dynamicDataSource = false;

	return expect;
}

void omxExpectationPrint(omxExpectation* ox, char* d) {
	if(ox->printFun != NULL) {
		ox->printFun(ox);
	} else {
		mxLog("(Expectation, type %s) ", (ox->expType==NULL?"Untyped":ox->expType));
	}
}
