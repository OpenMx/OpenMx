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

static void omxExpectationProcessDefinitionVars(omxExpectation* ox, SEXP rObj)
{
	if(!R_has_slot(rObj, Rf_install("definitionVars"))) {
		ox->numDefs = 0;
		ox->defVars = NULL;
		return;
	}

	if(OMX_DEBUG) {
		mxLog("Accessing definition variables structure.");
	}

	SEXP nextMatrix;
	ScopedProtect dv(nextMatrix, R_do_slot(rObj, Rf_install("definitionVars")));
	int numDefs = Rf_length(nextMatrix);
	ox->numDefs = numDefs;
	if(OMX_DEBUG) {
		mxLog("Number of definition variables is %d.", numDefs);
	}
	ox->defVars = (omxDefinitionVar *) R_alloc(numDefs, sizeof(omxDefinitionVar));
	for(int nextDef = 0; nextDef < numDefs; nextDef++) {
		SEXP dataSource, columnSource, depsSource; 
		int nextDataSource, numDeps;

		SEXP itemList;
		ScopedProtect p1(itemList, VECTOR_ELT(nextMatrix, nextDef));
		ScopedProtect p2(dataSource, VECTOR_ELT(itemList, 0));
		nextDataSource = INTEGER(dataSource)[0];
		if(OMX_DEBUG) {
			mxLog("Data source number is %d.", nextDataSource);
		}
		ox->defVars[nextDef].data = nextDataSource;
		ox->defVars[nextDef].source = ox->currentState->dataList[nextDataSource];
		ScopedProtect p3(columnSource, VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG) {
			mxLog("Data column number is %d.", INTEGER(columnSource)[0]);
		}
		ox->defVars[nextDef].column = INTEGER(columnSource)[0];
		ScopedProtect p4(depsSource, VECTOR_ELT(itemList, 2));
		numDeps = LENGTH(depsSource);
		ox->defVars[nextDef].numDeps = numDeps;
		ox->defVars[nextDef].deps = (int*) R_alloc(numDeps, sizeof(int));
		for(int i = 0; i < numDeps; i++) {
			ox->defVars[nextDef].deps[i] = INTEGER(depsSource)[i];
		}

		ox->defVars[nextDef].numLocations = Rf_length(itemList) - 3;
		ox->defVars[nextDef].matrices = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));
		ox->defVars[nextDef].rows = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));
		ox->defVars[nextDef].cols = (int *) R_alloc(Rf_length(itemList) - 3, sizeof(int));
		for(int index = 3; index < Rf_length(itemList); index++) {
			SEXP nextItem;
			ScopedProtect pi(nextItem, VECTOR_ELT(itemList, index));
			ox->defVars[nextDef].matrices[index-3] = INTEGER(nextItem)[0];
			ox->defVars[nextDef].rows[index-3] = INTEGER(nextItem)[1];
			ox->defVars[nextDef].cols[index-3] = INTEGER(nextItem)[2];
		}
	}
}

void omxExpectation::loadFakeData(double fake)
{
	// This is really just a quick hack that I wrote before I
	// understood the problem.  What we should really do is call
	// handleDefinitionVarList to load the first row of data.

	for (int dx=0; dx < numDefs; ++dx) {
		defVars[dx].loadFakeData(currentState, fake);
	}
}

void omxDefinitionVar::loadFakeData(omxState *state, double fake)
{
	for(int l = 0; l < numLocations; l++) {
		int matrixNumber = matrices[l];
		int matrow = rows[l];
		int matcol = cols[l];
		omxMatrix *matrix = state->matrixList[matrixNumber];
		if(OMX_DEBUG) {
			mxLog("Populating fake data (value %3.2f) into matrix '%s' at [%d,%d]",
			      fake, matrix->name, 1+matrow, 1+matcol);
		}
		omxSetMatrixElement(matrix, matrow, matcol, fake);
	}
}

omxExpectation* omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os) {

	SEXP ExpectationClass;
	const char *expType;
	{ScopedProtect p1(ExpectationClass, STRING_ELT(Rf_getAttrib(rObj, Rf_install("class")), 0));
		expType = CHAR(ExpectationClass);
	}

	omxExpectation* expect = omxNewInternalExpectation(expType, os);

	expect->rObj = rObj;
	expect->expNum = expNum;
	
	SEXP nextMatrix;
	{ScopedProtect p1(nextMatrix, R_do_slot(rObj, Rf_install("data")));
	expect->data = omxDataLookupFromState(nextMatrix, os);
	}

	if (rObj) {
		omxExpectationProcessDefinitionVars(expect, rObj);
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
		omxState *state = ox->currentState;
		std::string msg = string_snprintf("Expectation '%s' of type '%s' has"
						  " %d definition variables:\n", ox->name, ox->expType, ox->numDefs);
		for (int dx=0; dx < ox->numDefs; ++dx) {
			omxDefinitionVar *dv = &ox->defVars[dx];
			msg += string_snprintf("[%d] column '%s' ->", dx,
					       omxDataColumnName(dv->source, dv->column));
			for (int lx=0; lx < dv->numLocations; ++lx) {
				msg += string_snprintf(" %s[%d,%d]", state->matrixToName(~dv->matrices[lx]),
						       dv->rows[lx], dv->cols[lx]);
			}
			msg += "\n  dirty:";
			for (int mx=0; mx < dv->numDeps; ++mx) {
				msg += string_snprintf(" %s", state->matrixToName(dv->deps[mx]));
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
