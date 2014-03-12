/*
 *  Copyright 2007-2014 The OpenMx Project
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
	{"MxExpectationBA81", &omxInitExpectationBA81}
};

void omxFreeExpectationArgs(omxExpectation *ox) {
	if(ox==NULL) return;
    
	if (ox->destructFun) ox->destructFun(ox);
	Free(ox->submodels);
	Free(ox);
}

void omxExpectationRecompute(omxExpectation *ox) {
	if(ox->thresholds != NULL) {
		for(int i = 0; i < ox->numOrdinal; i++) {
			if (!ox->thresholds[i].matrix) continue;
			omxRecompute(ox->thresholds[i].matrix);
		}
	}

	omxExpectationCompute(ox, NULL);
}

void omxExpectationCompute(omxExpectation *ox, const char *what, const char *how)
{
	if (!ox) return;

	ox->computeFun(ox, what, how);
}

omxMatrix* omxGetExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component) {

	if(component == NULL) return NULL;

	/* Hard-wired expectation components */
	if(!strncmp("dataColumns", component, 11)) {
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

omxExpectation* omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os) {

	SEXP ExpectationClass;
	Rf_protect(ExpectationClass = STRING_ELT(Rf_getAttrib(rObj, Rf_install("class")), 0));
	const char* expType = CHAR(ExpectationClass);

	omxExpectation* expect = omxNewInternalExpectation(expType, os);

	expect->rObj = rObj;
	expect->expNum = expNum;
	
	return expect;
}

omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os)
{
	omxExpectation* ox = os->expectationList.at(expIndex);
	return ox;
}

void omxExpectationProcessDataStructures(omxExpectation* ox, SEXP rObj){

	int index, numDefs, nextDef, numCols, numOrdinal=0;
	SEXP nextMatrix, itemList, nextItem, threshMatrix; 
	
	if(rObj == NULL) return;

	if(OMX_DEBUG) { mxLog("Retrieving data."); }
	Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("data")));
	ox->data = omxDataLookupFromState(nextMatrix, ox->currentState);

	if(OMX_DEBUG) {
		mxLog("Accessing variable mapping structure.");
	}

	if (R_has_slot(rObj, Rf_install("dataColumns"))) {
		Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("dataColumns")));
		ox->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, ox->currentState, 0, 0);
		if(OMX_DEBUG) {
			omxPrint(ox->dataColumns, "Variable mapping");
		}
	
		numCols = ox->dataColumns->cols;

		if (R_has_slot(rObj, Rf_install("thresholds"))) {
			if(OMX_DEBUG) {
				mxLog("Accessing Threshold matrix.");
			}
			Rf_protect(threshMatrix = R_do_slot(rObj, Rf_install("thresholds")));

			if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
				if(OMX_DEBUG) {
					mxLog("Accessing Threshold Mappings.");
				}
        
				/* Process the data and threshold mapping structures */
				/* if (threshMatrix == NA_INTEGER), then we could ignore the slot "thresholdColumns"
				 * and fill all the thresholds with {NULL, 0, 0}.
				 * However the current path does not have a lot of overhead. */
				Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("thresholdColumns")));
				Rf_protect(itemList = R_do_slot(rObj, Rf_install("thresholdLevels")));
				int* thresholdColumn, *thresholdNumber;
				thresholdColumn = INTEGER(nextMatrix);
				thresholdNumber = INTEGER(itemList);
				ox->thresholds = (omxThresholdColumn *) R_alloc(numCols, sizeof(omxThresholdColumn));
				for(index = 0; index < numCols; index++) {
					if(thresholdColumn[index] == NA_INTEGER) {	// Continuous variable
						if(OMX_DEBUG) {
							mxLog("Column %d is continuous.", index);
						}
						ox->thresholds[index].matrix = NULL;
						ox->thresholds[index].column = 0;
						ox->thresholds[index].numThresholds = 0;
					} else {
						ox->thresholds[index].matrix = omxMatrixLookupFromState1(threshMatrix, 
												       ox->currentState);
						ox->thresholds[index].column = thresholdColumn[index];
						ox->thresholds[index].numThresholds = thresholdNumber[index];
						if(OMX_DEBUG) {
							mxLog("Column %d is ordinal with %d thresholds in threshold column %d.", 
								index, thresholdColumn[index], thresholdNumber[index]);
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
				ox->thresholds = NULL;
				ox->numOrdinal = 0;
			}
		}
	}

	if(!R_has_slot(rObj, Rf_install("definitionVars"))) {
		ox->numDefs = 0;
		ox->defVars = NULL;
	} else {	
		if(OMX_DEBUG) {
			mxLog("Accessing definition variables structure.");
		}
		Rf_protect(nextMatrix = R_do_slot(rObj, Rf_install("definitionVars")));
		numDefs = Rf_length(nextMatrix);
		ox->numDefs = numDefs;
		if(OMX_DEBUG) {
			mxLog("Number of definition variables is %d.", numDefs);
		}
		ox->defVars = (omxDefinitionVar *) R_alloc(numDefs, sizeof(omxDefinitionVar));
		for(nextDef = 0; nextDef < numDefs; nextDef++) {
			SEXP dataSource, columnSource, depsSource; 
			int nextDataSource, numDeps;

			Rf_protect(itemList = VECTOR_ELT(nextMatrix, nextDef));
			Rf_protect(dataSource = VECTOR_ELT(itemList, 0));
			nextDataSource = INTEGER(dataSource)[0];
			if(OMX_DEBUG) {
				mxLog("Data source number is %d.", nextDataSource);
			}
			ox->defVars[nextDef].data = nextDataSource;
			ox->defVars[nextDef].source = ox->currentState->dataList[nextDataSource];
			Rf_protect(columnSource = VECTOR_ELT(itemList, 1));
			if(OMX_DEBUG) {
				mxLog("Data column number is %d.", INTEGER(columnSource)[0]);
			}
			ox->defVars[nextDef].column = INTEGER(columnSource)[0];
			Rf_protect(depsSource = VECTOR_ELT(itemList, 2));
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
			for(index = 3; index < Rf_length(itemList); index++) {
				Rf_protect(nextItem = VECTOR_ELT(itemList, index));
				ox->defVars[nextDef].matrices[index-3] = INTEGER(nextItem)[0];
				ox->defVars[nextDef].rows[index-3] = INTEGER(nextItem)[1];
				ox->defVars[nextDef].cols[index-3] = INTEGER(nextItem)[2];
			}
		}
	}
	
}

void omxCompleteExpectation(omxExpectation *ox) {
	
	if(ox->isComplete) return;

	omxState* os = ox->currentState;

	if (ox->rObj) {
		SEXP slot;
		Rf_protect(slot = R_do_slot(ox->rObj, Rf_install("container")));
		if (Rf_length(slot) == 1) {
			int ex = INTEGER(slot)[0];
			ox->container = os->expectationList.at(ex);
		}

		Rf_protect(slot = R_do_slot(ox->rObj, Rf_install("submodels")));
		if (Rf_length(slot)) {
			ox->numSubmodels = Rf_length(slot);
			ox->submodels = Realloc(NULL, Rf_length(slot), omxExpectation*);
			int *submodel = INTEGER(slot);
			for (int ex=0; ex < ox->numSubmodels; ex++) {
				int sx = submodel[ex];
				ox->submodels[ex] = omxExpectationFromIndex(sx, os);
				omxCompleteExpectation(ox->submodels[ex]);
			}
		}

		omxExpectationProcessDataStructures(ox, ox->rObj);
	}

	ox->initFun(ox);

	if(ox->computeFun == NULL) {
		// Should never happen
		Rf_error("Could not initialize Expectation function %s", ox->expType);
	}

	ox->isComplete = TRUE;

}

static void defaultSetVarGroup(omxExpectation *ox, FreeVarGroup *fvg)
{
	if (ox->freeVarGroup && ox->freeVarGroup != fvg) {
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
		if(strncmp(expType, entry->name, MAX_STRING_LEN) == 0) {
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

	return expect;
}

void omxExpectationPrint(omxExpectation* ox, char* d) {
	if(ox->printFun != NULL) {
		ox->printFun(ox);
	} else {
		mxLog("(Expectation, type %s) ", (ox->expType==NULL?"Untyped":ox->expType));
	}
}
