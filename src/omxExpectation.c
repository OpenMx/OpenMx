/*
 *  Copyright 2007-2012 The OpenMx Project
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
	void (*initFun)(omxExpectation*, SEXP);
};

extern void omxInitNormalExpectation(omxExpectation *ox, SEXP rObj);
extern void omxInitLISRELExpectation(omxExpectation *ox, SEXP rObj);
extern void omxInitStateSpaceExpectation(omxExpectation *ox, SEXP rObj);
extern void omxInitRAMExpectation(omxExpectation *ox, SEXP rObj);

static const omxExpectationTableEntry omxExpectationSymbolTable[] = {
	{"MxExpectationLISREL",			&omxInitLISRELExpectation},
	{"MxExpectationStateSpace",			&omxInitStateSpaceExpectation},
	{"MxExpectationNormal", 		&omxInitNormalExpectation},
	{"MxExpectationRAM",			&omxInitRAMExpectation},
	{ "", 0 }
};

void omxInitEmptyExpectation(omxExpectation *ox) {
	/* Sets everything to NULL to avoid bad pointer calls */
	
  memset(ox, 0, sizeof(*ox));
}

void omxFreeExpectationArgs(omxExpectation *ox) {
	if(ox==NULL) return;
    
	/* Completely destroy the Expectation function tree */
	if(OMX_DEBUG) {Rprintf("Freeing %s Expectation object at 0x%x.\n", (ox->expType == NULL?"untyped":ox->expType), ox);}
	if(ox->destructFun != NULL) {
		if(OMX_DEBUG) {Rprintf("Calling Expectation destructor for 0x%x.\n", ox);}
		ox->destructFun(ox);
	}
}

void omxExpectationRecompute(omxExpectation *ox) {
	if(OMX_DEBUG_ALGEBRA) { 
	    Rprintf("Expectation recompute: 0x%0x\n", ox);
	}

	omxExpectationCompute(ox);
}

void omxExpectationCompute(omxExpectation *ox) {
	if(OMX_DEBUG_ALGEBRA) { 
	    Rprintf("Expectation compute: 0x%0x\n", ox);
	}

	ox->computeFun(ox);
}

omxMatrix* omxGetExpectationComponent(omxExpectation* ox, omxFitFunction* off, char* component) {

	if(component == NULL) return NULL;

	/* Hard-wired expectation components */
	if(!strncmp("dataColumns", component, 11)) {
		return ox->dataColumns;
	}

	if(ox->componentFun == NULL) return NULL;

	return(ox->componentFun(ox, off, component));
	
}

void omxSetExpectationComponent(omxExpectation* ox, omxFitFunction* off, char* component, omxMatrix* om) {
	if(!strcmp(ox->expType, "omxStateSpaceExpectation")) {
		ox->mutateFun(ox, off, component, om);
	}
}

omxExpectation* omxNewExpectationFromMxExpectation(SEXP mxobj, int expNum, omxState* os) {
	
	omxExpectation* expectation = omxNewIncompleteExpectation(mxobj, expNum, os);
	omxCompleteExpectation(expectation);
	return expectation;
}

omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState) {

	if(OMX_DEBUG) {Rprintf("Duplicating Expectation 0x%x\n", src);}

	// if(src == NULL) {
	// 	return NULL;
	// }
	// 
	// omxExpectation* tgt = (omxExpectation*) R_alloc(1, sizeof(omxExpectation));
	// omxInitEmptyExpectation(tgt);
	// 
	// tgt->initFun 					= src->initFun;
	// tgt->destructFun 				= src->destructFun;
	// tgt->repopulateFun 				= src->repopulateFun;
	// tgt->computeFun 				= src->computeFun;
	// tgt->componentFun				= src->componentFun;
	// tgt->populateAttrFun 			= src->populateAttrFun;
	// tgt->setFinalReturns 			= src->setFinalReturns;
	// tgt->sharedArgs					= src->sharedArgs;
	// tgt->currentState				= newState;
	// tgt->rObj						= src->rObj;
	// tgt->data						= src->data;
	// tgt->dataColumns				= omxLookupDuplicateElement(newState, src->dataColumns);
	// tgt->defVars					= src->defVars;
	// tgt->numDefs					= src->numDefs;
	// int numDefs = tgt->numDefs;
	// // for(int i = 0; i < numDefs; i++) {
	// // 	int thisCount = tgt->defVars[i].numLocations;
	// // 	for(int index = 0; index < thisCount; index++) {
	// // 		tgt->defVars[i].matrices[index] = omxLookupDuplicateElement(newState, src->defVars[i].matrices[index]);
	// // 	}
	// // }
	// 
	// tgt->numOrdinal					= src->numOrdinal;
	// tgt->thresholds					= src->thresholds;
	// int nCols = tgt->dataColumns->rows;
	// for(int i = 0; i < nCols; i++) {
	// 	if(tgt->thresholds[i].matrix != NULL) {
	// 		tgt->thresholds[i].matrix = omxLookupDuplicateElement(newState, src->thresholds[i].matrix);
	// 	}
	// }
	// 
	// tgt->expNum						= src->expNum;
	// 
	//     strncpy(tgt->expType, src->expType, MAX_STRING_LEN);
	// 
	// return tgt;

	return omxNewIncompleteExpectation(src->rObj, src->expNum, newState);

}

omxExpectation* omxNewIncompleteExpectation(SEXP rObj, int expNum, omxState* os) {

	SEXP ExpectationClass;
	const char* expType;
	omxExpectation* expect = (omxExpectation*) R_alloc(1, sizeof(omxExpectation));
	omxInitEmptyExpectation(expect);
	
	/* Get Expectation Type */
	PROTECT(ExpectationClass = STRING_ELT(getAttrib(rObj, install("class")), 0));
	expType = CHAR(ExpectationClass);

	/* Switch based on Expectation type. */ 
	const omxExpectationTableEntry *entry = omxExpectationSymbolTable;
	while (entry->initFun) {
		if(strncmp(expType, entry->name, MAX_STRING_LEN) == 0) {
		        expect->expType = entry->name;
			expect->initFun = entry->initFun;
			break;
		}
		entry += 1;
	}

	UNPROTECT(1);	/* ExpectationClass */

	if(!expect->initFun) {
		char newError[MAX_STRING_LEN];
		sprintf(newError, "Expectation function %s not implemented.\n", (expect->expType==NULL?"Untyped":expect->expType));
		omxRaiseError(os, -1, newError);
		return NULL;
	}

	expect->rObj = rObj;
	expect->expNum = expNum;
	expect->currentState = os;
	
	return expect;
}

omxExpectation* omxNewExpectationFromExpectationIndex(int expIndex, omxState* os) {

	omxExpectation* ox = os->expectationList[expIndex];
	
	if(!ox->isComplete) omxCompleteExpectation(ox);
	
	return ox;
}

void omxExpectationProcessDataStructures(omxExpectation* ox, SEXP rObj){

	int index, numDefs, nextDef, numCols, numOrdinal=0;
	SEXP nextMatrix, itemList, nextItem, threshMatrix; 
	
	if(rObj == NULL) return;

	if(OMX_DEBUG) { Rprintf("Retrieving data.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data")));
	ox->data = omxNewDataFromMxDataPtr(nextMatrix, ox->currentState);
	UNPROTECT(1); // nextMatrix

	if(OMX_DEBUG && ox->currentState->parentState == NULL) {
		Rprintf("Accessing variable mapping structure.\n");
	}

	if (R_has_slot(rObj, install("dataColumns"))) {
		PROTECT(nextMatrix = GET_SLOT(rObj, install("dataColumns")));
		ox->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, ox->currentState, 0, 0);
		if(OMX_DEBUG && ox->currentState->parentState == NULL) {
			omxPrint(ox->dataColumns, "Variable mapping");
		}
		UNPROTECT(1); // dataColumns
	
		numCols = ox->dataColumns->cols;

		if (R_has_slot(rObj, install("thresholds"))) {
			if(OMX_DEBUG && ox->currentState->parentState == NULL) {
				Rprintf("Accessing Threshold matrix.\n");
			}
			PROTECT(threshMatrix = GET_SLOT(rObj, install("thresholds")));

			if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
				if(OMX_DEBUG && ox->currentState->parentState == NULL) {
					Rprintf("Accessing Threshold Mappings.\n");
				}
        
				/* Process the data and threshold mapping structures */
				/* if (threshMatrix == NA_INTEGER), then we could ignore the slot "thresholdColumns"
				 * and fill all the thresholds with {NULL, 0, 0}.
				 * However the current path does not have a lot of overhead. */
				PROTECT(nextMatrix = GET_SLOT(rObj, install("thresholdColumns")));
				PROTECT(itemList = GET_SLOT(rObj, install("thresholdLevels")));
				int* thresholdColumn, *thresholdNumber;
				thresholdColumn = INTEGER(nextMatrix);
				thresholdNumber = INTEGER(itemList);
				ox->thresholds = (omxThresholdColumn *) R_alloc(numCols, sizeof(omxThresholdColumn));
				for(index = 0; index < numCols; index++) {
					if(thresholdColumn[index] == NA_INTEGER) {	// Continuous variable
						if(OMX_DEBUG && ox->currentState->parentState == NULL) {
							Rprintf("Column %d is continuous.\n", index);
						}
						ox->thresholds[index].matrix = NULL;
						ox->thresholds[index].column = 0;
						ox->thresholds[index].numThresholds = 0;
					} else {
						ox->thresholds[index].matrix = omxNewMatrixFromMxIndex(threshMatrix, 
												       ox->currentState);
						ox->thresholds[index].column = thresholdColumn[index];
						ox->thresholds[index].numThresholds = thresholdNumber[index];
						if(OMX_DEBUG && ox->currentState->parentState == NULL) {
							Rprintf("Column %d is ordinal with %d thresholds in threshold column %d.\n", 
								index, thresholdColumn[index], thresholdNumber[index]);
						}
						numOrdinal++;
					}
				}
				if(OMX_DEBUG && ox->currentState->parentState == NULL) {
					Rprintf("%d threshold columns processed.\n", numOrdinal);
				}
				UNPROTECT(2); /* nextMatrix and itemList ("thresholds" and "thresholdColumns") */
				ox->numOrdinal = numOrdinal;
			} else {
				if (OMX_DEBUG && ox->currentState->parentState == NULL) {
					Rprintf("No thresholds matrix; not processing thresholds.");
				}
				ox->thresholds = NULL;
				ox->numOrdinal = 0;
			}
			UNPROTECT(1); /* threshMatrix */
		}
	}

	if(!R_has_slot(rObj, install("definitionVars"))) {
		ox->numDefs = 0;
		ox->defVars = NULL;
	} else {	
		if(OMX_DEBUG && ox->currentState->parentState == NULL) {
			Rprintf("Accessing definition variables structure.\n");
		}
		PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
		numDefs = length(nextMatrix);
		ox->numDefs = numDefs;
		if(OMX_DEBUG && ox->currentState->parentState == NULL) {
			Rprintf("Number of definition variables is %d.\n", numDefs);
		}
		ox->defVars = (omxDefinitionVar *) R_alloc(numDefs, sizeof(omxDefinitionVar));
		for(nextDef = 0; nextDef < numDefs; nextDef++) {
			SEXP dataSource, columnSource, depsSource; 
			int nextDataSource, numDeps;

			PROTECT(itemList = VECTOR_ELT(nextMatrix, nextDef));
			PROTECT(dataSource = VECTOR_ELT(itemList, 0));
			nextDataSource = INTEGER(dataSource)[0];
			if(OMX_DEBUG && ox->currentState->parentState == NULL) {
				Rprintf("Data source number is %d.\n", nextDataSource);
			}
			ox->defVars[nextDef].data = nextDataSource;
			ox->defVars[nextDef].source = ox->currentState->dataList[nextDataSource];
			PROTECT(columnSource = VECTOR_ELT(itemList, 1));
			if(OMX_DEBUG && ox->currentState->parentState == NULL) {
				Rprintf("Data column number is %d.\n", INTEGER(columnSource)[0]);
			}
			ox->defVars[nextDef].column = INTEGER(columnSource)[0];
			PROTECT(depsSource = VECTOR_ELT(itemList, 2));
			numDeps = LENGTH(depsSource);
			ox->defVars[nextDef].numDeps = numDeps;
			ox->defVars[nextDef].deps = (int*) R_alloc(numDeps, sizeof(int));
			for(int i = 0; i < numDeps; i++) {
				ox->defVars[nextDef].deps[i] = INTEGER(depsSource)[i];
			}
			UNPROTECT(3); // unprotect dataSource, columnSource, and depsSource

			ox->defVars[nextDef].numLocations = length(itemList) - 3;
			ox->defVars[nextDef].matrices = (int *) R_alloc(length(itemList) - 3, sizeof(int));
			ox->defVars[nextDef].rows = (int *) R_alloc(length(itemList) - 3, sizeof(int));
			ox->defVars[nextDef].cols = (int *) R_alloc(length(itemList) - 3, sizeof(int));
			for(index = 3; index < length(itemList); index++) {
				PROTECT(nextItem = VECTOR_ELT(itemList, index));
				ox->defVars[nextDef].matrices[index-3] = INTEGER(nextItem)[0];
				ox->defVars[nextDef].rows[index-3] = INTEGER(nextItem)[1];
				ox->defVars[nextDef].cols[index-3] = INTEGER(nextItem)[2];
				UNPROTECT(1); // unprotect nextItem
			}
			UNPROTECT(1); // unprotect itemList
		}
		UNPROTECT(1); // unprotect nextMatrix
	}
	
}

void omxCompleteExpectation(omxExpectation *ox) {
	
	char errorCode[MAX_STRING_LEN];
	
	if(OMX_DEBUG) {Rprintf("Completing Expectation 0x%x, type %s.\n", 
		ox, ((ox==NULL || ox->expType==NULL)?"Untyped":ox->expType));}
		
	if(ox == NULL) {
		if(OMX_DEBUG) { Rprintf("Could not complete NULL expectation.\n Somebody passed NULL to omxCompleteExpectation.  And there's nothing I can do about it.  Hopefully it'll be caught later.");}
		return;
		// Nothing to do at this point.
	}
	
	if(ox->isComplete) return;	// Already done; nothing more to do.
	
	omxState* os = ox->currentState;

	if(ox->rObj == NULL || ox->initFun == NULL ) {
		char newError[MAX_STRING_LEN];
		sprintf(newError, "Could not complete expectation %s.\n", (ox->expType==NULL?"Untyped":ox->expType));
		omxRaiseError(os, -1, newError);
		return;
	}

	omxExpectationProcessDataStructures(ox, ox->rObj);

	ox->initFun(ox, ox->rObj);

	if(ox->computeFun == NULL) {// If initialization fails, error code goes in argStruct
		if(os->statusCode != 0) {
			strncpy(errorCode, os->statusMsg, 150); // Report a status error
		} else {
			// If no error code is reported, we report that.
  			strncpy(errorCode, "No error code reported.", 25);
		}
		if(ox->argStruct != NULL) {
			strncpy(errorCode, (char*)(ox->argStruct), 51);
		}
		char newError[MAX_STRING_LEN];
		sprintf(newError, "Could not initialize Expectation function %s.  Error: %s\n", 
				ox->expType, errorCode);
		omxRaiseError(os, -1, newError);
		return;
	}

	ox->isComplete = TRUE;

}

void omxExpectationPrint(omxExpectation* ox, char* d) {
	if(ox->printFun != NULL) {
		ox->printFun(ox);
	} else {
		Rprintf("(Expectation, type %s) ", (ox->expType==NULL?"Untyped":ox->expType));
	}
}
