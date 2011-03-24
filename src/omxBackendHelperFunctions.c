/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>

#include "omxDefines.h"
#include "omxState.h"

extern omxState* currentState;		

int matchCaseInsensitive(const char *source, int lenSource, const char *target) {
	int lenTarget = strlen(target);
	return((lenSource == lenTarget)	&& (strncasecmp(source, target, lenSource) == 0));
}

int omxProcessMxDataEntities(SEXP data) {
	int errOut = FALSE;
	SEXP nextLoc;
	if(OMX_DEBUG) { Rprintf("Processing %d data source(s).\n", length(data));}
	currentState->numData = length(data);
	currentState->dataList = (omxData**) R_alloc(length(data), sizeof(omxData*));

	for(int index = 0; index < length(data); index++) {
		PROTECT(nextLoc = VECTOR_ELT(data, index));			// Retrieve the data object
		currentState->dataList[index] = omxNewDataFromMxData(NULL, nextLoc, currentState);
		if(OMX_DEBUG) {
			Rprintf("Data initialized at 0x%0xd = (%d x %d).\n",
				currentState->dataList[index], currentState->dataList[index]->rows, currentState->dataList[index]->cols);
		}
		UNPROTECT(1); // nextMat
		if(currentState->statusCode < 0) {
			errOut = TRUE;
			currentState->numData = index+1;
			break;
		}
	}
	return(errOut);
}

int omxProcessMxMatrixEntities(SEXP matList) {
	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}
	int errOut = FALSE;
	SEXP nextLoc, nextMat;
	currentState->numMats = length(matList);
	currentState->matrixList = (omxMatrix**) R_alloc(length(matList), sizeof(omxMatrix*));

	for(int index = 0; index < length(matList); index++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, index));		// This is the matrix + populations
		PROTECT(nextMat = VECTOR_ELT(nextLoc, 0));		// The first element of the list is the matrix of values
		currentState->matrixList[index] = omxNewMatrixFromRPrimitive(nextMat, currentState);
		if(OMX_DEBUG) {
			Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n",
				currentState->matrixList[index], currentState->matrixList[index]->rows, currentState->matrixList[index]->cols);
		}
		UNPROTECT(2); // nextLoc and nextMat
		if(currentState->statusCode < 0) {
			if(OMX_DEBUG) { Rprintf("Initialization Error processing %dth matrix.\n", index+1);}
			errOut = TRUE;
			currentState->numMats = index+1;
			break;
		}
	}
	return(errOut);
}

int omxProcessMxAlgebraEntities(SEXP algList) {
	int errOut = FALSE;
	SEXP nextAlgTuple, nextAlg;
	currentState->numAlgs = length(algList);
	SEXP algListNames = getAttrib(algList, R_NamesSymbol);

	if(OMX_DEBUG) { Rprintf("Processing %d algebras.\n", currentState->numAlgs, length(algList)); }
	currentState->algebraList = (omxMatrix**) R_alloc(currentState->numAlgs, sizeof(omxMatrix*));

	for(int index = 0; index < currentState->numAlgs; index++) {
		currentState->algebraList[index] = omxInitMatrix(NULL, 0,0, TRUE, currentState);
	}

	for(int index = 0; index < currentState->numAlgs; index++) {
		PROTECT(nextAlgTuple = VECTOR_ELT(algList, index));		// The next algebra or objective to process
		if(OMX_DEBUG) { Rprintf("Initializing algebra %d at location 0x%0x.\n", index, currentState->algebraList + index); }
		if(IS_S4_OBJECT(nextAlgTuple)) {		// This is an objective object.
			omxFillMatrixFromMxObjective(currentState->algebraList[index], nextAlgTuple);
		} else {								// This is an algebra spec.
			PROTECT(nextAlg = VECTOR_ELT(nextAlgTuple, 0));
			omxFillMatrixFromRPrimitive(currentState->algebraList[index],
				nextAlg, currentState);
			UNPROTECT(1);	// nextAlg
			PROTECT(nextAlg = VECTOR_ELT(nextAlgTuple, 1));
			omxFillMatrixFromMxAlgebra(currentState->algebraList[index],
				nextAlg, CHAR(STRING_ELT(algListNames, index)));
			UNPROTECT(1);	// nextAlg
		}
		UNPROTECT(1);	// nextAlgTuple
		if(currentState->statusCode < 0) {
			if(OMX_DEBUG) { Rprintf("Initialization Error processing %dth algebra.\n", index+1);}
			errOut = TRUE;
			currentState->numAlgs = index+1;
			break;
		}
	}
	return(errOut);
}

