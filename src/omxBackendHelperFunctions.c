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

#include <sys/stat.h>

#include "omxDefines.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"

/* Outside R Functions */
static int isDir(const char *path);

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

int omxInitialMatrixAlgebraCompute() {
	int errOut = FALSE;
	if(OMX_DEBUG) {Rprintf("Completed Algebras and Matrices.  Beginning Initial Compute.\n");}
	omxStateNextEvaluation(currentState);

	for(int index = 0; index < currentState->numMats; index++) {
		omxRecompute(currentState->matrixList[index]);
	}

	for(int index = 0; index < currentState->numAlgs; index++) {
		omxRecompute(currentState->algebraList[index]);
	}
	return(errOut);
}

int omxProcessObjectiveFunction(SEXP objective, int *n) {
	int errOut = FALSE;
	if(!isNull(objective)) {
		if(OMX_DEBUG) { Rprintf("Processing objective function.\n"); }
		currentState->objectiveMatrix = omxNewMatrixFromMxIndex(objective, currentState);
	} else {
		currentState->objectiveMatrix = NULL;
		*n = 0;
		currentState->numFreeParams = *n;
	}
	return(errOut);
}

/*
checkpointList is a list().  Each element refers to one checkpointing request.
Each interval request is a list of length 5.
The first element is an integer that specifies type: 0 = file, 1 = socket, 2=R_connection
For a file, the next two are the directory(string)  and file name (string).
For a socket, they are server (string) and port number (int).
For a connection, the next one is the R_connection SEXP object.
After that is an integer <type> specifier.  0 means minutes, 1 means iterations.
The last element is an integer count, indicating the number of <type>s per checkpoint.
*/
void omxProcessCheckpointOptions(SEXP checkpointList) {
	if(OMX_VERBOSE) { Rprintf("Processing Checkpoint Requests.\n");}
	currentState->numCheckpoints = length(checkpointList);
	if(OMX_DEBUG) {Rprintf("Found %d checkpoints.\n", currentState->numCheckpoints); }
	currentState->checkpointList = (omxCheckpoint*) R_alloc(currentState->numCheckpoints, sizeof(omxCheckpoint));
	SEXP nextLoc;

	for(int index = 0; index < currentState->numCheckpoints; index++) {
		omxCheckpoint *oC = &(currentState->checkpointList[index]);

		/* Initialize Checkpoint object */
		oC->socket = -1;
		oC->file = NULL;
		oC->connection = NULL;
		oC->time = 0;
		oC->numIterations = 0;
		oC->lastCheckpoint = 0;

		const char *pathName, *fileName, *serverName;

		PROTECT(nextLoc = VECTOR_ELT(checkpointList, index));
		int next = 0;
		oC->type = INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		switch(oC->type) {
		case OMX_FILE_CHECKPOINT:
			pathName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));			// FIXME: Might need PROTECT()ion
			fileName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));
			char sep ='/';

			if(!isDir(pathName)) {
				error("Unable to open directory %s for checkpoint storage.\n", pathName);
			}

			char* fullname = Calloc(strlen(pathName) + strlen(fileName) + 5, char);
			sprintf(fullname, "%s%c%s", pathName, sep, fileName);
			if(OMX_VERBOSE) { Rprintf("Opening File: %s\n", fullname); }
			oC->file = fopen(fullname, "w");
			if(!oC->file) {
				error("Unable to open file %s for checkpoint storage: %s.\n", fullname, strerror(errno));
			}
			Free(fullname);
			oC->saveHessian = FALSE;	// TODO: Decide if this should be true.
			break;

		case OMX_SOCKET_CHECKPOINT:
			serverName = CHAR(VECTOR_ELT(nextLoc, next++));
			// int portno = INTEGER(AS_INTEGER(VECTOR_ELT(nextLoc, next++)))[0]; // Commented to suppress "unused variable" warning.
			Rprintf("Warning NYI: Socket checkpoints Not Yet Implemented.\n");
			oC->saveHessian = FALSE;
			break;

		case OMX_CONNECTION_CHECKPOINT:	// NYI :::DEBUG:::
			oC->connection = VECTOR_ELT(nextLoc, next++);
			Rprintf("Warning NYI: Socket checkpoints Not Yet Implemented.\n");
			oC->saveHessian = FALSE;
			break;
		}

		int isCount = INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		if(isCount) {
			oC->numIterations = INTEGER(AS_INTEGER(VECTOR_ELT(nextLoc, next++)))[0];
		} else {
			oC->time = REAL(AS_NUMERIC(VECTOR_ELT(nextLoc, next++)))[0] * 60;	// Constrained to seconds.
			if(oC->time < 1) oC->time = 1;										// Constrained to at least one.
		}

		UNPROTECT(1); /* nextLoc */

	}
}

/*
varList is a list().  Each element of this list corresponds to one free parameter.
Each free parameter is a list.  The first element of this list is the lower bound.
The second element of the list is the upper bound.  The remaining elements of this
list are 3-tuples.  These 3-tuples are (mxIndex, row, col).
*/
void omxProcessFreeVarList(SEXP varList, int n) {
	SEXP nextVar, nextLoc;
	if(OMX_VERBOSE) { Rprintf("Processing Free Parameters.\n"); }
	currentState->freeVarList = (omxFreeVar*) R_alloc (n, sizeof (omxFreeVar));			// Data for replacement of free vars
	for(int freeVarIndex = 0; freeVarIndex < n; freeVarIndex++) {
		PROTECT(nextVar = VECTOR_ELT(varList, freeVarIndex));
		int numLocs = length(nextVar) - 2;
		currentState->freeVarList[freeVarIndex].numLocations = numLocs;
		currentState->freeVarList[freeVarIndex].location = (double**) R_alloc(numLocs, sizeof(double*));
		currentState->freeVarList[freeVarIndex].matrices = (int*) R_alloc(numLocs, sizeof(int));
		currentState->freeVarList[freeVarIndex].name = CHAR(STRING_ELT(GET_NAMES(varList), freeVarIndex));

		/* Lower Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));							// Position 0 is lower bound.
		currentState->freeVarList[freeVarIndex].lbound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[freeVarIndex].lbound)) currentState->freeVarList[freeVarIndex].lbound = NEG_INF;
		if(currentState->freeVarList[freeVarIndex].lbound == 0.0) currentState->freeVarList[freeVarIndex].lbound = 0.0;
		UNPROTECT(1); // NextLoc
		/* Upper Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));							// Position 1 is upper bound.
		currentState->freeVarList[freeVarIndex].ubound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[freeVarIndex].ubound)) currentState->freeVarList[freeVarIndex].ubound = INF;
		if(currentState->freeVarList[freeVarIndex].ubound == 0.0) currentState->freeVarList[freeVarIndex].ubound = -0.0;
		UNPROTECT(1); // NextLoc

		if(OMX_DEBUG) { 
			Rprintf("Free parameter %d bounded (%f, %f): %d locations\n", freeVarIndex, 
				currentState->freeVarList[freeVarIndex].lbound, 
				currentState->freeVarList[freeVarIndex].ubound, numLocs);
		}
		for(int locIndex = 0; locIndex < currentState->freeVarList[freeVarIndex].numLocations; locIndex++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, locIndex+2));
			int* theVarList = INTEGER(nextLoc);			// These come through as integers.

			int theMat = theVarList[0];			// Matrix is zero-based indexed.
			int theRow = theVarList[1];			// Row is zero-based.
			int theCol = theVarList[2];			// Column is zero-based.

			double* locationMatrixElement = omxLocationOfMatrixElement(currentState->matrixList[theMat], theRow, theCol);
			currentState->freeVarList[freeVarIndex].location[locIndex] = locationMatrixElement;
			currentState->freeVarList[freeVarIndex].matrices[locIndex] = theMat;
			UNPROTECT(1); // nextLoc
		}
		UNPROTECT(1); // nextVar
	}
}

/*
*  Acknowledgement:
*  This function is duplicated from the function of the same name in the R source code.
*  The function appears in src/main/sysutils.c
*  Thanks to Michael Spiegel for finding it.
*/
static int isDir(const char *path)
{
    struct stat sb;
    int isdir = 0;
    if(!path) return 0;
    if(stat(path, &sb) == 0) {
        isdir = (sb.st_mode & S_IFDIR) > 0; /* is a directory */
        /* We want to know if the directory is writable by this user,
           which mode does not tell us */
        isdir &= (access(path, W_OK) == 0);
    }
    return isdir;
}

