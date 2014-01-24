/*
 *  Copyright 2007-2012 The OpenMx Project
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
#include <errno.h>

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
	SEXP matListNames = getAttrib(matList, R_NamesSymbol);

	for(int index = 0; index < length(matList); index++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, index));		// This is the matrix + populations
		PROTECT(nextMat = VECTOR_ELT(nextLoc, 0));		// The first element of the list is the matrix of values
		currentState->matrixList[index] = omxNewMatrixFromRPrimitive(
			nextMat, currentState, 1, -index - 1);
		currentState->matrixList[index]->name = CHAR(STRING_ELT(matListNames, index));
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
	SEXP nextAlgTuple;
	currentState->numAlgs = length(algList);
	SEXP algListNames = getAttrib(algList, R_NamesSymbol);

	if(OMX_DEBUG) { Rprintf("Processing %d algebras.\n", currentState->numAlgs, length(algList)); }
	currentState->algebraList = (omxMatrix**) R_alloc(currentState->numAlgs, sizeof(omxMatrix*));

	for(int index = 0; index < currentState->numAlgs; index++) {
		currentState->algebraList[index] = omxInitMatrix(NULL, 0, 0, TRUE, currentState);
	}

	for(int index = 0; index < currentState->numAlgs; index++) {
		PROTECT(nextAlgTuple = VECTOR_ELT(algList, index));		// The next algebra or objective to process
		if(OMX_DEBUG) { Rprintf("Initializing algebra %d at location 0x%0x.\n", index, currentState->algebraList + index); }
		if(IS_S4_OBJECT(nextAlgTuple)) {		// This is an objective object.
			omxFillMatrixFromMxObjective(currentState->algebraList[index], nextAlgTuple, 1, index);
		} else {								// This is an algebra spec.
			SEXP initialValue, formula, dependencies;
			PROTECT(initialValue = VECTOR_ELT(nextAlgTuple, 0));
			omxFillMatrixFromRPrimitive(currentState->algebraList[index],
				initialValue, currentState, 1, index);
			UNPROTECT(1);	// initialValue
			PROTECT(formula = VECTOR_ELT(nextAlgTuple, 1));
			omxFillMatrixFromMxAlgebra(currentState->algebraList[index],
				formula, CHAR(STRING_ELT(algListNames, index)));
			UNPROTECT(1);	// formula
			PROTECT(dependencies = VECTOR_ELT(nextAlgTuple, 2));
			UNPROTECT(1);	// dependencies
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
	int numMats = currentState->numMats;
	int numAlgs = currentState->numAlgs;

	if(OMX_DEBUG) {Rprintf("Completed Algebras and Matrices.  Beginning Initial Compute.\n");}
	omxStateNextEvaluation(currentState);

	currentState->markMatrices = (int*) R_alloc(numMats + numAlgs, sizeof(int));

	for(int index = 0; index < numMats; index++) {
		omxRecompute(currentState->matrixList[index]);
	}

	for(int index = 0; index < numAlgs; index++) {
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

		const char *pathName, *fileName;
		const char __attribute__((unused)) *serverName;

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
			int __attribute__((unused)) portno = INTEGER(AS_INTEGER(VECTOR_ELT(nextLoc, next++)))[0];
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
The second element of the list is the upper bound.  The third element of the list
is a vector of mxIndices specifying the dependencies of the free parameter. 
The remaining elements of the list are 3-tuples.  These 3-tuples are (mxIndex, row, col).
*/
void omxProcessFreeVarList(SEXP varList, int n) {
	SEXP nextVar, nextLoc;
	if(OMX_VERBOSE) { Rprintf("Processing Free Parameters.\n"); }
	currentState->freeVarList = (omxFreeVar*) R_alloc (n, sizeof (omxFreeVar));			// Data for replacement of free vars
	for(int freeVarIndex = 0; freeVarIndex < n; freeVarIndex++) {
		int numDeps;
		PROTECT(nextVar = VECTOR_ELT(varList, freeVarIndex));
		int numLocs = length(nextVar) - 3;
		currentState->freeVarList[freeVarIndex].numLocations = numLocs;
		currentState->freeVarList[freeVarIndex].matrices = (int*) R_alloc(numLocs, sizeof(int));
		currentState->freeVarList[freeVarIndex].row		 = (int*) R_alloc(numLocs, sizeof(int));
		currentState->freeVarList[freeVarIndex].col		 = (int*) R_alloc(numLocs, sizeof(int));
		currentState->freeVarList[freeVarIndex].name = CHAR(STRING_ELT(GET_NAMES(varList), freeVarIndex));

		/* Lower Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));							// Position 0 is lower bound.
		currentState->freeVarList[freeVarIndex].lbound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[freeVarIndex].lbound)) currentState->freeVarList[freeVarIndex].lbound = NEG_INF;
		if(currentState->freeVarList[freeVarIndex].lbound == 0.0) currentState->freeVarList[freeVarIndex].lbound = 0.0;
		UNPROTECT(1); // nextLoc
		/* Upper Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));							// Position 1 is upper bound.
		currentState->freeVarList[freeVarIndex].ubound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[freeVarIndex].ubound)) currentState->freeVarList[freeVarIndex].ubound = INF;
		if(currentState->freeVarList[freeVarIndex].ubound == 0.0) currentState->freeVarList[freeVarIndex].ubound = -0.0;
		UNPROTECT(1); // nextLoc

		PROTECT(nextLoc = VECTOR_ELT(nextVar, 2));							// Position 2 is a vector of dependencies.
		numDeps = LENGTH(nextLoc);
		currentState->freeVarList[freeVarIndex].numDeps = numDeps;
		currentState->freeVarList[freeVarIndex].deps = (int*) R_alloc(numDeps, sizeof(int));
		for(int i = 0; i < numDeps; i++) {
			currentState->freeVarList[freeVarIndex].deps[i] = INTEGER(nextLoc)[i];
		}
		UNPROTECT(1); // nextLoc


		if(OMX_DEBUG) { 
			Rprintf("Free parameter %d bounded (%f, %f): %d locations\n", freeVarIndex, 
				currentState->freeVarList[freeVarIndex].lbound, 
				currentState->freeVarList[freeVarIndex].ubound, numLocs);
		}
		for(int locIndex = 0; locIndex < currentState->freeVarList[freeVarIndex].numLocations; locIndex++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, locIndex+3));
			int* theVarList = INTEGER(nextLoc);			// These come through as integers.

			int theMat = theVarList[0];			// Matrix is zero-based indexed.
			int theRow = theVarList[1];			// Row is zero-based.
			int theCol = theVarList[2];			// Column is zero-based.

			currentState->freeVarList[freeVarIndex].matrices[locIndex] = theMat;
			currentState->freeVarList[freeVarIndex].row[locIndex] = theRow;
			currentState->freeVarList[freeVarIndex].col[locIndex] = theCol;
			UNPROTECT(1); // nextLoc
		}
		UNPROTECT(1); // nextVar
	}

	int numMats = currentState->numMats;

	for(int freeVarIndex = 0; freeVarIndex < n; freeVarIndex++) {
		omxFreeVar* freeVar = currentState->freeVarList + freeVarIndex;
		int *deps   = freeVar->deps;
		int numDeps = freeVar->numDeps;
		for (int index = 0; index < numDeps; index++) {
			currentState->markMatrices[deps[index] + numMats] = 1;
		}
	}

}

/*
	intervalList is a list().  Each element refers to one confidence interval request.
	Each interval request is a length 5 vector of REAL.
	The first three elements are the matrixPointer, Row, and Column of the element
	for which bounds are to be calculated, and are cast to ints here for speed.
	The last two are the upper and lower boundaries for the confidence space (respectively).
*/
void omxProcessConfidenceIntervals(SEXP intervalList)  {
	SEXP nextVar;
	if(OMX_VERBOSE) { Rprintf("Processing Confidence Interval Requests.\n");}
	currentState->numIntervals = length(intervalList);
	if(OMX_DEBUG) {Rprintf("Found %d requests.\n", currentState->numIntervals); }
	currentState->intervalList = (omxConfidenceInterval*) R_alloc(currentState->numIntervals, sizeof(omxConfidenceInterval));
	for(int index = 0; index < currentState->numIntervals; index++) {
		omxConfidenceInterval *oCI = &(currentState->intervalList[index]);
		PROTECT(nextVar = VECTOR_ELT(intervalList, index));
		double* intervalInfo = REAL(nextVar);
		oCI->matrix = omxNewMatrixFromMxIndex( nextVar, currentState);	// Expects an R object
		oCI->row = (int) intervalInfo[1];		// Cast to int in C to save memory/Protection ops
		oCI->col = (int) intervalInfo[2];		// Cast to int in C to save memory/Protection ops
		oCI->lbound = intervalInfo[3];
		oCI->ubound = intervalInfo[4];
		UNPROTECT(1);
		oCI->max = R_NaReal;					// NAs, in case something goes wrong
		oCI->min = R_NaReal;
	}
	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d intervals requested.\n", currentState->numIntervals); }
}

int omxProcessConstraints(SEXP constraints)  {
	int ncnln = 0; 
	if(OMX_VERBOSE) { Rprintf("Processing Constraints.\n");}
	omxMatrix *arg1, *arg2;
	SEXP nextVar, nextLoc;
	currentState->numConstraints = length(constraints);
	if(OMX_DEBUG) {Rprintf("Found %d constraints.\n", currentState->numConstraints); }
	currentState->conList = (omxConstraint*) R_alloc(currentState->numConstraints, sizeof(omxConstraint));
	ncnln = 0;
	for(int constraintIndex = 0; constraintIndex < currentState->numConstraints; constraintIndex++) {
		PROTECT(nextVar = VECTOR_ELT(constraints, constraintIndex));
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));
		arg1 = omxNewMatrixFromMxIndex(nextLoc, currentState);
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));
		arg2 = omxNewMatrixFromMxIndex(nextLoc, currentState);
		PROTECT(nextLoc = AS_INTEGER(VECTOR_ELT(nextVar, 2)));
		currentState->conList[constraintIndex].opCode = INTEGER(nextLoc)[0];
		UNPROTECT(4);
		omxMatrix *args[2] = {arg1, arg2};
		currentState->conList[constraintIndex].result = omxNewAlgebraFromOperatorAndArgs(10, args, 2, currentState); // 10 = binary subtract
		omxRecompute(currentState->conList[constraintIndex].result);
		int nrows = currentState->conList[constraintIndex].result->rows;
		int ncols = currentState->conList[constraintIndex].result->cols;
		currentState->conList[constraintIndex].size = nrows * ncols;
		ncnln += currentState->conList[constraintIndex].size;
	}
	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d effective constraints.\n", ncnln); }
	return(ncnln);
}

void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin) {
	/* Set min and max limits */
	for(int index = 0; index < n; index++) {
		bl[index] = currentState->freeVarList[index].lbound;				// -Infinity'd be -10^20.
		bu[index] = currentState->freeVarList[index].ubound;				// Infinity would be at 10^20.
	}

	for(int index = n; index < n+nclin; index++) {						// At present, nclin == 0.
		bl[index] = NEG_INF; 							// Linear constraints have no bounds.
		bu[index] = INF;								// Because there are no linear constraints.
	}												    // But if there were, they would go here.

	int index = n + nclin;
	for(int constraintIndex = 0; constraintIndex < currentState->numConstraints; constraintIndex++) {		// Nonlinear constraints:
		if(OMX_DEBUG) { Rprintf("Constraint %d: ", constraintIndex);}
		switch(currentState->conList[constraintIndex].opCode) {
		case 0:									// Less than: Must be strictly less than 0.
			if(OMX_DEBUG) { Rprintf("Bounded at (-INF, 0).\n");}
			for(int offset = 0; offset < currentState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = -0.0;
				index++;
			}
			break;
		case 1:									// Equal: Must be roughly equal to 0.
			if(OMX_DEBUG) { Rprintf("Bounded at (-0, 0).\n");}
			for(int offset = 0; offset < currentState->conList[constraintIndex].size; offset++) {
				bl[index] = -0.0;
				bu[index] = 0.0;
				index++;
			}
			break;
		case 2:									// Greater than: Must be strictly greater than 0.
			if(OMX_DEBUG) { Rprintf("Bounded at (0, INF).\n");}
			for(int offset = 0; offset < currentState->conList[constraintIndex].size; offset++) {
				if(OMX_DEBUG) { Rprintf("\tBounds set for constraint %d.%d.\n", constraintIndex, offset);}
				bl[index] = 0.0;
				bu[index] = INF;
				index++;
			}
			break;
		default:
			if(OMX_DEBUG) { Rprintf("Bounded at (-INF, INF).\n");}
			for(int offset = 0; offset < currentState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = INF;
				index++;
			}
			break;
		}
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

