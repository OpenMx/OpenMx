/*
 *  Copyright 2007-2013 The OpenMx Project
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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <sys/stat.h>
#include <errno.h>

#include "omxDefines.h"
#include "glue.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "Compute.h"

/* Outside R Functions */
static int isDir(const char *path);

int matchCaseInsensitive(const char *source, const char *target) {
	return strcasecmp(source, target) == 0;
}

void omxProcessMxDataEntities(SEXP data) {
	SEXP nextLoc;
	if(OMX_DEBUG) { mxLog("Processing %d data source(s).", length(data));}

	for(int index = 0; index < length(data); index++) {
		PROTECT(nextLoc = VECTOR_ELT(data, index));			// Retrieve the data object
		omxNewDataFromMxData(nextLoc, globalState);
		if(OMX_DEBUG) {
			mxLog("Data initialized at %p = (%d x %d).",
				globalState->dataList[index], globalState->dataList[index]->rows, globalState->dataList[index]->cols);
		}
	}
}

void omxProcessMxMatrixEntities(SEXP matList) {
	if(OMX_DEBUG) { mxLog("Processing %d matrix(ces).", length(matList));}
	SEXP nextLoc, nextMat;
	globalState->matrixList.clear();
	SEXP matListNames = getAttrib(matList, R_NamesSymbol);

	for(int index = 0; index < length(matList); index++) {
		omxManageProtectInsanity protectManager;
		PROTECT(nextLoc = VECTOR_ELT(matList, index));		// This is the matrix + populations
		PROTECT(nextMat = VECTOR_ELT(nextLoc, 0));		// The first element of the list is the matrix of values
		omxMatrix *mat = omxNewMatrixFromRPrimitive(nextMat, globalState, 1, -index - 1);
		globalState->matrixList.push_back(mat);
		globalState->matrixList[index]->name = CHAR(STRING_ELT(matListNames, index));
		if(OMX_DEBUG) {
			mxLog("Matrix initialized at %p = (%d x %d).",
				globalState->matrixList[index], globalState->matrixList[index]->rows, globalState->matrixList[index]->cols);
		}
		if (isErrorRaised(globalState)) return;
	}
}

void omxProcessMxAlgebraEntities(SEXP algList) {
	SEXP nextAlgTuple;
	SEXP algListNames = getAttrib(algList, R_NamesSymbol);

	if(OMX_DEBUG) { mxLog("Processing %d algebras.", length(algList)); }

	for(int index = 0; index < length(algList); index++) {
		globalState->algebraList.push_back(omxInitMatrix(NULL, 0, 0, TRUE, globalState));
	}

	for(int index = 0; index < length(algList); index++) {
		omxManageProtectInsanity protectManager;
		PROTECT(nextAlgTuple = VECTOR_ELT(algList, index));		// The next algebra or fit function to process
		if(IS_S4_OBJECT(nextAlgTuple)) {
			// delay until expectations are ready
		} else {								// This is an algebra spec.
			SEXP initialValue, formula;
			PROTECT(initialValue = VECTOR_ELT(nextAlgTuple, 0));
			omxFillMatrixFromRPrimitive(globalState->algebraList[index],
				initialValue, globalState, 1, index);
			PROTECT(formula = VECTOR_ELT(nextAlgTuple, 1));
			omxFillMatrixFromMxAlgebra(globalState->algebraList[index],
				formula, CHAR(STRING_ELT(algListNames, index)));
		}
		if (isErrorRaised(globalState)) return;
	}
}

void omxProcessMxFitFunction(SEXP algList)
{
	SEXP nextAlgTuple;
	SEXP algListNames = getAttrib(algList, R_NamesSymbol);

	for(int index = 0; index < length(algList); index++) {
		PROTECT(nextAlgTuple = VECTOR_ELT(algList, index));		// The next algebra or fit function to process
		if(IS_S4_OBJECT(nextAlgTuple)) {
			SEXP fitFunctionClass;
			PROTECT(fitFunctionClass = STRING_ELT(getAttrib(nextAlgTuple, install("class")), 0));
			const char *fitType = CHAR(fitFunctionClass);
			omxMatrix *fm = globalState->algebraList[index];
			omxFillMatrixFromMxFitFunction(fm, fitType, index);
			fm->fitFunction->rObj = nextAlgTuple;
			fm->name = CHAR(STRING_ELT(algListNames, index));
			UNPROTECT(1);	// fitFunctionClass
		}
		if (isErrorRaised(globalState)) return;
		UNPROTECT(1); //nextAlgTuple
	}
}

void omxCompleteMxFitFunction(SEXP algList)
{
	SEXP nextAlgTuple;

	for(int index = 0; index < length(algList); index++) {
		PROTECT(nextAlgTuple = VECTOR_ELT(algList, index));             // The next algebra or fit function to process
		if(IS_S4_OBJECT(nextAlgTuple)) {
			omxMatrix *fm = globalState->algebraList[index];
			setFreeVarGroup(fm->fitFunction, Global->freeGroup[0]);
			omxCompleteFitFunction(fm);
		}
		UNPROTECT(1);
	}
}

void omxProcessMxExpectationEntities(SEXP expList) {
	if(OMX_DEBUG) { mxLog("Initializing %d Model Expectation(s).", length(expList));}
	SEXP nextExp;
	SEXP eNames = getAttrib(expList, R_NamesSymbol);

	for(int index = 0; index < length(expList); index++) {
		PROTECT(nextExp = VECTOR_ELT(expList, index));
		omxExpectation *ex = omxNewIncompleteExpectation(nextExp, index, globalState);
		ex->name = CHAR(STRING_ELT(eNames, index));
		globalState->expectationList.push_back(ex);
		if(OMX_DEBUG) {
			mxLog("%s incomplete expectation set up at %p.",
				(globalState->expectationList[index]->expType
					== NULL ? "Untyped" : globalState->expectationList[index]->expType),
					 globalState->expectationList[index]);
		}
		if (isErrorRaised(globalState)) return;
	}
}


void omxCompleteMxExpectationEntities() {
	if(OMX_DEBUG) { mxLog("Completing %lu Model Expectation(s).", globalState->expectationList.size());}
	
	for(size_t index = 0; index < globalState->expectationList.size(); index++) {
		omxCompleteExpectation(globalState->expectationList[index]);
		if(OMX_DEBUG) {
			mxLog("%s expectation completed at %p.",
				(globalState->expectationList[index]->expType
					== NULL ? "Untyped" : globalState->expectationList[index]->expType),
					 globalState->expectationList[index]);
		}
		if (isErrorRaised(globalState)) return;
	}
}

void omxProcessMxComputeEntities(SEXP computeList)
{
	SEXP rObj, s4class;

	for(int index = 0; index < length(computeList); index++) {
		PROTECT(rObj = VECTOR_ELT(computeList, index));
		PROTECT(s4class = STRING_ELT(getAttrib(rObj, install("class")), 0));
		omxCompute *compute = omxNewCompute(globalState, CHAR(s4class));
		compute->initFromFrontend(rObj);
		Global->computeList.push_back(compute);
	}
}

void omxInitialMatrixAlgebraCompute() {
	size_t numMats = globalState->matrixList.size();
	int numAlgs = globalState->algebraList.size();

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}
	omxStateNextEvaluation(globalState);

	for(size_t index = 0; index < numMats; index++) {
		omxRecompute(globalState->matrixList[index]);
	}

	for(int index = 0; index < numAlgs; index++) {
		omxMatrix *matrix = globalState->algebraList[index];
		omxInitialCompute(matrix);
	}
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
	if(OMX_VERBOSE) { mxLog("Processing Checkpoint Requests.");}
	globalState->numCheckpoints = length(checkpointList);
	if(OMX_DEBUG) {mxLog("Found %d checkpoints.", globalState->numCheckpoints); }
	globalState->checkpointList = (omxCheckpoint*) R_alloc(globalState->numCheckpoints, sizeof(omxCheckpoint));
	SEXP nextLoc;

	for(int index = 0; index < globalState->numCheckpoints; index++) {
		omxCheckpoint *oC = &(globalState->checkpointList[index]);

		/* Initialize Checkpoint object */
		oC->file = NULL;
		oC->connection = NULL;
		oC->time = 0;
		oC->numIterations = 0;
		oC->lastCheckpoint = 0;

		const char *pathName, *fileName;
		const char __attribute__((unused)) *serverName;

		PROTECT(nextLoc = VECTOR_ELT(checkpointList, index));
		int next = 0;
		oC->type = (omxCheckpointType) INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		switch(oC->type) {
		case OMX_FILE_CHECKPOINT:{
			pathName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));			// FIXME: Might need PROTECT()ion
			fileName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));
			char sep ='/';

			if(!isDir(pathName)) {
				error("Unable to open directory %s for checkpoint storage.\n", pathName);
			}

			char* fullname = Calloc(strlen(pathName) + strlen(fileName) + 5, char);
			sprintf(fullname, "%s%c%s", pathName, sep, fileName);
			if(OMX_VERBOSE) { mxLog("Opening File: %s", fullname); }
			oC->file = fopen(fullname, "w");
			if(!oC->file) {
				error("Unable to open file %s for checkpoint storage: %s.\n", fullname, strerror(errno));
			}
			Free(fullname);
			oC->saveHessian = FALSE;	// TODO: Decide if this should be true.
			break;}

		case OMX_CONNECTION_CHECKPOINT:{	// NYI :::DEBUG:::
			oC->connection = VECTOR_ELT(nextLoc, next++);
			error("Warning NYI: Socket checkpoints Not Yet Implemented.\n");
			oC->saveHessian = FALSE;
			break;}
		}

		int isCount = INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		if(isCount) {
			oC->numIterations = INTEGER(AS_INTEGER(VECTOR_ELT(nextLoc, next++)))[0];
		} else {
			oC->time = REAL(AS_NUMERIC(VECTOR_ELT(nextLoc, next++)))[0] * 60;	// Constrained to seconds.
			if(oC->time < 1) oC->time = 1;										// Constrained to at least one.
		}
	}
}

void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues)
{
	if(OMX_VERBOSE) { mxLog("Processing Free Parameters."); }

	{
		FreeVarGroup *fvg = new FreeVarGroup;
		fvg->id = FREEVARGROUP_ALL;   // all variables
		Global->freeGroup.push_back(fvg);

		fvg = new FreeVarGroup;
		fvg->id = FREEVARGROUP_NONE;  // no variables
		Global->freeGroup.push_back(fvg);
	}

	SEXP nextVar, nextLoc;
	int numVars = length(varList);
	startingValues->resize(numVars);
	for (int fx = 0; fx < numVars; fx++) {
		omxManageProtectInsanity mpi;

		omxFreeVar *fv = new omxFreeVar;
		// default group has free all variables
		Global->freeGroup[0]->vars.push_back(fv);

		fv->name = CHAR(STRING_ELT(GET_NAMES(varList), fx));
		PROTECT(nextVar = VECTOR_ELT(varList, fx));

		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));
		fv->lbound = REAL(nextLoc)[0];
		if (ISNA(fv->lbound)) fv->lbound = NEG_INF;
		if (fv->lbound == 0.0) fv->lbound = 0.0;

		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));
		fv->ubound = REAL(nextLoc)[0];
		if (ISNA(fv->ubound)) fv->ubound = INF;
		if (fv->ubound == 0.0) fv->ubound = -0.0;

		PROTECT(nextLoc = VECTOR_ELT(nextVar, 2));
		int groupCount = length(nextLoc);
		for (int gx=0; gx < groupCount; ++gx) {
			int group = INTEGER(nextLoc)[gx];
			if (group == 0) continue;
			Global->findOrCreateVarGroup(group)->vars.push_back(fv);
		}

		PROTECT(nextLoc = VECTOR_ELT(nextVar, 3));
		int numDeps = LENGTH(nextLoc);
		fv->numDeps = numDeps;
		fv->deps = (int*) R_alloc(numDeps, sizeof(int));
		for (int i = 0; i < numDeps; i++) {
			fv->deps[i] = INTEGER(nextLoc)[i];
		}

		int numLocs = length(nextVar) - 5;
		if(OMX_DEBUG) { 
			mxLog("Free parameter %d bounded (%f, %f): %d locations", fx, 
			      fv->lbound, fv->ubound, numLocs);
		}
		for(int locIndex = 0; locIndex < numLocs; locIndex++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, locIndex+4));
			int* theVarList = INTEGER(nextLoc);

			omxFreeVarLocation loc;
			loc.matrix = theVarList[0];
			loc.row = theVarList[1];
			loc.col = theVarList[2];

			fv->locations.push_back(loc);
		}
		PROTECT(nextLoc = VECTOR_ELT(nextVar, length(nextVar)-1));
		double sv = REAL(nextLoc)[0];
		/*if (sv < fv->lbound) {
			warning("Moving starting value of parameter '%s' within bounds %f -> %f",
				fv->name, sv, fv->lbound);
			sv = fv->lbound;
		} else if (sv > fv->ubound) {
			warning("Moving starting value of parameter '%s' within bounds %f -> %f",
				fv->name, sv, fv->ubound);
			sv = fv->ubound;
		}*/
		(*startingValues)[fx] = sv;
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
	if(OMX_VERBOSE) { mxLog("Processing Confidence Interval Requests.");}
	Global->numIntervals = length(intervalList);
	if(OMX_DEBUG) {mxLog("Found %d requests.", Global->numIntervals); }
	Global->intervalList = (omxConfidenceInterval*) R_alloc(Global->numIntervals, sizeof(omxConfidenceInterval));
	for(int index = 0; index < Global->numIntervals; index++) {
		omxConfidenceInterval *oCI = &(Global->intervalList[index]);
		PROTECT(nextVar = VECTOR_ELT(intervalList, index));
		double* intervalInfo = REAL(nextVar);
		oCI->matrix = omxMatrixLookupFromState1( nextVar, globalState);	// Expects an R object
		oCI->row = (int) intervalInfo[1];		// Cast to int in C to save memory/Protection ops
		oCI->col = (int) intervalInfo[2];		// Cast to int in C to save memory/Protection ops
		oCI->lbound = intervalInfo[3];
		oCI->ubound = intervalInfo[4];
		oCI->max = R_NaReal;					// NAs, in case something goes wrong
		oCI->min = R_NaReal;
	}
	if(OMX_VERBOSE) { mxLog("Processed."); }
	if(OMX_DEBUG) { mxLog("%d intervals requested.", Global->numIntervals); }
}

void omxProcessConstraints(SEXP constraints)  {
	int ncnln = 0; 
	if(OMX_VERBOSE) { mxLog("Processing Constraints.");}
	omxMatrix *arg1, *arg2;
	SEXP nextVar, nextLoc;
	globalState->numConstraints = length(constraints);
	if(OMX_DEBUG) {mxLog("Found %d constraints.", globalState->numConstraints); }
	globalState->conList = (omxConstraint*) R_alloc(globalState->numConstraints, sizeof(omxConstraint));
	ncnln = 0;
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {
		PROTECT(nextVar = VECTOR_ELT(constraints, constraintIndex));
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));
		arg1 = omxMatrixLookupFromState1(nextLoc, globalState);
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));
		arg2 = omxMatrixLookupFromState1(nextLoc, globalState);
		PROTECT(nextLoc = AS_INTEGER(VECTOR_ELT(nextVar, 2)));
		globalState->conList[constraintIndex].opCode = INTEGER(nextLoc)[0];
		omxMatrix *args[2] = {arg1, arg2};
		globalState->conList[constraintIndex].result = omxNewAlgebraFromOperatorAndArgs(10, args, 2, globalState); // 10 = binary subtract
		omxRecompute(globalState->conList[constraintIndex].result);
		int nrows = globalState->conList[constraintIndex].result->rows;
		int ncols = globalState->conList[constraintIndex].result->cols;
		globalState->conList[constraintIndex].size = nrows * ncols;
		ncnln += globalState->conList[constraintIndex].size;
	}
	if(OMX_VERBOSE) { mxLog("Processed."); }
	if(OMX_DEBUG) { mxLog("%d effective constraints.", ncnln); }
	globalState->ncnln = ncnln;
}

/*
*  Acknowledgement:
*  This function is duplicated from the function of the same name in the R source code.
*  The function appears in src/main/sysutils.c
*  Thanks to Michael Spiegel for finding it.
*  This code is licensed under the terms of the GNU General Public License.
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

