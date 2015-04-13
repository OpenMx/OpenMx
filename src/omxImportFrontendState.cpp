/*
 *  Copyright 2007-2015 The OpenMx Project
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

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <sys/stat.h>
#include <errno.h>

#include "omxDefines.h"
#include "glue.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "Compute.h"

int matchCaseInsensitive(const char *source, const char *target) {
	return strcasecmp(source, target) == 0;
}

void omxState::omxProcessMxDataEntities(SEXP data)
{
	SEXP nextLoc;
	if(OMX_DEBUG) { mxLog("Processing %d data source(s).", Rf_length(data));}

	SEXP listNames = Rf_getAttrib(data, R_NamesSymbol);

	for(int index = 0; index < Rf_length(data); index++) {
		ScopedProtect p1(nextLoc, VECTOR_ELT(data, index));			// Retrieve the data object
		omxNewDataFromMxData(nextLoc, CHAR(STRING_ELT(listNames, index)));
	}
}

void omxState::omxProcessMxMatrixEntities(SEXP matList)
{
	if(OMX_DEBUG) { mxLog("Processing %d matrix(ces).", Rf_length(matList));}
	SEXP nextLoc, nextMat;
	matrixList.clear();
	SEXP matListNames = Rf_getAttrib(matList, R_NamesSymbol);

	for(int index = 0; index < Rf_length(matList); index++) {
		omxManageProtectInsanity protectManager;
		Rf_protect(nextLoc = VECTOR_ELT(matList, index));		// This is the matrix + populations
		Rf_protect(nextMat = VECTOR_ELT(nextLoc, 0));		// The first element of the list is the matrix of values
		omxMatrix *mat = omxNewMatrixFromRPrimitive(nextMat, this, 1, -index - 1);
		mat->name = CHAR(STRING_ELT(matListNames, index));
		matrixList.push_back(mat);

		if(OMX_DEBUG) { omxPrintMatrix(mat, "Imported"); }

		if (isErrorRaised()) return;
	}
}

void omxState::omxProcessMxAlgebraEntities(SEXP algList)
{
	SEXP nextAlgTuple;
	SEXP algListNames = Rf_getAttrib(algList, R_NamesSymbol);

	if(OMX_DEBUG) { mxLog("Processing %d algebras.", Rf_length(algList)); }

	for(int index = 0; index < Rf_length(algList); index++) {
		algebraList.push_back(omxInitMatrix(0, 0, TRUE, this));
	}

	for(int index = 0; index < Rf_length(algList); index++) {
		omxManageProtectInsanity protectManager;
		Rf_protect(nextAlgTuple = VECTOR_ELT(algList, index));
		if(IS_S4_OBJECT(nextAlgTuple)) {
			SEXP fitFunctionClass;
			ScopedProtect p1(fitFunctionClass, STRING_ELT(Rf_getAttrib(nextAlgTuple, Rf_install("class")), 0));
			const char *fitType = CHAR(fitFunctionClass);
			omxMatrix *fm = algebraList[index];
			omxFillMatrixFromMxFitFunction(fm, fitType, index, nextAlgTuple);
			fm->name = CHAR(STRING_ELT(algListNames, index));
		} else {								// This is an algebra spec.
			SEXP dimnames, formula;
			omxMatrix *amat = algebraList[index];
			Rf_protect(dimnames = VECTOR_ELT(nextAlgTuple, 0));
			omxFillMatrixFromRPrimitive(amat, NULL, this, 1, index);
			Rf_protect(formula = VECTOR_ELT(nextAlgTuple, 1));
			omxFillMatrixFromMxAlgebra(amat, formula, CHAR(STRING_ELT(algListNames, index)), dimnames);
		}
		if (isErrorRaised()) return;
	}
}

void omxState::omxCompleteMxFitFunction(SEXP algList)
{
	SEXP nextAlgTuple;

	for(int index = 0; index < Rf_length(algList); index++) {
		bool s4;
		{
			ScopedProtect p1(nextAlgTuple, VECTOR_ELT(algList, index));
			s4 = IS_S4_OBJECT(nextAlgTuple);
		}
		if(!s4) continue;
		omxMatrix *fm = algebraList[index];
		if (!fm->fitFunction->freeVarGroup) {
			setFreeVarGroup(fm->fitFunction, Global->findVarGroup(FREEVARGROUP_ALL));
		}
		omxCompleteFitFunction(fm);
	}
}

void omxState::omxProcessMxExpectationEntities(SEXP expList)
{
	if(OMX_DEBUG) { mxLog("Initializing %d Model Expectation(s).", Rf_length(expList));}
	SEXP nextExp;
	SEXP eNames = Rf_getAttrib(expList, R_NamesSymbol);

	for(int index = 0; index < Rf_length(expList); index++) {
		if (isErrorRaised()) return;
		Rf_protect(nextExp = VECTOR_ELT(expList, index));
		omxExpectation *ex = omxNewIncompleteExpectation(nextExp, index, this);
		ex->name = CHAR(STRING_ELT(eNames, index));
		expectationList.push_back(ex);
	}
}


void omxState::omxCompleteMxExpectationEntities()
{
	if(OMX_DEBUG) { mxLog("Completing %d Model Expectation(s).", (int) expectationList.size());}
	
	for(size_t index = 0; index < expectationList.size(); index++) {
		if (isErrorRaised()) return;
		omxCompleteExpectation(expectationList[index]);
	}
}

void omxGlobal::omxProcessMxComputeEntities(SEXP rObj, omxState *currentState)
{
	if (Rf_isNull(rObj)) return;

	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(rObj, Rf_install("class")), 0));
	omxCompute *compute = omxNewCompute(currentState, CHAR(s4class));
	compute->initFromFrontend(currentState, rObj);
	computeList.push_back(compute);
}

void omxInitialMatrixAlgebraCompute(omxState *state, FitContext *fc)
{
	// We use FF_COMPUTE_INITIAL_FIT because an expectation
	// could depend on the value of an algebra. However, we
	// don't mark anything clean because an algebra could
	// depend on an expectation (via a fit function).

	size_t numMats = state->matrixList.size();
	int numAlgs = state->algebraList.size();

	state->setWantStage(FF_COMPUTE_DIMS);

	for (int ax=0; ax < numAlgs; ++ax) {
		omxMatrix *matrix = state->algebraList[ax];
		omxRecompute(matrix, fc);
	}

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}

	state->setWantStage(FF_COMPUTE_INITIAL_FIT);

	// Need something finite for definition variables to avoid exceptions

	for (int ex = 0; ex < (int) state->expectationList.size(); ++ex) {
		state->expectationList[ex]->loadFakeData(1);
	}

	// This is required because FF_COMPUTE_DIMS cannot compute
	// dims for all algebra ops.

	for(size_t index = 0; index < numMats; index++) {
		omxRecompute(state->matrixList[index], fc);
	}

	// This is required because we have chosen to compute algebras
	// in models without a fitfunction. This is the only place
	// that we loop over _all_ algebras and compute them.

	for(int index = 0; index < numAlgs; index++) {
		omxMatrix *matrix = state->algebraList[index];
		omxRecompute(matrix, fc);
	}
}

/*
checkpointList is a list().  Each element refers to one checkpointing request.
Each interval request is a list of Rf_length 5.
The first element is an integer that specifies type: 0 = file, 1 = socket, 2=R_connection
For a file, the next two are the directory(string)  and file name (string).
For a socket, they are server (string) and port number (int).
For a connection, the next one is the R_connection SEXP object.
After that is an integer <type> specifier.  0 means minutes, 1 means iterations.
The last element is an integer count, indicating the number of <type>s per checkpoint.
*/
void omxProcessCheckpointOptions(SEXP checkpointList)
{
	if(OMX_DEBUG) {mxLog("Found %d checkpoints", Rf_length(checkpointList)); }

	SEXP nextLoc;

	for(int index = 0; index < Rf_length(checkpointList); ++index) {
		omxCheckpoint *oC = new omxCheckpoint;

		Rf_protect(nextLoc = VECTOR_ELT(checkpointList, index));
		int next = 0;
		oC->type = (omxCheckpointType) INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		switch(oC->type) {
		case OMX_FILE_CHECKPOINT:{
			const char *fullname = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));

			if(OMX_DEBUG) { mxLog("Opening File: %s", fullname); }
			oC->file = fopen(fullname, "w");
			if(!oC->file) {
				Rf_error("Unable to open file %s for checkpoint storage: %s",
					 fullname, strerror(errno));
			}
			break;}

		case OMX_CONNECTION_CHECKPOINT:{
			Rf_error("Warning NYI: Socket checkpoints Not Yet Implemented.\n");
			break;}
		}

		const char *units = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));
		if (strEQ(units, "iterations")) {
			oC->iterPerCheckpoint = Rf_asInteger(VECTOR_ELT(nextLoc, next++));
		} else if (strEQ(units, "minutes")) {
			oC->timePerCheckpoint = Rf_asReal(VECTOR_ELT(nextLoc, next++)) * 60.0; // Constrained to seconds.
			if(oC->timePerCheckpoint < 1) oC->timePerCheckpoint = 1;
		} else if (strEQ(units, "evaluations")) {
			oC->evalsPerCheckpoint = Rf_asInteger(VECTOR_ELT(nextLoc, next++));
		} else {
			Rf_error("In 'Checkpoint Units' model option, '%s' not recognized", units);
		}
		Global->checkpointList.push_back(oC);
	}
}

void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues)
{
	if(OMX_DEBUG) { mxLog("Processing Free Parameters."); }

	SEXP nextVar, nextLoc;
	int numVars = Rf_length(varList);
	startingValues->resize(numVars);
	for (int fx = 0; fx < numVars; fx++) {
		omxManageProtectInsanity mpi;

		omxFreeVar *fv = new omxFreeVar;
		// default group has free all variables
		Global->findVarGroup(FREEVARGROUP_ALL)->vars.push_back(fv);

		fv->id = fx;
		fv->name = CHAR(Rf_asChar(STRING_ELT(Rf_getAttrib(varList, R_NamesSymbol), fx)));
		Rf_protect(nextVar = VECTOR_ELT(varList, fx));

		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 0));
		fv->lbound = REAL(nextLoc)[0];
		if (ISNA(fv->lbound)) fv->lbound = NEG_INF;
		if (fv->lbound == 0.0) fv->lbound = 0.0;

		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 1));
		fv->ubound = REAL(nextLoc)[0];
		if (ISNA(fv->ubound)) fv->ubound = INF;
		if (fv->ubound == 0.0) fv->ubound = -0.0;

		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 2));
		int groupCount = Rf_length(nextLoc);
		for (int gx=0; gx < groupCount; ++gx) {
			int group = INTEGER(nextLoc)[gx];
			if (group == 0) continue;
			Global->findOrCreateVarGroup(group)->vars.push_back(fv);
		}

		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 3));
		int numDeps = LENGTH(nextLoc);
		fv->numDeps = numDeps;
		fv->deps = (int*) R_alloc(numDeps, sizeof(int));
		for (int i = 0; i < numDeps; i++) {
			fv->deps[i] = INTEGER(nextLoc)[i];
		}

		int numLocs = Rf_length(nextVar) - 5;
		for(int locIndex = 0; locIndex < numLocs; locIndex++) {
			Rf_protect(nextLoc = VECTOR_ELT(nextVar, locIndex+4));
			int* theVarList = INTEGER(nextLoc);

			omxFreeVarLocation loc;
			loc.matrix = theVarList[0];
			loc.row = theVarList[1];
			loc.col = theVarList[2];

			fv->locations.push_back(loc);
		}
		Rf_protect(nextLoc = VECTOR_ELT(nextVar, Rf_length(nextVar)-1));
		double sv = REAL(nextLoc)[0];
		(*startingValues)[fx] = sv;
		if(OMX_DEBUG) {
			mxLog("Free var %d %s, bounds (%.3g, %.3g), %d loc, starting %f", fx, fv->name,
			      fv->lbound, fv->ubound, numLocs, sv);
		}
	}

	Global->deduplicateVarGroups();
}

/*
	intervalList is a list().  Each element refers to one confidence interval request.
	Each interval request is a Rf_length 5 vector of REAL.
	The first three elements are the matrixPointer, Row, and Column of the element
	for which bounds are to be calculated, and are cast to ints here for speed.
	The last two are the upper and lower boundaries for the confidence space (respectively).
*/
void omxGlobal::omxProcessConfidenceIntervals(SEXP intervalList, omxState *currentState)
{
	SEXP names = Rf_getAttrib(intervalList, R_NamesSymbol);
	SEXP nextVar;
	int numIntervals = Rf_length(intervalList);
	if(OMX_DEBUG) {mxLog("Found %d Confidence Interval requests.", numIntervals); }
	Global->intervalList.reserve(numIntervals);
	for(int index = 0; index < numIntervals; index++) {
		omxConfidenceInterval *oCI = new omxConfidenceInterval;
		Rf_protect(nextVar = VECTOR_ELT(intervalList, index));
		double* intervalInfo = REAL(nextVar);
		oCI->name = CHAR(Rf_asChar(STRING_ELT(names, index)));
		oCI->matrix = omxMatrixLookupFromState1( nextVar, currentState);	// Expects an R object
		oCI->row = (int) intervalInfo[1];		// Cast to int in C to save memory/Protection ops
		oCI->col = (int) intervalInfo[2];		// Cast to int in C to save memory/Protection ops
		oCI->lbound = intervalInfo[3];
		oCI->ubound = intervalInfo[4];
		oCI->max = R_NaReal;					// NAs, in case something goes wrong
		oCI->min = R_NaReal;
		Global->intervalList.push_back(oCI);
	}
}

void omxState::omxProcessConstraints(SEXP constraints, FitContext *fc)
{
	SEXP names = Rf_getAttrib(constraints, R_NamesSymbol);

	if(OMX_DEBUG) { mxLog("Processing Constraints.");}
	SEXP nextVar, nextLoc;
	int numConstraints = Rf_length(constraints);
	if(OMX_DEBUG) {mxLog("Found %d constraints.", numConstraints); }
	conList.reserve(numConstraints + 1);  // reserve 1 extra for confidence intervals
	for(int ci = 0; ci < numConstraints; ci++) {
		Rf_protect(nextVar = VECTOR_ELT(constraints, ci));
		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 0));
		omxMatrix *arg1 = omxMatrixLookupFromState1(nextLoc, this);
		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 1));
		omxMatrix *arg2 = omxMatrixLookupFromState1(nextLoc, this);
		omxConstraint *constr = new omxConstraint(fc, CHAR(Rf_asChar(STRING_ELT(names, ci))), arg1, arg2);
		constr->opCode = (omxConstraint::Type) Rf_asInteger(VECTOR_ELT(nextVar, 2));
		conList.push_back(constr);
	}
	if(OMX_DEBUG) {
		int equality, inequality;
		countNonlinearConstraints(equality, inequality);
		mxLog("Found %d equality and %d inequality constraints", equality, inequality);
	}
}
