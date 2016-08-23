/*
 *  Copyright 2007-2016 The OpenMx Project
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

#include <sys/stat.h>
#include <errno.h>

#include "omxDefines.h"
#include "glue.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "Compute.h"

#pragma GCC diagnostic warning "-Wshadow"

int matchCaseInsensitive(const char *source, const char *target) {
	return strcasecmp(source, target) == 0;
}

void omxState::omxProcessMxDataEntities(SEXP data, SEXP defvars)
{
	SEXP nextLoc;
	if(OMX_DEBUG) { mxLog("Processing %d data source(s).", Rf_length(data));}

	SEXP listNames = Rf_getAttrib(data, R_NamesSymbol);

	for(int index = 0; index < Rf_length(data); index++) {
		ScopedProtect p1(nextLoc, VECTOR_ELT(data, index));			// Retrieve the data object
		omxNewDataFromMxData(nextLoc, CHAR(STRING_ELT(listNames, index)));
	}

	int numDefs = Rf_length(defvars);
	for(int nextDef = 0; nextDef < numDefs; nextDef++) {
		omxDefinitionVar dvar;
		
		SEXP itemList;
		ScopedProtect p1(itemList, VECTOR_ELT(defvars, nextDef));
		int *ilist = INTEGER(itemList);
		omxData *od = dataList[ ilist[0] ];
		dvar.column = ilist[1];
		dvar.matrix = ilist[2];
		dvar.row = ilist[3];
		dvar.col = ilist[4];

		int numDeps = Rf_length(itemList) - 5;
		dvar.numDeps = numDeps;
		dvar.deps = (int*) R_alloc(numDeps, sizeof(int));
		for(int i = 0; i < numDeps; i++) {
			dvar.deps[i] = ilist[5+i];
		}

		od->defVars.push_back(dvar);
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
		mat->nameStr = CHAR(STRING_ELT(matListNames, index));
		matrixList.push_back(mat);

		if(OMX_DEBUG) { omxPrintMatrix(mat, NULL); }

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
			ScopedProtect p1(fitFunctionClass, STRING_ELT(Rf_getAttrib(nextAlgTuple, R_ClassSymbol), 0));
			const char *fitType = CHAR(fitFunctionClass);
			omxMatrix *fm = algebraList[index];
			omxFillMatrixFromMxFitFunction(fm, fitType, index, nextAlgTuple);
			fm->nameStr = CHAR(STRING_ELT(algListNames, index));
		} else {								// This is an algebra spec.
			SEXP dimnames, formula;
			omxMatrix *amat = algebraList[index];
			Rf_protect(dimnames = VECTOR_ELT(nextAlgTuple, 0));
			omxFillMatrixFromRPrimitive(amat, NULL, this, 1, index);
			amat->setJoinInfo(VECTOR_ELT(nextAlgTuple, 1), VECTOR_ELT(nextAlgTuple, 2));
			Rf_protect(formula = VECTOR_ELT(nextAlgTuple, 3));
			std::string name = CHAR(STRING_ELT(algListNames, index));
			omxFillMatrixFromMxAlgebra(amat, formula, name, dimnames);
		}
		if (isErrorRaised()) return;
	}
}

void omxState::omxCompleteMxFitFunction(SEXP algList, FitContext *fc)
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
		omxFitFunctionComputeAuto(fm->fitFunction, FF_COMPUTE_INITIAL_FIT, fc);
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
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
	omxCompute *compute = omxNewCompute(currentState, CHAR(s4class));
	compute->initFromFrontend(currentState, rObj);
	computeList.push_back(compute);
}

// This is called at initialization and when we copy
// state to facilitate multiple threads. When we copy
// state, it would be more efficient to deeply copy the
// state instead of recomputing it. However, it is less
// error prone from a programming standpoint to recompute.
void omxState::omxInitialMatrixAlgebraCompute(FitContext *fc)
{
	// Note: some overlap with cacheDependencies
	FreeVarGroup *fvg = Global->findVarGroup(FREEVARGROUP_ALL);
	for (size_t vx=0; vx < fvg->vars.size(); ++vx) {
		omxFreeVar &fv = *fvg->vars[vx];
		for (size_t lx = 0; lx < fv.locations.size(); ++lx) {
			omxMatrix *mat = matrixList[ fv.locations[lx].matrix ];
			if (OMX_DEBUG && !mat->dependsOnParameters()) {
				mxLog("Matrix '%s' depends on free parameter '%s'", mat->name(), fv.name);
			}
			mat->setDependsOnParameters();
		}
	}

	for (size_t dx=0; dx < dataList.size(); ++dx) {
		std::vector<omxDefinitionVar> &defVars = dataList[dx]->defVars;
		for (size_t vx=0; vx < defVars.size(); ++vx) {
			omxDefinitionVar &dv = defVars[vx];
			int matrixNumber = dv.matrix;
			omxMatrix *matrix = matrixList[matrixNumber];
			matrix->setDependsOnDefinitionVariables();
		}
	}

	// We use FF_COMPUTE_INITIAL_FIT because an expectation
	// could depend on the value of an algebra. However, we
	// don't mark anything clean because an algebra could
	// depend on an expectation (via a fit function).

	size_t numMats = matrixList.size();
	int numAlgs = algebraList.size();

	if(OMX_DEBUG) mxLog("omxInitialMatrixAlgebraCompute(state[%d], ...)", getId());

	setWantStage(FF_COMPUTE_INITIAL_FIT);

	// Need something finite for definition variables to avoid exceptions

	for (int ex = 0; ex < (int) dataList.size(); ++ex) {
		// It is necessary to load some number (like 1) instead
		// of NAs because algebra can use definition variables
		// for indexing. We will load real data later.
		dataList[ex]->loadFakeData(this, 1);
	}

	for(size_t index = 0; index < numMats; index++) {
		omxRecompute(matrixList[index], fc);
	}

	// This is required because we have chosen to compute algebras
	// in models without a fitfunction. This is the only place
	// that we loop over _all_ algebras and compute them.

	for(int index = 0; index < numAlgs; index++) {
		omxMatrix *matrix = algebraList[index];
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

void omxState::omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues)
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

ConfidenceInterval::ConfidenceInterval()
{
	row = -1;
	col = -1;
	varIndex = -1;
	bound.setConstant(R_NaReal);
	val.setConstant(R_NaReal);
	code.setConstant(INFORM_UNINITIALIZED);
}

/*
	intervalList is a list().  Each element refers to one confidence interval request.
	Each interval request is a Rf_length 5 vector of REAL.
	The first three elements are the matrixPointer, Row, and Column of the element
	for which bounds are to be calculated, and are cast to ints here for speed.
	The last two are the upper and lower boundaries for the confidence space (respectively).
*/
void omxGlobal::omxProcessConfidenceIntervals(SEXP iList, omxState *currentState)
{
	SEXP names = Rf_getAttrib(iList, R_NamesSymbol);
	SEXP nextVar;
	int numIntervals = Rf_length(iList);
	if(OMX_DEBUG) {mxLog("Found %d Confidence Interval requests.", numIntervals); }
	Global->intervalList.reserve(numIntervals);
	for(int index = 0; index < numIntervals; index++) {
		ConfidenceInterval *oCI = new ConfidenceInterval;
		Rf_protect(nextVar = VECTOR_ELT(iList, index));
		double* intervalInfo = REAL(nextVar);
		oCI->name = CHAR(Rf_asChar(STRING_ELT(names, index)));
		oCI->matrixNumber = Rf_asInteger(nextVar);
		oCI->row = (int) intervalInfo[1];		// Cast to int in C to save memory/Protection ops
		oCI->col = (int) intervalInfo[2];		// Cast to int in C to save memory/Protection ops
		oCI->bound[ConfidenceInterval::Lower] = intervalInfo[3];
		oCI->bound[ConfidenceInterval::Upper] = intervalInfo[4];
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
		omxConstraint *constr = new UserConstraint(fc, CHAR(Rf_asChar(STRING_ELT(names, ci))), arg1, arg2);
		constr->opCode = (omxConstraint::Type) Rf_asInteger(VECTOR_ELT(nextVar, 2));
		if (OMX_DEBUG) mxLog("constraint '%s' is type %d", constr->name, constr->opCode);
		conList.push_back(constr);
	}
	if(OMX_DEBUG) {
		int equality, inequality;
		countNonlinearConstraints(equality, inequality);
		mxLog("Found %d equality and %d inequality constraints", equality, inequality);
	}
}
