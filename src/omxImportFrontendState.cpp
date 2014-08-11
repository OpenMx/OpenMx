/*
 *  Copyright 2007-2014 The OpenMx Project
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

/* Outside R Functions */
static int isDir(const char *path);

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
			setFreeVarGroup(fm->fitFunction, Global->freeGroup[0]);
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

	for (int ax=0; ax < numAlgs; ++ax) {
		omxMatrix *matrix = state->algebraList[ax];
		omxRecompute(matrix, FF_COMPUTE_DIMS, fc);
	}

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}

	// This is required because FF_COMPUTE_DIMS cannot compute
	// dims for all algebra ops.

	for(size_t index = 0; index < numMats; index++) {
		omxRecompute(state->matrixList[index], FF_COMPUTE_INITIAL_FIT, fc);
	}

	// This is required because we have chosen to compute algebras
	// in models without a fitfunction. This is the only place
	// that we loop over _all_ algebras and compute them.

	for(int index = 0; index < numAlgs; index++) {
		omxMatrix *matrix = state->algebraList[index];
		omxRecompute(matrix, FF_COMPUTE_INITIAL_FIT, fc);
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

		const char *pathName, *fileName;

		Rf_protect(nextLoc = VECTOR_ELT(checkpointList, index));
		int next = 0;
		oC->type = (omxCheckpointType) INTEGER(VECTOR_ELT(nextLoc, next++))[0];
		switch(oC->type) {
		case OMX_FILE_CHECKPOINT:{
			pathName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));
			fileName = CHAR(STRING_ELT(VECTOR_ELT(nextLoc, next++), 0));

			if(!isDir(pathName)) {
				Rf_error("Unable to open directory %s for checkpoint storage.\n", pathName);
			}

			std::string fullname = string_snprintf("%s/%s", pathName, fileName);
			if(OMX_DEBUG) { mxLog("Opening File: %s", fullname.c_str()); }
			oC->file = fopen(fullname.c_str(), "w");
			if(!oC->file) {
				Rf_error("Unable to open file %s for checkpoint storage: %s.\n",
					 fullname.c_str(), strerror(errno));
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

	{
		FreeVarGroup *fvg = new FreeVarGroup;
		fvg->id.push_back(FREEVARGROUP_ALL);   // all variables
		Global->freeGroup.push_back(fvg);

		fvg = new FreeVarGroup;
		fvg->id.push_back(FREEVARGROUP_NONE);  // no variables
		Global->freeGroup.push_back(fvg);
	}

	SEXP nextVar, nextLoc;
	int numVars = Rf_length(varList);
	startingValues->resize(numVars);
	for (int fx = 0; fx < numVars; fx++) {
		omxManageProtectInsanity mpi;

		omxFreeVar *fv = new omxFreeVar;
		// default group has free all variables
		Global->freeGroup[FREEVARGROUP_ALL]->vars.push_back(fv);

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
		/*if (sv < fv->lbound) {
			Rf_warning("Moving starting value of parameter '%s' within bounds %f -> %f",
				fv->name, sv, fv->lbound);
			sv = fv->lbound;
		} else if (sv > fv->ubound) {
			Rf_warning("Moving starting value of parameter '%s' within bounds %f -> %f",
				fv->name, sv, fv->ubound);
			sv = fv->ubound;
		}*/
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

	int ncnln = 0; 
	if(OMX_DEBUG) { mxLog("Processing Constraints.");}
	omxMatrix *arg1, *arg2;
	SEXP nextVar, nextLoc;
	numConstraints = Rf_length(constraints);
	if(OMX_DEBUG) {mxLog("Found %d constraints.", numConstraints); }
	conList = (omxConstraint*) R_alloc(numConstraints, sizeof(omxConstraint));
	ncnln = 0;
	for(int ci = 0; ci < numConstraints; ci++) {
		omxConstraint &constr = conList[ci];
		constr.name = CHAR(Rf_asChar(STRING_ELT(names, ci)));
		Rf_protect(nextVar = VECTOR_ELT(constraints, ci));
		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 0));
		arg1 = omxMatrixLookupFromState1(nextLoc, this);
		Rf_protect(nextLoc = VECTOR_ELT(nextVar, 1));
		arg2 = omxMatrixLookupFromState1(nextLoc, this);
		constr.opCode = Rf_asInteger(VECTOR_ELT(nextVar, 2));
		omxMatrix *args[2] = {arg1, arg2};
		constr.result = omxNewAlgebraFromOperatorAndArgs(10, args, 2, this); // 10 = binary subtract
		omxRecompute(constr.result, FF_COMPUTE_DIMS, fc);
		omxRecompute(constr.result, FF_COMPUTE_INITIAL_FIT, fc);
		int nrows = constr.result->rows;
		int ncols = constr.result->cols;
		constr.size = nrows * ncols;
		if (constr.size == 0) {
			Rf_warning("Constraint '%s' evaluated to a 0x0 matrix and will have no effect",
				   constr.name);
		}
		ncnln += constr.size;
	}
	if(OMX_DEBUG) { mxLog("%d effective constraints.", ncnln); }
	this->ncnln = ncnln;
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

// NOTE: All non-linear constraints are applied regardless of free TODO
// variable group.  This is probably wrong. TODO
void omxSetupBoundsAndConstraints(FitContext *fc, double * bl, double * bu)
{
	FreeVarGroup *freeVarGroup = fc->varGroup;
	omxState *globalState = fc->state;
	size_t n = freeVarGroup->vars.size();

	/* Set min and max limits */
	for(size_t index = 0; index < n; index++) {
		bl[index] = freeVarGroup->vars[index]->lbound;
		bu[index] = freeVarGroup->vars[index]->ubound;
	}

	int index = n;
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {		// Nonlinear constraints:
		if(OMX_DEBUG) { mxLog("Constraint %d: ", constraintIndex);}
		switch(globalState->conList[constraintIndex].opCode) {
		case 0:									// Less than: Must be strictly less than 0.
			if(OMX_DEBUG) { mxLog("Bounded at (-INF, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = -0.0;
				index++;
			}
			break;
		case 1:									// Equal: Must be roughly equal to 0.
			if(OMX_DEBUG) { mxLog("Bounded at (-0, 0).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = -0.0;
				bu[index] = 0.0;
				index++;
			}
			break;
		case 2:									// Greater than: Must be strictly greater than 0.
			if(OMX_DEBUG) { mxLog("Bounded at (0, INF).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				if(OMX_DEBUG) { mxLog("\tBounds set for constraint %d.%d.", constraintIndex, offset);}
				bl[index] = 0.0;
				bu[index] = INF;
				index++;
			}
			break;
		default:
			if(OMX_DEBUG) { mxLog("Bounded at (-INF, INF).");}
			for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
				bl[index] = NEG_INF;
				bu[index] = INF;
				index++;
			}
			break;
		}
	}
}

