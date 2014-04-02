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

#include <stdarg.h>
#include <errno.h>

#include "omxState.h"
#include "Compute.h"
#include "omxOpenmpWrap.h"

struct omxGlobal *Global = NULL;

FreeVarGroup *omxGlobal::findVarGroup(int id)
{
	size_t numGroups = Global->freeGroup.size();
	for (size_t vx=0; vx < numGroups; ++vx) {
		std::vector<int> &ids = Global->freeGroup[vx]->id;
		for (size_t ix=0; ix < ids.size(); ++ix) {
			if (ids[ix] == id) return Global->freeGroup[vx];
		}
	}
	return NULL;
}

FreeVarGroup *omxGlobal::findOrCreateVarGroup(int id)
{
	FreeVarGroup *old = findVarGroup(id);
	if (old) return old;

	FreeVarGroup *fvg = new FreeVarGroup;
	fvg->id.push_back(id);
	Global->freeGroup.push_back(fvg);
	return fvg;
}

bool FreeVarGroup::hasSameVars(FreeVarGroup *g2)
{
	if (vars.size() != g2->vars.size()) return false;

	for (size_t vx=0; vx < vars.size(); ++vx) {
		if (vars[vx] != g2->vars[vx]) return false;
	}
	return true;
}

int FreeVarGroup::lookupVar(const char *name)
{
	for (size_t vx=0; vx < vars.size(); ++vx) {
		if (strcmp(name, vars[vx]->name) == 0) return vx;
	}
	return -1;
}

void FreeVarGroup::cacheDependencies()
{
	omxState *os = globalState;
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	dependencies.assign(numMats + numAlgs, false);
	locations.assign(numMats, false);

	for (size_t vx = 0; vx < vars.size(); vx++) {
		omxFreeVar *fv = vars[vx];
		int *deps   = fv->deps;
		int numDeps = fv->numDeps;
		for (int index = 0; index < numDeps; index++) {
			dependencies[deps[index] + numMats] = true;
		}
		for (size_t lx=0; lx < fv->locations.size(); ++lx) {
			locations[fv->locations[lx].matrix] = true;
		}
	}

	// Everything is set up. This is a good place to log.
	if (OMX_DEBUG) { log(); }
}

void FreeVarGroup::markDirty(omxState *os)
{
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	for(size_t i = 0; i < numMats; i++) {
		if (!locations[i]) continue;
		omxMarkClean(os->matrixList[i]);
	}

	for(size_t i = 0; i < numMats; i++) {
		if (dependencies[i]) {
			int offset = ~(i - numMats);
			omxMarkDirty(os->matrixList[offset]);
		}
	}

	for(size_t i = 0; i < numAlgs; i++) {
		if (dependencies[i + numMats]) {
			omxMarkDirty(os->algebraList[i]);
		}
	}
}

void FreeVarGroup::log()
{
	omxState *os = globalState;
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();
	std::string str;

	str += string_snprintf("FreeVarGroup(id=%d", id[0]);
	for (size_t ix=1; ix < id.size(); ++ix) {
		str += string_snprintf(",%d", id[ix]);
	}
	str += string_snprintf(") with %d variables:", (int) vars.size());

	for (size_t vx=0; vx < vars.size(); ++vx) {
		str += " ";
		str += vars[vx]->name;
	}
	if (vars.size()) str += "\nwill dirty:";

	for(size_t i = 0; i < numMats; i++) {
		if (dependencies[i]) {
			int offset = ~(i - numMats);
			str += " ";
			str += os->matrixList[offset]->name;
		}
	}

	for(size_t i = 0; i < numAlgs; i++) {
		if (dependencies[i + numMats]) {
			str += " ";
			str += os->algebraList[i]->name;
		}
	}
	str += "\n";

	mxLogBig(str);
}

omxGlobal::omxGlobal()
{
	ciMaxIterations = 5;
	numThreads = 1;
	analyticGradients = 0;
	numChildren = 0;
	llScale = -2.0;
}

void omxGlobal::deduplicateVarGroups()
{
	for (size_t g1=0; g1 < freeGroup.size(); ++g1) {
		for (size_t g2=freeGroup.size()-1; g2 > g1; --g2) {
			if (freeGroup[g1]->hasSameVars(freeGroup[g2])) {
				freeGroup[g1]->id.insert(freeGroup[g1]->id.end(),
							 freeGroup[g2]->id.begin(), freeGroup[g2]->id.end());
				delete freeGroup[g2];
				freeGroup.erase(freeGroup.begin() + g2);
			}
		}
	}
}

/* Initialize and Destroy */
	void omxInitState(omxState* state) {
		state->stale = FALSE;
		state->numConstraints = 0;
		state->conList = NULL;

		state->majorIteration = 0;
		state->minorIteration = 0;
		state->startTime = 0;
		state->endTime = 0;
		state->numCheckpoints = 0;
		state->checkpointList = NULL;
		state->chkptText1 = NULL;
		state->chkptText2 = NULL;

		state->computeCount = 0;
		state->currentRow = -1;
	}

	void omxSetMajorIteration(omxState *state, int value) {
		state->majorIteration = value;
		if (!state->childList.size()) return;
		for(int i = 0; i < Global->numChildren; i++) {
			omxSetMajorIteration(state->childList[i], value);
		}
	}

	void omxSetMinorIteration(omxState *state, int value) {
		state->minorIteration = value;
		if (!state->childList.size()) return;
		for(int i = 0; i < Global->numChildren; i++) {
			omxSetMinorIteration(state->childList[i], value);
		}
	}
	
	void omxDuplicateState(omxState* tgt, omxState* src) {
		tgt->dataList			= src->dataList;
		
		for(size_t mx = 0; mx < src->matrixList.size(); mx++) {
			// TODO: Smarter inference for which matrices to duplicate
			tgt->matrixList.push_back(omxDuplicateMatrix(src->matrixList[mx], tgt));
		}

		tgt->numConstraints     = src->numConstraints;
		tgt->conList			= (omxConstraint*) R_alloc(tgt->numConstraints, sizeof(omxConstraint));
		for(int j = 0; j < tgt->numConstraints; j++) {
			tgt->conList[j].size   = src->conList[j].size;
			tgt->conList[j].opCode = src->conList[j].opCode;
			tgt->conList[j].lbound = src->conList[j].lbound;
			tgt->conList[j].ubound = src->conList[j].ubound;
			tgt->conList[j].result = omxDuplicateMatrix(src->conList[j].result, tgt);
		}

		for(size_t j = 0; j < src->algebraList.size(); j++) {
			// TODO: Smarter inference for which algebras to duplicate
			tgt->algebraList.push_back(omxDuplicateMatrix(src->algebraList[j], tgt));
		}

		for(size_t j = 0; j < src->expectationList.size(); j++) {
			// TODO: Smarter inference for which expectations to duplicate
			tgt->expectationList.push_back(omxDuplicateExpectation(src->expectationList[j], tgt));
		}

		for(size_t j = 0; j < tgt->algebraList.size(); j++) {
			omxDuplicateAlgebra(tgt->algebraList[j], src->algebraList[j], tgt);
		}

		for(size_t j = 0; j < src->expectationList.size(); j++) {
			// TODO: Smarter inference for which expectations to duplicate
			omxCompleteExpectation(tgt->expectationList[j]);
		}

		tgt->majorIteration 	= 0;
		tgt->minorIteration 	= 0;
		tgt->startTime 			= src->startTime;
		tgt->endTime			= 0;
		
		// TODO: adjust checkpointing based on parallelization method
		tgt->numCheckpoints     = 0;
		tgt->checkpointList 	= NULL;
		tgt->chkptText1 		= NULL;
		tgt->chkptText2 		= NULL;
                                  
		tgt->computeCount 		= src->computeCount;
		tgt->currentRow 		= src->currentRow;
	}

	omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element) {
		if(os == NULL || element == NULL) return NULL;

		if (element->hasMatrixNumber) {
			int matrixNumber = element->matrixNumber;
			if (matrixNumber >= 0) {
				return(os->algebraList[matrixNumber]);
			} else {
				return(os->matrixList[-matrixNumber - 1]);
			}
		}

		omxConstraint* parentConList = globalState->conList;

		for(int i = 0; i < os->numConstraints; i++) {
			if(parentConList[i].result == element) {
				if(os->conList[i].result != NULL) {   // Not sure of proper failure behavior here.
            	return(os->conList[i].result);
				} else {
                    omxRaiseError(os, -2, "Initialization Copy Error: Constraint required but not yet processed.");
            }
			}
		}

		return NULL;
	}
	
void omxFreeChildStates(omxState *state)
{
	if (state->childList.size() == 0) return;

	for(int k = 0; k < Global->numChildren; k++) {
		// Data are not modified and not copied. The same memory
		// is shared across all instances of state. We only need
		// to free the data once, so let the parent do it.
		state->childList[k]->dataList.clear();

		omxFreeState(state->childList[k]);
	}
	state->childList.clear();
	Global->numChildren = 0;
}

	void omxFreeState(omxState *state) {
		omxFreeChildStates(state);

		if(OMX_DEBUG) { mxLog("Freeing %d Constraints.", (int) state->numConstraints);}
		for(int k = 0; k < state->numConstraints; k++) {
			omxFreeMatrix(state->conList[k].result);
		}

		for(size_t ax = 0; ax < state->algebraList.size(); ax++) {
			// free argument tree
			omxFreeMatrix(state->algebraList[ax]);
		}

		for(size_t ax = 0; ax < state->algebraList.size(); ax++) {
			state->algebraList[ax]->hasMatrixNumber = false;
			omxFreeMatrix(state->algebraList[ax]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Matrices.", (int) state->matrixList.size());}
		for(size_t mk = 0; mk < state->matrixList.size(); mk++) {
			state->matrixList[mk]->hasMatrixNumber = false;
			omxFreeMatrix(state->matrixList[mk]);
		}
		
		if(OMX_DEBUG) { mxLog("Freeing %d Model Expectations.", (int) state->expectationList.size());}
		for(size_t ex = 0; ex < state->expectationList.size(); ex++) {
			omxFreeExpectationArgs(state->expectationList[ex]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Data Sets.", (int) state->dataList.size());}
		for(size_t dx = 0; dx < state->dataList.size(); dx++) {
			omxFreeData(state->dataList[dx]);
		}

		if(OMX_DEBUG) { mxLog("Freeing %d Checkpoints.", state->numCheckpoints);}
		for(int k = 0; k < state->numCheckpoints; k++) {
			omxCheckpoint oC = state->checkpointList[k];
			switch(oC.type) {
				case OMX_FILE_CHECKPOINT:
					fclose(oC.file);
					break;
				case OMX_CONNECTION_CHECKPOINT:	// NYI :::DEBUG:::
					// Do nothing: this should be handled by R upon return.
					break;
			}
			if(state->chkptText1 != NULL) {
				Free(state->chkptText1);
			}
			if(state->chkptText2 != NULL) {
				Free(state->chkptText2);
			}
			// Checkpoint list itself is freed by R.
		}

		delete state;

		if(OMX_DEBUG) { mxLog("State Freed.");}
	}

omxGlobal::~omxGlobal()
{
	for (size_t cx=0; cx < computeList.size(); ++cx) {
		delete computeList[cx];
	}
	for (size_t cx=0; cx < algebraList.size(); ++cx) {
		delete algebraList[cx];
	}
	if (freeGroup.size()) {
		std::vector< omxFreeVar* > &vars = freeGroup[0]->vars;  // has all vars
		for (size_t vx=0; vx < vars.size(); ++vx) {
			delete vars[vx];
		}
	}
	for (size_t gx=0; gx < freeGroup.size(); ++gx) {
		delete freeGroup[gx];
	}
}

void omxResetStatus(omxState *)
{
	Global->bads.clear();
}

std::string string_vsnprintf(const char *fmt, va_list ap)
{
    int size = 100;
    std::string str;
    while (1) {
        str.resize(size);
        int n = vsnprintf((char *)str.c_str(), size, fmt, ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

std::string string_snprintf(const char *fmt, ...)
{
	va_list ap;
        va_start(ap, fmt);
	std::string str = string_vsnprintf(fmt, ap);
        va_end(ap);
	return str;
}

void mxLogBig(const std::string str)   // thread-safe
{
	ssize_t len = ssize_t(str.size());
	ssize_t wrote = 0;
	int maxRetries = 20;
	ssize_t got;
#pragma omp critical(stderp)
	{
		while (--maxRetries > 0) {
			got = write(2, str.data() + wrote, len - wrote);
			if (got == -EINTR) continue;
			if (got <= 0) break;
			wrote += got;
			if (wrote == len) break;
		}
	}
	if (got <= 0) Rf_error("mxLogBig failed with errno=%d", got);

}

void mxLog(const char* msg, ...)   // thread-safe
{
	const int maxLen = 240;
	char buf1[maxLen];
	char buf2[maxLen];

	va_list ap;
	va_start(ap, msg);
	vsnprintf(buf1, maxLen, msg, ap);
	va_end(ap);

	int len = snprintf(buf2, maxLen, "[%d] %s\n", omx_absolute_thread_num(), buf1);

	int maxRetries = 20;
	ssize_t wrote = 0;
	ssize_t got;
#pragma omp critical(stderp)
	{
		while (--maxRetries > 0) {
			got = write(2, buf2 + wrote, len - wrote);
			if (got == -EINTR) continue;
			if (got <= 0) break;
			wrote += got;
			if (wrote == len) break;
		}
	}
	if (got <= 0) Rf_error("mxLog failed with errno=%d", got);
}

void _omxRaiseError()
{
	// keep for debugger breakpoints
}

void omxRaiseErrorf(omxState *, const char* msg, ...)
{
	va_list ap;
	va_start(ap, msg);
	std::string str = string_vsnprintf(msg, ap);
	va_end(ap);
	_omxRaiseError();

	if(OMX_DEBUG) {
		mxLog("Error raised: %s", str.c_str());
	}

#pragma omp critical(bads)
	Global->bads.push_back(str);
}

const char *omxGlobal::getBads()
{
	if (bads.size() == 0) return NULL;

	std::string str;
	for (size_t mx=0; mx < bads.size(); ++mx) {
		if (bads.size() > 1) str += string_snprintf("%d:", (int)mx+1);
		str += bads[mx];
		if (mx < bads.size() - 1) str += "\n";
	}

	size_t sz = str.size();
	char *mem = R_alloc(sz+1, 1);  // use R's memory
	memcpy(mem, str.c_str(), sz);
	mem[sz] = 0;
	return mem;
}

void omxRaiseError(omxState *, int, const char* msg) { // DEPRECATED
	omxRaiseErrorf(NULL, "%s", msg);
}

	void omxStateNextRow(omxState *state) {
		state->currentRow++;
	};

	void omxStateNextEvaluation(omxState *state) {
		state->currentRow = -1;
		state->computeCount++;
	};

static void omxWriteCheckpointHeader(omxState *os, omxCheckpoint* oC) {
	// rewrite with std::string TODO
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		os->chkptText1 = (char*) Calloc((24 + 15 * numParam), char);
		os->chkptText2 = (char*) Calloc(1.0 + 15.0 * numParam*
			(numParam + 1.0) / 2.0, char);
		if (oC->type == OMX_FILE_CHECKPOINT) {
			fprintf(oC->file, "iterations\ttimestamp\tobjective\t");
			for(size_t j = 0; j < numParam; j++) {
				if(strcmp(vars[j]->name, CHAR(NA_STRING)) == 0) {
					fprintf(oC->file, "%s", vars[j]->name);
				} else {
					fprintf(oC->file, "\"%s\"", vars[j]->name);
				}
				if (j != numParam - 1) fprintf(oC->file, "\t");
			}
			fprintf(oC->file, "\n");
			fflush(oC->file);
		}
	}
 
void omxWriteCheckpointMessage(char *msg) {
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		omxState *os = globalState;
		for(int i = 0; i < os->numCheckpoints; i++) {
			omxCheckpoint* oC = &(os->checkpointList[i]);
			if(os->chkptText1 == NULL) {    // First one: set up output
				omxWriteCheckpointHeader(os, oC);
			}
			if (oC->type == OMX_FILE_CHECKPOINT) {
				fprintf(oC->file, "%d \"%s\" NA ", os->majorIteration, msg);
				for(size_t j = 0; j < numParam; j++) {
					fprintf(oC->file, "NA ");
				}
				fprintf(oC->file, "\n");
			}
		}
	}

/*
 * Next time we rewrite this code, Mike Neale suggested that the
 * checkpoint be taken before evaluating the fit function and then
 * again after obtaining the fit statistic. This will help with
 * debugging fit functions that fail to evaluate in some spots.
 */
void omxSaveCheckpoint(double* x, double f, int force) {
	// rewrite with std::string TODO
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

		omxState *os = globalState;
		time_t now = time(NULL);
		int soFar = now - os->startTime;		// Translated into minutes
		int n;
		for(int i = 0; i < os->numCheckpoints; i++) {
			n = 0;
			omxCheckpoint* oC = &(os->checkpointList[i]);
			// Check based on time            
			if((oC->time > 0 && (soFar - oC->lastCheckpoint) >= oC->time) || force) {
				oC->lastCheckpoint = soFar;
				n = 1;
			}
			// Or iterations
			if((oC->numIterations > 0 && (os->majorIteration - oC->lastCheckpoint) >= oC->numIterations) || force) {
				oC->lastCheckpoint = os->majorIteration;
				n = 1;
			}

			if(n) {		//In either case, save a checkpoint.
				if(os->chkptText1 == NULL) {	// First one: set up output
					omxWriteCheckpointHeader(os, oC);
				}
				char tempstring[25];
				sprintf(tempstring, "%d", os->majorIteration);

				if(strncmp(os->chkptText1, tempstring, strlen(tempstring))) {	// Returns zero if they're the same.
					struct tm * nowTime = localtime(&now);						// So this only happens if the text is out of date.
					strftime(tempstring, 25, "%b %d %Y %I:%M:%S %p", nowTime);
					sprintf(os->chkptText1, "%d \"%s\" %9.5f", os->majorIteration, tempstring, f);
					for(size_t j = 0; j < numParam; j++) {
						sprintf(tempstring, " %9.5f", x[j]);
						strncat(os->chkptText1, tempstring, 14);
					}
				}

				if(oC->type == OMX_FILE_CHECKPOINT) {
					fprintf(oC->file, "%s", os->chkptText1);
					if(oC->saveHessian)
						fprintf(oC->file, "%s", os->chkptText2);
					fprintf(oC->file, "\n");
					fflush(oC->file);
				} else if(oC->type == OMX_CONNECTION_CHECKPOINT) {
					Rf_warning("NYI: R_connections are not yet implemented.");
					oC->numIterations = 0;
					oC->time = 0;
				}
			}
		}
	}

void omxExamineFitOutput(omxState *state, omxMatrix *fitMatrix, int *mode)
{
	if (!R_FINITE(fitMatrix->data[0])) {
		omxRaiseErrorf(state, "Fit function returned %g", fitMatrix->data[0]);
		*mode = -1;
	}
}

omxFreeVarLocation *omxFreeVar::getLocation(int matrix)
{
	for (size_t lx=0; lx < locations.size(); lx++) {
		omxFreeVarLocation *loc = &locations[lx];
		if (~loc->matrix == matrix) return loc;
	}
	return NULL;
}
