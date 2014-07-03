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

	for(size_t i = 0; i < numMats; i++) {
		if (!locations[i]) continue;
		os->matrixList[i]->setNotConstant();
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
	llScale = -2.0;
	computeCount = 0;
	anonAlgebra = 0;
	rowLikelihoodsWarning = false;
	unpackedConfidenceIntervals = false;
}

void omxGlobal::unpackConfidenceIntervals()
{
	if (unpackedConfidenceIntervals) return;
	unpackedConfidenceIntervals = true;

	// take care to preserve order
	std::vector<omxConfidenceInterval*> tmp;
	std::swap(tmp, intervalList);

	for (int ix=0; ix < (int) tmp.size(); ++ix) {
		omxConfidenceInterval *ci = tmp[ix];
		if (!ci->isWholeAlgebra()) {
			intervalList.push_back(ci);
			continue;
		}
		omxMatrix *mat = ci->matrix;
		for (int cx=0; cx < mat->cols; ++cx) {
			for (int rx=0; rx < mat->rows; ++rx) {
				omxConfidenceInterval *cell = new omxConfidenceInterval(*ci);
				std::string name = string_snprintf("%s[%d,%d]", ci->name, 1+rx, 1+cx);
				cell->name = CHAR(Rf_mkChar(name.c_str()));
				cell->row = rx;
				cell->col = cx;
				intervalList.push_back(cell);
			}
		}
		delete ci;
	}
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
		state->currentRow = -1;
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
                    omxRaiseError("Initialization Copy Error: Constraint required but not yet processed.");
            }
			}
		}

		return NULL;
	}
	
void omxFreeChildStates(omxState *state)
{
	if (state->childList.size() == 0) return;

	for(size_t k = 0; k < state->childList.size(); k++) {
		// Data are not modified and not copied. The same memory
		// is shared across all instances of state. We only need
		// to free the data once, so let the parent do it.
		state->childList[k]->dataList.clear();

		omxFreeState(state->childList[k]);
	}
	state->childList.clear();
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

		delete state;

		if(OMX_DEBUG) { mxLog("State Freed.");}
	}

omxGlobal::~omxGlobal()
{
	for (size_t cx=0; cx < intervalList.size(); ++cx) {
		delete intervalList[cx];
	}
	for (size_t cx=0; cx < computeList.size(); ++cx) {
		delete computeList[cx];
	}
	for (size_t cx=0; cx < algebraList.size(); ++cx) {
		delete algebraList[cx];
	}
	for (size_t cx=0; cx < checkpointList.size(); ++cx) {
		delete checkpointList[cx];
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

std::string string_vsnprintf(const char *fmt, va_list orig_ap)
{
    int size = 100;
    std::string str;
    while (1) {
        str.resize(size);
	va_list ap;
	va_copy(ap, orig_ap);
        int n = vsnprintf((char *)str.c_str(), size, fmt, ap);
	va_end(ap);
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
	if (len == 0) Rf_error("Attempt to log 0 characters with mxLogBig");
	ssize_t wrote = 0;
	int maxRetries = 20;
	ssize_t got;
	const char *outBuf = str.c_str();
#pragma omp critical(stderp)
	{
		while (--maxRetries > 0) {
			got = write(2, outBuf + wrote, len - wrote);
			if (got == -EINTR) continue;
			if (got < 0) break;
			wrote += got;
			if (wrote == len) break;
		}
	}
	if (wrote != len) Rf_error("mxLogBig only wrote %d/%d, errno %d", wrote, len, errno);

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

void omxRaiseErrorf(const char* msg, ...)
{
	va_list ap;
	va_start(ap, msg);
	std::string str = string_vsnprintf(msg, ap);
	va_end(ap);
	_omxRaiseError();

	if(OMX_DEBUG) {
		mxLog("Error raised: %s", str.c_str());
	}

	bool overflow = false;
#pragma omp critical(bads)
        {
		if (Global->bads.size() > 100) {
			overflow = true;
		} else {
			Global->bads.push_back(str);
		}
	}

        // mxLog takes a lock too, so call it outside of critical section
        if (overflow) mxLog("Too many errors: %s", str.c_str());
}

const char *omxGlobal::getBads()
{
	if (bads.size() == 0) return NULL;

	std::string str;
	for (size_t mx=0; mx < bads.size(); ++mx) {
		if (bads.size() > 1) str += string_snprintf("%d:", (int)mx+1);
		str += bads[mx];
		if (str.size() > (1<<14)) break;
		if (mx < bads.size() - 1) str += "\n";
	}

	size_t sz = str.size();
	char *mem = R_alloc(sz+1, 1);  // use R's memory
	memcpy(mem, str.c_str(), sz);
	mem[sz] = 0;
	return mem;
}

void omxRaiseError(const char* msg) { // DEPRECATED
	omxRaiseErrorf("%s", msg);
}

	void omxStateNextRow(omxState *state) {
		state->currentRow++;
	};

void omxGlobal::checkpointMessage(FitContext *fc, double *est, const char *fmt, ...)
{
	va_list ap;
        va_start(ap, fmt);
	std::string str = string_vsnprintf(fmt, ap);
        va_end(ap);

	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->message(fc, est, str.c_str());
	}
}

void omxGlobal::checkpointPrefit(FitContext *fc, double *est, bool force)
{
	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->prefit(fc, est, force);
	}
}

void omxGlobal::checkpointPostfit(FitContext *fc)
{
	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->postfit(fc);
	}
}

omxCheckpoint::omxCheckpoint() : wroteHeader(false), lastCheckpoint(0), lastIterations(0),
				 lastEvaluation(0), fitPending(false),
				 timePerCheckpoint(0), iterPerCheckpoint(0), evalsPerCheckpoint(0), file(NULL)
{}

omxCheckpoint::~omxCheckpoint()
{
	if (file) fclose(file);
}

/* We need to re-design checkpointing when it is possible to run
   more than 1 optimization in parallel. */
void omxCheckpoint::omxWriteCheckpointHeader()
{
	if (wroteHeader) return;
	std::vector< omxFreeVar* > &vars = Global->freeGroup[0]->vars;
	size_t numParam = vars.size();

	// New columns should use the OpenMx prefit to avoid clashing with
	// free parameter names.
	fprintf(file, "OpenMxContext\tOpenMxNumFree\tOpenMxEvals\titerations\ttimestamp");
	for(size_t j = 0; j < numParam; j++) {
		fprintf(file, "\t\"%s\"", vars[j]->name);
	}
	fprintf(file, "\tobjective\n");
	fflush(file);
	wroteHeader = true;
}
 
void omxCheckpoint::message(FitContext *fc, double *est, const char *msg)
{
	_prefit(fc, est, true, msg);
	postfit(fc);
}

void omxCheckpoint::_prefit(FitContext *fc, double *est, bool force, const char *context)
{
	const int timeBufSize = 32;
	char timeBuf[timeBufSize];
	time_t now = time(NULL); // avoid checking unless we need it

	bool doit = force;
	if ((timePerCheckpoint && timePerCheckpoint <= now - lastCheckpoint) ||
	    (iterPerCheckpoint && iterPerCheckpoint <= fc->iterations - lastIterations) ||
	    (evalsPerCheckpoint && evalsPerCheckpoint <= Global->computeCount - lastEvaluation)) {
		doit = true;
	}
	if (!doit) return;

	omxWriteCheckpointHeader();

	std::vector< omxFreeVar* > &vars = fc->varGroup->vars;
	struct tm *nowTime = localtime(&now);
	strftime(timeBuf, timeBufSize, "%b %d %Y %I:%M:%S %p", nowTime);
	fprintf(file, "%s\t%d\t%d\t%d\t%s", context, int(vars.size()), lastEvaluation, lastIterations, timeBuf);

	size_t lx=0;
	size_t numParam = Global->freeGroup[0]->vars.size();
	for (size_t px=0; px < numParam; ++px) {
		if (lx < vars.size() && vars[lx]->id == (int)px) {
			fprintf(file, "\t%.10g", est[lx]);
			++lx;
		} else {
			fprintf(file, "\tNA");
		}
	}
	fflush(file);
	if (fitPending) Rf_error("Checkpoint not reentrant");
	fitPending = true;
	lastCheckpoint = now;
	lastIterations = fc->iterations;
	lastEvaluation = Global->computeCount;
}

void omxCheckpoint::prefit(FitContext *fc, double *est, bool force)
{
	_prefit(fc, est, force, "opt");
}

void omxCheckpoint::postfit(FitContext *fc)
{
	if (!fitPending) return;
	fprintf(file, "\t%.10g\n", fc->fit);
	fflush(file);
	fitPending = false;
}

omxFreeVarLocation *omxFreeVar::getLocation(int matrix)
{
	for (size_t lx=0; lx < locations.size(); lx++) {
		omxFreeVarLocation *loc = &locations[lx];
		if (~loc->matrix == matrix) return loc;
	}
	return NULL;
}
