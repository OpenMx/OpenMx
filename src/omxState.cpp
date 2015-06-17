/*
 *  Copyright 2007-2015 The OpenMx Project
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
#include <set>

#include "omxState.h"
#include "Compute.h"
#include "omxImportFrontendState.h"

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

int FreeVarGroup::lookupVar(int matrix, int row, int col)
{
	for (size_t vx=0; vx < vars.size(); ++vx) {
		std::vector<omxFreeVarLocation> &locations = vars[vx]->locations;
		for (size_t lx=0; lx < locations.size(); lx++) {
			const omxFreeVarLocation &loc = locations[lx];
			if (loc.matrix != matrix) continue;
			if (loc.row == row && loc.col == col) return vx;
		}
	}
	return -1;
}

int FreeVarGroup::lookupVar(omxMatrix *matrix, int row, int col)
{
	return lookupVar(~matrix->matrixNumber, row, col);
}

int FreeVarGroup::lookupVar(const char *name)
{
	for (size_t vx=0; vx < vars.size(); ++vx) {
		if (strcmp(name, vars[vx]->name) == 0) return vx;
	}
	return -1;
}

void FreeVarGroup::cacheDependencies(omxState *os)
{
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
	if (OMX_DEBUG) { log(os); }
}

static int freeVarComp(omxFreeVar *fv1, omxFreeVar *fv2)
{ return fv1->id < fv2->id; }

// NOTE: This assumes that free variables are sorted.
bool FreeVarGroup::isDisjoint(FreeVarGroup *other)
{
	std::vector< omxFreeVar* > overlap(std::max(vars.size(), other->vars.size()));
	std::vector< omxFreeVar* >::iterator it =
		std::set_intersection(vars.begin(), vars.end(),
				      other->vars.begin(), other->vars.end(),
				      overlap.begin(), freeVarComp);
	return it - overlap.begin() == 0;
}

void FreeVarGroup::markDirty(omxState *os)
{
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

        for(size_t i = 0; i < numMats; i++) {
		if (!locations[i]) continue;

		// The point of this is to increment the version numbers
		// on matrices holding free parameters. Otherwise there
		// is no way to tell which free parameters changes when
		// freeSet is used to partition parameters.
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

void FreeVarGroup::log(omxState *os)
{
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
			str += os->matrixList[offset]->name();
		}
	}

	for(size_t i = 0; i < numAlgs; i++) {
		if (dependencies[i + numMats]) {
			str += " ";
			str += os->algebraList[i]->name();
		}
	}
	str += "\n";

	mxLogBig(str);
}

omxGlobal::omxGlobal()
{
	numThreads = 1;
	analyticGradients = 0;
	llScale = -2.0;
	debugProtectStack = OMX_DEBUG;
	computeCount = 0;
	anonAlgebra = 0;
	rowLikelihoodsWarning = false;
	unpackedConfidenceIntervals = false;
	fc = NULL;
	intervals = true;
	gradientTolerance = 1e-6;

	FreeVarGroup *fvg = new FreeVarGroup;
	fvg->id.push_back(FREEVARGROUP_ALL);   // all variables
	freeGroup.push_back(fvg);

	fvg = new FreeVarGroup;
	fvg->id.push_back(FREEVARGROUP_NONE);  // no variables
	freeGroup.push_back(fvg);
}

const char *omxState::matrixToName(int matnum)
{
	return matnum<0? matrixList[~matnum]->name() : algebraList[matnum]->name();
}

void omxState::setWantStage(int stage)
{
	wantStage = stage;
	if (OMX_DEBUG) mxLog("wantStage set to 0x%x", stage);
}

struct ciCmp {
	bool operator() (const omxConfidenceInterval* x, const omxConfidenceInterval* y) const
	{
		if (x->matrix->matrixNumber != y->matrix->matrixNumber) {
			return x->matrix->matrixNumber < y->matrix->matrixNumber;
		} else if (x->row != y->row) {
			return x->row < y->row;
		} else {
			return x->col < y->col;
		}
	}
};

void omxGlobal::unpackConfidenceIntervals()
{
	if (unpackedConfidenceIntervals) return;
	unpackedConfidenceIntervals = true;

	// take care to preserve order
	std::vector<omxConfidenceInterval*> tmp;
	std::swap(tmp, intervalList);
	std::set<omxConfidenceInterval*, ciCmp> uniqueCIs;

	for (int ix=0; ix < (int) tmp.size(); ++ix) {
		omxConfidenceInterval *ci = tmp[ix];
		if (!ci->isWholeAlgebra()) {
			if (uniqueCIs.count(ci) == 0) {
				uniqueCIs.insert(ci);
				intervalList.push_back(ci);
			}
			continue;
		}
		omxMatrix *mat = ci->matrix;
		for (int cx=0; cx < mat->cols; ++cx) {
			for (int rx=0; rx < mat->rows; ++rx) {
				omxConfidenceInterval *cell = new omxConfidenceInterval(*ci);
				cell->name = string_snprintf("%s[%d,%d]", ci->name.c_str(), 1+rx, 1+cx);
				cell->row = rx;
				cell->col = cx;
				if (uniqueCIs.count(cell) == 0) {
					uniqueCIs.insert(cell);
					intervalList.push_back(cell);
				} else {
					delete cell;
				}
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

void omxState::init()
{
	wantStage = 0;
	currentRow = -1;
}

void omxState::loadDefinitionVariables(bool start)
{
	for(int ex = 0; ex < int(expectationList.size()); ++ex) {
		omxExpectation *e1 = expectationList[ex];
		if (e1->defVars.size() == 0) continue;
		int row = 0;
		if (start) {
			if (e1->data->rows != 1) {
				e1->loadFakeData(NA_REAL);
				continue;
			}
		} else {
			int obs = omxDataNumObs(e1->data);
			for (int dx=0; dx < obs; ++dx) {
				if (omxDataIndex(e1->data, dx) == 0) {
					row = dx;
					break;
				}
			}
		}
		Eigen::VectorXd oldDefs(e1->defVars.size());
		oldDefs.setConstant(NA_REAL);
		e1->handleDefinitionVarList(this, row, oldDefs.data());
	}
}

omxState::omxState(omxState *src)
{
	init();

	dataList			= src->dataList;
		
	for(size_t mx = 0; mx < src->matrixList.size(); mx++) {
		// TODO: Smarter inference for which matrices to duplicate
		matrixList.push_back(omxDuplicateMatrix(src->matrixList[mx], this));
	}

	for(size_t j = 0; j < src->expectationList.size(); j++) {
		// TODO: Smarter inference for which expectations to duplicate
		expectationList.push_back(omxDuplicateExpectation(src->expectationList[j], this));
	}

	for(size_t j = 0; j < src->algebraList.size(); j++) {
		// TODO: Smarter inference for which algebras to duplicate
		algebraList.push_back(omxDuplicateMatrix(src->algebraList[j], this));
	}

	for(size_t j = 0; j < algebraList.size(); j++) {
		omxDuplicateAlgebra(algebraList[j], src->algebraList[j], this);
	}

	omxInitialMatrixAlgebraCompute(this, NULL); // pass in FC TODO

	for(size_t j = 0; j < src->expectationList.size(); j++) {
		// TODO: Smarter inference for which expectations to duplicate
		omxCompleteExpectation(expectationList[j]);
	}

	for (int ax=0; ax < (int) algebraList.size(); ++ax) {
		omxMatrix *matrix = algebraList[ax];
		if (!matrix->fitFunction) continue;
		omxCompleteFitFunction(matrix);
	}

	currentRow 		= src->currentRow;
}

omxState::~omxState()
{
	if(OMX_DEBUG) { mxLog("Freeing %d Constraints.", (int) conList.size());}
	for(int k = 0; k < (int) conList.size(); k++) {
		delete conList[k];
	}

	for(size_t ax = 0; ax < algebraList.size(); ax++) {
		// free argument tree
		omxFreeMatrix(algebraList[ax]);
	}

	for(size_t ax = 0; ax < algebraList.size(); ax++) {
		algebraList[ax]->hasMatrixNumber = false;
		omxFreeMatrix(algebraList[ax]);
	}

	if(OMX_DEBUG) { mxLog("Freeing %d Matrices.", (int) matrixList.size());}
	for(size_t mk = 0; mk < matrixList.size(); mk++) {
		matrixList[mk]->hasMatrixNumber = false;
		omxFreeMatrix(matrixList[mk]);
	}
		
	if(OMX_DEBUG) { mxLog("Freeing %d Model Expectations.", (int) expectationList.size());}
	for(size_t ex = 0; ex < expectationList.size(); ex++) {
		omxFreeExpectationArgs(expectationList[ex]);
	}
}

omxGlobal::~omxGlobal()
{
	if (fc) {
		omxState *state = fc->state;
		delete fc;
		delete state;
	}
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

void mxLogBig(const std::string &str)   // thread-safe
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
	ssize_t got = 0;
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
	if (got <= 0) Rf_error("mxLog(%s) failed with errno=%d", buf2, got);
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

	if (OMX_DEBUG) mxLog("checkpointMessage: %s", str.c_str());

	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->message(fc, est, str.c_str());
	}
}

void omxGlobal::checkpointPrefit(const char *callerName, FitContext *fc, double *est, bool force)
{
	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->prefit(callerName, fc, est, force);
	}
}

void omxGlobal::checkpointPostfit(FitContext *fc)
{
	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->postfit(fc);
	}
}

UserConstraint::UserConstraint(FitContext *fc, const char *name, omxMatrix *arg1, omxMatrix *arg2) :
	super(name)
{
	omxState *state = fc->state;
	omxMatrix *args[2] = {arg1, arg2};
	pad = omxNewAlgebraFromOperatorAndArgs(10, args, 2, state); // 10 = binary subtract
	state->setWantStage(FF_COMPUTE_INITIAL_FIT);
	refresh(fc);
	int nrows = pad->rows;
	int ncols = pad->cols;
	size = nrows * ncols;
	if (size == 0) {
		Rf_warning("Constraint '%s' evaluated to a 0x0 matrix and will have no effect", name);
	}
}

UserConstraint::~UserConstraint()
{
	omxFreeMatrix(pad);
}

void UserConstraint::refresh(FitContext *fc)
{
	omxRecompute(pad, fc);
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
	std::vector< omxFreeVar* > &vars = Global->findVarGroup(FREEVARGROUP_ALL)->vars;
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
	size_t numParam = Global->findVarGroup(FREEVARGROUP_ALL)->vars.size();
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

void omxCheckpoint::prefit(const char *callerName, FitContext *fc, double *est, bool force)
{
	_prefit(fc, est, force, callerName);
}

void omxCheckpoint::postfit(FitContext *fc)
{
	if (!fitPending) return;
	fprintf(file, "\t%.10g\n", fc->fit);
	fflush(file);
	fitPending = false;
}

const omxFreeVarLocation *omxFreeVar::getLocation(int matrix) const
{
	for (size_t lx=0; lx < locations.size(); lx++) {
		const omxFreeVarLocation &loc = locations[lx];
		if (loc.matrix == matrix) return &loc;
	}
	return NULL;
}

const omxFreeVarLocation *omxFreeVar::getLocation(omxMatrix *mat) const
{ return getLocation(~mat->matrixNumber); }
