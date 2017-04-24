/*
 *  Copyright 2007-2017 The OpenMx Project
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

#include <time.h>
#include <stdarg.h>
#include <errno.h>
#include <set>

#include "omxState.h"
#include "Compute.h"
#include "omxImportFrontendState.h"
#include "EnableWarnings.h"

struct omxGlobal *Global = NULL;
static bool mxLogEnabled = false;

SEXP enableMxLog()
{
	mxLogEnabled = true;
	return Rf_ScalarLogical(1);
}

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
		std::vector<omxFreeVarLocation> &locVec = vars[vx]->locations;
		for (size_t lx=0; lx < locVec.size(); lx++) {
			const omxFreeVarLocation &loc = locVec[lx];
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

/* might be useful?
int FreeVarGroup::lookupVar(int id)
{
	std::vector<int>::iterator low =
		std::lower_bound(vars.begin(), vars.end(), id);
	if (low == vars.end()) return -1;
	int got = low - vars.begin();
	if (vars[got]->id == id) return got;
	return -1;
}
*/

void FreeVarGroup::cacheDependencies(omxState *os)
{
	size_t numMats = os->matrixList.size();
	size_t numAlgs = os->algebraList.size();

	dependencies.assign(numMats + numAlgs, false);
	locations.assign(numMats, false);

	for (size_t vx = 0; vx < vars.size(); vx++) {
		omxFreeVar *fv = vars[vx];
		auto deps = fv->getDeps();
		for (int index = 0; index < deps.size(); index++) {
			dependencies[deps[index] + numMats] = true;
		}
		for (size_t lx=0; lx < fv->locations.size(); ++lx) {
			locations[fv->locations[lx].matrix] = true;
		}
	}

	for(size_t i = 0; i < numMats; i++) {
		if (!locations[i]) continue;
		os->matrixList[i]->unshareMemroyWithR();
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
	silent = true;
	lastProgressReport = time(0);
	previousReportLength = 0;
	previousReportFit = 0;
	previousComputeCount = 0;
	mxLogSetCurrentRow(-1);
	numThreads = 1;
	analyticGradients = 0;
	llScale = -2.0;
	debugProtectStack = OMX_DEBUG;
	anonAlgebra = 0;
	rowLikelihoodsWarning = false;
	unpackedConfidenceIntervals = false;
	topFc = NULL;
	intervals = true;
	gradientTolerance = 1e-6;
	boundsUpdated = false;
	dataTypeWarningCount = 0;

	RAMInverseOpt = true;
	RAMMaxDepth = 30;

	FreeVarGroup *fvg = new FreeVarGroup;
	fvg->id.push_back(FREEVARGROUP_ALL);   // all variables
	freeGroup.push_back(fvg);

	fvg = new FreeVarGroup;
	fvg->id.push_back(FREEVARGROUP_NONE);  // no variables
	freeGroup.push_back(fvg);
}

omxMatrix *omxState::getMatrixFromIndex(int matnum) const
{
	return matnum<0? matrixList[~matnum] : algebraList[matnum];
}

omxMatrix *omxState::lookupDuplicate(omxMatrix *element) const
{
	if (!element->hasMatrixNumber) Rf_error("lookupDuplicate without matrix number");
	return getMatrixFromIndex(element->matrixNumber);
}

void omxState::setWantStage(int stage)
{
	wantStage = stage;
	if (OMX_DEBUG) mxLog("wantStage set to 0x%x", stage);
}

omxMatrix *ConfidenceInterval::getMatrix(omxState *st) const
{
	return st->getMatrixFromIndex(matrixNumber);
}

struct ciCmp {
	bool operator() (const ConfidenceInterval* x, const ConfidenceInterval* y) const
	{
		if (x->matrixNumber != y->matrixNumber) {
			return x->matrixNumber < y->matrixNumber;
		} else if (x->row != y->row) {
			return x->row < y->row;
		} else {
			return x->col < y->col;
		}
	}
};

void omxGlobal::unpackConfidenceIntervals(omxState *currentState)
{
	if (unpackedConfidenceIntervals) return;
	unpackedConfidenceIntervals = true;

	// take care to preserve order
	std::vector<ConfidenceInterval*> tmp;
	std::swap(tmp, intervalList);
	std::set<ConfidenceInterval*, ciCmp> uniqueCIs;

	for (int ix=0; ix < (int) tmp.size(); ++ix) {
		ConfidenceInterval *ci = tmp[ix];
		if (!ci->isWholeAlgebra()) {
			auto iter = uniqueCIs.find(ci);
			if (iter == uniqueCIs.end()) {
				uniqueCIs.insert(ci);
				intervalList.push_back(ci);
			} else if (ci->cmpBoundAndType(**iter)) {
				Rf_warning("Different confidence intervals '%s' and '%s' refer to the same thing",
					   ci->name.c_str(), (*iter)->name.c_str());
			}
			continue;
		}
		omxMatrix *mat = ci->getMatrix(currentState);
		for (int cx=0; cx < mat->cols; ++cx) {
			for (int rx=0; rx < mat->rows; ++rx) {
				ConfidenceInterval *cell = new ConfidenceInterval(*ci);
				cell->name = string_snprintf("%s[%d,%d]", ci->name.c_str(), 1+rx, 1+cx);
				cell->row = rx;
				cell->col = cx;
				auto iter = uniqueCIs.find(cell);
				if (iter == uniqueCIs.end()) {
					uniqueCIs.insert(cell);
					intervalList.push_back(cell);
				} else {
					if (cell->cmpBoundAndType(**iter)) {
						Rf_warning("Different confidence intervals '%s' and '%s' refer to the same thing",
							   cell->name.c_str(), (*iter)->name.c_str());
					}
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

int omxState::nextId = 0;

void omxState::init()
{
	stateId = ++nextId;
	setWantStage(FF_COMPUTE_INITIAL_FIT);
}

void omxState::loadDefinitionVariables(bool start)
{
	if (OMX_DEBUG) {
		mxLog("omxState[%d]::loadDefinitionVariables(start=%d)", getId(), start);
	}
	for(int ex = 0; ex < int(dataList.size()); ++ex) {
		omxData *e1 = dataList[ex];
		if (e1->defVars.size() == 0) continue;
		if (start && e1->rows != 1) {
			e1->loadFakeData(this, NA_REAL);
			continue;
		}
		e1->loadDefVars(this, 0);
	}
}

omxState::omxState(omxState *src) : clone(true)
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

	for(size_t mx = 0; mx < src->matrixList.size(); mx++) {
		// TODO: Smarter inference for which matrices to duplicate
		matrixList[mx]->copyAttr(src->matrixList[mx]);
	}

	for (size_t xx=0; xx < src->conListX.size(); ++xx) {
		conListX.push_back(src->conListX[xx]->duplicate(this));
	}
}

void omxState::initialRecalc(FitContext *fc)
{
	omxInitialMatrixAlgebraCompute(fc);

	for(size_t j = 0; j < expectationList.size(); j++) {
		// TODO: Smarter inference for which expectations to duplicate
		omxCompleteExpectation(expectationList[j]);
	}

	for (int ax=0; ax < (int) algebraList.size(); ++ax) {
		omxMatrix *matrix = algebraList[ax];
		if (!matrix->fitFunction) continue;
		omxCompleteFitFunction(matrix);
		omxFitFunctionCompute(matrix->fitFunction, FF_COMPUTE_INITIAL_FIT, fc);
	}

	for (size_t xx=0; xx < conListX.size(); ++xx) {
		conListX[xx]->prep(fc);
	}
}

void omxState::invalidateCache()
{
	for(size_t ex = 0; ex < expectationList.size(); ex++) {
		expectationList[ex]->invalidateCache();
	}
	for (int ax=0; ax < (int) algebraList.size(); ++ax) {
		omxMatrix *matrix = algebraList[ax];
		if (!matrix->fitFunction) continue;
		matrix->fitFunction->invalidateCache();
	}
}

omxState::~omxState()
{
	if(OMX_DEBUG) { mxLog("Freeing %d Constraints.", (int) conListX.size());}
	for(int k = 0; k < (int) conListX.size(); k++) {
		delete conListX[k];
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
	if (previousReportLength) reportProgressStr("");
	if (topFc) {
		omxState *state = topFc->state;
		delete topFc;
		delete state;
	}
	for (size_t cx=0; cx < intervalList.size(); ++cx) {
		delete intervalList[cx];
	}
	for (size_t cx=0; cx < computeList.size(); ++cx) {
		delete computeList[cx];
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

#if OMX_DEBUG
// __thread is a gcc extension. I'm not sure about portability to
// other compilers. This is only here for debugging so you can remove
// it if you have trouble getting it through the compiler on your
// system.

static __thread int mxLogCurrentRow;

void mxLogSetCurrentRow(int row)
{
	mxLogCurrentRow = row;
}
#else
static const int mxLogCurrentRow = -1;
void mxLogSetCurrentRow(int row) {}
#endif

static ssize_t mxLogWriteSynchronous(const char *outBuf, int len)
{
	if (!mxLogEnabled) return len;

	int maxRetries = 20;
	ssize_t wrote = 0;
	ssize_t got;
#pragma omp critical
	{
		while (--maxRetries > 0) {
			got = write(2, outBuf + wrote, len - wrote);
			if (got == -EINTR) continue;
			if (got < 0) break;
			wrote += got;
			if (wrote == len) break;
		}
	}
	return wrote;
}

void mxLogBig(const std::string &str)   // thread-safe
{
	ssize_t len = ssize_t(str.size());
	if (len == 0) Rf_error("Attempt to log 0 characters with mxLogBig");

	std::string fullstr;
	if (mxLogCurrentRow == -1) {
		fullstr = string_snprintf("[%d] ", omx_absolute_thread_num());
	} else {
		fullstr = string_snprintf("[%d@%d] ", omx_absolute_thread_num(), mxLogCurrentRow);
	}
	fullstr += str;
	len = ssize_t(fullstr.size());
	
	const char *outBuf = fullstr.c_str();
	ssize_t wrote = mxLogWriteSynchronous(outBuf, len);
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

	int len;
	if (mxLogCurrentRow == -1) {
		len = snprintf(buf2, maxLen, "[%d] %s\n", omx_absolute_thread_num(), buf1);
	} else {
		len = snprintf(buf2, maxLen, "[%d@%d] %s\n", omx_absolute_thread_num(), mxLogCurrentRow, buf1);
	}

	ssize_t wrote = mxLogWriteSynchronous(buf2, len);
	if (wrote != len) Rf_error("mxLog only wrote %d/%d, errno=%d", wrote, len, errno);
}

void omxGlobal::reportProgressStr(const char *msg)
{
	ProtectedSEXP theCall(Rf_allocVector(LANGSXP, 3));
	SETCAR(theCall, Rf_install("imxReportProgress"));
	ProtectedSEXP Rmsg(Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(Rmsg, 0, Rf_mkChar(msg));
	SETCADR(theCall, Rmsg);
	SETCADDR(theCall, Rf_ScalarInteger(previousReportLength));
	Rf_eval(theCall, R_GlobalEnv);
}

void omxGlobal::reportProgress(const char *context, FitContext *fc)
{
	if (omx_absolute_thread_num() != 0) {
		mxLog("omxGlobal::reportProgress called in a thread context (report this bug to developers)");
		return;
	}

	R_CheckUserInterrupt();

	time_t now = time(0);
	if (silent || now - lastProgressReport < 1 || fc->getGlobalComputeCount() == previousComputeCount) return;

	lastProgressReport = now;

	std::string str;
	if (previousReportFit == 0.0 || previousReportFit == fc->fit) {
		str = string_snprintf("%s %d %.6g",
				      context, fc->getGlobalComputeCount(), fc->fit);
	} else {
		str = string_snprintf("%s %d %.6g %.4g",
				      context, fc->getGlobalComputeCount(), fc->fit, fc->fit - previousReportFit);
	}

	reportProgressStr(str.c_str());
	previousReportLength = str.size();
	previousReportFit = fc->fit;
	previousComputeCount = fc->getGlobalComputeCount();
}

void diagParallel(int verbose, const char* msg, ...)
{
	if (!verbose && !Global->parallelDiag) return;

	const int maxLen = 240;
	char buf1[maxLen];

	va_list ap;
	va_start(ap, msg);
	vsnprintf(buf1, maxLen, msg, ap);
	va_end(ap);

	if (verbose) {
		mxLog("%s", buf1);
	} else if (Global->parallelDiag) {
		ProtectedSEXP theCall(Rf_allocVector(LANGSXP, 2));
		SETCAR(theCall, Rf_install("message"));
		ProtectedSEXP Rmsg(Rf_allocVector(STRSXP, 1));
		SET_STRING_ELT(Rmsg, 0, Rf_mkChar(buf1));
		SETCADR(theCall, Rmsg);
		Rf_eval(theCall, R_GlobalEnv);
	}
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

void omxGlobal::checkpointPostfit(const char *callerName, FitContext *fc, double *est, bool force)
{
	for(size_t i = 0; i < checkpointList.size(); i++) {
		checkpointList[i]->postfit(callerName, fc, est, force);
	}
}

void UserConstraint::prep(FitContext *fc)
{
	fc->state->setWantStage(FF_COMPUTE_INITIAL_FIT);
	refresh(fc);
	nrows = pad->rows;
	ncols = pad->cols;
	size = nrows * ncols;
	if (size == 0) {
		Rf_warning("Constraint '%s' evaluated to a 0x0 matrix and will have no effect", name);
	}
	omxAlgebraPreeval(pad, fc);
	if(jacobian){
		jacMap.resize(jacobian->cols);
		std::vector<const char*> *jacColNames = &jacobian->colnames;
		for (size_t nx=0; nx < jacColNames->size(); ++nx) {
			int to = fc->varGroup->lookupVar((*jacColNames)[nx]);
			jacMap[nx] = to;
		}
	}
}

UserConstraint::UserConstraint(FitContext *fc, const char *_name, omxMatrix *arg1, omxMatrix *arg2, omxMatrix *jac, int lin) :
	super(_name)
{
	omxState *state = fc->state;
	omxMatrix *args[2] = {arg1, arg2};
	pad = omxNewAlgebraFromOperatorAndArgs(10, args, 2, state); // 10 = binary subtract
	jacobian = jac;
	linear = lin;
}

omxConstraint *UserConstraint::duplicate(omxState *dest)
{
	omxMatrix *args[2] = {
		dest->lookupDuplicate(pad->algebra->algArgs[0]),
		dest->lookupDuplicate(pad->algebra->algArgs[1])
	};

	UserConstraint *uc = new UserConstraint(name);
	uc->opCode = opCode;
	uc->pad = omxNewAlgebraFromOperatorAndArgs(10, args, 2, dest); // 10 = binary subtract
	uc->jacobian = jacobian;
	uc->linear = linear;
	return uc;
}

void UserConstraint::refreshAndGrab(FitContext *fc, Type ineqType, double *out)
{
	fc->incrComputeCount();
	refresh(fc);

	for(int k = 0; k < size; k++) {
		double got = pad->data[k];
		if (opCode != ineqType) got = -got;
		out[k] = got;
	}
}

UserConstraint::~UserConstraint()
{
	omxFreeMatrix(pad);
}

void UserConstraint::refresh(FitContext *fc)
{
	omxRecompute(pad, fc);
	//omxRecompute(jacobian, fc); //<--Not sure if Jacobian needs to be recomputed every time constraint function does.
}

omxCheckpoint::omxCheckpoint() : wroteHeader(false), lastCheckpoint(0), lastIterations(0),
				 lastEvaluation(0),
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
	postfit(msg, fc, est, true);
}

void omxCheckpoint::postfit(const char *context, FitContext *fc, double *est, bool force)
{
	const int timeBufSize = 32;
	char timeBuf[timeBufSize];
	time_t now = time(NULL); // avoid checking unless we need it
	int curEval = fc->getGlobalComputeCount();

	bool doit = force;
	if ((timePerCheckpoint && timePerCheckpoint <= now - lastCheckpoint) ||
	    (iterPerCheckpoint && iterPerCheckpoint <= fc->iterations - lastIterations) ||
	    (evalsPerCheckpoint && evalsPerCheckpoint <= curEval - lastEvaluation)) {
		doit = true;
	}
	if (!doit) return;

#pragma omp critical
	{
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
		fprintf(file, "\t%.10g\n", fc->fit);
		fflush(file);
		lastCheckpoint = now;
		lastIterations = fc->iterations;
		lastEvaluation = curEval;
	}
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

const omxFreeVarLocation *omxFreeVar::getOnlyOneLocation(int matrix, bool &moreThanOne) const
{
	moreThanOne = false;
	const omxFreeVarLocation *result = NULL;
	for (size_t lx=0; lx < locations.size(); lx++) {
		const omxFreeVarLocation &loc = locations[lx];
		if (loc.matrix == matrix) {
			if (result) { moreThanOne = true; return NULL; }
			result = &loc;
		}
	}
	return result;
}

const omxFreeVarLocation *omxFreeVar::getOnlyOneLocation(omxMatrix *mat, bool &moreThanOne) const
{ return getOnlyOneLocation(~mat->matrixNumber, moreThanOne); }

void omxFreeVar::markDirty(omxState *os)
{
	auto deps = getDeps();
	for (int dx=0; dx < deps.size(); ++dx) {
		int dep = deps[dx];
		if (dep < 0) {
			omxMarkDirty(os->matrixList[~dep]);
		} else {
			omxMarkDirty(os->algebraList[dep]);
		}
	}

	for (int lx=0; lx < int(locations.size()); ++lx) {
		omxMarkClean(os->matrixList[locations[lx].matrix]);
	}
}

void omxFreeVar::copyToState(omxState *os, double val)
{
	for(size_t l = 0; l < locations.size(); l++) {
		omxFreeVarLocation *loc = &locations[l];
		omxMatrix *matrix = os->matrixList[loc->matrix];
		int row = loc->row;
		int col = loc->col;
		omxSetMatrixElement(matrix, row, col, val);
		if (OMX_DEBUG) {
			mxLog("free var %s, matrix %s[%d, %d] = %f",
			      name, matrix->name(), row, col, val);
		}
	}
}
