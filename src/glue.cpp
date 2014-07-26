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

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxDefines.h"
#include "glue.h"
#include "omxOpenmpWrap.h"
#include "omxState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxFitFunction.h"
#include "omxExpectation.h"
#include "omxNPSOLSpecific.h"
#include "omxImportFrontendState.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "dmvnorm.h"
#include "npsolswitch.h"
#include "matrix.h"

static SEXP do_logm_eigen(SEXP x)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;

    if (!Rf_isNumeric(x) || !Rf_isMatrix(x)) Rf_error("invalid argument");

    dims = Rf_getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m) Rf_error("non-square matrix");
    if (n == 0) return(Rf_allocVector(REALSXP, 0));

    PROTECT(z = Rf_allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    logm_eigen(n, rx, rz);

    Rf_setAttrib(z, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);

    return z;
}

static SEXP do_expm_eigen(SEXP x)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;

    if (!Rf_isNumeric(x) || !Rf_isMatrix(x)) Rf_error("invalid argument");

    dims = Rf_getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m) Rf_error("non-square matrix");
    if (n == 0) return(Rf_allocVector(REALSXP, 0));

    PROTECT(z = Rf_allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    expm_eigen(n, rx, rz);

    Rf_setAttrib(z, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);

    return z;
}

static SEXP has_NPSOL()
{ return Rf_ScalarLogical(HAS_NPSOL); }

static SEXP has_openmp()
{
#if defined(_OPENMP)
	bool opm = true;
#else
	bool opm = false;
#endif
	return Rf_ScalarLogical(opm);
}

static R_CallMethodDef callMethods[] = {
	{"backend", (DL_FUNC) omxBackend, 10},
	{"callAlgebra", (DL_FUNC) omxCallAlgebra, 3},
	{"findIdenticalRowsData", (DL_FUNC) findIdenticalRowsData, 5},
	{"Dmvnorm_wrapper", (DL_FUNC) dmvnorm_wrapper, 3},
	{"hasNPSOL_wrapper", (DL_FUNC) has_NPSOL, 0},
	{"sparseInvert_wrapper", (DL_FUNC) sparseInvert_wrapper, 1},
	{"hasOpenMP_wrapper", (DL_FUNC) has_openmp, 0},
	{"do_logm_eigen", (DL_FUNC) &do_logm_eigen, 1},
	{"do_expm_eigen", (DL_FUNC) &do_expm_eigen, 1},
	{NULL, NULL, 0}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_OpenMx(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	// There is no code that will change behavior whether openmp
	// is set for nested or not. I'm just keeping this in case it
	// makes a difference with older versions of openmp. 2012-12-24 JNP
#if defined(_OPENMP) && _OPENMP <= 200505
	omp_set_nested(0);
#endif
}

void R_unload_OpenMx(DllInfo *) {
	// keep this stub in case we need it
}

#ifdef  __cplusplus
}
#endif

void string_to_try_Rf_error( const std::string& str )
{
	Rf_error("%s", str.c_str());
}

void exception_to_try_Rf_error( const std::exception& ex )
{
	string_to_try_Rf_error(ex.what());
}

SEXP MxRList::asR()
{
	// detect duplicate keys? TODO
	SEXP names, ans;
	int len = size();
	Rf_protect(names = Rf_allocVector(STRSXP, len));
	Rf_protect(ans = Rf_allocVector(VECSXP, len));
	for (int lx=0; lx < len; ++lx) {
		const char *p1 = (*this)[lx].first;
		SEXP p2 = (*this)[lx].second;
		if (!p1 || !p2) Rf_error("Attempt to return NULL pointer to R");
		SET_STRING_ELT(names, lx, Rf_mkChar(p1));
		SET_VECTOR_ELT(ans,   lx, p2);
	}
	Rf_namesgets(ans, names);
	return ans;
}

static void
friendlyStringToLogical(const char *key, const char *str, int *out)
{
	int understood = FALSE;
	int newVal;
	if (matchCaseInsensitive(str, "Yes")) {
		understood = TRUE;
		newVal = 1;
	} else if (matchCaseInsensitive(str, "No")) {
		understood = TRUE;
		newVal = 0;
	} else if (isdigit(str[0]) && (atoi(str) == 1 || atoi(str) == 0)) {
		understood = TRUE;
		newVal = atoi(str);
	}
	if (!understood) {
		Rf_warning("Expecting 'Yes' or 'No' for '%s' but got '%s', ignoring", key, str);
		return;
	}
	if(OMX_DEBUG) { mxLog("%s=%d", key, newVal); }
	*out = newVal;
}

// TODO: make member of omxGlobal class
static void readOpts(SEXP options, int *ciMaxIterations, int *numThreads,
		     int *analyticGradients)
{
		int numOptions = Rf_length(options);
		SEXP optionNames;
		Rf_protect(optionNames = Rf_getAttrib(options, R_NamesSymbol));
		for(int i = 0; i < numOptions; i++) {
			const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
			const char *nextOptionValue = CHAR(Rf_asChar(VECTOR_ELT(options, i)));
			if (matchCaseInsensitive(nextOptionName, "CI Max Iterations")) {
				int newvalue = atoi(nextOptionValue);
				if (newvalue > 0) *ciMaxIterations = newvalue;
			} else if(matchCaseInsensitive(nextOptionName, "Analytic Gradients")) {
				friendlyStringToLogical(nextOptionName, nextOptionValue, analyticGradients);
			} else if(matchCaseInsensitive(nextOptionName, "loglikelihoodScale")) {
				Global->llScale = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "Number of Threads")) {
#ifdef _OPENMP
				*numThreads = atoi(nextOptionValue);
				if (*numThreads < 1) {
					Rf_warning("Computation will be too slow with %d threads; using 1 thread instead", *numThreads);
					*numThreads = 1;
				}
#endif
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsA")) {
				Global->maxptsa = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsB")) {
				Global->maxptsb = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsC")) {
				Global->maxptsc = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnAbsEps")) {
				Global->absEps = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnRelEps")) {
				Global->relEps = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "maxStackDepth")) {
				Global->maxStackDepth = atoi(nextOptionValue);
			} else {
				// ignore
			}
		}
		Rf_unprotect(1); // optionNames
}

/* Main functions */
SEXP omxCallAlgebra2(SEXP matList, SEXP algNum, SEXP options) {

	omxManageProtectInsanity protectManager;

	if(OMX_DEBUG) { mxLog("-----------------------------------------------------------------------");}
	if(OMX_DEBUG) { mxLog("Explicit call to algebra %d.", INTEGER(algNum)[0]);}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;

	FitContext::setRFitFunction(NULL);
	Global = new omxGlobal;

	globalState = new omxState;
	omxInitState(globalState);

	readOpts(options, &Global->ciMaxIterations, &Global->numThreads, 
			&Global->analyticGradients);

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { mxLog("Processing %d matrix(ces).", Rf_length(matList));}

	std::vector<omxMatrix *> args(Rf_length(matList));
	for(k = 0; k < Rf_length(matList); k++) {
		Rf_protect(nextMat = VECTOR_ELT(matList, k));	// This is the matrix + populations
		args[k] = omxNewMatrixFromRPrimitive(nextMat, globalState, 1, - k - 1);
		globalState->matrixList.push_back(args[k]);
		if(OMX_DEBUG) {
			mxLog("Matrix[%d] initialized (%d x %d)",
				k, globalState->matrixList[k]->rows, globalState->matrixList[k]->cols);
		}
	}

	algebra = omxNewAlgebraFromOperatorAndArgs(algebraNum, args.data(), Rf_length(matList), globalState);

	if(algebra==NULL) {
		Rf_error("Failed to build algebra");
	}

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}

	omxRecompute(algebra, FF_COMPUTE_FIT, NULL);

	Rf_protect(ans = Rf_allocMatrix(REALSXP, algebra->rows, algebra->cols));
	for(l = 0; l < algebra->rows; l++)
		for(j = 0; j < algebra->cols; j++)
			REAL(ans)[j * algebra->rows + l] =
				omxMatrixElement(algebra, l, j);

	if(OMX_DEBUG) { mxLog("All Algebras complete."); }

	const char *bads = Global->getBads();

	omxFreeMatrix(algebra);
	omxFreeState(globalState);
	delete Global;

	if (bads) Rf_error(bads);

	return ans;
}

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options)
{
	try {
		return omxCallAlgebra2(matList, algNum, options);
	} catch( std::exception& __ex__ ) {
		exception_to_try_Rf_error( __ex__ );
	} catch(...) {
		string_to_try_Rf_error( "c++ exception (unknown reason)" );
	}
}

SEXP omxBackend2(SEXP constraints, SEXP matList,
		 SEXP varList, SEXP algList, SEXP expectList, SEXP computeList,
		 SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options)
{
	SEXP nextLoc;

	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these.  For now, all checking is done on the front-end. */
//	if(!isVector(matList)) Rf_error ("matList must be a list");
//	if(!isVector(algList)) Rf_error ("algList must be a list");

	omxManageProtectInsanity protectManager;

	FitContext::setRFitFunction(NULL);
	Global = new omxGlobal;

	/* Create new omxState for current state storage and initialize it. */
	globalState = new omxState;
	omxInitState(globalState);

	readOpts(options, &Global->ciMaxIterations, &Global->numThreads, 
			&Global->analyticGradients);
#if HAS_NPSOL
	omxSetNPSOLOpts(options);
#endif

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxDataEntities(data);
    
	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxExpectationEntities(expectList);

	for (int dx=0; dx < (int) globalState->dataList.size(); ++dx) {
		globalState->dataList[dx]->connectDynamicData();
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxMatrixEntities(matList);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	std::vector<double> startingValues;
	omxProcessFreeVarList(varList, &startingValues);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxAlgebraEntities(algList);

	/* Process Matrix and Algebra Population Function */
	/*
	  Each matrix is a list containing a matrix and the other matrices/algebras that are
	  populated into it at each iteration.  The first element is already processed, above.
	  The rest of the list will be processed here.
	*/
	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	for(int j = 0; j < Rf_length(matList); j++) {
		Rf_protect(nextLoc = VECTOR_ELT(matList, j));		// This is the matrix + populations
		globalState->matrixList[j]->omxProcessMatrixPopulationList(nextLoc);
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxInitialMatrixAlgebraCompute(globalState, NULL);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxComputeEntities(computeList);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxCompleteMxExpectationEntities();

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxCompleteMxFitFunction(algList);

	FitContext fc(startingValues);

	// Nothing depend on constraints so we can process them last.
	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessConstraints(constraints, &fc);

	if (isErrorRaised()) {
		Rf_error(Global->getBads());
	}

	omxCompute *topCompute = NULL;
	if (Global->computeList.size()) topCompute = Global->computeList[0];

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessConfidenceIntervals(intervalList);

	omxProcessCheckpointOptions(checkpointList);

	for (size_t vg=0; vg < Global->freeGroup.size(); ++vg) {
		Global->freeGroup[vg]->cacheDependencies();
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	if (protectManager.getDepth() > Global->maxStackDepth) {
		Rf_error("Protection stack too large; report this problem to the OpenMx forum");
	}

	if (topCompute && !isErrorRaised()) {
		topCompute->compute(&fc);
	}

	SEXP evaluations;
	Rf_protect(evaluations = Rf_allocVector(REALSXP,1));

	REAL(evaluations)[0] = Global->computeCount;

	if (topCompute && !isErrorRaised() && globalState->stale) {
		fc.copyParamToModel(globalState);
	}

	MxRList result;

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxExportResults(globalState, &result);

	if (topCompute && !isErrorRaised()) {
		LocalComputeResult cResult;
		topCompute->collectResults(&fc, &cResult, &result);

		if (cResult.size()) {
			SEXP computes;
			Rf_protect(computes = Rf_allocVector(VECSXP, cResult.size() * 2));
			for (size_t cx=0; cx < cResult.size(); ++cx) {
				std::pair<int, MxRList*> &c1 = cResult[cx];
				SET_VECTOR_ELT(computes, cx*2, Rf_ScalarInteger(c1.first));
				SET_VECTOR_ELT(computes, cx*2+1, c1.second->asR());
				delete c1.second;
			}
			result.add("computes", computes);
		}

		if (fc.wanted & FF_COMPUTE_FIT) {
			result.add("fit", Rf_ScalarReal(fc.fit));
			result.add("Minus2LogLikelihood", Rf_ScalarReal(fc.fit));
		}
		if (fc.wanted & FF_COMPUTE_BESTFIT) {
			result.add("minimum", Rf_ScalarReal(fc.fit));
		}

		size_t numFree = Global->freeGroup[FREEVARGROUP_ALL]->vars.size();
		if (numFree) {
			// move other global reporting here TODO

			SEXP estimate;
			Rf_protect(estimate = Rf_allocVector(REALSXP, numFree));
			memcpy(REAL(estimate), fc.est, sizeof(double)*numFree);
			result.add("estimate", estimate);

			if (fc.stderrs) {
				SEXP stdErrors;
				Rf_protect(stdErrors = Rf_allocMatrix(REALSXP, numFree, 1));
				memcpy(REAL(stdErrors), fc.stderrs, sizeof(double) * numFree);
				result.add("standardErrors", stdErrors);
			}
			if (fc.wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
				result.add("infoDefinite", Rf_ScalarLogical(fc.infoDefinite));
				result.add("conditionNumber", Rf_ScalarReal(fc.infoCondNum));
			}
		}
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	MxRList backwardCompatStatus;
	backwardCompatStatus.add("code", Rf_ScalarInteger(fc.inform));
	backwardCompatStatus.add("status", Rf_ScalarInteger(-isErrorRaised()));

	if (isErrorRaised()) {
		SEXP msg;
		Rf_protect(msg = Rf_allocVector(STRSXP, 1));
		SET_STRING_ELT(msg, 0, Rf_mkChar(Global->getBads()));
		result.add("error", msg);
		backwardCompatStatus.add("statusMsg", msg);
	}

	result.add("status", backwardCompatStatus.asR());
	result.add("iterations", Rf_ScalarInteger(fc.iterations));
	result.add("evaluations", evaluations);

	omxFreeState(globalState);
	delete Global;

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	return result.asR();
}

SEXP omxBackend(SEXP constraints, SEXP matList,
		SEXP varList, SEXP algList, SEXP expectList, SEXP computeList,
		SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options)
{
	try {
		return omxBackend2(constraints, matList,
				   varList, algList, expectList, computeList,
				   data, intervalList, checkpointList, options);
	} catch( std::exception& __ex__ ) {
		exception_to_try_Rf_error( __ex__ );
	} catch(...) {
		string_to_try_Rf_error( "c++ exception (unknown reason)" );
	}
}

