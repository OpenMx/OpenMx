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

static SEXP has_NPSOL()
{ return Rf_ScalarLogical(HAS_NPSOL); }

static R_CallMethodDef callMethods[] = {
	{"backend", (DL_FUNC) omxBackend, 10},
	{"callAlgebra", (DL_FUNC) omxCallAlgebra, 3},
	{"findIdenticalRowsData", (DL_FUNC) findIdenticalRowsData, 5},
	{"Dmvnorm_wrapper", (DL_FUNC) dmvnorm_wrapper, 3},
	{"hasNPSOL_wrapper", (DL_FUNC) has_NPSOL, 0},
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

/* Main functions */
SEXP omxCallAlgebra2(SEXP matList, SEXP algNum, SEXP) {

	omxManageProtectInsanity protectManager;

	if(OMX_DEBUG) { mxLog("-----------------------------------------------------------------------");}
	if(OMX_DEBUG) { mxLog("Explicit call to algebra %d.", INTEGER(algNum)[0]);}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;
	char output[MAX_STRING_LEN];

	FitContext::setRFitFunction(NULL);
	Global = new omxGlobal;

	globalState = new omxState;
	omxInitState(globalState);

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { mxLog("Processing %d matrix(ces).", Rf_length(matList));}

	omxMatrix *args[Rf_length(matList)];
	for(k = 0; k < Rf_length(matList); k++) {
		Rf_protect(nextMat = VECTOR_ELT(matList, k));	// This is the matrix + populations
		args[k] = omxNewMatrixFromRPrimitive(nextMat, globalState, 1, - k - 1);
		globalState->matrixList.push_back(args[k]);
		if(OMX_DEBUG) {
			mxLog("Matrix[%d] initialized (%d x %d)",
				k, globalState->matrixList[k]->rows, globalState->matrixList[k]->cols);
		}
	}

	algebra = omxNewAlgebraFromOperatorAndArgs(algebraNum, args, Rf_length(matList), globalState);

	if(algebra==NULL) {
		Rf_error(globalState->statusMsg);
	}

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}
	omxStateNextEvaluation(globalState);

	omxRecompute(algebra);

	Rf_protect(ans = Rf_allocMatrix(REALSXP, algebra->rows, algebra->cols));
	for(l = 0; l < algebra->rows; l++)
		for(j = 0; j < algebra->cols; j++)
			REAL(ans)[j * algebra->rows + l] =
				omxMatrixElement(algebra, l, j);

	if(OMX_DEBUG) { mxLog("All Algebras complete."); }

	output[0] = 0;
	if (isErrorRaised(globalState)) {
		strncpy(output, globalState->statusMsg, MAX_STRING_LEN);
	}

	omxFreeAllMatrixData(algebra);
	omxFreeState(globalState);
	delete Global;

	if(output[0]) Rf_error(output);

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
				*numThreads = atoi(nextOptionValue);
				if (*numThreads < 1) {
					Rf_warning("Computation will be too slow with %d threads; using 1 thread instead", *numThreads);
					*numThreads = 1;
				}
			} else {
				// ignore
			}
		}
		Rf_unprotect(1); // optionNames
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
	omxProcessMxExpectationEntities(expectList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxDataEntities(data);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);
    
	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxMatrixEntities(matList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	std::vector<double> startingValues;
	omxProcessFreeVarList(varList, &startingValues);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxAlgebraEntities(algList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxFitFunction(algList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessMxComputeEntities(computeList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxCompleteMxExpectationEntities();
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxCompleteMxFitFunction(algList);
	if (isErrorRaised(globalState)) Rf_error(globalState->statusMsg);

	// This is the chance to check for matrix
	// conformability, etc.  Any Rf_errors encountered should
	// be reported using R's Rf_error() function, not
	// omxRaiseErrorf.

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxInitialMatrixAlgebraCompute();
	omxResetStatus(globalState);

	/*
	// Fit functions may not have computed anything because want=0.
	// We shouldn't leave them marked clean because, for example,
	// ComputeNumericDeriv will see an invalid reference fit value.
        //
	// Wait, does dirty/clean make sense for fitfunctions?
        //
	for(size_t index = 0; index < globalState->algebraList.size(); ++index) {
		omxMatrix *mat = globalState->algebraList[index];
		if (mat->fitFunction) omxMarkDirty(mat);
	}
	*/

	omxCompute *topCompute = NULL;
	if (Global->computeList.size()) topCompute = Global->computeList[0];

	/* Process Matrix and Algebra Population Function */
	/*
	  Each matrix is a list containing a matrix and the other matrices/algebras that are
	  populated into it at each iteration.  The first element is already processed, above.
	  The rest of the list will be processed here.
	*/
	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	for(int j = 0; j < Rf_length(matList); j++) {
		Rf_protect(nextLoc = VECTOR_ELT(matList, j));		// This is the matrix + populations
		omxProcessMatrixPopulationList(globalState->matrixList[j], nextLoc);
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessConstraints(constraints);

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxProcessConfidenceIntervals(intervalList);

	omxProcessCheckpointOptions(checkpointList);

	for (size_t vg=0; vg < Global->freeGroup.size(); ++vg) {
		Global->freeGroup[vg]->cacheDependencies();
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	FitContext fc(startingValues);

	if (topCompute && !isErrorRaised(globalState)) {
		topCompute->compute(&fc);
	}

	SEXP evaluations;
	Rf_protect(evaluations = Rf_allocVector(REALSXP,2));

	REAL(evaluations)[0] = globalState->computeCount;

	if (topCompute && !isErrorRaised(globalState) && globalState->stale) {
		fc.copyParamToModel(globalState);
	}

	MxRList result;

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	omxExportResults(globalState, &result); 

	REAL(evaluations)[1] = globalState->computeCount;

	double optStatus = NA_REAL;
	if (topCompute && !isErrorRaised(globalState)) {
		LocalComputeResult cResult;
		topCompute->collectResults(&fc, &cResult, &result);
		optStatus = topCompute->getOptimizerStatus();

		if (cResult.size()) {
			SEXP computes;
			Rf_protect(computes = Rf_allocVector(VECSXP, cResult.size() * 2));
			for (size_t cx=0; cx < cResult.size(); ++cx) {
				std::pair<int, MxRList*> &c1 = cResult[cx];
				SET_VECTOR_ELT(computes, cx*2, Rf_ScalarInteger(c1.first));
				SET_VECTOR_ELT(computes, cx*2+1, c1.second->asR());
				delete c1.second;
			}
			result.push_back(std::make_pair("computes", computes));
		}

		if (fc.wanted & FF_COMPUTE_FIT) {
			result.push_back(std::make_pair("fit", Rf_ScalarReal(fc.fit)));
			result.push_back(std::make_pair("Minus2LogLikelihood", Rf_ScalarReal(fc.fit)));
		}
		if (fc.wanted & FF_COMPUTE_BESTFIT) {
			result.push_back(std::make_pair("minimum", Rf_ScalarReal(fc.fit)));
		}

		size_t numFree = Global->freeGroup[FREEVARGROUP_ALL]->vars.size();
		if (numFree) {
			// move other global reporting here TODO

			SEXP estimate;
			Rf_protect(estimate = Rf_allocVector(REALSXP, numFree));
			memcpy(REAL(estimate), fc.est, sizeof(double)*numFree);
			result.push_back(std::make_pair("estimate", estimate));

			if (fc.stderrs) {
				SEXP stdErrors;
				Rf_protect(stdErrors = Rf_allocMatrix(REALSXP, numFree, 1));
				memcpy(REAL(stdErrors), fc.stderrs, sizeof(double) * numFree);
				result.push_back(std::make_pair("standardErrors", stdErrors));
			}
			if (fc.wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
				result.push_back(std::make_pair("infoDefinite",
								Rf_ScalarLogical(fc.infoDefinite)));
				result.push_back(std::make_pair("conditionNumber",
								Rf_ScalarReal(fc.infoCondNum)));
			}
		}
	}

	if(OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	MxRList backwardCompatStatus;
	backwardCompatStatus.push_back(std::make_pair("code", Rf_ScalarReal(optStatus)));
	backwardCompatStatus.push_back(std::make_pair("status", Rf_ScalarInteger(-isErrorRaised(globalState))));

	if (isErrorRaised(globalState)) {
		SEXP msg;
		Rf_protect(msg = Rf_allocVector(STRSXP, 1));
		SET_STRING_ELT(msg, 0, Rf_mkChar(globalState->statusMsg));
		result.push_back(std::make_pair("error", msg));
		backwardCompatStatus.push_back(std::make_pair("statusMsg", msg));
	}

	result.push_back(std::make_pair("status", backwardCompatStatus.asR()));
	result.push_back(std::make_pair("evaluations", evaluations));

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

