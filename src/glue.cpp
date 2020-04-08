/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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

#include "omxDefines.h"
#include <R_ext/Rdynload.h>

#include "glue.h"
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
#include "omxCsolnp.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "EnableWarnings.h"
#include "omxSadmvnWrapper.h"

void loadCharVecFromR(const char *context, SEXP names, std::vector<const char *> &dest)
{
	if (!Rf_isNull(names) && !Rf_isString(names)) {
		Rf_warning("%s: found type '%s' instead of a character vector (ignored)",
			   context, Rf_type2char(TYPEOF(names)));
	} else {
		int nlen = Rf_length(names);
		dest.resize(nlen);
		for (int nx=0; nx < nlen; ++nx) {
			dest[nx] = CHAR(STRING_ELT(names, nx));
		}
	}
}

void markAsDataFrame(SEXP list, int rows)
{
	if (rows >= 0) {
		// Use the special form c(NA, #rows) to avoid creating actual
		// rownames.
		ProtectedSEXP rownames(Rf_allocVector(INTSXP, 2));
		INTEGER(rownames)[0] = NA_INTEGER;
		INTEGER(rownames)[1] = rows;
		Rf_setAttrib(list, R_RowNamesSymbol, rownames);
	}

	ProtectedSEXP classes(Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("data.frame"));
	Rf_setAttrib(list, R_ClassSymbol, classes);
}

SEXP makeFactor(SEXP vec, int levels, const char **labels)
{
	bool debug = false;
	Rf_protect(vec);
	if (debug) {
		int len = Rf_length(vec);
		int *data = INTEGER(vec);
		for (int xx=0; xx < len; ++xx) data[xx] = 1;
	}

	SEXP classes;
	Rf_protect(classes = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(classes, 0, Rf_mkChar("factor"));
	Rf_setAttrib(vec, R_ClassSymbol, classes);

	SEXP Rlev;
	Rf_protect(Rlev = Rf_allocVector(STRSXP, levels));
	for (int lx=0; lx < levels; ++lx) {
		SET_STRING_ELT(Rlev, lx, Rf_mkChar(labels[lx]));
	}

	Rf_setAttrib(vec, Rf_install("levels"), Rlev);
	return vec;
}

static SEXP do_logm_eigen(SEXP x)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;

    if (!Rf_isNumeric(x) || !Rf_isMatrix(x)) mxThrow("invalid argument");

    dims = Rf_getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m) mxThrow("non-square matrix");
    if (n == 0) return(Rf_allocVector(REALSXP, 0));

    ScopedProtect p1(z, Rf_allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    logm_eigen(n, rx, rz);

    Rf_setAttrib(z, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));

    return z;
}

static SEXP do_expm_eigen(SEXP x)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;

    if (!Rf_isNumeric(x) || !Rf_isMatrix(x)) mxThrow("invalid argument");

    dims = Rf_getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m) mxThrow("non-square matrix");
    if (n == 0) return(Rf_allocVector(REALSXP, 0));

    ScopedProtect(z, Rf_allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    expm_eigen(n, rx, rz);

    Rf_setAttrib(z, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));

    return z;
}

SEXP dtmvnorm_marginal2(SEXP Rxq, SEXP Rxr, SEXP Rq, SEXP Rr,
			SEXP Rsigma, SEXP Rlower, SEXP Rupper)
{
	using Eigen::Map;
	using Eigen::VectorXd;
	using Eigen::MatrixXd;
	using Rcpp::as;
	const Map<VectorXd> xq(as<Map<VectorXd> >(Rxq));
	const Map<VectorXd> xr(as<Map<VectorXd> >(Rxr));
	int qq = Rf_asInteger(Rq) - 1;
	int rr = Rf_asInteger(Rr) - 1;
	const Map<MatrixXd> sigma(as<Map<MatrixXd> >(Rsigma));
	const Map<VectorXd> lower(as<Map<VectorXd> >(Rlower));
	const Map<VectorXd> upper(as<Map<VectorXd> >(Rupper));
	VectorXd density(4);

	omxState *globalState = new omxState;
	FitContext *fc = new FitContext(globalState);
	_dtmvnorm_marginal2(fc, NA_REAL, xq, xr, qq, rr, sigma, lower, upper, density);
	delete fc;
	delete globalState;

	return Rcpp::wrap(density);
}

SEXP dtmvnorm_marginal(SEXP Rxn, SEXP Rn, SEXP Rsigma, SEXP Rlower, SEXP Rupper)
{
	using Eigen::Map;
	using Eigen::VectorXd;
	using Eigen::MatrixXd;
	using Rcpp::as;
	const Map<VectorXd> xn(as<Map<VectorXd> >(Rxn));
	int nn = Rf_asInteger(Rn) - 1;
	const Map<MatrixXd> sigma(as<Map<MatrixXd> >(Rsigma));
	const Map<VectorXd> lower(as<Map<VectorXd> >(Rlower));
	const Map<VectorXd> upper(as<Map<VectorXd> >(Rupper));
	VectorXd density(2);

	omxState *globalState = new omxState;
	FitContext *fc = new FitContext(globalState);
	_dtmvnorm_marginal(fc, NA_REAL, xn, nn, sigma, lower, upper, density);
	delete fc;
	delete globalState;

	return Rcpp::wrap(density);
}

SEXP mtmvnorm(SEXP Rsigma, SEXP Rlower, SEXP Rupper)
{
	using Eigen::Map;
	using Eigen::VectorXd;
	using Eigen::MatrixXd;
	using Rcpp::as;
	const Map<MatrixXd> sigma(as<Map<MatrixXd> >(Rsigma));
	const Map<VectorXd> fullLower(as<Map<VectorXd> >(Rlower));
	const Map<VectorXd> fullUpper(as<Map<VectorXd> >(Rupper));

	VectorXd tmean;
	MatrixXd tcov;
	omxState *globalState = new omxState;
	FitContext *fc = new FitContext(globalState);
	_mtmvnorm(fc, NA_REAL, sigma, fullLower, fullUpper, tmean, tcov);
	delete fc;
	delete globalState;

	ProtectAutoBalanceDoodad mpi;
	MxRList result;
	result.add("tmean", Rcpp::wrap(tmean));
	result.add("tvar", Rcpp::wrap(tcov));
	return result.asR();
}

void omxSadmvnWrapper(FitContext *fc, int numVars, 
	double *corList, double *lThresh, double *uThresh, int *Infin, double *likelihood, int *inform)
{
	// Eigen::Map< Eigen::ArrayXd > elt(lThresh, numVars);
	// Eigen::Map< Eigen::ArrayXd > eut(uThresh, numVars);
	// mxPrintMat("lower", elt);
	// mxPrintMat("upper", eut);

    // SADMVN calls Alan Genz's sadmvn.f--see appropriate file for licensing info.
   	// TODO: Check with Genz: should we be using sadmvn or sadmvn?
   	// Parameters are:
   	// 	N 		int			# of vars
   	//	Lower	double*		Array of lower bounds
   	//	Upper	double*		Array of upper bounds
   	//	Infin	int*		Array of flags: 0 = (-Inf, upper] 1 = [lower, Inf), 2 = [lower, upper]
   	//	Correl	double*		Array of correlation coeffs: in row-major lower triangular order
   	//	MaxPts	int			Maximum # of function values (use 1000*N or 1000*N*N)
   	//	Abseps	double		Absolute stop tolerance.  Yick.
   	//	Releps	double		Relative stop tolerance.  Use EPSILON.
   	//	Error	&double		On return: absolute real stop, 99% confidence
   	//	Value	&double		On return: evaluated value
   	//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
   	double Error;
	double absEps = 0.0;  // Inappropriate to use for our application
	double relEps = Global->relEps;
	int MaxPts = Global->calcNumIntegrationPoints(numVars);
	int fortranThreadId = omx_absolute_thread_num() + 1;
   	/* FOR DEBUGGING PURPOSES */
    /*	numVars = 2;
   	lThresh[0] = -2;
   	uThresh[0] = -1.636364;
   	Infin[0] = 2;
   	lThresh[1] = 0;
   	uThresh[1] = 0;
   	Infin[1] = 0;
   	smallCor[0] = 1.0; smallCor[1] = 0; smallCor[2] = 1.0; */
   	F77_CALL(sadmvn)(&numVars, lThresh, uThresh, Infin, corList, &MaxPts, 
		&absEps, &relEps, &Error, likelihood, inform, &fortranThreadId);

	if (fc) fc->recordOrdinalRelativeError(Error / *likelihood);

   	if (0) {
   		char infinCodes[3][20];
   		strcpy(infinCodes[0], "(-INF, upper]");
   		strcpy(infinCodes[1], "[lower, INF)");
   		strcpy(infinCodes[2], "[lower, upper]");
   		mxLog("Input to sadmvn is (%d rows):", numVars); //:::DEBUG:::
		for(int i = 0; i < numVars; i++) {
			mxLog("Row %d: %f, %f, %d(%s)", i, lThresh[i], uThresh[i], Infin[i], infinCodes[Infin[i]]);
		}

		mxLog("Cor: (Lower %d x %d):", numVars, numVars); //:::DEBUG:::
		for(int i = 0; i < numVars*(numVars-1)/2; i++) {
			// mxLog("Row %d of Cor: ", i);
			// for(int j = 0; j < i; j++)
			mxLog(" %f", corList[i]); // (i*(i-1)/2) + j]);
			// mxLog("");
		}
	}
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

static SEXP testMxLog(SEXP Rstr) {
	mxLog("%s", CHAR(Rf_asChar(Rstr)));
	return Rf_ScalarLogical(1);
}

static int untitledCounter = 0;

static SEXP untitledNumberReset() {
	untitledCounter = 0;
	return Rf_ScalarLogical(1);
}

static SEXP loadedHook() {
	void ComputeLoadDataLoadedHook();
	untitledCounter = 0;
	ComputeLoadDataLoadedHook();
	return Rf_ScalarLogical(1);
}

static SEXP untitledNumber() {
	return Rf_ScalarInteger(++untitledCounter);
}

void string_to_Rf_error( const std::string& str )
{
	Rf_error("%s", str.c_str());
}

void exception_to_Rf_error( const std::exception& ex )
{
	string_to_Rf_error(ex.what());
}

SEXP MxRList::asR()
{
	// detect duplicate keys? TODO
	SEXP names, ans;
	int len = size();
	Rf_protect(names = Rf_allocVector(STRSXP, len));
	Rf_protect(ans = Rf_allocVector(VECSXP, len));
	for (int lx=0; lx < len; ++lx) {
		SEXP p1 = (*this)[lx].first;
		SEXP p2 = (*this)[lx].second;
		if (!p1 || !p2) mxThrow("Attempt to return NULL pointer to R");
		SET_STRING_ELT(names, lx, p1);
		SET_VECTOR_ELT(ans,   lx, p2);
	}
	Rf_namesgets(ans, names);
	return ans;
}

void friendlyStringToLogical(const char *key, SEXP rawValue, int *out)
{
	if (TYPEOF(rawValue) == LGLSXP) {
		*out = Rf_asLogical(rawValue);
		return;
	}

	const char *str = CHAR(Rf_asChar(rawValue));
	if (TYPEOF(rawValue) != STRSXP) {
		Rf_warning("Not sure how to interpret '%s' (type %s) for mxOption '%s'",
			   str, Rf_type2char(TYPEOF(rawValue)), key);
		return;
	}

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
static void readOpts(SEXP options, int *numThreads, int *analyticGradients)
{
		int numOptions = Rf_length(options);
		SEXP optionNames;
		Rf_protect(optionNames = Rf_getAttrib(options, R_NamesSymbol));
		for(int i = 0; i < numOptions; i++) {
			const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
			SEXP rawValue = VECTOR_ELT(options, i);
			const char *nextOptionValue = CHAR(Rf_asChar(rawValue));
			if(matchCaseInsensitive(nextOptionName, "Analytic Gradients")) {
				friendlyStringToLogical(nextOptionName, rawValue, analyticGradients);
			} else if(matchCaseInsensitive(nextOptionName, "loglikelihoodScale")) {
				Global->llScale = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "debug protect stack")) {
				friendlyStringToLogical(nextOptionName, rawValue, &Global->debugProtectStack);
			} else if(matchCaseInsensitive(nextOptionName, "Number of Threads")) {
#ifdef _OPENMP
				*numThreads = atoi(nextOptionValue);
				if (*numThreads < 1) {
					Rf_warning("Computation will be too slow with %d threads; using 1 thread instead", *numThreads);
					*numThreads = 1;
				}
				char *ont = getenv("OMP_NUM_THREADS");
				if (ont && *numThreads > atoi(ont)) {
					mxThrow("I'm confused! %d threads requested. "
									"Either request fewer threads, i.e., %s, in the mxOption() "
									"statement, or submit your batch job with OMP_NUM_THREADS "
									"environment varible set to %d.  This env variable may be "
									"controlled by PBS’s -ncpus argument, "
									"or similar on other batch systems.",
									*numThreads, ont, *numThreads);
				}
#endif
			} else if(matchCaseInsensitive(nextOptionName, "Parallel diagnostics")) {
				friendlyStringToLogical(nextOptionName, rawValue, &Global->parallelDiag);
			} else if(matchCaseInsensitive(nextOptionName, "maxOrdinalPerBlock")) {
				Global->maxOrdinalPerBlock = atoi(nextOptionValue);
				if (Global->maxOrdinalPerBlock < 1) mxThrow("maxOrdinalPerBlock must be strictly positive");
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsA")) {
				Global->maxptsa = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsB")) {
				Global->maxptsb = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsC")) {
				Global->maxptsc = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsD")) {
				Global->maxptsd = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnMaxPointsE")) {
				Global->maxptse = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "mvnAbsEps")) {
				double absEps = atof(nextOptionValue);
				if (absEps != 0.0) Rf_warning("mxOption mvnAbsEps ignored");
			} else if(matchCaseInsensitive(nextOptionName, "mvnRelEps")) {
				Global->relEps = atof(nextOptionValue);
			} else if(matchCaseInsensitive(nextOptionName, "maxStackDepth")) {
				Global->maxStackDepth = atoi(nextOptionValue);
			} else if (matchCaseInsensitive(nextOptionName, "Feasibility tolerance")) {
				Global->feasibilityTolerance = atof(nextOptionValue);
			} else if (matchCaseInsensitive(nextOptionName, "max minutes")) {
				Global->maxSeconds = nearbyint(atof(nextOptionValue) * 60);
			} else if (matchCaseInsensitive(nextOptionName, "Optimality tolerance")) {
				Global->optimalityTolerance = atof(nextOptionValue);
			} else if (matchCaseInsensitive(nextOptionName, "Major iterations")) {
				Global->majorIterations = atoi(nextOptionValue);
			} else if (matchCaseInsensitive(nextOptionName, "Intervals")) {
				Global->intervals = Rf_asLogical(VECTOR_ELT(options, i));
			} else if (matchCaseInsensitive(nextOptionName, "RAM Inverse Optimization")) {
				friendlyStringToLogical(nextOptionName, rawValue, &Global->RAMInverseOpt);
				//mxLog("inv opt = %s %d", nextOptionValue, Global->RAMInverseOpt);
			} else if (matchCaseInsensitive(nextOptionName, "RAM Max Depth")) {
				if (strEQ(nextOptionValue, "NA")) {
					Global->RAMMaxDepth = NA_INTEGER;
				} else {
					Global->RAMMaxDepth = atoi(nextOptionValue);
				}
				//mxLog("ram max depth = %s %d", nextOptionValue, Global->RAMMaxDepth);
			} else {
				// ignore
			}
		}
}

/* Main functions */
SEXP omxCallAlgebra2(SEXP matList, SEXP algNum, SEXP options) {

	ProtectAutoBalanceDoodad protectManager;

	if(OMX_DEBUG) { mxLog("-----------------------------------------------------------------------");}
	if(OMX_DEBUG) { mxLog("Explicit call to algebra %d.", INTEGER(algNum)[0]);}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;

	FitContext::setRFitFunction(NULL);
	Global = new omxGlobal;

	omxState *globalState = new omxState;

	readOpts(options, &Global->numThreads, &Global->analyticGradients);

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
		mxThrow("Failed to build algebra");
	}

	if(OMX_DEBUG) {mxLog("Completed Algebras and Matrices.  Beginning Initial Compute.");}

	omxRecompute(algebra, NULL);

	Rf_protect(ans = Rf_allocMatrix(REALSXP, algebra->rows, algebra->cols));
	for(l = 0; l < algebra->rows; l++)
		for(j = 0; j < algebra->cols; j++)
			REAL(ans)[j * algebra->rows + l] =
				omxMatrixElement(algebra, l, j);

	if(OMX_DEBUG) { mxLog("All Algebras complete."); }

	const char *bads = Global->getBads();

	omxFreeMatrix(algebra);
	delete globalState;
	delete Global;

	if (bads) mxThrow("%s", bads);

	return ans;
}

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options)
{
	try {
		return omxCallAlgebra2(matList, algNum, options);
	} catch( std::exception& __ex__ ) {
		exception_to_Rf_error( __ex__ );
	} catch(...) {
		string_to_Rf_error( "c++ exception (unknown reason)" );
	}
}

static double internalToUserBound(double val, double inf)
{
	if (val == inf) return NA_REAL;
	return val;
}

SEXP omxBackend2(SEXP constraints, SEXP matList,
		 SEXP varList, SEXP algList, SEXP expectList, SEXP computeList,
		 SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options,
		 SEXP defvars, bool silent)
{
	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these.  For now, all checking is done on the front-end. */
//	if(!isVector(matList)) stop ("matList must be a list");
//	if(!isVector(algList)) stop ("algList must be a list");

	ProtectAutoBalanceDoodad protectManager;

	FitContext::setRFitFunction(NULL);
	Global = new omxGlobal;
	Global->silent = silent;
	Global->mpi = &protectManager;

	/* Create new omxState for current state storage and initialize it. */
	omxState *globalState = new omxState;

	readOpts(options, &Global->numThreads, &Global->analyticGradients);
#if HAS_NPSOL
	omxSetNPSOLOpts(options);
#endif

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessMxDataEntities(data, defvars);
    
	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessMxExpectationEntities(expectList);

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessMxMatrixEntities(matList);

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessFreeVarList(varList);
	FitContext *fc = new FitContext(globalState);
	Global->topFc = fc;
	fc->copyParamToModelClean();

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessMxAlgebraEntities(algList);

	/* Process Matrix and Algebra Population Function */
	/*
	  Each matrix is a list containing a matrix and the other matrices/algebras that are
	  populated into it at each iteration.  The first element is already processed, above.
	  The rest of the list will be processed here.
	*/
	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	for(int j = 0; j < Rf_length(matList); j++) {
		ProtectedSEXP nextLoc(VECTOR_ELT(matList, j));		// This is the matrix + populations
		globalState->matrixList[j]->omxProcessMatrixPopulationList(nextLoc);
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxInitialMatrixAlgebraCompute(fc);
	if (isErrorRaised()) {
		mxThrow("%s", Global->getBads());
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxCompleteMxExpectationEntities();

	for (int dx=0; dx < (int) globalState->dataList.size(); ++dx) {
		globalState->dataList[dx]->connectDynamicData(globalState);
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxCompleteMxFitFunction(algList, fc);

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	Global->omxProcessMxComputeEntities(computeList, globalState);

	// Nothing depend on constraints so we can process them last.
	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	globalState->omxProcessConstraints(constraints, fc);

	if (isErrorRaised()) {
		mxThrow("%s", Global->getBads());
	}

	globalState->loadDefinitionVariables(true);

	omxCompute *topCompute = NULL;
	if (Global->computeList.size()) {
		topCompute = Global->computeList[0];
		Global->ComputePersist = topCompute->isPersist();
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	Global->omxProcessConfidenceIntervals(intervalList, globalState);

	omxProcessCheckpointOptions(checkpointList);

	Global->cacheDependencies(globalState);

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	if (protectManager.getDepth() > Global->maxStackDepth) {
		mxThrow("Protection stack too large; report this problem to the OpenMx forum");
	}

	if (globalState->getWantStage() != FF_COMPUTE_INITIAL_FIT) {
		mxThrow("globalState->getWantStage() != FF_COMPUTE_INITIAL_FIT");
	}
	globalState->setWantStage(FF_COMPUTE_FIT);

	if (topCompute && !isErrorRaised()) {
		topCompute->compute(fc);
		if (Global->computeLoopContext.size() != 0) mxThrow("computeLoopContext imbalance");

		if (fc->wanted & FF_COMPUTE_FIT) {
			if (!std::isfinite(fc->fit)) {
				std::string diag = fc->getIterationError();
				if (diag.size()) {
					omxRaiseErrorf("fit is not finite (%s)", diag.c_str());
				} else if (fc->getInform() == INFORM_CONVERGED_OPTIMUM ||
					   fc->getInform() == INFORM_UNCONVERGED_OPTIMUM) {
					omxRaiseErrorf("fit is not finite");
				}
			} else if (fc->skippedRows) {
				Rf_warning("%d rows obtained probability of exactly zero; "
					   "You may wish to try again with better starting values.", fc->skippedRows);
			}
		}
	}

	SEXP evaluations;
	Rf_protect(evaluations = Rf_allocVector(REALSXP,1));

	REAL(evaluations)[0] = fc->getGlobalComputeCount();

	MxRList result;

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	if (!isErrorRaisedIgnTime()) globalState->omxExportResults(&result, fc);

	if (topCompute && !isErrorRaisedIgnTime()) {
		LocalComputeResult cResult;
		topCompute->collectResults(fc, &cResult, &result);

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

		if (fc->wanted & FF_COMPUTE_FIT) {
			result.add("fit", Rf_ScalarReal(fc->fit));
			if (fc->fitUnits) {
				SEXP units;
				Rf_protect(units = Rf_allocVector(STRSXP, 1));
				SET_STRING_ELT(units, 0, Rf_mkChar(fitUnitsToName(fc->fitUnits)));
				result.add("fitUnits", units);
			}
			result.add("Minus2LogLikelihood", Rf_ScalarReal(fc->fit));
			result.add("maxRelativeOrdinalError",
				   Rf_ScalarReal(fc->getOrdinalRelativeError()));
		}
		if (fc->wanted & FF_COMPUTE_BESTFIT) {
			result.add("minimum", Rf_ScalarReal(fc->fit));
		}

		FreeVarGroup *varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
		int numFree = int(varGroup->vars.size());
		if (numFree) {
			SEXP estimate;
			Rf_protect(estimate = Rf_allocVector(REALSXP, numFree));
			memcpy(REAL(estimate), fc->est, sizeof(double)*numFree);
			result.add("estimate", estimate);

			if (Global->boundsUpdated) {
				MxRList bret;
				SEXP Rlb = Rf_allocVector(REALSXP, numFree);
				bret.add("l", Rlb);
				SEXP Rub = Rf_allocVector(REALSXP, numFree);
				bret.add("u", Rub);
				double *lb = REAL(Rlb);
				double *ub = REAL(Rub);
				for(int px = 0; px < numFree; px++) {
					lb[px] = internalToUserBound(varGroup->vars[px]->lbound, NEG_INF);
					ub[px] = internalToUserBound(varGroup->vars[px]->ubound, INF);
				}
				result.add("bounds", bret.asR());
			}
			if (fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
				result.add("infoDefinite", Rf_ScalarLogical(fc->infoDefinite));
				result.add("conditionNumber", Rf_ScalarReal(fc->infoCondNum));
			}
		}
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	MxRList backwardCompatStatus;
	backwardCompatStatus.add("code", Rf_ScalarInteger(fc->getInform()));
	backwardCompatStatus.add("status", Rf_ScalarInteger(-isErrorRaisedIgnTime()));

	if (Global->getBads()) {
		SEXP msg;
		Rf_protect(msg = Rf_allocVector(STRSXP, 1));
		SET_STRING_ELT(msg, 0, Rf_mkChar(Global->getBads()));
		result.add("error", msg);
		backwardCompatStatus.add("statusMsg", msg);
	}

	result.add("status", backwardCompatStatus.asR());
	result.add("iterations", Rf_ScalarInteger(fc->iterations));
	result.add("evaluations", evaluations);

	// Data are not modified and not copied. The same memory
	// is shared across all instances of state.
	// NOTE: This may need to change for MxDataDynamic
	for(size_t dx = 0; dx < globalState->dataList.size(); dx++) {
		omxFreeData(globalState->dataList[dx]);
	}

	if (Global->debugProtectStack) mxLog("Protect depth at line %d: %d", __LINE__, protectManager.getDepth());
	delete Global;

	return result.asR();
}

static SEXP omxBackend(SEXP constraints, SEXP matList,
		SEXP varList, SEXP algList, SEXP expectList, SEXP computeList,
		SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options,
		       SEXP defvars, SEXP Rsilent)
{
	try {
		return omxBackend2(constraints, matList,
				   varList, algList, expectList, computeList,
				   data, intervalList, checkpointList, options, defvars,
				   Rf_asLogical(Rsilent));
	} catch( std::exception& __ex__ ) {
		exception_to_Rf_error( __ex__ );
	} catch(...) {
		string_to_Rf_error( "c++ exception (unknown reason)" );
	}
}

static SEXP testEigenDebug()
{
	Eigen::VectorXd v1(2);
	if (!std::isfinite(v1[0]) && !std::isfinite(v1[1])) {
		Eigen::VectorXd v2(3);
		Eigen::VectorXd v(3);
		v = v1+v2;  // will call abort() unless EIGEN_NO_DEBUG defined
	}
	return Rf_ScalarLogical(false);
}

static R_CallMethodDef callMethods[] = {
	{"backend", (DL_FUNC) omxBackend, 12},
	{"callAlgebra", (DL_FUNC) omxCallAlgebra, 3},
	{"Dmvnorm_wrapper", (DL_FUNC) dmvnorm_wrapper, 3},
	{"hasNPSOL_wrapper", (DL_FUNC) has_NPSOL, 0},
	{"sparseInvert_wrapper", (DL_FUNC) sparseInvert_wrapper, 1},
	{"hasOpenMP_wrapper", (DL_FUNC) has_openmp, 0},
	{"do_logm_eigen", (DL_FUNC) &do_logm_eigen, 1},
	{"do_expm_eigen", (DL_FUNC) &do_expm_eigen, 1},
	{"Log_wrapper", (DL_FUNC) &testMxLog, 1},
	{"untitledNumberReset", (DL_FUNC) &untitledNumberReset, 0},
	{"untitledNumber", (DL_FUNC) &untitledNumber, 0},
	{".EigenDebuggingEnabled", (DL_FUNC) testEigenDebug, 0},
	{".dtmvnorm.marginal", (DL_FUNC) dtmvnorm_marginal, 5},
	{".dtmvnorm.marginal2", (DL_FUNC) dtmvnorm_marginal2, 7},
	{".mtmvnorm", (DL_FUNC) mtmvnorm, 3},
	{".enableMxLog", (DL_FUNC) &enableMxLog, 0},
	{".OpenMxLoaded", (DL_FUNC) &loadedHook, 0},
	{NULL, NULL, 0}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_OpenMx(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_RegisterCCallable("OpenMx", "AddLoadDataProvider", (DL_FUNC) AddLoadDataProvider);
	R_forceSymbols(info, TRUE);

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

