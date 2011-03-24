/*
 *  Copyright 2007-2009 The OpenMx Project
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

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include "npsolWrap.h"

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>
#include "omxState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxObjective.h"
#include "omxNPSOLSpecific.h"
#include "omxBackendHelperFunctions.h"
//#include "omxSymbolTable.h"

/* NPSOL-related functions */
extern void F77_SUB(npoptn)(char* string, int length);
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,
							double *bl,	double *bu, void* funcon, void* funobj, int *inform, int *iter, 
							int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,
							double *x, int *iw,	int *leniw, double *w, int *lenw);

/* Objective Function */
void F77_SUB(objectiveFunction)	( int* mode, int* n, double* x, double* f, double* g, int* nstate );				// For general computation
void F77_SUB(limitObjectiveFunction)(	int* mode, int* n, double* x, double* f, double* g, int* nstate );			// For limit computations

/* Constraint Function */
void F77_SUB(constraintFunction)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate);	// Constraints in the style of old Mx

/* Helper functions */
void handleFreeVarList(omxState *os, double* x, int numVars);	// Locates and inserts elements from the optimizer into matrices.
SEXP getListElement(SEXP list, const char *str); 				// Gets the value named str from SEXP list.  From "Writing R Extensions"
SEXP getVar(SEXP str, SEXP env);								// Gets the object named str from environment env.  From "Writing R Extensions"
unsigned short omxEstimateHessian(omxMatrix** matList, int numHessians, double functionPrecision, 
										int r, omxState* currentState, double optimum);

/* Globals for function evaluation */
SEXP RObjFun, RConFun;			// Pointers to the functions NPSOL calls
SEXP env;						// Environment for evaluation and object hunting

/* Made global for objective functions that want them */
int n, nclin, ncnln;			// Number of Free Params, linear, and nonlinear constraints
double f;						// Objective Function Value
double *g;						// Gradient Pointer
double *R, *cJac;				// Hessian (Approx) and Jacobian
int *istate;					// Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))

omxState* currentState;			// Current State of optimization

/* Main functions */
SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options) {

	if(OMX_DEBUG) { Rprintf("-----------------------------------------------------------------------\n");}
	if(OMX_DEBUG) { Rprintf("Explicit call to algebra %d.\n", INTEGER(algNum));}

	int j,k,l;
	omxMatrix* algebra;
	int algebraNum = INTEGER(algNum)[0];
	SEXP ans, nextMat;
	char output[250];
	int errOut = 0;

	/* Create new omxState for current state storage and initialize it. */
	currentState = (omxState*) R_alloc(1, sizeof(omxState));
	omxInitState(currentState);
	currentState->numFreeParams = n;
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", currentState);}

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}
	currentState->numMats = length(matList);
	currentState->matrixList = (omxMatrix**) R_alloc(length(matList), sizeof(omxMatrix*));

	for(k = 0; k < length(matList); k++) {
		PROTECT(nextMat = VECTOR_ELT(matList, k));	// This is the matrix + populations
		currentState->matrixList[k] = omxNewMatrixFromRPrimitive(nextMat, currentState);
		if(OMX_DEBUG) {
			Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n",
				currentState->matrixList[k], currentState->matrixList[k]->rows, currentState->matrixList[k]->cols);
		}
		UNPROTECT(1); // nextMat
	}

	algebra = omxNewAlgebraFromOperatorAndArgs(algebraNum, currentState->matrixList, currentState->numMats, currentState);

	if(algebra==NULL) {
		error(currentState->statusMsg);
	}

	if(OMX_DEBUG) {Rprintf("Completed Algebras and Matrices.  Beginning Initial Compute.\n");}
	omxStateNextEvaluation(currentState);

	omxRecompute(algebra);

	PROTECT(ans = allocMatrix(REALSXP, algebra->rows, algebra->cols));
	for(l = 0; l < algebra->rows; l++)
		for(j = 0; j < algebra->cols; j++)
			REAL(ans)[j * algebra->rows + l] =
				omxMatrixElement(algebra, l, j);

	UNPROTECT(1);	/* algebra */

	if(OMX_DEBUG) { Rprintf("All Algebras complete.\n"); }

	if(currentState->statusCode != 0) {
		errOut = currentState->statusCode;
		strncpy(output, currentState->statusMsg, 250);
	}

	omxFreeState(currentState);

	if(errOut != 0) {
		error(output);
	}

	return ans;
}

SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state) {

	/* NPSOL Arguments */
	void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);

	int ldA, ldJ, ldR, inform, iter, leniw, lenw; // nclin, ncnln
	int *iw = NULL; // , istate;

//	double f;
	double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *x = NULL, *w=NULL; //  *g, *R, *cJac,
	double *est, *grad, *hess;

	/* Helpful variables */

	char optionCharArray[250] = "";			// For setting options

	int j, k, l, m;					// Index Vars

	int nctotl, nlinwid, nlnwid;	// Helpful side variables.
	
    int errOut = 0;                 // Error state: Clear

	SEXP nextLoc, nextMat, nextVar;

	omxMatrix** calculateHessians = NULL;
	int calculateStdErrors = FALSE;
	int numHessians = 0;
	int ciMaxIterations = 5;
	int disableOptimizer = 0;

	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these.  For now, all checking is done on the front-end. */
//	if(!isVector(startVals)) error ("startVals must be a vector");
//	if(!isVector(matList)) error ("matList must be a list");
//	if(!isVector(algList)) error ("algList must be a list");

	n = length(startVals);

	/* Create new omxState for current state storage and initialize it. */
	currentState = (omxState*) R_alloc(1, sizeof(omxState));
	omxInitState(currentState);
	currentState->numFreeParams = n;
	if(OMX_DEBUG) { Rprintf("Created state object at 0x%x.\n", currentState);}

	/* Retrieve Data Objects */
	if(!errOut) errOut = omxProcessMxDataEntities(data);
    
	/* Retrieve All Matrices From the MatList */
	if(!errOut) errOut = omxProcessMxMatrixEntities(matList);

	/* Process Algebras Here */
	if(!errOut) errOut = omxProcessMxAlgebraEntities(algList);

	/* Initial Matrix and Algebra Calculations */
	if(!errOut) errOut = omxInitialMatrixAlgebraCompute();
	
	/* Process Objective Function */
	if(!errOut) errOut = omxProcessObjectiveFunction(objective, &n);

	// TODO: Make calculateHessians an option instead.

  if(!errOut) {	// In the event of an initialization error, skip all this.

	/* Process Matrix and Algebra Population Function */
	/*
	 Each matrix is a list containing a matrix and the other matrices/algebras that are
	 populated into it at each iteration.  The first element is already processed, above.
	 The rest of the list will be processed here.
	*/
	for(int j = 0; j < currentState->numMats; j++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, j));		// This is the matrix + populations
		omxProcessMatrixPopulationList(currentState->matrixList[j], nextLoc);
		UNPROTECT(1);
	}

	/* Process Free Var List */
	omxProcessFreeVarList(varList, n);

	/* Processing Constraints */
	if(OMX_VERBOSE) { Rprintf("Processing Constraints.\n");}
	omxMatrix *arg1, *arg2;
	currentState->numConstraints = length(constraints);
	if(OMX_DEBUG) {Rprintf("Found %d constraints.\n", currentState->numConstraints); }
	currentState->conList = (omxConstraint*) R_alloc(currentState->numConstraints, sizeof(omxConstraint));
	ncnln = 0;
	for(k = 0; k < currentState->numConstraints; k++) {
		PROTECT(nextVar = VECTOR_ELT(constraints, k));
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));
		arg1 = omxNewMatrixFromMxIndex(nextLoc, currentState);
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));
		arg2 = omxNewMatrixFromMxIndex(nextLoc, currentState);
		PROTECT(nextLoc = AS_INTEGER(VECTOR_ELT(nextVar, 2)));
		currentState->conList[k].opCode = INTEGER(nextLoc)[0];
		UNPROTECT(4);
		omxMatrix *args[2] = {arg1, arg2};
		currentState->conList[k].result = omxNewAlgebraFromOperatorAndArgs(10, args, 2, currentState); // 10 = binary subtract
		omxRecompute(currentState->conList[k].result);
		currentState->conList[k].size = currentState->conList[k].result->rows * currentState->conList[k].result->cols;
		ncnln += currentState->conList[k].size;
	}
	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d effective constraints.\n", ncnln); }
	funcon = F77_SUB(constraintFunction);

	/* Process Confidence Interval List */
	/*
	 intervalList is a list().  Each element refers to one confidence interval request.
	 Each interval request is a length 5 vector of REAL.
	 The first three elements are the matrixPointer, Row, and Column of the element
	 for which bounds are to be calculated, and are cast to ints here for speed.
	 The last two are the upper and lower boundaries for the confidence space (respectively).
    */
	if(OMX_VERBOSE) { Rprintf("Processing Confidence Interval Requests.\n");}
	currentState->numIntervals = length(intervalList);
	if(OMX_DEBUG) {Rprintf("Found %d requests.\n", currentState->numIntervals); }
	currentState->intervalList = (omxConfidenceInterval*) R_alloc(currentState->numIntervals, sizeof(omxConfidenceInterval));
	for(k = 0; k < currentState->numIntervals; k++) {
		omxConfidenceInterval *oCI = &(currentState->intervalList[k]);
		PROTECT(nextVar = VECTOR_ELT(intervalList, k));
		double* intervalInfo = REAL(nextVar);
		oCI->matrix = omxNewMatrixFromMxIndex( nextVar, currentState);	// Expects an R object
		oCI->row = (int) intervalInfo[1];		// Cast to int in C to save memory/Protection ops
		oCI->col = (int) intervalInfo[2];		// Cast to int in C to save memory/Protection ops
		oCI->lbound = intervalInfo[3];
		oCI->ubound = intervalInfo[4];
		UNPROTECT(1);
		oCI->max = R_NaReal;					// NAs, in case something goes wrong
		oCI->min = R_NaReal;
	}
	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d intervals requested.\n", currentState->numIntervals); }

	/* Process Checkpoint List */
	omxProcessCheckpointOptions(checkpointList);

  } // End if(errOut)

	/* Set up Optimization Memory Allocations */
	if(n == 0) {			// Special Case for the evaluation-only condition

		if(OMX_DEBUG) { Rprintf("No free parameters.  Avoiding Optimizer Entirely.\n"); }
		int mode = 0, nstate = -1;
		f = 0;
		double* x = NULL, *g = NULL;

		if(currentState->objectiveMatrix != NULL) {
			F77_SUB(objectiveFunction)(&mode, &n, x, &f, g, &nstate);
		};
		numHessians = 0;					// No hessian if there's no free params
		currentState->numIntervals = 0; 	// No intervals if there's no free params
		inform = 0;
		iter = 0;

		omxStateNextEvaluation(currentState);	// Advance for a final evaluation.

	} else {
		/* Initialize Scalar Variables. */
		nclin = 0;						// No linear constraints.

		/* Set boundaries and widths. */
		if(nclin <= 0) {
			nclin = 0;					// This should never matter--nclin should always be non-negative.
			nlinwid = 1;				// For memory allocation purposes, nlinwid > 0
		} else {						// nlinwid is  used to calculate ldA, and eventually the size of A.
			nlinwid = nclin;
		}

		if(ncnln <= 0) {
			ncnln = 0;					// This should never matter--ncnln should always be non-negative.
			nlnwid = 1;					// For memory allocation purposes nlnwid > 0
		} else {						// nlnwid is used to calculate ldJ, and eventually the size of J.
			nlnwid = ncnln;
		}

		nctotl = n + nlinwid + nlnwid;

		leniw = 3 * n + nclin + 2 * ncnln;
		lenw = 2 * n * n + n * nclin + 2 * n * ncnln + 20 * n + 11 * nclin + 21 * ncnln;

		ldA = nlinwid;  		// NPSOL specifies this should be either 1 or nclin, whichever is greater
		ldJ = nlnwid; 			// NPSOL specifies this should be either 1 or nclin, whichever is greater
		ldR = n;				// TODO: Test alternative versions of the size of R to see what's best.

	/* Allocate arrays */
		A		= (double*) R_alloc (ldA * n, sizeof ( double )  );
		bl		= (double*) R_alloc ( nctotl, sizeof ( double ) );
		bu		= (double*) R_alloc (nctotl, sizeof ( double ) );
		c		= (double*) R_alloc (nlnwid, sizeof ( double ));
		cJac	= (double*) R_alloc (ldJ * n, sizeof ( double ) );
		clambda = (double*) R_alloc (nctotl, sizeof ( double )  );
		g		= (double*) R_alloc (n, sizeof ( double ) );
		R		= (double*) R_alloc (ldR * n, sizeof ( double ));
		x		= (double*) R_alloc ((n+1), sizeof ( double ));
		w		= (double*) R_alloc (lenw, sizeof ( double ));

		istate	= (int*) R_alloc (nctotl, sizeof ( int ) );
		iw		= (int*) R_alloc (leniw, sizeof ( int ));

	/* Set up actual run */

		/* Set min and max limits */
		for(k = 0; k < n; k++) {
			bl[k] = currentState->freeVarList[k].lbound;				// -Infinity'd be -10^20.
			bu[k] = currentState->freeVarList[k].ubound;				// Infinity would be at 10^20.
		}

		for(; k < n+nclin; k++) {						// At present, nclin == 0.
			bl[k] = NEG_INF; 							// Linear constraints have no bounds.
			bu[k] = INF;								// Because there are no linear constraints.
		}												// But if there were, they would go here.

		for(l = 0; l < currentState->numConstraints; l++) {					// Nonlinear constraints:
			if(OMX_DEBUG) { Rprintf("Constraint %d: ", l);}
			switch(currentState->conList[l].opCode) {
				case 0:									// Less than: Must be strictly less than 0.
					if(OMX_DEBUG) { Rprintf("Bounded at (-INF, 0).\n");}
					for(m = 0; m < currentState->conList[l].size; m++) {
						bl[k] = NEG_INF;
						bu[k] = -0.0;
						k++;
					}
					break;
				case 1:									// Equal: Must be roughly equal to 0.
					if(OMX_DEBUG) { Rprintf("Bounded at (-0, 0).\n");}
					for(m = 0; m < currentState->conList[l].size; m++) {
						bl[k] = -0.0;
						bu[k] = 0.0;
						k++;
					}
					break;
				case 2:									// Greater than: Must be strictly greater than 0.
					if(OMX_DEBUG) { Rprintf("Bounded at (0, INF).\n");}
					for(m = 0; m < currentState->conList[l].size; m++) {
						if(OMX_DEBUG) { Rprintf("\tBounds set for constraint %d.%d.\n", l, m);}
						bl[k] = 0.0;
						bu[k] = INF;
						k++;
					}
					break;
				default:
					if(OMX_DEBUG) { Rprintf("Bounded at (-INF, INF).\n");}
					for(m = 0; m < currentState->conList[l].size; m++) {
						bl[k] = NEG_INF;
						bu[k] = INF;
						k++;
					}
					break;
			}
		}


		/* Initialize Starting Values */
		if(OMX_VERBOSE) {
			Rprintf("--------------------------\n");
			Rprintf("Starting Values (%d) are:\n", n);
		}
		for(k = 0; k < n; k++) {
			x[k] = REAL(startVals)[k];
			if(x[k] == 0.0) {
				x[k] += .1;
			}
			if(OMX_VERBOSE) { Rprintf("%d: %f\n", k, x[k]); }
		}
		if(OMX_DEBUG) {
			Rprintf("--------------------------\n");
			Rprintf("Setting up optimizer...");
		}

		/* 	Set NPSOL options  (Maybe separate this into a different function?) */
		/* Options That Change The Optimizer */
		int numOptions = length(options);
		SEXP optionNames;
		PROTECT(optionNames = GET_NAMES(options));
		for(int i = 0; i < numOptions; i++) {
			const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
			const char *nextOptionValue = STRING_VALUE(VECTOR_ELT(options, i));
			int lenName = strlen(nextOptionName);
			int lenValue = strlen(nextOptionValue);
			if(matchCaseInsensitive(nextOptionName, lenName, "Calculate Hessian")) {
				if(OMX_DEBUG) { Rprintf("Found hessian option... Value: %s. ", nextOptionValue);};
				if(!matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Enabling explicit hessian calculation.\n");}
					calculateHessians = (omxMatrix**) R_alloc(1, sizeof(omxMatrix*));
					calculateHessians[0] = currentState->objectiveMatrix;		// TODO: move calculateHessians default
					numHessians = 1;
				}
			} else if(matchCaseInsensitive(nextOptionName, lenName, "Standard Errors")) {
				if(OMX_DEBUG) { Rprintf("Found standard error option...Value: %s. ", nextOptionValue);};
				if(!matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Enabling explicit standard error calculation.\n");}
					calculateStdErrors = TRUE;
					if(calculateHessians == NULL) {
						calculateHessians = (omxMatrix**) R_alloc(1, sizeof(omxMatrix*));
						calculateHessians[0] = currentState->objectiveMatrix;
						numHessians = 1;
					}
				}
			} else if(matchCaseInsensitive(nextOptionName, lenName, "CI Max Iterations")) { 
				int newvalue = atoi(nextOptionValue);
				if (newvalue > 0) ciMaxIterations = newvalue;
			} else if(matchCaseInsensitive(nextOptionName, lenName, "useOptimizer")) {
				if(OMX_DEBUG) { Rprintf("Found useOptimizer option...");};	
				if(matchCaseInsensitive(nextOptionValue, lenValue, "No")) {
					if(OMX_DEBUG) { Rprintf("Disabling optimization.\n");}
					disableOptimizer = 1;
				}
			} else {
				sprintf(optionCharArray, "%s %s", nextOptionName, nextOptionValue);
				F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
				if(OMX_DEBUG) { Rprintf("Option %s \n", optionCharArray); }
			}
		}
		UNPROTECT(1); // optionNames

	/*  F77_CALL(npsol)
		(	int *n,					-- Number of variables
			int *nclin,				-- Number of linear constraints
			int *ncnln,				-- Number of nonlinear constraints
			int *ldA,				-- Row dimension of A (Linear Constraints)
			int *ldJ,				-- Row dimension of cJac (Jacobian)
			int *ldR,				-- Row dimension of R (Hessian)
			double *A,				-- Linear Constraints Array A (in Column-major order)
			double *bl,				-- Lower Bounds Array (at least n + nclin + ncnln long)
			double *bu,				-- Upper Bounds Array (at least n + nclin + ncnln long)
			function funcon,		-- Nonlinear constraint function
			function funobj,		-- Objective function
			int *inform,			-- Used to report state.  Need not be initialized.
			int *iter,				-- Used to report number of major iterations performed.  Need not be initialized.
			int *istate,			-- Initial State.  Need not be initialized unless using Warm Start.
			double *c,				-- Array of length ncnln.  Need not be initialized.  Reports nonlinear constraints at final iteration.
			double *cJac,			-- Array of Row-length ldJ.  Unused if ncnln = 0. Generally need not be initialized.
			double *clambda,		-- Array of length n+nclin+ncnln.  Need not be initialized unless using Warm Start. Reports final QP multipliers.
			double *f,				-- Used to report final objective value.  Need not be initialized.
			double *g,				-- Array of length n. Used to report final objective gradient.  Need not be initialized.
			double *R,				-- Array of length ldR.  Need not be intialized unless using Warm Start.
			double *x,				-- Array of length n.  Contains initial solution estimate.
			int *iw,				-- Array of length leniw. Need not be initialized.  Provides workspace.
			int *leniw,				-- Length of iw.  Must be at least 3n + nclin + ncnln.
			double *w,				-- Array of length lenw. Need not be initialized.  Provides workspace.
			int *lenw				-- Length of w.  Must be at least 2n^2 + n*nclin + 2*n*ncnln + 20*n + 11*nclin +21*ncnln
		)

		bl, bu, istate, and clambda are all length n+nclin+ncnln.
			First n elements refer to the vars, in order.
			Next nclin elements refer to bounds on Ax
			Last ncnln elements refer to bounds on c(x)

		All arrays must be in column-major order.

		*/

		if(OMX_DEBUG) {
			Rprintf("Set.\n");
		}

		if (disableOptimizer) {
			int mode = 0, nstate = -1;		
			if(currentState->objectiveMatrix != NULL) {
				F77_SUB(objectiveFunction)(&mode, &n, x, &f, g, &nstate);
			};

			inform = 0;
			iter = 0;

			omxStateNextEvaluation(currentState);	// Advance for a final evaluation.		
		} else {
			F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
							(void*) F77_SUB(objectiveFunction), &inform, &iter, istate, c, cJac,
							clambda, &f, g, R, x, iw, &leniw, w, &lenw);
		}
		if(OMX_DEBUG) { Rprintf("Final Objective Value is: %f.\n", f); }
		
		handleFreeVarList(currentState, x, n);
		
	} // END OF PERFORM OPTIMIZATION CASE

	SEXP minimum, estimate, gradient, hessian, code, status, statusMsg, iterations;
	SEXP evaluations, ans=NULL, names=NULL, algebras, algebra, matrices, optimizer;
	SEXP intervals, NAmat, intervalCodes, calculatedHessian, stdErrors;

	int numReturns = 13;

	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(status = allocVector(VECSXP, 3));
	PROTECT(iterations = NEW_NUMERIC(1));
	PROTECT(evaluations = NEW_NUMERIC(2));
	PROTECT(matrices = NEW_LIST(currentState->numMats));
	PROTECT(algebras = NEW_LIST(currentState->numAlgs));

	/* N-dependent SEXPs */
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(optimizer = allocVector(VECSXP, 2));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));
	PROTECT(calculatedHessian = allocMatrix(REALSXP, n, n));
	PROTECT(stdErrors = allocMatrix(REALSXP, n, 1)); // for optimizer
	PROTECT(names = allocVector(STRSXP, 2)); // for optimizer
	PROTECT(intervals = allocMatrix(REALSXP, currentState->numIntervals, 2)); // for optimizer
	PROTECT(intervalCodes = allocMatrix(INTSXP, currentState->numIntervals, 2)); // for optimizer
	PROTECT(NAmat = allocMatrix(REALSXP, 1, 1)); // In case of missingness
	REAL(NAmat)[0] = R_NaReal;

	/* Store outputs for return */
	if(objective != NULL) {
		REAL(minimum)[0] = f;
		currentState->optimum = f;
	} else {
		REAL(minimum)[0] = R_NaReal;
	}

	est = REAL(estimate);  // Aliases to avoid repeated function calls.
	grad = REAL(gradient);
	hess = REAL(hessian);

	for(k = 0; k < n; k++) {		// Do these with memcpy instead, probably.
		est[k] = x[k];
		grad[k] = g[k];
		for(l = 0; l < n; l++) {			// Save the NPSOL hessian, in case somebody wants it
			hess[k*n + l] = R[k*n + l];		// This is ok, because they're both in Col-Major already.
		}
	}

	omxSaveState(currentState, x, f);		// Keep the current values for the currentState.

	/* Fill in details from the optimizer */
	SET_VECTOR_ELT(optimizer, 0, gradient);
	SET_VECTOR_ELT(optimizer, 1, hessian);

	SET_STRING_ELT(names, 0, mkChar("minimum"));
	SET_STRING_ELT(names, 1, mkChar("estimate"));
	namesgets(optimizer, names);

	REAL(code)[0] = inform;
	REAL(iterations)[0] = iter;
	REAL(evaluations)[0] = currentState->computeCount;

	/* Fill Status code. */
	SET_VECTOR_ELT(status, 0, code);
	PROTECT(code = NEW_NUMERIC(1));
	REAL(code)[0] = currentState->statusCode;
	SET_VECTOR_ELT(status, 1, code);
	PROTECT(statusMsg = allocVector(STRSXP, 1));
	SET_STRING_ELT(statusMsg, 0, mkChar(currentState->statusMsg));
	SET_VECTOR_ELT(status, 2, statusMsg);

	if(numHessians && currentState->optimumStatus >= 0) {		// No hessians or standard errors if the optimum is invalid
		if(currentState->numConstraints == 0) {
			if(OMX_DEBUG) { Rprintf("Calculating Hessian for Objective Function.\n");}
			int gotHessians = omxEstimateHessian(calculateHessians, numHessians, .0001, /*2.0,*/ 4, currentState, currentState->optimum);
			if(gotHessians) {
				if(calculateStdErrors) {
					for(int j = 0; j < numHessians; j++) {		//TODO: Fix Hessian calculation to allow more if requested
						if(OMX_DEBUG) { Rprintf("Calculating Standard Errors for Objective Function.\n");}
						omxObjective* oo = calculateHessians[j]->objective;
						if(oo->getStandardErrorFun != NULL) {
							oo->getStandardErrorFun(oo);
						} else {
							omxCalculateStdErrorFromHessian(2.0, oo);
						}
					}
				}
			} else {
				numHessians = 0;
			}
		} else {
			numHessians = 0;
		}
	}

	/* Likelihood-based Confidence Interval Calculation */
	if(currentState->numIntervals) {
		if(inform == 0 || inform == 1 || inform == 6) {
			if(OMX_DEBUG) { Rprintf("Calculating likelihood-based confidence intervals.\n"); }
			currentState->optimizerState = (omxOptimizerState*) R_alloc(1, sizeof(omxOptimizerState));
			for(int i = 0; i < currentState->numIntervals; i++) {

				memcpy(x, currentState->optimalValues, n * sizeof(double)); // Reset to previous optimum
				currentState->currentInterval = i;
				omxConfidenceInterval *currentCI = &(currentState->intervalList[i]);
				currentCI->lbound += currentState->optimum;			// Convert from offsets to targets
				currentCI->ubound += currentState->optimum;			// Convert from offsets to targets

				/* Set up for the lower bound */
				inform = -1;
				// Number of times to keep trying.
				int cycles = ciMaxIterations;
				double value = INF;
				double objDiff = 1.e-4;		// TODO : Use function precision to determine CI jitter?
				while(inform != 0 && cycles > 0) {
					/* Find lower limit */
					currentCI->calcLower = TRUE;
					F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
									(void*) F77_SUB(limitObjectiveFunction), &inform, &iter, istate, c, cJac,
									clambda, &f, g, R, x, iw, &leniw, w, &lenw);

					currentCI->lCode = inform;
					if(f < value) {
						currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
						value = f;
					}

					if(inform != 0 && OMX_DEBUG) {
						Rprintf("Calculation of lower interval %d failed: Bad inform value of %d\n",
							i, inform);
					}
					cycles--;
					if(inform != 0) {
						unsigned int jitter = TRUE;
						for(int j = 0; j < n; j++) {
							if(fabs(x[j] - currentState->optimalValues[j]) > objDiff){
								jitter = FALSE;
								break;
							}
						}
						if(jitter) {
							for(int j = 0; j < n; j++) {
								x[j] = currentState->optimalValues[j] + objDiff;
							}
						}
					}
				}

				if(OMX_DEBUG) { Rprintf("Found lower bound %d.  Seeking upper.\n", i); }
				// TODO: Repopulate original optimizer state in between CI calculations

				memcpy(x, currentState->optimalValues, n * sizeof(double));

				/* Reset for the upper bound */
				value = INF;
				inform = -1;
				cycles = ciMaxIterations;

				while(inform != 0 && cycles >= 0) {
					/* Find upper limit */
					currentCI->calcLower = FALSE;
					F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
									(void*) F77_SUB(limitObjectiveFunction), &inform, &iter, istate, c, cJac,
									clambda, &f, g, R, x, iw, &leniw, w, &lenw);

					currentCI->uCode = inform;
					if(f < value) {
						currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
					}

					if(inform != 0 && OMX_DEBUG) {
						Rprintf("Calculation of upper interval %d failed: Bad inform value of %d\n",
							i, inform);
					}
					cycles--;
					if(inform != 0) {
						unsigned int jitter = TRUE;
						for(int j = 0; j < n; j++) {
							if(fabs(x[j] - currentState->optimalValues[j]) > objDiff){
								jitter = FALSE;
								break;
							}
						}
						if(jitter) {
							for(int j = 0; j < n; j++) {
								x[j] = currentState->optimalValues[j] + objDiff;
							}
						}
					}
				}
				if(OMX_DEBUG) {Rprintf("Found Upper bound %d.\n", i);}
			}
		} else {					// Improper code. No intervals calculated.
									// TODO: Throw a warning, allow force()
			warning("Not calculating confidence intervals because of error status.");
			if(OMX_DEBUG) {
				Rprintf("Calculation of all intervals failed: Bad inform value of %d", inform);
			}
		}
	}

	handleFreeVarList(currentState, currentState->optimalValues, n);  // Restore to optima for final compute
	if(!errOut) {
		for(k = 0; k < currentState->numMats; k++) {
			if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Matrix %d.\n", k); }
			omxMatrix* nextMatrix = currentState->matrixList[k];
			omxRecompute(nextMatrix);
			PROTECT(nextMat = allocMatrix(REALSXP, nextMatrix->rows, nextMatrix->cols));
			for(l = 0; l < nextMatrix->rows; l++)
				for(j = 0; j < nextMatrix->cols; j++)
					REAL(nextMat)[j * nextMatrix->rows + l] =
						omxMatrixElement(nextMatrix, l, j);
			SET_VECTOR_ELT(matrices, k, nextMat);

			UNPROTECT(1);	/* nextMat */
		}

		for(k = 0; k < currentState->numAlgs; k++) {
			if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Algebra %d.\n", k); }
			omxMatrix* nextAlgebra = currentState->algebraList[k];
			omxRecompute(nextAlgebra);
			PROTECT(algebra = allocMatrix(REALSXP, nextAlgebra->rows, nextAlgebra->cols));
			if (nextAlgebra->objective != NULL && nextAlgebra->objective->populateAttrFun != NULL) {
				nextAlgebra->objective->populateAttrFun(nextAlgebra->objective, algebra);
			}
			for(l = 0; l < nextAlgebra->rows; l++)
				for(j = 0; j < nextAlgebra->cols; j++)
					REAL(algebra)[j * nextAlgebra->rows + l] =
						omxMatrixElement(nextAlgebra, l, j);
			SET_VECTOR_ELT(algebras, k, algebra);

			UNPROTECT(1);	/* algebra */
		}
		if(OMX_DEBUG) { Rprintf("All Algebras complete.\n"); }
	}

	omxMatrix* om = currentState->objectiveMatrix;
	if(om != NULL) {					// In the event of a no-objective run.
		omxObjective* oo = om->objective;
		if(OMX_DEBUG) { Rprintf("Checking for additional objective info.\n"); }

		if(oo != NULL && oo->setFinalReturns != NULL) {
			if(OMX_DEBUG) { Rprintf("Expecting objective Info....");}
			int numEls;
			SEXP oElement;
			omxRListElement* orle = oo->setFinalReturns(oo, &numEls);
			PROTECT(ans = allocVector(VECSXP, numReturns + numEls));
			PROTECT(names = allocVector(STRSXP, numReturns + numEls));
			if(numEls != 0) {
				if(OMX_DEBUG) { Rprintf("Adding %d sets of objective Info....", numEls);}
				for(int i = 0; i < numEls; i++) {
					PROTECT(oElement = allocVector(REALSXP, orle[i].numValues));
					for(int j = 0; j < orle[i].numValues; j++)
						REAL(oElement)[j] = orle[i].values[j];
					SET_STRING_ELT(names, i+numReturns, mkChar(orle[i].label));
					SET_VECTOR_ELT(ans, i+numReturns, oElement);
					UNPROTECT(1); // oElement
				}
			}
		} else {
			PROTECT(ans = allocVector(VECSXP, numReturns));
			PROTECT(names = allocVector(STRSXP, numReturns));
		}
		if(OMX_DEBUG) { Rprintf("Done.\n");}
	} else {
		PROTECT(ans = allocVector(VECSXP, numReturns));
		PROTECT(names = allocVector(STRSXP, numReturns));
	}

	if(numHessians) {
		if(OMX_DEBUG) { Rprintf("Populating hessians for %d objectives.\n", numHessians); }
		for(int j = 0; j < numHessians; j++) {		//TODO: Fix Hessian calculation to allow more if requested
			omxObjective* oo = calculateHessians[j]->objective;
			if(oo->hessian == NULL) {
				if(OMX_DEBUG) { Rprintf("Objective %d has no hessian. Aborting.\n", j);}
				continue;
			}

			if(OMX_DEBUG) { Rprintf("Objective %d has hessian at 0x%x.\n", j, oo->hessian);}

			double* hessian  = REAL(calculatedHessian);
			double* stdError = REAL(stdErrors);
			for(int k = 0; k < n * n; k++) {
				if(OMX_DEBUG) {Rprintf("Populating hessian at %d.\n", k);}
				hessian[k] = oo->hessian[k];		// For expediency, ignore majority for symmetric matrices.
			}
			if(OMX_DEBUG) {Rprintf("Done.\n", k);}
			if(calculateStdErrors) {
				if(oo->stdError == NULL) {
					for(int k = 0; k < n; k++) {
						if(OMX_DEBUG) {Rprintf("Populating NA standard error at %d.\n", k);}
						stdError[k] = R_NaReal;
					}
				} else {
					for(int k = 0; k < n; k++) {
						if(OMX_DEBUG) {Rprintf("Populating standard error at %d.\n", k);}
						stdError[k] = oo->stdError[k];
					}
				}
			}
		}
	}

	if(currentState->numIntervals) {	// Populate CIs
		int numInts = currentState->numIntervals;
		if(OMX_DEBUG) { Rprintf("Populating CIs for %d objectives.\n", numInts); }
		double* interval = REAL(intervals);
		int* intervalCode = INTEGER(intervalCodes);
		for(int j = 0; j < numInts; j++) {
			omxConfidenceInterval *oCI = &(currentState->intervalList[j]);
			interval[j] = oCI->min;
			interval[j + numInts] = oCI->max;
			intervalCode[j] = oCI->lCode;
			intervalCode[j + numInts] = oCI->uCode;
		}
	}
	
	REAL(evaluations)[1] = currentState->computeCount;

	int nextEl = 0;

	SET_STRING_ELT(names, nextEl++, mkChar("minimum"));
	SET_STRING_ELT(names, nextEl++, mkChar("estimate"));
	SET_STRING_ELT(names, nextEl++, mkChar("gradient"));
	SET_STRING_ELT(names, nextEl++, mkChar("hessianCholesky"));
	SET_STRING_ELT(names, nextEl++, mkChar("status"));
	SET_STRING_ELT(names, nextEl++, mkChar("iterations"));
	SET_STRING_ELT(names, nextEl++, mkChar("evaluations"));
	SET_STRING_ELT(names, nextEl++, mkChar("matrices"));
	SET_STRING_ELT(names, nextEl++, mkChar("algebras"));
	SET_STRING_ELT(names, nextEl++, mkChar("confidenceIntervals"));
	SET_STRING_ELT(names, nextEl++, mkChar("confidenceIntervalCodes"));
	SET_STRING_ELT(names, nextEl++, mkChar("calculatedHessian"));
	SET_STRING_ELT(names, nextEl++, mkChar("standardErrors"));

	nextEl = 0;

	SET_VECTOR_ELT(ans, nextEl++, minimum);
	SET_VECTOR_ELT(ans, nextEl++, estimate);
	SET_VECTOR_ELT(ans, nextEl++, gradient);
	SET_VECTOR_ELT(ans, nextEl++, hessian);
	SET_VECTOR_ELT(ans, nextEl++, status);
	SET_VECTOR_ELT(ans, nextEl++, iterations);
	SET_VECTOR_ELT(ans, nextEl++, evaluations);
	SET_VECTOR_ELT(ans, nextEl++, matrices);
	SET_VECTOR_ELT(ans, nextEl++, algebras);
	SET_VECTOR_ELT(ans, nextEl++, intervals);
	SET_VECTOR_ELT(ans, nextEl++, intervalCodes);
	if(numHessians == 0) {
		SET_VECTOR_ELT(ans, nextEl++, NAmat);
	} else {
		SET_VECTOR_ELT(ans, nextEl++, calculatedHessian);
	}
	if(!calculateStdErrors) {
		SET_VECTOR_ELT(ans, nextEl++, NAmat);
	} else {
		SET_VECTOR_ELT(ans, nextEl++, stdErrors);
	}
	namesgets(ans, names);

	if(OMX_VERBOSE) {
		Rprintf("Inform Value: %d\n", currentState->optimumStatus);
		Rprintf("--------------------------\n");
	}

	/* Free data memory */
	omxFreeState(currentState);

	UNPROTECT(numReturns);						// Unprotect Output Parameters
	UNPROTECT(8);								// Unprotect internals

	if(OMX_DEBUG) {Rprintf("All vectors freed.\n");}

	return(ans);

}

/****** Objective Function *********/
void F77_SUB(objectiveFunction)
	(	int* mode, int* n, double* x,
		double* f, double* g, int* nstate )
{
	unsigned short int checkpointNow = FALSE;

	if(OMX_DEBUG) {Rprintf("Starting Objective Run.\n");}

	if(*mode == 1) {
		currentState->majorIteration++;
		currentState->minorIteration = 0;
		checkpointNow = TRUE;					// Only checkpoint at major iterations.
	} else currentState->minorIteration++;

	omxMatrix* objectiveMatrix = currentState->objectiveMatrix;
	char errMsg[5] = "";
	omxRaiseError(currentState, 0, errMsg);						// Clear Error State
	/* Interruptible? */
	R_CheckUserInterrupt();
    /* This allows for abitrary repopulation of the free parameters.
     * Typically, the default is for repopulateFun to be NULL,
     * and then handleFreeVarList is invoked */
	if (objectiveMatrix->objective->repopulateFun != NULL) {
		objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, x, *n);
	} else {
		handleFreeVarList(currentState, x, *n);
	}
	omxRecompute(objectiveMatrix);

	/* Derivative Calculation Goes Here. */

	if(isnan(objectiveMatrix->data[0])) {
		if(OMX_DEBUG) {
			Rprintf("Objective value is NaN.\n");
		}
		omxRaiseError(currentState, -1, "Objective function returned a value of NaN.");
		*mode = -1;
	}

	if(isinf(objectiveMatrix->data[0])) {
		if(OMX_DEBUG) {
			Rprintf("Objective Value is infinite.\n");
		}
		omxRaiseError(currentState, -1, "Objective function returned an infinite value.");
		*mode = -1;
	}

	if(currentState->statusCode <= -1) {		// At some point, we'll add others
		if(OMX_DEBUG) {
			Rprintf("Error status reported.\n");
		}
		*mode = -1;
	}

	*f = objectiveMatrix->data[0];
	if(OMX_VERBOSE) {
		Rprintf("Objective Value is: %f, Mode is %d.\n", objectiveMatrix->data[0], *mode);
	}

	if(OMX_DEBUG) { Rprintf("-======================================================-\n"); }

	if(checkpointNow && currentState->numCheckpoints != 0) {	// If it's a new major iteration
		omxSaveCheckpoint(currentState, x, f);		// Check about saving a checkpoint
	}

}

/* (Non)Linear Constraint Functions */
void F77_SUB(constraintFunction)
	(	int *mode, int *ncnln, int *n,
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	if(OMX_DEBUG) { Rprintf("Constraint function called.\n");}

	if(*mode==1) {
		if(OMX_DEBUG) {
			Rprintf("But only gradients requested.  Returning.\n");
			Rprintf("-=====================================================-\n");
		}
		return;
	}

	int j, k, l = 0, size;

	handleFreeVarList(currentState, x, *n);

	for(j = 0; j < currentState->numConstraints; j++) {
		omxRecompute(currentState->conList[j].result);
		size = currentState->conList[j].result->rows * currentState->conList[j].result->cols;
		if(OMX_VERBOSE) { omxPrint(currentState->conList[j].result, "Constraint evaluates as:"); }
		for(k = 0; k < currentState->conList[j].size; k++){
			c[l++] = currentState->conList[j].result->data[k];
		}
	}

	if(OMX_DEBUG) { Rprintf("-=======================================================-\n"); }

	return;

}

/****** HELPER FUNCTIONS ******* (integrate into omxOptimizer) */
/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(omxState* os, double* x, int numVars) {

	if(OMX_DEBUG) {
		Rprintf("Processing Free Parameter Estimates.\n");
		Rprintf("Number of free parameters is %d.\n", numVars);
	}

	if(numVars == 0) return;

	omxFreeVar* freeVarList = os->freeVarList;
	omxMatrix** matrixList = os->matrixList;

	os->computeCount++;

	if(OMX_VERBOSE) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d.%d (%d)\n", os->majorIteration, os->minorIteration, os->computeCount);
		Rprintf("Estimates: [");
		for(int k = 0; k < numVars; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
		Rprintf("--------------------------\n");
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < numVars; k++) {
		// if(OMX_DEBUG) { Rprintf("%d: %f - %d\n", k,  x[k], freeVarList[k].numLocations); }
		for(int l = 0; l < freeVarList[k].numLocations; l++) {
			*(freeVarList[k].location[l]) = x[k];
			if(OMX_DEBUG) {
				Rprintf("Setting location:%ld to value %f for var %d\n",
					freeVarList[k].location[l], x[k], k);
			}
			omxMarkDirty(matrixList[freeVarList[k].matrices[l]]);
		}
	}
}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str) {
/* Attribution: modified from the code given in Writing R Extensions */
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

SEXP getVar(SEXP str, SEXP env) {
/* Attribution: modified from the code given in Writing R Extensions */
   SEXP ans;
   if(!isString(str) || length(str) != 1)
        error("getVar: variable name is not a single string");
   if(!isEnvironment(env))
	error("getVar: env should be an environment");
   ans = findVar(install(CHAR(STRING_ELT(str, 0))), env);
   return(ans);
}

/* Objective function for confidence interval limit finding. 
 * Replaces the standard objective function when finding confidence intervals. */
void F77_SUB(limitObjectiveFunction)
	(	int* mode, int* n, double* x, double* f, double* g, int* nstate ) {
		
		if(OMX_VERBOSE) Rprintf("Calculating interval %d, %s boundary:", currentState->currentInterval, (currentState->intervalList[currentState->currentInterval].calcLower?"lower":"upper"));

		F77_CALL(objectiveFunction)(mode, n, x, f, g, nstate);	// Standard objective function call

		omxConfidenceInterval *oCI = &(currentState->intervalList[currentState->currentInterval]);
		
		omxRecompute(oCI->matrix);
		
		double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);

		if(OMX_DEBUG) {
			Rprintf("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.\n",
				oCI->lbound, oCI->ubound, *f, CIElement);
		}

		/* Catch boundary-passing condition */
		if(isnan(CIElement) || isinf(CIElement)) {
			omxRaiseError(currentState, -1, 
				"Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			*mode = -1;
			return;
		}

		if(oCI->calcLower) {
			double diff = oCI->lbound - *f;		// Offset - likelihood
			*f = diff*diff + CIElement;
				// Minimize element for lower bound.
		} else {
			double diff = oCI->ubound - *f;			// Offset - likelihood
			*f = diff*diff - CIElement;
				// Maximize element for upper bound.
		}

		if(OMX_DEBUG) {
			Rprintf("Interval Objective in previous iteration was calculated to be %f.\n", *f);
		}
}


/*************************************************************************************
 *
 *   omxEstimateHessian
 *
 *  author: tbrick, 2010-02-04
 *
 *  Based on code in the numDeriv library for R <citation needed> // TODO: Find NumDeriv Citation
 *
 *  @params matList				list of objectives to evaluate
 *  @params gradient 			an nx1 structure to be filled with the gradient
 *  @params functionPrecision	functional precision for the calculation
 *  @params v					divisor for Richardson approximation
 *  @params r					number of repetitions for Richardson approximation
 *  @params currentState		the current omxState
 *  @params optimum				the current optimum value of the objective function
 *
 ************************************************************************************/
unsigned short omxEstimateHessian(omxMatrix** matList, int numHessians, double functionPrecision, /*double v,*/ int r, omxState* currentState, double optimum) {

	// TODO: Check for nonlinear constraints and adjust algorithm accordingly.
	// TODO: Allow more than one hessian value for calculation

	double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
	double eps = 1E-4;	// Kept here for access purposes.

	if(numHessians > 1) {
		error("NYI: Cannot yet calculate more than a single hessian per optimization.\n");
	}

	if(numHessians == 0) return FALSE;

	omxObjective* oo = matList[0]->objective;
	int numParams = currentState->numFreeParams;
	double* freeParams = (double*) Calloc(numParams, double);
	double* optima = currentState->optimalValues;

	if(oo->hessian == NULL) {
		oo->hessian = (double*) R_alloc(numParams * numParams, sizeof(double));
		if(OMX_DEBUG) {Rprintf("Generated hessian memory, (%d x %d), at 0x%x.\n", numParams, numParams, oo->hessian);}
	}

	if(oo->gradient == NULL) {
		oo->gradient = (double*) R_alloc(numParams, sizeof(double));
		if(OMX_DEBUG) {Rprintf("Generated gradient memory, (%d), at 0x%x.\n", numParams, oo->gradient);}
	}

	double* gradient = oo->gradient;
	double* hessian = oo->hessian;

	double* Haprox = (double*) Calloc(r, double);						// Hessian Workspace
	double* Gaprox = (double*) Calloc(r, double);						// Gradient Workspace

	for(int i = 0; i < numParams; i++) {
		freeParams[i] = optima[i];
	}

	omxMatrix* objectiveMatrix = oo->matrix;

	if (objectiveMatrix->objective->repopulateFun != NULL) {	//  Just in case
		objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, freeParams, numParams);
	} else {
		handleFreeVarList(currentState, freeParams, numParams);
	}
	omxRecompute(objectiveMatrix);								// Initial recompute in case it matters.
	double f0 = omxMatrixElement(objectiveMatrix, 0, 0);

	for(int i = 0; i < numParams; i++) {			// Which parameter?

		/* Part the first: Gradient and diagonal */
		double iOffset = fabs(functionPrecision*optima[i]);
		if(fabs(iOffset) < eps) iOffset += eps;
		if(OMX_DEBUG) {Rprintf("Hessian estimation: iOffset: %f.\n", iOffset);}
		for(int k = 0; k < r; k++) {			// Decreasing step size, starting at k == 0
			if(OMX_DEBUG) {Rprintf("Hessian estimation: Parameter %d at refinement level %d (%f). One Step Forward.\n", i, k, iOffset);}
			freeParams[i] = optima[i] + iOffset;
			if (objectiveMatrix->objective->repopulateFun != NULL) {	//  Just in case
				objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, freeParams, numParams);
			} else {
				handleFreeVarList(currentState, freeParams, numParams);
			}
			omxRecompute(objectiveMatrix);
			double f1 = omxMatrixElement(objectiveMatrix, 0, 0);

			if(OMX_DEBUG) {Rprintf("Hessian estimation: One Step Back.\n");}

			freeParams[i] = optima[i] - iOffset;
			if (objectiveMatrix->objective->repopulateFun != NULL) {	// Just in case
				objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, freeParams, numParams);
			} else {
				handleFreeVarList(currentState, freeParams, numParams);
			}
			omxRecompute(objectiveMatrix);
			double f2 = omxMatrixElement(objectiveMatrix, 0, 0);

			Gaprox[k] = (f1 - f2) / (2.0*iOffset); 						// This is for the gradient
			Haprox[k] = (f1 - 2.0*f0 + f2) / (iOffset * iOffset);		// This is second derivative
			freeParams[i] = optima[i];									// Reset parameter value
			iOffset /= v;
			if(OMX_DEBUG) {Rprintf("Hessian estimation: (%d, %d)--Calculating F1: %f F2: %f, Haprox: %f...\n", i, i, f1, f2, Haprox[k]);}
		}

		for(int m = 1; m < r; m++) {						// Richardson Step
			for(int k = 0; k < (r - m); k++) {
				Gaprox[k] = (Gaprox[k+1] * pow(4.0, m) - Gaprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
				Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k])/(pow(4.0, m)-1); // NumDeriv Hard-wires 4s for r here. Why?
			}
		}

		if(OMX_DEBUG) { Rprintf("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d) with value %f...\n", hessian, i, i, i*numParams+i, Haprox[0]); }
		gradient[i] = Gaprox[0];						// NPSOL reports a gradient that's fine.  Why report two?
		hessian[i*numParams + i] = Haprox[0];
	// }

	/* Part the Second: Off-diagonals */  			// These go here because they depend on their associated diagonals
	// for(int i = 0; i < numParams; i++) {			// Which parameter (Again).
		for(int l = (i-1); l >= 0; l--) {				// Crossed with which parameter? (Start with diagonal)
			double iOffset = fabs(functionPrecision*optima[i]);
			if(fabs(iOffset) < eps) iOffset += eps;
			double lOffset = fabs(functionPrecision*optima[l]);
			if(fabs(lOffset) < eps) lOffset += eps;

			for(int k = 0; k < r; k++) {
				freeParams[i] = optima[i] + iOffset;
				freeParams[l] = optima[l] + lOffset;
				if (objectiveMatrix->objective->repopulateFun != NULL) {	//  Just in case
					objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, freeParams, numParams);
				} else {
					handleFreeVarList(currentState, freeParams, numParams);
				}
				omxRecompute(objectiveMatrix);
				double f1 = omxMatrixElement(objectiveMatrix, 0, 0);

				if(OMX_DEBUG) {Rprintf("Hessian estimation: One Step Back.\n");}

				freeParams[i] = optima[i] - iOffset;
				freeParams[l] = optima[l] - lOffset;
				if (objectiveMatrix->objective->repopulateFun != NULL) {	// Just in case
					objectiveMatrix->objective->repopulateFun(objectiveMatrix->objective, freeParams, numParams);
				} else {
					handleFreeVarList(currentState, freeParams, numParams);
				}
				omxRecompute(objectiveMatrix);
				double f2 = omxMatrixElement(objectiveMatrix, 0, 0);

				Haprox[k] = (f1 - 2.0*f0 + f2 - hessian[i*numParams+i]*iOffset*iOffset -
							hessian[l*numParams+l]*lOffset*lOffset)/(2.0*iOffset*lOffset);
				if(OMX_DEBUG) {Rprintf("Hessian first off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f).\n", Haprox[k], iOffset, lOffset, f1, hessian[i*numParams+i], hessian[l*numParams+l], v, k, pow(v, k), functionPrecision*optima[i], functionPrecision*optima[l]);}

				freeParams[i] = optima[i];				// Reset parameter values
				freeParams[l] = optima[l];

				iOffset = iOffset / v;					//  And shrink step
				lOffset = lOffset / v;
			}

			for(int m = 1; m < r; m++) {						// Richardson Step
				for(int k = 0; k < (r - m); k++) {
					if(OMX_DEBUG) {Rprintf("Hessian off-diagonal calculation: Haprox = %f, iOffset = %f, lOffset=%f from params %f, %f and %f, %f and %d (also: %f, %f and %f, and %f).\n", Haprox[k], iOffset, lOffset, functionPrecision, optima[i], optima[l], v, m, pow(4.0, m), functionPrecision*optima[i], functionPrecision*optima[l], k);}
					Haprox[k] = (Haprox[k+1] * pow(4.0, m) - Haprox[k]) / (pow(4.0, m)-1);
				}
			}

			if(OMX_DEBUG) {Rprintf("Hessian estimation: Populating Hessian (0x%x) at ([%d, %d] = %d and %d) with value %f...", hessian, i, l, i*numParams+l, l*numParams+i, Haprox[0]);}
			hessian[i*numParams+l] = Haprox[0];
			hessian[l*numParams+i] = Haprox[0];
		}

		if(OMX_DEBUG) {Rprintf("Done with parameter %d.\n", i);}
	}

	if(OMX_DEBUG) {Rprintf("Hessian Computation complete.\n");}

	Free(freeParams);
	Free(Haprox);
	Free(Gaprox);

	return TRUE;

}
