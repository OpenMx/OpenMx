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

#include <stdio.h>
#include "omxState.h"
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxObjective.h"
//#include "omxSymbolTable.h"

/* NPSOL-related functions */
extern void F77_SUB(npoptn)(char* string, int length);
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,	double *bl,	double *bu, void* funcon, void* funobj,
						int *inform, int *iter, int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,
						double *x, int *iw,	int *leniw, double *w, int *lenw);

/* Objective Function */
void F77_SUB(objectiveFunction)	( int* mode, int* n, double* x, double* f, double* g, int* nstate );

/* Constraint Function */
void F77_SUB(constraintFunction)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate);	// Constraints in the style of old Mx

/* Helper functions */
void handleFreeVarList(omxState *os, double* x, int numVars);	// Locates and inserts elements from the optimizer into matrices.
SEXP getListElement(SEXP list, const char *str); 				// Gets the value named str from SEXP list.  From "Writing R Extensions"
SEXP getVar(SEXP str, SEXP env);								// Gets the object named str from environment env.  From "Writing R Extensions"

/* NPSOL-specific globals */
const double NPSOL_BIGBND = 1e20;
const double NEG_INF = -2e20;
const double INF = 2e20;


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


/* Functions for Export */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList,
	SEXP data, SEXP options, SEXP state);  // Calls NPSOL.  Duh.

/* Set up R .Call info */
R_CallMethodDef callMethods[] = {
{"callNPSOL", (void*(*)())&callNPSOL, 9},
{NULL, NULL, 0}
};

void
R_init_mylib(DllInfo *info)
{
/* Register routines, allocate resources. */
R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void
R_unload_mylib(DllInfo *info)
{
/* Release resources. */
}




/* Main functions */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList,
	SEXP data, SEXP options, SEXP state) {
	// For now, assume no constraints.

	/* NPSOL Arguments */
	void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);

	int ldA, ldJ, ldR, inform, iter, leniw, lenw; // nclin, ncnln,
	int *iw; // , istate;

//	double f;
	double *A, *bl, *bu, *c, *clambda, *x = NULL, *w; //  *g, *R, *cJac,
	double *est, *grad, *hess;

	/* Helpful variables */

	char optionCharArray[250] = "";			// For setting options

	int j, k, l, m;					// Index Vars

	int nctotl, nlinwid, nlnwid;	// Helpful side variables.

	SEXP nextLoc, nextMat, nextAlg, nextVar;

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
	if(OMX_DEBUG) { Rprintf("Processing %d data source(s).\n", length(data));}
	currentState->numData = length(data);
	currentState->dataList = (omxData**) R_alloc(length(data), sizeof(omxData*));

	for(k = 0; k < length(data); k++) {
		PROTECT(nextLoc = VECTOR_ELT(data, k));			// Retrieve the data object
		currentState->dataList[k] = omxNewDataFromMxData(NULL, nextLoc, currentState);
		if(OMX_DEBUG) {
			Rprintf("Data initialized at 0x%0xd = (%d x %d).\n",
				currentState->dataList[k], currentState->dataList[k]->rows, currentState->dataList[k]->cols);
		}
		UNPROTECT(1); // nextMat
	}

	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}
	currentState->numMats = length(matList);
	currentState->matrixList = (omxMatrix**) R_alloc(length(matList), sizeof(omxMatrix*));

	for(k = 0; k < length(matList); k++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, k));		// This is the matrix + populations
		PROTECT(nextMat = VECTOR_ELT(nextLoc, 0));		// The first element of the list is the matrix of values
		currentState->matrixList[k] = omxNewMatrixFromMxMatrix(nextMat, currentState);
		if(OMX_DEBUG) {
			Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n",
				currentState->matrixList[k], currentState->matrixList[k]->rows, currentState->matrixList[k]->cols);
		}
		UNPROTECT(2); // nextLoc and nextMat
	}

	/* Process Algebras Here */
	currentState->numAlgs = length(algList);
	SEXP algListNames = getAttrib(algList, R_NamesSymbol);

	l = 1;
	if(OMX_DEBUG) { Rprintf("Processing %d algebras.\n", currentState->numAlgs, length(algList)); }
	currentState->algebraList = (omxMatrix**) R_alloc(currentState->numAlgs, sizeof(omxMatrix*));

	for(int j = 0; j < currentState->numAlgs; j++) {
		currentState->algebraList[j] = omxInitMatrix(NULL, 0,0, TRUE, currentState);
	}

	for(int j = 0; j < currentState->numAlgs; j++) {
		PROTECT(nextAlg = VECTOR_ELT(algList, j));		// The next algebra or objective to process
		if(OMX_DEBUG) { Rprintf("Intializing algebra %d at location 0x%0x.\n", j, currentState->algebraList + j); }
		if(IS_S4_OBJECT(nextAlg)) {											// This is an objective object.
			omxFillMatrixFromMxObjective(currentState->algebraList[j], nextAlg);
		} else {															// This is an algebra spec.
			omxFillMatrixFromMxAlgebra(currentState->algebraList[j], 
				nextAlg, CHAR(STRING_ELT(algListNames, j)));
		}
		UNPROTECT(1);	// nextAlg
	}

	/* Process Objective Function */
	if(!isNull(objective)) {
		if(OMX_DEBUG) { Rprintf("Processing objective function.\n"); }
		currentState->objectiveMatrix = omxNewMatrixFromMxIndex(objective, currentState);
	} else {
		currentState->objectiveMatrix = NULL;
		n = 0;
		currentState->numFreeParams = n;
	}

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
	/*
	 varList is a list().  Each element of this list corresponds to one free parameter.
	 Each free parameter is a list.  The first element of this list is the lower bound.
	 The second element of the list is the upper bound.  The remaining elements of this
	 list are 3-tuples.  These 3-tuples are (mxIndex, row, col).
    */
	if(OMX_VERBOSE) { Rprintf("Processing Free Parameters.\n"); }
	currentState->freeVarList = (omxFreeVar*) R_alloc (n, sizeof (omxFreeVar));				// Data for replacement of free vars
	for(k = 0; k < n; k++) {
		PROTECT(nextVar = VECTOR_ELT(varList, k));
		int numLocs = length(nextVar) - 2;
		currentState->freeVarList[k].numLocations = numLocs;
		currentState->freeVarList[k].location = (double**) R_alloc(numLocs, sizeof(double*));
		currentState->freeVarList[k].matrices = (int*) R_alloc(numLocs, sizeof(int));

		/* Lower Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));							// Position 0 is lower bound.
		currentState->freeVarList[k].lbound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[k].lbound)) currentState->freeVarList[k].lbound = NEG_INF;
		if(currentState->freeVarList[k].lbound == 0.0) currentState->freeVarList[k].lbound = 0.0;
		UNPROTECT(1); // NextLoc
		/* Upper Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));							// Position 1 is upper bound.
		currentState->freeVarList[k].ubound = REAL(nextLoc)[0];
		if(ISNA(currentState->freeVarList[k].ubound)) currentState->freeVarList[k].ubound = INF;
		if(currentState->freeVarList[k].ubound == 0.0) currentState->freeVarList[k].ubound = -0.0;
		UNPROTECT(1); // NextLoc

		if(OMX_DEBUG) { Rprintf("Free parameter %d bounded (%f, %f): %d locations\n", k, currentState->freeVarList[k].lbound, currentState->freeVarList[k].ubound, numLocs); }
		for(l = 0; l < currentState->freeVarList[k].numLocations; l++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, l+2));
			int* theVarList = INTEGER(nextLoc);			// These come through as integers.

			int theMat = theVarList[0];			// Matrix is zero-based indexed.
			int theRow = theVarList[1];			// Row is zero-based.
			int theCol = theVarList[2];			// Column is zero-based.

			currentState->freeVarList[k].location[l] = omxLocationOfMatrixElement(currentState->matrixList[theMat], theRow, theCol);
			currentState->freeVarList[k].matrices[l] = theMat;
			UNPROTECT(1); // nextLoc
		}
		UNPROTECT(1); // nextVar
	}

	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }

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
		currentState->conList[k].result = omxNewAlgebraFromOperatorAndArgs(10, arg1, arg2, currentState); // 10 = binary subtract
		omxRecompute(currentState->conList[k].result);
		currentState->conList[k].size = currentState->conList[k].result->rows * currentState->conList[k].result->cols;
		ncnln += currentState->conList[k].size;
	}
	if(OMX_VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d effective constraints.\n", ncnln); }
	funcon = F77_SUB(constraintFunction);

	/* Set up Optimization Memory Allocations */

	if(n == 0) {			// Special Case for the evaluation-only condition

		if(OMX_DEBUG) { Rprintf("No free parameters.  Avoiding Optimizer Entirely.\n"); }
		int mode = 0, nstate = -1;
		f = 0;
		double* x = NULL, *g = NULL;

		if(currentState->objectiveMatrix != NULL) {
			F77_SUB(objectiveFunction)(&mode, &n, x, &f, g, &nstate);
		};

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
			sprintf(optionCharArray, "%s %s", CHAR(STRING_ELT(optionNames, i)), STRING_VALUE(VECTOR_ELT(options, i)));
			F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
			if(OMX_DEBUG) { Rprintf("Option %s \n", optionCharArray); }
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

		F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
						(void*) F77_SUB(objectiveFunction), &inform, &iter, istate, c, cJac,
						clambda, &f, g, R, x, iw, &leniw, w, &lenw);
		if(OMX_DEBUG) { Rprintf("Final Objective Value is: %f.\n", f); }
		handleFreeVarList(currentState, x, n);
	}

	SEXP minimum, estimate, gradient, hessian, code, status, statusMsg, iterations, ans, names, algebras, algebra, matrices;

	int numReturns = 8;

	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(status = allocVector(VECSXP, 3));
	PROTECT(iterations = NEW_NUMERIC(1));
	PROTECT(matrices = NEW_LIST(currentState->numMats));
	PROTECT(algebras = NEW_LIST(currentState->numAlgs));

	/* N-dependent SEXPs */
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));

	/* Store outputs for return */
	if(objective != NULL) {
		REAL(minimum)[0] = f;
	} else {
		REAL(minimum)[0] = R_NaReal;
	}

	est = REAL(estimate);  // Aliases to avoid repeated function calls.
	grad = REAL(gradient);
	hess = REAL(hessian);

	for(k = 0; k < n; k++) {		// Do these with memcpy instead, probably.
		est[k] = x[k];
		grad[k] = g[k];
		for(l = 0; l < n; l++) {
			hess[k*n + l] = R[k*n + l];  // This is ok, because they're both in Col-Major already.
		}
	}

	for(k = 0; k < currentState->numMats; k++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Matrix %d.\n", k); }
		omxRecompute(currentState->matrixList[k]);
		PROTECT(nextMat = allocMatrix(REALSXP, currentState->matrixList[k]->rows, currentState->matrixList[k]->cols));
		for(l = 0; l < currentState->matrixList[k]->rows; l++)
			for(j = 0; j < currentState->matrixList[k]->cols; j++)
				REAL(nextMat)[j * currentState->matrixList[k]->rows + l] =
					omxMatrixElement(currentState->matrixList[k], l, j);
		SET_VECTOR_ELT(matrices, k, nextMat);

		UNPROTECT(1);	/* nextMat */
	}

	for(k = 0; k < currentState->numAlgs; k++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Algebra %d.\n", k); }
		omxRecompute(currentState->algebraList[k]);
		PROTECT(algebra = allocMatrix(REALSXP, currentState->algebraList[k]->rows, currentState->algebraList[k]->cols));
		for(l = 0; l < currentState->algebraList[k]->rows; l++)
			for(j = 0; j < currentState->algebraList[k]->cols; j++)
				REAL(algebra)[j * currentState->algebraList[k]->rows + l] =
					omxMatrixElement(currentState->algebraList[k], l, j);
		SET_VECTOR_ELT(algebras, k, algebra);

		UNPROTECT(1);	/* algebra */
	}
	if(OMX_DEBUG) { Rprintf("All Algebras complete.\n"); }

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
		}
		if(OMX_DEBUG) { Rprintf("Done.\n");}
	} else {
		PROTECT(ans = allocVector(VECSXP, numReturns));
		PROTECT(names = allocVector(STRSXP, numReturns));
	}


	REAL(code)[0] = inform;
	REAL(iterations)[0] = iter;

	/* Fill Status code. */
	SET_VECTOR_ELT(status, 0, code);
	PROTECT(code = NEW_NUMERIC(1));
	REAL(code)[0] = currentState->statusCode;
	SET_VECTOR_ELT(status, 1, code);
	PROTECT(statusMsg = allocVector(STRSXP, 1));
	SET_STRING_ELT(statusMsg, 0, mkChar(currentState->statusMsg));
	SET_VECTOR_ELT(status, 2, statusMsg);
	SET_STRING_ELT(names, 0, mkChar("minimum"));
	SET_STRING_ELT(names, 1, mkChar("estimate"));
	SET_STRING_ELT(names, 2, mkChar("gradient"));
	SET_STRING_ELT(names, 3, mkChar("hessianCholesky"));
	SET_STRING_ELT(names, 4, mkChar("status"));
	SET_STRING_ELT(names, 5, mkChar("iterations"));
	SET_STRING_ELT(names, 6, mkChar("matrices"));
	SET_STRING_ELT(names, 7, mkChar("algebras"));

	SET_VECTOR_ELT(ans, 0, minimum);
	SET_VECTOR_ELT(ans, 1, estimate);
	SET_VECTOR_ELT(ans, 2, gradient);
	SET_VECTOR_ELT(ans, 3, hessian);
	SET_VECTOR_ELT(ans, 4, status);
	SET_VECTOR_ELT(ans, 5, iterations);
	SET_VECTOR_ELT(ans, 6, matrices);
	SET_VECTOR_ELT(ans, 7, algebras);
	namesgets(ans, names);

	if(OMX_VERBOSE) {
		Rprintf("Inform Value: %d\n", inform);
		Rprintf("--------------------------\n");
	}

	/* Free data memory */
	omxFreeState(currentState);

	UNPROTECT(13);						// Unprotect Output Parameters

	return(ans);

}


/****** Objective Function *********/
void F77_SUB(objectiveFunction)
	(	int* mode, int* n, double* x,
		double* f, double* g, int* nstate )
{
	if(OMX_DEBUG) {Rprintf("Starting Objective Run.\n");}

	omxMatrix* objectiveMatrix = currentState->objectiveMatrix;
	char errMsg[5] = "";
	omxRaiseError(currentState, 0, errMsg);						//Clear Error State
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

	if(isnan(objectiveMatrix->data[0]) || isinf(objectiveMatrix->data[0])) {
		if(OMX_DEBUG) {
			Rprintf("Objective Value is incorrect.\n", objectiveMatrix->data[0]);
		}
		*mode = -1;
	}

	if(currentState->statusCode == -1) {		// At some point, we'll add others
		*mode = -1;
	}

	*f = objectiveMatrix->data[0];
	if(OMX_VERBOSE) {
		Rprintf("Objective Value is: %f.\n", objectiveMatrix->data[0]);
	}

	if(OMX_DEBUG) { Rprintf("-======================================================-\n"); }

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

	omxFreeVar* freeVarList = os->freeVarList;
	omxMatrix** matrixList = os->matrixList;

	if(OMX_DEBUG) {Rprintf("Processing Free Parameter Estimates.\n");}
	os->computeCount++;

	if(OMX_VERBOSE) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d\n", os->computeCount);
		Rprintf("Estimates: [");
		for(int k = 0; k < os->numFreeParams; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
		Rprintf("--------------------------\n");
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < os->numFreeParams; k++) {
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

