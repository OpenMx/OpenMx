#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h> 

#include <stdio.h>
#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxObjective.h"
//#include "omxSymbolTable.h"

#define M(y,z) m[z][y]
#define EPSILON 1e-16
#define TRUE 1
#define FALSE 0

#ifdef DEBUGMX
#define OMX_DEBUG 1
#define VERBOSE 1
#else
#define OMX_DEBUG 0
#define VERBOSE 0
#endif /* DEBUGMX */

/* Structure definitions for object evaluation */
typedef struct omxFreeVar {			// Free Variables
	double lbound, ubound;	// Bounds
	int numLocations;
	double** location;		// And where they go.
	int* matrices;			// Matrix numbers for dirtying.
} omxFreeVar;

typedef struct omxConstraint {			// Free Variables
	int size;
	int opCode;
	omxMatrix* result;
} omxConstraint;


/* NPSOL-related functions */
extern void F77_SUB(npoptn)(char* string, int length);
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,	double *bl,	double *bu, void* funcon, void* funobj, 
						int *inform, int *iter, int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,				
						double *x, int *iw,	int *leniw, double *w, int *lenw);

/* Objective Function */
void F77_SUB(objectiveFunction)	( int* mode, int* n, double* x, double* f, double* g, int* nstate );

/* Constraint Function Variants */
void F77_SUB(noConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); 		// No constraints
void F77_SUB(callRConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); 	// Call R for constraints
void F77_SUB(AlgConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); 		// Algebra constraints
void F77_SUB(oldMxConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); 	
// Constraints in the style of old Mx

/* Helper functions */
void handleFreeVarList(double* x, int numVars);					// Locates and inserts elements from the optimizer.  Should handle Algebras, as well.
SEXP getListElement(SEXP list, const char *str); 				// Gets the value named str from SEXP list.  From "Writing R Extensions"
SEXP getVar(SEXP str, SEXP env);								// Gets the object named str from environment env.  From "Writing R Extensions"

/* NPSOL-specific globals */
const double NPSOL_BIGBND = 1e20;
const double NEG_INF = -2e20;
const double INF = 2e20;


/* Globals for function evaluation */
SEXP RObjFun, RConFun;			// Pointers to the functions NPSOL calls
SEXP env;						// Environment for evaluation and object hunting
int ncalls, nminor;				// For debugging--how many calls?
omxMatrix** matrixList;			// Data matrices and their data.
omxMatrix** algebraList;		// List of the matrices of all algebras.
omxConstraint* conList;			// List of constraints
omxFreeVar* freeVarList;		// List of Free Variables and where they go.
omxMatrix *objMatrix;			// Objective Function Matrix

/* Made global for objective functions that want them */
int n, nclin, ncnln;			// Number of Free Params, linear, and nonlinear constraints
int numCons;					// Number of constraint algebras
double f;						// Objective Function
double *g;						// Gradient Pointer
double *R, *cJac;				// Hessian (Approx) and Jacobian
int *istate;					// Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))


/* Functions for Export */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints, SEXP matList, SEXP varList, SEXP algList, SEXP data, SEXP state);  // Calls NPSOL.  Duh.

/* Set up R .Call info */
R_CallMethodDef callMethods[] = { 
{"callNPSOL", (void*(*)())&callNPSOL, 8}, 
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
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints, SEXP matList, SEXP varList, SEXP algList, SEXP data, SEXP state) {
	// For now, assume no constraints.
	
	ncalls = 0;
	
	int N; // n, 
	SEXP nameString;

	/* NPSOL Arguments */
	void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);
	
	int ldA, ldJ, ldR, inform, iter, leniw, lenw; // nclin, ncnln, 
	int *iw; // , istate;
	
//	double f;
	double *A, *bl, *bu, *c, *clambda, *x, *w; //  *g, *R, *cJac,
	double *est, *grad, *hess;
	
	/* Helpful variables */
	
	char option[250] = "";			// For setting options
	
	int j, k, l, m;					// Index Vars
	
	int nctotl, nlinwid, nlnwid;	// Helpful side variables
	
	int numAlgs, numMats, numObjMats;// Number of algebras/matrices.  Used for memory management.
	
	SEXP nextLoc, nextVar;
	
	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these. */
//	if(!isVector(startVals)) error ("startVals must be a vector");
//	if(!isVector(matList)) error ("matList must be a list");
//	if(!isVector(algList)) error ("algList must be a list");

	n = length(startVals);	
	
	/* Store Data from MxMatrices */
	
	/* Retrieve All Matrices From the MatList */

	if(OMX_DEBUG) { Rprintf("Processing %d matrix(ces).\n", length(matList));}
	numMats = length(matList);
	matrixList = (omxMatrix**) R_alloc(sizeof(omxMatrix*), length(matList));
	
	for(k = 0; k < length(matList); k++) {
		PROTECT(nextLoc = VECTOR_ELT(matList, k));				// TODO: Find out if this duplicates the matrix.
		matrixList[k] = omxNewMatrixFromMxMatrix(nextLoc);
		if(OMX_DEBUG) { Rprintf("Matrix initialized at 0x%0xd = (%d x %d).\n", matrixList[k], matrixList[k]->rows, matrixList[k]->cols); }
		UNPROTECT(1); // nextLoc
	}

	/* Process Algebras Here */
	numAlgs = length(algList);
	l = 1;
	if(OMX_DEBUG) { Rprintf("Processing %d algebras.\n", numAlgs, length(algList)); }
	algebraList = (omxMatrix**) R_alloc(sizeof(omxMatrix*), numAlgs);

	for(int j = 0; j < numAlgs; j++) {
		algebraList[j] = omxInitMatrix(NULL, 0,0, TRUE);
	}
	
	for(int j = 0; j < numAlgs; j++) {
		PROTECT(nextLoc = VECTOR_ELT(algList, j));
		if(OMX_DEBUG) { Rprintf("Intializing algebra %d at location 0x%0x.\n", j, algebraList + j); }
		if(IS_S4_OBJECT(nextLoc)) {											// This is an objective object.
			omxFillMatrixFromMxObjective(algebraList[j], nextLoc, data);
		} else {															// This is an algebra spec.
			omxFillMatrixFromMxAlgebra(algebraList[j], nextLoc);
		}
		UNPROTECT(1);	// nextLoc
	}
	
	/* Process Objective Function */
	if(!isNull(objective)) {
		objMatrix = omxNewMatrixFromMxMatrixPtr(objective);
	} else {
		objMatrix = NULL;
		n = 0;
	}
	
/*	for(int j = 0; j < length(algList); j++) {
		if(OMX_DEBUG) { Rprintf("Computing Algebra %d.\n", j); }
		algebraList[j].compute();
	} */

	/* Process Free Var List */
	if(VERBOSE) { Rprintf("Processing Free Parameters.\n"); }
	omxMatrix dm;
	freeVarList = (omxFreeVar*) R_alloc (sizeof ( omxFreeVar ), n);				// Data for replacement of free vars
	for(k = 0; k < n; k++) {
		PROTECT(nextVar = VECTOR_ELT(varList, k));
		int numLocs = length(nextVar) - 2;
		freeVarList[k].numLocations = numLocs;
		freeVarList[k].location = (double**) R_alloc(sizeof(double*), numLocs);
		freeVarList[k].matrices = (int*) R_alloc(sizeof(int), numLocs);
		
		/* Lower Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));							// Position 0 is lower bound.
		freeVarList[k].lbound = REAL(nextLoc)[0];
		if(ISNA(freeVarList[k].lbound)) freeVarList[k].lbound = NEG_INF;
		if(freeVarList[k].lbound == 0.0) freeVarList[k].lbound = EPSILON;
		UNPROTECT(1); // NextLoc
		/* Upper Bound */
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));							// Position 1 is upper bound.
		freeVarList[k].ubound = REAL(nextLoc)[0];
		if(ISNA(freeVarList[k].ubound)) freeVarList[k].ubound = INF;
		if(freeVarList[k].ubound == 0.0) freeVarList[k].ubound = -EPSILON;
		UNPROTECT(1); // NextLoc
		
		if(OMX_DEBUG) { Rprintf("Free parameter %d bounded (%f, %f): %d locations\n", k, freeVarList[k].lbound, freeVarList[k].ubound, numLocs); }
		for(l = 0; l < freeVarList[k].numLocations; l++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, l+2));
			double* theVarList = REAL(nextLoc);			// These come through as doubles.
			int theMat = (int)theVarList[0];			// Matrix is zero-based indexed.
			int theRow = (int)theVarList[1];			// Row is zero-based.
			int theCol = (int)theVarList[2];			// Column is zero-based.
			freeVarList[k].location[l] = omxLocationOfMatrixElement(matrixList[theMat], theRow, theCol);
			freeVarList[k].matrices[l] = theMat;
			UNPROTECT(1); // nextLoc
		}
		UNPROTECT(1); // nextVar
	}

	if(VERBOSE) { Rprintf("Processed.\n"); }
	
	/* Processing Constraints */
	if(VERBOSE) { Rprintf("Processing Constraints.\n");}
	omxMatrix *arg1, *arg2;
	numCons = length(constraints);
	if(OMX_DEBUG) {Rprintf("Found %d constraints.\n", numCons); }
	conList = (omxConstraint*) R_alloc(sizeof(omxConstraint), numCons);
	ncnln = 0;
	for(k = 0; k < numCons; k++) {
		PROTECT(nextVar = VECTOR_ELT(constraints, k));
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 0));
		arg1 = omxNewMatrixFromMxMatrixPtr(nextLoc);
		PROTECT(nextLoc = VECTOR_ELT(nextVar, 1));
		arg2 = omxNewMatrixFromMxMatrixPtr(nextLoc);
		PROTECT(nextLoc = AS_INTEGER(VECTOR_ELT(nextVar, 2)));
		conList[k].opCode = INTEGER(nextLoc)[0];
		UNPROTECT(4);
		conList[k].result = omxNewAlgebraFromOperatorAndArgs(10, arg1, arg2); // 10 = binary subtract
		omxRecomputeMatrix(conList[k].result);
		conList[k].size = conList[k].result->rows * conList[k].result->cols;
		ncnln += conList[k].size;
	}
	if(VERBOSE) { Rprintf("Processed.\n"); }
	if(OMX_DEBUG) { Rprintf("%d effective constraints.\n", ncnln); }
	funcon = F77_SUB(oldMxConFun);
	
	/* Set up Optimization Memory Allocations */
	
	SEXP minimum, estimate, gradient, hessian, code, iterations, ans, names, algebras, algebra;
	
	PROTECT(ans = allocVector(VECSXP, 7));
	PROTECT(names = allocVector(STRSXP, 7));
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(iterations = NEW_NUMERIC(1));
	PROTECT(algebras = NEW_LIST(numAlgs));
	
	if(n == 0) {			// Special Case for the evaluation-only condition
		
		if(OMX_DEBUG) { Rprintf("No free parameters.  Avoiding Optimizer Entirely.\n"); }
		int mode = 0, nstate = -1;
		f = 0;
		double* x = NULL, *g = NULL;
		
		if(objMatrix != NULL) {
			F77_SUB(objectiveFunction)(&mode, &n, x, &f, g, &nstate);
		};
		
		inform = 0;
		iter = 0;
	
		// Allocate vectors & matrices of length 0,
		// because the front-end will read these values
		// and expects 0 elements, not 1 element that is NA.
		PROTECT(estimate = allocVector(REALSXP, 0));
		PROTECT(gradient = allocVector(REALSXP, 0));
		PROTECT(hessian = allocMatrix(REALSXP, 0, 0));
		
	} else {
		
		/* N-dependent SEXPs */
		PROTECT(estimate = allocVector(REALSXP, n));
		PROTECT(gradient = allocVector(REALSXP, n));
		PROTECT(hessian = allocMatrix(REALSXP, n, n));
	
		/* Initialize Scalar Variables. */	
		nclin = 0;						// No linear constraints.
		
		/* Set boundaries and widths. */
		if(nclin <= 0) {
			nlinwid = 1;				// For memory allocation purposes, nlinwid > 0
		} else {
			nlinwid = nclin;
		}
	
		if(ncnln <= 0) {
			nlnwid = 1;					// For memory allocation purposes nlnwid > 0
		} else {
			nlnwid = ncnln;
		}
		
		nctotl = n + nlinwid + nlnwid;
		
		leniw = 3 * n + nclin + 2 * ncnln;
		lenw = 2 * n * n + n * nclin + 2 * n * ncnln + 20 * n + 11 * nclin + 21 * ncnln;
		
		ldA = nlinwid;  		// nclin;
		ldJ = nlnwid; 			// ncnln;
		ldR = n;				// Need to check on size of this one.

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
			bl[k] = freeVarList[k].lbound; 				// No negative covariances.  -Infinity'd be -10^20.
			bu[k] = freeVarList[k].ubound;				// Infinity would be at 10^20.
		}
		
		for(; k < n+nclin; k++) {
			bl[k] = NEG_INF; 							// Linear constraints have no bounds.
			bu[k] = INF;								// Because there are no linear constraints.
		}
		
		for(l = 0; l < numCons; l++) {						// Nonlinear constraints:
			switch(conList[l].opCode) {
				case 0:									// Less than: Must be strictly less than 0.
					for(m = 0; m < conList[l].size; m++) {
						bl[k] = NEG_INF;
						bu[k] = EPSILON;
					}
					break;
				case 1:									// Equal: Must be roughly equal to 0.
					for(m = 0; m < conList[l].size; m++) {
						bl[k] = -EPSILON;
						bu[k] = EPSILON;
					}
					break;
				case 2:									// Greater than: Must be strictly greater than 0.
					for(m = 0; m < conList[l].size; m++) {
						bl[k] = -EPSILON;
						bu[k] = INF;
					}
					break;
				default:
					for(m = 0; m < conList[l].size; m++) {
						bl[k] = NEG_INF;
						bu[k] = INF;
					}
					break;
			}
			k++;
		}
	
	
		/* Initialize Starting Values */
		if(VERBOSE) {
			Rprintf("--------------------------\n");	
			Rprintf("Starting Values (%d) are:\n", n);
		}
		for(k = 0; k < n; k++) {
			x[k] = REAL(startVals)[k];
			if(x[k] == 0.0) {
				x[k] += .1;
			}
			if(VERBOSE) { Rprintf("%d: %f\n", k, x[k]); }
		}
		if(OMX_DEBUG) {
			Rprintf("--------------------------\n");	
			Rprintf("Setting up optimizer...");
		}
		
	/* 	Set NPSOL options  (Maybe separate this into a different function?)
	
		F77_CALL(npoptn)(char* string, int *length)
	
		String is of the form:
			'Option = Value'
	
	*/
	
		/* Output Options */
		strcpy(option, "Nolist"); 						// Suppress that annoying output file
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Print level 0");  				// 0 = No Output, 20=Verbose
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Minor print level 0");			// 0 = No Output, 20=Verbose
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Print file 0");
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Summary file 0");
		F77_CALL(npoptn)(option, strlen(option));

		/* Options That Change The Optimizer */
		sprintf(option, "Function Precision 1e-14");		// Set epsilon
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Verify level 3");
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Line search tolerance .2");		// Set accuracy to Merit Function
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Step Limit 0");					// At 0, defaults to 2.0
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Derivative Level 0");			// Always estimate gradient and hessian
		F77_CALL(npoptn)(option, strlen(option));	
		strcpy(option, "Hessian Yes");					// Evaluate Hessian
		F77_CALL(npoptn)(option, strlen(option));
														// Iteration limit is max(50, 3(numFreeParams + numLinearConstraints) + 10numNonlinear)
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
						(void*) F77_SUB(objectiveFunction), &inform, &iter, istate, c, cJac, clambda, &f, g, R,
						x, iw, &leniw, w, &lenw);
		if(OMX_DEBUG) {Rprintf("Objective Value is: %f.\n", f); }
	}
	
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
	
	for(k = 0; k < numAlgs; k++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Algebra %d.\n", k); }
		omxRecomputeMatrix(algebraList[k]);
		PROTECT(algebra = allocMatrix(REALSXP, algebraList[k]->rows, algebraList[k]->cols));
		for(j = 0; j < algebraList[k]->cols; j++) 
			for(l = 0; l <algebraList[k]->rows; l++)
				REAL(algebra)[j * algebraList[k]->rows + l] =
					omxMatrixElement(algebraList[k], j, l);
		SET_VECTOR_ELT(algebras, k, algebra);
		UNPROTECT(1);	/* algebra */
	}
	
	REAL(code)[0] = inform;
	REAL(iterations)[0] = iter;
	
	SET_STRING_ELT(names, 0, mkChar("minimum"));
	SET_STRING_ELT(names, 1, mkChar("estimate"));
	SET_STRING_ELT(names, 2, mkChar("gradient"));
	SET_STRING_ELT(names, 3, mkChar("hessian"));
	SET_STRING_ELT(names, 4, mkChar("code"));
	SET_STRING_ELT(names, 5, mkChar("iterations"));
	SET_STRING_ELT(names, 6, mkChar("algebras"));
	
	SET_VECTOR_ELT(ans, 0, minimum);
	SET_VECTOR_ELT(ans, 1, estimate);
	SET_VECTOR_ELT(ans, 2, gradient);
	SET_VECTOR_ELT(ans, 3, hessian);
	SET_VECTOR_ELT(ans, 4, code);
	SET_VECTOR_ELT(ans, 5, iterations);
	SET_VECTOR_ELT(ans, 6, algebras);
	namesgets(ans, names);
		
	if(VERBOSE) { 
		Rprintf("Inform Value: %d\n", inform); 
		Rprintf("--------------------------\n");
	}

	/* Free data memory */
	if(OMX_DEBUG) { Rprintf("Freeing Algebras.\n");}
	for(k = 0; k < numAlgs; k++) {
		if(OMX_DEBUG) { Rprintf("Freeing Algebra %d.\n", k); }
		omxFreeAllMatrixData(algebraList[k]);
	}
	
	if(OMX_DEBUG) { Rprintf("Freeing Matrices.\n");}
	for(k = 0; k < numMats; k++) {
		omxFreeAllMatrixData(matrixList[k]);
	}
	
	UNPROTECT(9);						// Unprotect NPSOL Parameters
	
	return(ans);

}


/****** Objective Function *********/
void F77_SUB(objectiveFunction)
	(	int* mode, int* n, double* x, 
		double* f, double* g, int* nstate )
{

	/* Interruptible? */
	R_CheckUserInterrupt();

	handleFreeVarList(x, *n);
	
	omxRecomputeMatrix(objMatrix);
	
	/* Derivative Calculation Goes Here. */
	
	*f = objMatrix->data[0];
	if(VERBOSE) {
		Rprintf("Objective Value is: %f.\n", objMatrix->data[0]);
	}

}

/* (Non)Linear Constraint Functions */
void F77_SUB(oldMxConFun)
	(	int *mode, int *ncnln, int *n, 
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	if(OMX_DEBUG) {Rprintf("Constraint function called.\n");}

	ncalls--;
	nminor++;
	int j, k, l = 0, size;
	
	handleFreeVarList(x, *n);

	for(j = 0; j < numCons; j++) {
		omxRecomputeMatrix(conList[j].result);
		size = conList[j].result->rows * conList[j].result->cols;
		if(VERBOSE) { omxPrintMatrix(conList[j].result, "Constraint evaluates as:"); }
		for(k = 0; k < conList[j].size; k++){
			c[l++] = conList[j].result->data[k];
		}
	}

return;

}


void F77_SUB(noConFun)
	(	int *mode, int *ncnln, int *n, 
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	Rprintf("-=====================-\n");
	Rprintf("Funcon called.\n");
	Rprintf("Constraint functions not yet implemented.\n");
	Rprintf("-=====================-\n");	
/* Defines the nonlinear constraints for the run of npsol. */

return;

}

void F77_SUB(callRConFun)
/* Calls R to determine constraints.  Not yet functional. */
	(	int *mode, int *ncnln, int *n, 
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	Rprintf("-=====================================================-\n");
	Rprintf("Funcon called.\n");
	Rprintf("R should be called next.  But won't be.\n");
	Rprintf("Constraint function implementation will follow shortly.\n");
	Rprintf("-=====================================================-\n");	
/* Defines the nonlinear constraints for the run of npsol. */

return;

}

void F77_SUB(AlgConFun)
/* Constraints based on algebra statements */
	(	int *mode, int *ncnln, int *n, 
		int *ldJ, int *needc, double *x,
		double *c, double *cJac, int *nstate)
{

	Rprintf("-=========================================-\n");
	Rprintf("Algebraic constraint function called.\n");
	Rprintf("Constraint functions not yet implemented.\n");
	Rprintf("-=========================================-\n");	
/* Defines the nonlinear constraints for the run of npsol. */

return;

}

/****** HELPER FUNCTIONS ******* (put in separate file) */
/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(double* x, int numVars) {
	
	if(OMX_DEBUG) {Rprintf("Processing Free Parameter Estimates.\n");}
	ncalls++;
	if(VERBOSE) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d (Minor Calls: %d)\n", ncalls, nminor);
		Rprintf("Estimates: [");
		for(int k = 0; k < numVars; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
		Rprintf("--------------------------\n");
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < numVars; k++) {
		if(OMX_DEBUG) { Rprintf("%d: %f - %d\n", k,  x[k], freeVarList[k].numLocations); }
		for(int l = 0; l < freeVarList[k].numLocations; l++) {
			*(freeVarList[k].location[l]) = x[k];
			if(OMX_DEBUG) { Rprintf("Setting location:%ld to value %f for var %d\n", freeVarList[k].location[l], x[k], k); }
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

