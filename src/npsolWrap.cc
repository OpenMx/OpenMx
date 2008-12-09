#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h> 

#include <stdio.h>
#include "omxDataMatrix.h"
//#include "omxSymbolTable.h"

#define M(y,z) m[z][y]
#define EPSILON 0.00000000001
#define TRUE 1
#define FALSE 0
#define DEBUG 1
#define VERBOSE 1

/* Structure definitions for object evaluation */
typedef struct {			// Free Variables
	int numLocations;		
	double** location;		// And where they go.
} freeVar;

/* NPSOL-related functions */
extern "C" {
extern void F77_SUB(npoptn)(char* string, int length);
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,	double *bl,	double *bu, void* funcon, void* funobj, 
						int *inform, int *iter, int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,				
						double *x, int *iw,	int *leniw, double *w, int *lenw);
}
/* Objective Function Variants */
void F77_SUB(callRObjFun)(int* mode, int* n, double* x, double* f, double* g, int* nstate); // Call an R function
void F77_SUB(covObjFun)(int* mode, int* n, double* x, double* f, double* g, int* nstate);	// Covariance Fit
void F77_SUB(FIMLObjFun)(int* mode, int* n, double* x, double* f, double* g, int* nstate);	// FIML Fit

/* Constraint Function Variants */
void F77_SUB(noConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); // No constraints
void F77_SUB(callRConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); // Call R for constraints

/* Helper functions */
void handleFreeVarList(double* x, int n);						// Locates and inserts elements from the optimizer.  Should handle Algebras, as well.
SEXP getListElement(SEXP list, const char *str); 				// Gets the value named str from SEXP list.  From "Writing R Extensions"
SEXP getVar(SEXP str, SEXP env);								// Gets the object named str from environment env.  From "Writing R Extensions"

/* Globals for function evaluation */
SEXP RObjFun, RConFun;			// Pointers to the functions NPSOL calls
SEXP env;						// Environment for evaluation and object hunting
int ncalls;						// For debugging--how many calls?
omxDataMatrix* matrixList;		// Data matrices and their data.
freeVar* freeVarList;			// List of Free Variables and where they go.
omxDataMatrix dataRows;			// All the data, for now.
omxDataMatrix cov;				// The Covariance Matrix, probably a link to the matrixList.
omxDataMatrix means;			// Vector of means, probably a link to the matrixList.

/* Globals for Covariance Evaluation */
omxDataMatrix I; //, Z, C, Y;

extern "C" {
	/* Functions for Export */
	SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP bounds, SEXP matList, SEXP varList, SEXP data, SEXP state);  // Calls NPSOL.  Duh.
	
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
}	/* Extern "C" */



/* Main functions */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP bounds, SEXP matList, SEXP varList, SEXP data, SEXP state) {
	// For now, assume no constraints.
	
	ncalls = 0;
	
	int n, N;
	SEXP nameString;

	/* NPSOL Arguments */
	void (*funobj)(int*, int*, double*, double*, double*, int*);
	void (*funcon)(int*, int*, int*, int*, int*, double*, double*, double*, int*);
	
	int nclin, ncnln, ldA, ldJ, ldR, inform, iter, leniw, lenw;
	int *istate, *iw;
	
	double f;
	double *A, *bl, *bu, *c, *cJac, *clambda, *g, *R, *x, *w;
	double *est, *grad, *hess;
	
	/* Helpful variables */
	
	char option[250] = "";			// For setting options
	
	int k, l, m;					// Index Vars
	
	int nctotl, nlinwid, nlnwid;	// Helpful side variables
	
	SEXP nextLoc, nextVar;
	
	
	/* Sanity Check and Parse Inputs */
	/* TODO: Need to find a way to account for nullness in these. */
//	if(!isVector(startVals)) error ("startVals must be a vector");
//	if(!isVector(matList)) error ("matList must be a list");

	n = length(startVals);	

	/* Parse Mx Objective Structure. */
	PROTECT(nameString = getAttrib(objective, install("class")));
	/* Type Checking removed; we either need a common objective class or a comprehensive list. */
//	if(strcmp(CHAR(STRING_ELT(nameString, 0)), "MxObjectiveFunction") != 0) error("objective must be an MxObjectiveFunction."); 

	/* Determine Type of Optimization, and Set Up Structures */
	/* This should be done with a hash-table lookup for extensibility. */
	if(strncmp(CHAR(STRING_ELT(nameString, 0)), "MxCovObjective", 21) == 0) { // Covariance-style optimization.
		/* In Covariance Optimization, matList contains A, S, F, and the observed Covariance Matrix. */
		if(DEBUG) { Rprintf("Using covariance optimization.\n"); }
		funobj = F77_SUB(covObjFun);
		int *dimList;

		/* Store Data from MxMatrices */
		/* Retrieve All Matrices From the MatList */
		matrixList = (omxDataMatrix*) R_alloc(sizeof(omxDataMatrix), length(matList));		// Stores links to data/covariance matrices
		if(DEBUG) {Rprintf("Processing %d matrices.", length(matList));}
		for(k = 0; k < length(matList) - 1; k++) {										// Last is Covariance Matrix, not an MxMatrix
			PROTECT(nextLoc = VECTOR_ELT(matList, k));									// TODO: Find out if even this one duplicates the matrix.
			matrixList[k].fillFromMatrix(nextLoc);
			UNPROTECT(1); // nextLoc
		}

		/* Last Matrix is the Covariance Matrix */
		PROTECT(nextLoc = VECTOR_ELT(matList, k));
		matrixList[k].fillFromMatrix(nextLoc);
		UNPROTECT(1); // nextLoc

		/* Identity Matrix, Size Of A */
		I.resize(matrixList[0].rows, matrixList[0].cols, FALSE);
		for(k =0; k < I.rows; k++) {
			for(l = 0; l < I.cols; l++) {
				if(l == k) {
					I.setElement(k, l, 1);
				} else {
					I.setElement(k, l, 0);
				}
			}
		}
	} else if(strncmp(CHAR(STRING_ELT(nameString, 0)), "MxFIMLObjective", 15) == 0) {	// FIML optimization.
		if(DEBUG) { Rprintf("Using covariance optimization.\n"); }
		
		funobj = F77_SUB(FIMLObjFun);								// FIML Objective Function
		funcon = F77_SUB(noConFun);									// Check for constraint functions once they're implemented.
		
		matrixList = (omxDataMatrix*) R_alloc(sizeof(omxDataMatrix), length(matList));		// Stores links to data matrices
		if(DEBUG) {Rprintf("Processing %d matrices.\n", length(matList));}
		for(k = 0; k < length(matList) - 1; k++) {
			PROTECT(nextLoc = VECTOR_ELT(matList, k));
			matrixList[k].fillFromMatrix(nextLoc);
			UNPROTECT(1);
		}
		
		/* Process Algebras Here */
		// TODO: Process Algebras.
		if(DEBUG) { Rprintf("Processing %d algebras.\n", 0);}

		/* In FIML optimization, the FIMLObjective Structure contains means and covariances. */
		if(DEBUG) {Rprintf("Processing means.\n");}
		SEXP meanStruct, covStruct;
		PROTECT(meanStruct = AS_NUMERIC(GET_SLOT(objective, install("means"))));		// TODO: Needs sanity check
		means.fillFromMatrix(meanStruct);
		UNPROTECT(1);
		
		if(DEBUG) {Rprintf("Processing covariances.\n");}
		PROTECT(covStruct = GET_SLOT(objective, install("covariance")));	// TODO: Needs sanity check 
		cov.fillFromMatrix(covStruct);										// Until MxMatrices are processable.
		UNPROTECT(1);
		
		/* Process Data Into Data Matrix */
		dataRows.fillFromMatrix(data);
		
	} else if(strncmp(CHAR(STRING_ELT(nameString, 0)), "MxRObjective", 12) == 0) {		// For all others, assume 'R' function
	/* In R-style optimization, an R function is passed in that evaluates the objective function at each step. */
		SEXP objfun, confun;
		PROTECT(objfun = GET_SLOT(objective, install("objective")));
		if(!isFunction(objfun)) {
			error("For 'R' style optimization, MxRObjective must contain a valid R function.");
		} else {
			RObjFun = objfun;
		}
		funobj = F77_SUB(callRObjFun);
		PROTECT(confun = GET_SLOT(objective, install("constraint")));
		if(isFunction(confun)) {	  		  		// If confun is a function, call R.
				funcon = F77_SUB(callRConFun);
				RConFun = confun;
		} else {									// Otherwise, assume no constraints.
			funcon = F77_SUB(noConFun);
			if(confun != R_NilValue) Rprintf("Constraint function is not a function.  Assuming no constraint function.");
		} 											// Someday, we may want to add other options.
		
		UNPROTECT(2);
	} else {
		error("Optimization function type %s not implemented on this kernel.", STRING_ELT(nameString, 1));
	}

	/* Process Free Var List */
	if(VERBOSE) { Rprintf("Processing Free Parameters.\n"); }
	omxDataMatrix dm;
	freeVarList = (freeVar*) R_alloc (sizeof ( freeVar ), n);				// Data for replacement of free vars
	for(k = 0; k < n; k++) {
		PROTECT(nextVar = VECTOR_ELT(varList, k));
		freeVarList[k].numLocations = length(nextVar);
		freeVarList[k].location = (double**) R_alloc(sizeof(double*), length(nextVar));
		if(DEBUG) { Rprintf("Free parameter %d: %d locations\n", k, length(nextVar)); }
		for(l = 0; l < freeVarList[k].numLocations; l++) {
			PROTECT(nextLoc = VECTOR_ELT(nextVar, l));
			double* theVarList = REAL(nextLoc);			// These come through as doubles.
			int theMat = (int)theVarList[0];			// Matrix is zero-based indexed.
			int theRow = (int)theVarList[1] - 1;		// Row is one-based.
			int theCol = (int)theVarList[2] - 1;		// Column is one-based.
			freeVarList[k].location[l] = matrixList[theMat].locationOfElement(theRow, theCol);
			UNPROTECT(1);
		}
		UNPROTECT(1);
	}
	UNPROTECT(1); // nameString

	if(VERBOSE) { Rprintf("Processed.\n"); }
	
	/* Set up Optimization Memory Allocations */
	
	SEXP minimum, estimate, gradient, hessian, code, iterations, ans, names;
	
	PROTECT(ans = allocVector(VECSXP, 6));
	PROTECT(names = allocVector(STRSXP, 6));
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(iterations = NEW_NUMERIC(1));
	
	/* Initialize Scalar Variables. */	
	nclin = 0;						// No constraints.
	ncnln = 0;						// None.  At all.
	
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

	if(n == 0) {			// Special Case for the evaluation-only condition
		if(DEBUG) { Rprintf("No free parameters.  Avoiding Optimizer Entirely.\n"); }
		int mode = 0, nstate = 1;
		f = 0;
		double* x = NULL, *g = NULL;
		funobj(&mode, &n, x, &f, g, &nstate);
		inform = 0;
		iter = 0;
	} else {
	/* Set up actual run */
	
		/* Set min and max limits */
		for(k = 0; k < nctotl; k++) {
			bl[k] = EPSILON; 				// No negative covariances.  -Infinity'd be -10^20.
			bu[k] = 2;						// Infinity would be at 10^20.
		}
	
	
		/* Initialize Starting Values */
		Rprintf("--------------------------\n");	
		Rprintf("Starting Values (%d) are:\n", n);
		for(k = 0; k < n; k++) {
			x[k] = REAL(startVals)[k];
			Rprintf("%d: %f\n", k, x[k]);
		}
		if(DEBUG) {
			Rprintf("--------------------------\n");	
			Rprintf("Setting up optimizer...");
		}
		
	/* 	Set NPSOL options  (Maybe separate this into a different function?)
	
		F77_CALL(npoptn)(char* string, int *length)
	
		String is of the form:
			'Option = Value'
	
	*/
	
		/* Options That Change The Optimizer */
		strcpy(option, "Step Limit 0");					// At 0, defaults to 2.0
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Derivative Level 0");
		F77_CALL(npoptn)(option, strlen(option));	
		strcpy(option, "Hessian Yes");	
		F77_CALL(npoptn)(option, strlen(option));
	
		/* Output Options */
		strcpy(option, "Print level 0");  				// 0 = No Output, 20=Verbose
		F77_CALL(npoptn)(option, strlen(option));
		strcpy(option, "Minor print level 0");			// 0 = No Output, 20=Verbose
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Print file 0");
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Summary file 0");
		F77_CALL(npoptn)(option, strlen(option));
		sprintf(option, "Verify level 3");
		F77_CALL(npoptn)(option, strlen(option));
		
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
	
		if(DEBUG) {
			Rprintf("Set.\n");
		}
	
		F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, (void*)funcon,
		 				(void*)funobj, &inform, &iter, istate, c, cJac, clambda, &f, g, R,
		 				x, iw, &leniw, w, &lenw);
	}
	
	/* Store outputs for return */
	
	REAL(minimum)[0] = f;	
	
	est = REAL(estimate);  // Aliases to avoid repeated function calls.
	grad = REAL(gradient);
	hess = REAL(hessian); 
	
	for(k = 0; k < n; k++) {		// Do these with memcpy instead, probably.
		REAL(estimate)[k] = x[k];
		REAL(gradient)[k] = g[k];
		for(l = 0; l < n; l++) {
			REAL(hessian)[k*n + l] = R[k*n + l];  // This is ok, because they're both in Col-Major already.
		}
	}
	
	REAL(code)[0] = inform;
	REAL(iterations)[0] = iter;
	
	SET_STRING_ELT(names, 0, mkChar("minimum"));
	SET_STRING_ELT(names, 1, mkChar("estimate"));
	SET_STRING_ELT(names, 2, mkChar("gradient"));
	SET_STRING_ELT(names, 3, mkChar("hessian"));
	SET_STRING_ELT(names, 4, mkChar("code"));
	SET_STRING_ELT(names, 5, mkChar("iterations"));
	
	SET_VECTOR_ELT(ans, 0, minimum);
	SET_VECTOR_ELT(ans, 1, estimate);
	SET_VECTOR_ELT(ans, 2, gradient);
	SET_VECTOR_ELT(ans, 3, hessian);
	SET_VECTOR_ELT(ans, 4, code);
	SET_VECTOR_ELT(ans, 5, iterations);
	namesgets(ans, names);
		
	Rprintf("Inform Value: %d\n", inform);

	UNPROTECT(8);						// Unprotect NPSOL Parameters
	
	Rprintf("--------------------------\n");
	
	return(ans);

}


/****** Objective Functions ******** (put in separate files) */
void F77_SUB(callRObjFun) 
	(	int* mode, int* n, double* x, 
		double* f, double* g, int* nstate )
{
/* 
	Generic Objective Function.
	For now, call R. 
*/

	SEXP theCall, theVars, theReturn;
	int k;
	double *vars;

	ncalls++;
	//if(ncalls > 5) *mode = -1;
	//fprintf(stderr, "Call: %d\n", ncalls);
	Rprintf("--------------------------\n");
	Rprintf("Call: %d\n", ncalls);
	Rprintf("Estimates:");
	for(k = 0; k < *n; k++) {
		Rprintf("%d: %3.17f", k, x[k]);
	}
	Rprintf("\n");

	PROTECT(theCall = allocList(2));
	PROTECT(theVars = allocVector(REALSXP, *n));
	vars = REAL(theVars);
	SET_TYPEOF(theCall, LANGSXP);
	SETCAR(theCall, RObjFun);

	for(k = 0; k < *n; k++) {
		vars[k] = x[k];
	}
	SETCADR(theCall, theVars);
	
	PROTECT(theReturn = eval(theCall, R_GlobalEnv));
	*f = REAL(AS_NUMERIC(theReturn))[0];
	UNPROTECT(3);

	Rprintf("Obj. Function Value: %3.17f\n", *f);

	return;
}

//	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1, A, 2, S, 2, 1, F, 2);	

void F77_SUB(covObjFun)
	(	int* mode, int* n, double* x, 
		double* f, double* g, int* nstate )
{
/* Objective function.
	For covariance-style optimization
*/

	int k, l;
	omxDataMatrix *A, *F, *S, *R;

	handleFreeVarList(x, *n);

	/* Aliases for ease of use and lack of dereference adds */
	A = &(matrixList[0]);
	S = &(matrixList[1]);
	F = &(matrixList[2]);
	R = &(matrixList[3]);
	
	A->recompute();
	S->recompute();
	F->recompute();
	R->recompute();
	
	if(DEBUG) { Rprintf("A is at %d\n", A->data); A->print("A"); }
	if(DEBUG) { Rprintf("S is at %d\n", S->data); S->print("S"); }
	if(DEBUG) { Rprintf("F is at %d\n", F->data); F->print("F"); }
	if(DEBUG) { Rprintf("R is at %d\n", R->data); R->print("R"); }


	// Z will be (I-A)^-1, C will be our running total.
	omxDataMatrix Z(*A);
	omxDataMatrix Y(I);
	omxDataMatrix C(I);

	const char NoTrans = 'n';
	const char Trans = 'T';
	double MinusOne = -1.0;
	double Zero = 0.0;
	double One = 1.0;
	double Two = 2.0;	
	int ipiv[I.rows];
	double work[5];
	int lwork = 5;
	
	/* Z = (I-A)^-1 */
	if(DEBUG) { Rprintf("Beginning Objective Calculation.\n"); }
	
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(I.cols), &(I.rows), &(Z.rows), &One, I.data, &(I.cols), I.data, &(I.cols), &MinusOne, Z.data, &(Z.cols));
	F77_CALL(dgetrf)(&(Z.cols), &(Z.rows), Z.data, &(Z.cols), ipiv, &l);
	F77_CALL(dgetri)(&(Z.rows), Z.data, &(Z.cols), ipiv, work, &lwork, &l);

	/* C = FZSZ'F' */ // There MUST be an easier way to do this.  I'm thinking matrix class.
	F77_CALL(dgemm)(&Trans, &Trans, &(Z.cols), &(Z.rows), &(F->cols), &One, Z.data, &(Z.cols), F->data, &(F->cols), &Zero, Y.data, &(Y.cols)); 		// C = ...Z'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(S->cols), &(S->rows), &(Y.rows), &One, S->data, &(S->cols), Y.data, &(Y.cols), &Zero, C.data, &(C.cols)); 	// C = ..SZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(Z.cols), &(Z.rows), &(C.rows), &One, Z.data, &(Z.cols), C.data, &(C.cols), &Zero, Y.data, &(Y.cols)); 	// C = .ZSZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(F->cols), &(F->rows), &(Y.rows), &One, F->data, &(F->cols), Y.data, &(Y.cols), &Zero, C.data, &(C.cols));	// C = FZSZ'F'
	
	/* Val = sum(diag(tempCov %*% solve(PredictedCov))) + log(det(PredictedCov)) */
	/* Alternately, Val = sum (tempCov .* PredictedCov^-1) + log(det(PredictedCov)) */	
	
	F77_CALL(dgetrf)(&(C.cols), &(C.rows), C.data, &(C.cols), ipiv, &l);

	if(DEBUG) { Rprintf("Info on Invert: %d\n", l); }
	double det = 1.0;
	double sum = 0;

	for(k = 0; k < C.cols; k++) { 		// |A| is the sum of the diagonal elements of U from the LU factorization.
		det *= C.data[k+C.rows*k];  	// Normally, we'd need to worry about transformations made during LU, but
	}									// we're safe here because the determinant of a covariance matrix > 0.

	if(DEBUG) { Rprintf("Determinant of F(I-A)^-1*S*(I-A)^1'*F': %f\n", det); }
	det = log(det);
	
	F77_CALL(dgetri)(&(C.rows), C.data, &(C.cols), ipiv, work, &lwork, &l);
	
	for(k = 0; k < (C.cols * C.rows); k++) { 
		sum += C.data[k] * R->data[k];
	}

	*f = (sum + det);

	if(VERBOSE) {
		Rprintf("Calculated Obj. Function Value: %3.17f\n", *f);
	}
	
	return;
}

void F77_SUB(FIMLObjFun)
	(int* mode, int* n, double* x, 
		double* f, double* g, int* nstate )
{
	if(DEBUG) { Rprintf("Beginning FIML Evaluation.\n");}
	// Requires: Data, means, covariances.

	SEXP matrixDims;
	int *dimList;
	double sum;
	char u = 'U';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	int onei = 1;
	int mainDist = 0;
	double Q = 0.0;
	double logDet = 0;
	int	*toRemove;
	int nextRow, nextCol, numCols, numRemoves;

	handleFreeVarList(x, *n);

	omxDataMatrix *mainMatrix, *expected, smallRow(1, cov.cols, TRUE), smallCov(cov.rows, cov.cols, TRUE), RCX(1, dataRows.cols, TRUE);

	/* These will be filled in somwhere else.  Check to see that we're not breaking stuff here.*/
	expected = &cov;
	mainMatrix = &dataRows;
	
	smallCov.alias(cov);

//	mainDist = mainMatrix->rows;

	sum = 0.0; // = log(2 * M_PI) * mainMatrix->cols * mainMatrix->rows;

	for(int row = 0; row < mainMatrix->rows; row++) {
		logDet = 0.0;
		Q = 0.0;

		/** HANDLE MISSINGNESS HERE **/
		// Note:  This really aught to be done using a matrix multiply.  Why isn't it?
		numCols = 0;
		numRemoves = 0;
		for(int j = 0; j < mainMatrix->cols; j++) {
			if(ISNA(mainMatrix->element(row, j))) {
				numRemoves++;
			}
		}
		smallRow.resize(1, cov.cols - numRemoves, TRUE);
		toRemove = (int*)malloc(cov.cols * sizeof(int));
		for(int j = 0; j < mainMatrix->cols; j++) {
			if(ISNA(mainMatrix->element(row, j))) {
				toRemove[j] = 1;
				continue;
			} else {
				smallRow.setElement(0, numCols, mainMatrix->element(row, j));
				toRemove[j] = 0;
				numCols++;
			}
		}
		if(numCols==0) continue;
		smallCov.reset();
		if(DEBUG) { Rprintf("Reducing by (%d, %d): %d, %d, %d\n", numCols, numRemoves, toRemove[0], toRemove[1], toRemove[2]); } 
		smallCov.removeRowsAndColumns(numRemoves, numRemoves, toRemove, toRemove);
		/** MISSINGNESS HANDLED **/
		// Rprintf("Row: %d -- sum is %f, Q is %f, logDet is %f", row, sum, Q, logDet);
		if(DEBUG) { smallRow.print("Next Row:"); smallCov.print("Covariance is:"); }

		F77_CALL(dpotrf)(&u, &(smallCov.rows), smallCov.data, &(smallCov.cols), &info);
		if(info != 0) error("Covariance Matrix is not positive-definite.");
		for(int diag = 0; diag < (smallCov.rows); diag++) {
			logDet += log(fabs(smallCov.data[diag + (diag * smallCov.rows)]));
		}
		logDet *= 2.0;
		F77_CALL(dpotri)(&u, &(smallCov.rows), smallCov.data, &(smallCov.cols), &info);
		if(info != 0) error("Cannot invert covariance matrix.");
		F77_CALL(dsymv)(&u, &(smallCov.rows), &oned, smallCov.data, &(smallCov.cols), smallRow.data, &onei, &zerod, RCX.data, &onei);
		for(int col = 0; col < smallRow.cols; col++) {
			Q += RCX.data[col] * smallRow.data[col];
		}
		sum += logDet + Q + (log(2 * M_PI) * smallRow.cols);
//		if(DEBUG) { Rprintf("Row %d: Q is %f, logDet is %f, Sum is now %f\n", row+1, Q, logDet, sum); }
	}
	
	*f = sum;
	return;

}

/* (Non)Linear Constraint Functions */

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

/****** HELPER FUNCTIONS ******* (put in separate file) */
/* Sub Free Vars Into Appropriate Slots */
void handleFreeVarList(double* x, int n) {
	
	if(DEBUG) {Rprintf("Processing Free Parameter Estimates.\n");}
	ncalls++;
	if(VERBOSE) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d\n", ncalls);
		Rprintf("Estimates: [");
		for(int k = 0; k < n; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
	}

	/* Fill in Free Var Estimates */
	for(int k = 0; k < n; k++) {
		if(DEBUG) { Rprintf("%d: %f - %d\n", k,  x[k], freeVarList[k].numLocations); }
		for(int l = 0; l < freeVarList[k].numLocations; l++) {
			*(freeVarList[k].location[l]) = x[k];
			if(DEBUG) { Rprintf("Setting location:%ld to value %f for var %d\n", freeVarList[k].location[l], x[k], k); }
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
