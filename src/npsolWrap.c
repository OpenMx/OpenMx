#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <stdio.h> 

#define M(y,z) m[z][y]
#define EPSILON 0.00000000001
#define TRUE 1
#define FALSE 0

/* Structure definitions for object evaluation */
typedef struct {			// Free Variables
	int numLocations;		
	double** location;		// And where they go.
} freeVar;

typedef struct {			// A matrix
	int rows, cols;			// and it size  (specifically, its leading edge)
	double* data;			// We'll assume column-major order.
} dataMatrix;

/* NPSOL-related functions */
extern void F77_SUB(npoptn)(char* string, int length);
extern void F77_SUB(npsol)(int *n, int *nclin, int *ncnln, int *ldA, int *ldJ, int *ldR, double *A,	double *bl,	double *bu, void* funcon, void* funobj, 
						int *inform, int *iter, int *istate, double *c, double *cJac, double *clambda, double *f, double *g, double *R,				
						double *x, int *iw,	int *leniw, double *w, int *lenw);
						
/* Objective Function Variants */
void F77_SUB(callRObjFun)(int* mode, int* n, double* x, double* f, double* g, int* nstate); // Call an R function
void F77_SUB(covObjFun)(int* mode, int* n, double* x, double* f, double* g, int* nstate);	// Covariance Fit

/* Constraint Function Variants */
void F77_SUB(noConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); // No constraints
void F77_SUB(callRConFun)(int *mode, int *ncnln, int *n, int *ldJ, int *needc, double *x, double *c, double *cJac, int *nstate); // Call R for constraints

/* Functions for Export */
SEXP callNPSOL(SEXP opType, SEXP startVals, SEXP bounds, SEXP matList, SEXP varList, SEXP objfun, SEXP confun, SEXP rho);  // Calls NPSOL.  Duh.

/* Helper functions */
SEXP getListElement(SEXP list, const char *str); 				// Gets the value named str from SEXP list.  From "Writing R Extensions"
SEXP getVar(SEXP str, SEXP env);								// Gets the object named str from environment env.  From "Writing R Extensions"
void printMatrix(char* d, dataMatrix* m);						// Pretty-print a (small) matrix
dataMatrix* allocDataMatrix(int ncols, int nrows);				// Allocate a data matrix of size (nrows x ncols)
void fillDataMatrixFromMatrix(dataMatrix *dM, SEXP matrix); // Populate a data matrix object with the values of an R matrix
void fillDataMatrixFromMxMatrix(dataMatrix *dM, SEXP mxMatrix); // Populate a data matrix to represent the $values field of an MxMatrix object

/* Globals for function evaluation */
SEXP RObjFun, RConFun;		// Pointers to the functions NPSOL calls
SEXP env;					// Environment for evaluation and object hunting
int ncalls;					// For debugging--how many calls?
dataMatrix* matrixList;		// Data matrices and their data.
freeVar* freeVarList;		// List of Free Variables and where they go.

/* Globals for Covariance Evaluation */
dataMatrix I, *Z = NULL, *C = NULL, *Y=NULL;

/* Set up R .Call info */
R_CallMethodDef callMethods[] = { 
{"callNPSOL", &callNPSOL, 8}, 
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
SEXP callNPSOL(SEXP opType, SEXP startVals, SEXP bounds, SEXP matList, SEXP varList, SEXP objfun, SEXP confun, SEXP rho) {
	// For now, assume no constraints.
	
	ncalls = 0;
	
	if(!isEnvironment(rho)) error("'rho' must be an environment.");
	if(!isString(opType)) error("'opType must be 'cov', 'R', or 'FIML'");
	if(!isVector(matList)) error ("matList must be a list");

	int n;
	n = length(startVals);
	
	SEXP minimum, estimate, gradient, hessian, code, iterations, ans, names;
	
	PROTECT(ans = allocVector(VECSXP, 6));
	PROTECT(names = allocVector(STRSXP, 6));
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, n));
	PROTECT(gradient = allocVector(REALSXP, n));
	PROTECT(hessian = allocMatrix(REALSXP, n, n));
	PROTECT(code = NEW_NUMERIC(1));
	PROTECT(iterations = NEW_NUMERIC(1));
	
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
	
/* Store globals for function evaluation  (Do Elsewhere!) */

	if(strncmp(CHAR(STRING_ELT(opType, 0)), "cov", length(opType)) == 0) { // Covariance-style optimization.  This should be a table lookup.
		Rprintf("Using covariance optimization.\n");
		funobj = F77_SUB(covObjFun);
		int *dimList;
		
		/* Store Data from MxMatrices */
	  	// This part is just for S3 DataMatrix Stuff
#ifdef USEMXMATRIX
		matrixList = (dataMatrix*) R_alloc(sizeof(dataMatrix), length(matList));		// Stores links to data/covariance matrices
		for(k = 0; k < length(matList)-1; k++) {										// Last is Covariance Matrix, not an MxMatrix
			PROTECT(nextLoc = VECTOR_ELT(matList, k));									// TODO: Find out if even this one duplicates the matrix.
			fillDataMatrixFromMxMatrix(&(matrixList[k]), nextLoc);
			UNPROTECT(1);
		} 
#endif
		/* Retrieve All Matrices From the MatList */
		matrixList = (dataMatrix*) R_alloc(sizeof(dataMatrix), length(matList));		// Stores links to data/covariance matrices
		for(k = 0; k < length(matList)-1; k++) {										// Last is Covariance Matrix, not an MxMatrix
			PROTECT(nextLoc = VECTOR_ELT(matList, k));									// TODO: Find out if even this one duplicates the matrix.
			fillDataMatrixFromMatrix(&(matrixList[k]), nextLoc);
			UNPROTECT(1);
		}
		
		/* Last Matrix is the Covariance Matrix */
		PROTECT(nextLoc = VECTOR_ELT(matList, k));
		fillDataMatrixFromMatrix(&(matrixList[k]), nextLoc);
		UNPROTECT(1);
		
		/* Identity Matrix, Size Of A */
		I.rows = matrixList[0].rows;
		I.cols = matrixList[0].cols;
		I.data = (double*) R_alloc(I.rows * I.cols, sizeof( double ));
		for(k =0; k < I.rows; k++) {
			for(l = 0; l < I.cols; l++) {
				if(l == k) {
					I.data[l*I.rows + k] = 1;
				} else {
					I.data[l*I.rows + k] = 0;
				}
			}
		}
		
		/* Data matrices for Covariance Calculation */
		Z = allocDataMatrix(I.rows, I.cols);
		C = allocDataMatrix(I.rows, I.cols);
		Y = allocDataMatrix(I.rows, I.cols);
		
		Rprintf("Processing Free Parameters.\n"); //:::
		
		/* Store Free Vars */
//		dataMatrix dm;
		freeVarList = (freeVar*) R_alloc (sizeof ( freeVar ), n);				// Data for replacement of free vars
		for(k = 0; k < n; k++) {
			PROTECT(nextVar = VECTOR_ELT(varList, k));
			freeVarList[k].numLocations = length(nextVar);
			freeVarList[k].location = (double**) R_alloc(sizeof(double*), length(nextVar));
			Rprintf("Top: length %d\n", length(nextVar)); //:::
			for(l = 0; l < freeVarList[k].numLocations; l++) {
				PROTECT(nextLoc = VECTOR_ELT(nextVar, l));
				double* theVarList = REAL(nextLoc);			// These come through as doubles.
				int theMat = (int)theVarList[0];			// Matrix is zero-based indexed.
				int theRow = (int)theVarList[1] - 1;		// Row is one-based.
				int theCol = (int)theVarList[2] - 1;		// Column is one-based.
				int theIndex = theRow + matrixList[theMat].rows * theCol;
				freeVarList[k].location[l] = &(matrixList[theMat].data[theIndex]);
				Rprintf("FreeVar (%ld, %ld, %ld) at location %ld, A:%ld, S:%ld, F:%ld, R:%ld\n", theMat, theRow, theCol, freeVarList[k].location[l], matrixList[0].data, &(matrixList[1].data[0]), &(matrixList[2].data[0]), &(matrixList[3].data[0])); //:::
				Rprintf("Double: %ld, Float: %ld, Int: %ld\n", sizeof(double), sizeof(float), sizeof(int));
				UNPROTECT(1);
			}
			UNPROTECT(1);
		}
		funobj = F77_SUB(covObjFun);
		Rprintf("Processed.\n");
	} else {													// For all others, assume 'R' function
		if(!isFunction(objfun)) {
			error("For 'R' style optimization, 'objfun' must be a function.  For others, change opType.");
		} else {
			RObjFun = objfun;
			env = rho;
		}
		funobj = F77_SUB(callRObjFun);
	}	
	
	if(isFunction(confun)) {	  		  		// If confun is a function, call R.
			funcon = F77_SUB(callRConFun);
			RConFun = confun;
	} else {									// Otherwise, assume no constraints.
		funcon = F77_SUB(noConFun); 	
	} 											// Someday, we may want to add other options.
		
	
/* Set up actual run */

	/* Set min and max limits */
	for(k = 0; k < nctotl; k++) {
		bl[k] = EPSILON; 				// No negative covariances.  -Infinity'd be -10^20.
		bu[k] = 2;						// Infinity would be at 10^20.
	}


	/* Initialize Starting Values */
	Rprintf("--------------------------\n");	
	Rprintf("Starting Values are:\n");
	for(k = 0; k < n; k++) {
		x[k] = REAL(startVals)[k];
		Rprintf("%d: %f\n", k, x[k]);
	}
	
	
/* 	Set NPSOL options  (Maybe separate this into a different function?)

	F77_CALL(npoptn)(char* string, int *length)

	String is of the form:
		'Option = Value'

*/

	strcpy(option, "Print level 20");
	F77_CALL(npoptn)(option, strlen(option));
	strcpy(option, "Minor print level 20");
	F77_CALL(npoptn)(option, strlen(option));
	strcpy(option, "Step Limit .08");
	F77_CALL(npoptn)(option, strlen(option));
	strcpy(option, "Derivative Level 0");
	F77_CALL(npoptn)(option, strlen(option));	
	strcpy(option, "Hessian Yes");
	F77_CALL(npoptn)(option, strlen(option));
	sprintf(option, "Verify level 3");
	F77_CALL(npoptn)(option, strlen(option));
//	sprintf(option, "Print file 2");
//	F77_CALL(npoptn)(option, strlen(option));
	
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
	
	F77_CALL(npsol)(&n, &nclin, &ncnln, &ldA, &ldJ, &ldR, A, bl, bu, funcon,
	 				funobj, &inform, &iter, istate, c, cJac, clambda, &f, g, R,
	 				x, iw, &leniw, w, &lenw);
	
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

	if(strncmp(CHAR(STRING_ELT(opType, 0)), "cov", length(opType)) == 0) {
	//	UNPROTECT(length(matList)); 	// Unprotect Data Matrices
	}
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
	dataMatrix *A, *F, *S, *R;
	const int VERBOSE = 1;

	ncalls++;
	if(VERBOSE) {
		Rprintf("--------------------------\n");
		Rprintf("Call: %d\n", ncalls);
		Rprintf("Estimates: [");
		for(k = 0; k < *n; k++) {
			Rprintf(" %f", x[k]);
		}
		Rprintf("] \n");
	}
	
	/* Fill in Free Var Estimates */
	for(k = 0; k < *n; k++) {
		Rprintf("%d: %f - %d", k,  x[k], freeVarList[k].numLocations);
		for(l = 0; l < freeVarList[k].numLocations; l++) {
			*(freeVarList[k].location[l]) = x[k];
//			Rprintf("Setting location:%ld to value %f for var %d\n", freeVarList[k].location[l], x[k], k);
		}
	}

	/* Aliases for ease of use and lack of dereference adds */
	A = &(matrixList[0]);
	S = &(matrixList[1]);
	F = &(matrixList[2]);
	R = &(matrixList[3]);
	
	printMatrix("A", A); //:::
	printMatrix("S", S); //:::
	printMatrix("F", F); //:::
	printMatrix("R", R); //:::
	
	// Z will be (I-A)^-1, C will be our running total.
	memcpy(Z->data, A->data, A->rows * A->cols * sizeof(double));
	memcpy(C->data, I.data, I.rows * I.cols * sizeof(double));

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
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(I.cols), &(I.rows), &(Z->rows), &One, I.data, &(I.cols), I.data, &(I.cols), &MinusOne, Z->data, &(Z->cols));
	F77_CALL(dgetrf)(&(Z->cols), &(Z->rows), Z->data, &(Z->cols), ipiv, &l);
	F77_CALL(dgetri)(&(Z->rows), Z->data, &(Z->cols), ipiv, work, &lwork, &l);

	/* C = FZSZ'F' */ // There MUST be an easier way to do this.  I'm thinking matrix class.
	F77_CALL(dgemm)(&Trans, &Trans, &(Z->cols), &(Z->rows), &(F->cols), &One, Z->data, &(Z->cols), F->data, &(F->cols), &Zero, Y->data, &(Y->cols)); 		// C = ...Z'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(S->cols), &(S->rows), &(Y->rows), &One, S->data, &(S->cols), Y->data, &(Y->cols), &Zero, C->data, &(C->cols)); 	// C = ..SZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(Z->cols), &(Z->rows), &(C->rows), &One, Z->data, &(Z->cols), C->data, &(C->cols), &Zero, Y->data, &(Y->cols)); 	// C = .ZSZ'F'
	F77_CALL(dgemm)(&NoTrans, &NoTrans, &(F->cols), &(F->rows), &(Y->rows), &One, F->data, &(F->cols), Y->data, &(Y->cols), &Zero, C->data, &(C->cols));	// C = FZSZ'F'
	
	/* Val = sum(diag(tempCov %*% solve(PredictedCov))) + log(det(PredictedCov)) */
	/* Alternately, Val = sum (tempCov .* PredictedCov^-1) + log(det(PredictedCov)) */	
	
	F77_CALL(dgetrf)(&(C->cols), &(C->rows), C->data, &(C->cols), ipiv, &l);

//	Rprintf("Info: %d\n", l); 
	double det = 1.0;
	double sum = 0;

	for(k = 0; k < C->cols; k++) { 		// |A| is the sum of the diagonal elements of U from the LU factorization.
		det *= C->data[k+C->rows*k];  	// Normally, we'd need to worry about transformations made during LU, but
	}									// we're safe here because the determinant of a covariance matrix > 0.

//	Rprintf("Determinant: %f\n", det); //:::
	det = log(det);
	
	F77_CALL(dgetri)(&(C->rows), C->data, &(C->cols), ipiv, work, &lwork, &l);
	
	for(k = 0; k < (C->cols * C->rows); k++) { 
		sum += C->data[k] * R->data[k];
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
/** Not Yet Implemented **/
	
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

/****** HELPER FUNCTIONS ******* (put in separate file */
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

void printMatrix(char* data, dataMatrix *mat) {
	int j, k, rows, cols;
	
	rows = mat->rows;
	cols = mat->cols;
	Rprintf("%s: (%d x %d)\n", data, rows, cols);
	
	for(j = 0; j < rows; j++) {
		for(k = 0; k < cols; k++) {
			Rprintf("\t%3.6f", mat->data[k*rows+j]);
		}
		Rprintf("\n");
	}
}

dataMatrix* allocDataMatrix(int nrows, int ncols) {
	dataMatrix* mat = (dataMatrix*) R_alloc(1, sizeof(dataMatrix));
	mat->rows = nrows;
	mat->cols = ncols;
	mat->data = (double*) R_alloc(nrows * ncols, sizeof( double ));
	return mat;
}
/* This one doesn't work.  Sorry. 
dataMatrix* allocDataMatrixFromMxMatrix(SEXP mxMatrix) {
	dataMatrix* mat = (dataMatrix*) R_alloc(1, sizeof(dataMatrix));
	mat->rows = nrows;
	mat->cols = ncols;
	mat->data = (double*) R_alloc(nrows * ncols, sizeof( double ));
	return mat;
}
*/
void fillDataMatrixFromMxMatrix(dataMatrix *dM, SEXP mxMatrix) {
/* Populates the fields of a dataMatrix with details from an mxMatrix. */ 

	Rprintf("fillDataMatrixFromMxMatrix() Should never be called.\n");

	SEXP objectEnv, valueMatrix, nameString, matrixDims;
	int* dimList;
	double* ans; 

	PROTECT(nameString = getAttrib(mxMatrix, install("class")));
	
	if(STRING_ELT(nameString,0) == R_NilValue) error("Not an MxMatrix object.\n");  // TODO: Need better class-checking here.
	
    SET_STRING_ELT(nameString, 0, mkChar(".env"));				
	PROTECT(objectEnv = getAttrib(mxMatrix, nameString));							// May duplicate.  TODO: Check.

	SET_STRING_ELT(nameString, 0, mkChar(".values"));
	PROTECT(valueMatrix = findVar(install(CHAR(STRING_ELT(nameString, 0))), objectEnv)); // May duplicate.  TODO: Check.
	dM->data = REAL(valueMatrix);													// Probably duplicates.  TODO: Need to fix.

	PROTECT(matrixDims = getAttrib(valueMatrix, R_DimSymbol));
	dimList = INTEGER(matrixDims);
	dM->rows = dimList[0];
	dM->cols = dimList[1];

	UNPROTECT(4);
	return;
}

void fillDataMatrixFromMatrix(dataMatrix *dM, SEXP matrix) {
/* Populates the fields of a dataMatrix with details from an R Matrix. */ 


	SEXP matrixDims;
	int* dimList;
								// TODO: Class-check first?
	dM->data = REAL(matrix);								// Probably duplicates.  TODO: Need to fix.

	PROTECT(matrixDims = getAttrib(matrix, R_DimSymbol));
	dimList = INTEGER(matrixDims);
	dM->rows = dimList[0];
	dM->cols = dimList[1];

	UNPROTECT(1);
	return;
}
/* Not implemented yet.
DataMatrix* dataMatrixMissingRowsAndColumns(dataMatrix *dM, int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[])
{
	if(numRowsRemoved ==0 && numColsRemoved == 0) return dM;
	/** NOT YET IMPLEMENTED **/
	
	
	
//}
