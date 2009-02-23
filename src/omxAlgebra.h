/***********************************************************
* 
*  omxMatrix.h
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Contains header information for the omxMatrix class
*   omxDataMatrices hold necessary information to simplify
* 	dealings between the OpenMX back end and BLAS.
*
**********************************************************/

#ifndef _OMXALGEBRA_H_
#define _OMXALGEBRA_H_

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifdef DEBUGMX
#define OMX_DEBUG 1
#else
#define OMX_DEBUG 0
#endif /* DEBUGMX */

typedef struct omxAlgebra omxAlgebra;
typedef struct omxMatrix omxMatrix;

#include "omxSymbolTable.h"

struct omxAlgebra {						// A matrix
										//TODO: Improve encapsulation
/* Fields unique to Algebras */
	void* funWrapper;					// Wrapper for the algebra function itself
	omxMatrix** args;					// Arguments to the above function
	int numArgs;						// Length of args

	omxMatrix* myMatrix;				// The matrix populated by this algebra

};

/* Initialize and Destroy */
	omxMatrix* omxInitAlgebra(omxAlgebra *oa);							// Constructor 
	void omxInitAlgebraWithMatrix(omxAlgebra *oa, omxMatrix* om);		// Constructor (with matrix)
	void omxFreeAlgebraArgs(omxAlgebra* algebra);						// Frees all args
	omxMatrix* omxNewMatrixFromMxAlgebra(SEXP mxmat); 					// Create an Algebra from an R mxMatrix
	void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP mxmat);			// Populate an Algebra from an R mxMatrix
	omxMatrix* omxNewMatrixFromMxMatrixPtr(SEXP matrix);				// Create a matrix/algebra from a matrix pointer
	omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix* arg1, omxMatrix* arg2); // For constraints.

/* Other Functions */
	void omxFillAlgebraFromTableEntry(omxAlgebra *algebra, const omxAlgebraTableEntry* oate);
	 																	// Adjust an algebra for a table entry
	void omxAlgebraCopyAlgebra(omxAlgebra *dest, omxAlgebra *src);		// Copy across another element.  
																		// NOTE: Duplicates.
	omxMatrix* omxAlgebraParseHelper(SEXP algebraArg);

/* Algebra-specific implementations of matrix functions */
	void omxAlgebraRecompute(omxAlgebra *oa);
	void omxAlgebraCompute(omxAlgebra *oa);
	unsigned short int omxAlgebraNeedsUpdate(omxAlgebra *oa);

	void omxAlgebraPrint(omxAlgebra *source, char* d);					// Pretty-print a (small) matrix

#endif /* _OMXALGEBRA_H_ */


