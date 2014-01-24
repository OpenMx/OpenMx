/*
 *  Copyright 2007-2012 The OpenMx Project
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
#include "omxDefines.h"

typedef struct omxAlgebra omxAlgebra;

#include "omxSymbolTable.h"
#include "omxState.h"

struct omxAlgebra {						// A matrix
										//TODO: Improve encapsulation
/* Fields unique to Algebras */
	void* funWrapper;					// Wrapper for the algebra function itself
	omxMatrix** args;					// Arguments to the above function
	int numArgs;						// Length of args

	omxMatrix* matrix;				// The matrix populated by this algebra
    SEXP sexpAlgebra;               // The SEXP MxAlgebra object
};

/* Initialize and Destroy */
	omxMatrix* omxInitAlgebra(omxAlgebra *oa, omxState* os);			// Constructor 
	void omxInitAlgebraWithMatrix(omxAlgebra *oa, omxMatrix* om);		// Constructor (with matrix)
	void omxFreeAlgebraArgs(omxAlgebra* algebra);						// Frees all args
	omxMatrix* omxNewMatrixFromMxAlgebra(SEXP mxmat, omxState* os, const char *name);	// Create an Algebra from an R mxMatrix
	void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP mxmat, const char* name);	// Populate an Algebra from an R mxMatrix
	omxMatrix* omxNewMatrixFromMxIndex(SEXP matrix, omxState* os);		// Create a matrix/algebra from a matrix pointer
	omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix** args, int numArgs, omxState* os); // For constraints.

/* Other Functions */
	void omxFillAlgebraFromTableEntry(omxAlgebra *algebra, const omxAlgebraTableEntry* oate);
	 																	// Adjust an algebra for a table entry
	void omxAlgebraCopyAlgebra(omxAlgebra *dest, omxAlgebra *src);		// Copy across another element.  
																		// NOTE: Duplicates.
	omxMatrix* omxAlgebraParseHelper(SEXP algebraArg, omxState* os, const char *name);

/* Algebra-specific implementations of matrix functions */
	void omxAlgebraRecompute(omxAlgebra *oa);
	void omxAlgebraCompute(omxAlgebra *oa);
	int omxAlgebraNeedsUpdate(omxAlgebra *oa);
	void omxDuplicateAlgebra(omxMatrix *tgt, omxMatrix* src, omxState* tgtState);

	void omxAlgebraPrint(omxAlgebra *source, char* d);					// Pretty-print a (small) matrix

#endif /* _OMXALGEBRA_H_ */


