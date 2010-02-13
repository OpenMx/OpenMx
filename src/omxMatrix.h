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
 *
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

#ifndef _OMXMATRIX_H_
#define _OMXMATRIX_H_

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"

typedef struct omxMatrix omxMatrix;

#include "omxAlgebra.h"
#include "omxObjective.h"
#include "omxState.h"


struct omxMatrix {						// A matrix
										//TODO: Improve encapsulation
/* Actually Useful Members */
	int rows, cols;						// Matrix size  (specifically, its leading edge)
	double* data;						// Actual Data Pointer
	unsigned short colMajor;			// and column-majority.

/* For Memory Administrivia */
	unsigned short localData;			// If data has been malloc'd, and must be freed.

/* For aliased matrices */				// Maybe this should be a subclass, as well.
	omxMatrix* aliasedPtr;				// For now, assumes outside data if aliased.
	unsigned short originalColMajor;	// Saved for reset of aliased matrix.
	unsigned short originalRows;		// Saved for reset of aliased matrix.
	unsigned short originalCols;		// Saved for reset of aliased matrix.

/* For BLAS Multiplication Speedup */ 	// TODO: Replace some of these with inlines or macros.
	const char* majority;				// Filled by compute(), included for speed
	const char* minority;				// Filled by compute(), included for speed
	int leading;						// Leading edge; depends on original majority
	int lagging;						// Non-leading edge.

/* Curent State */
	omxState* currentState;				// Optimizer State
	unsigned short isDirty;				// Retained, for historical purposes.
	int lastCompute;					// Compute Count Number at last computation
	int lastRow;						// Compute Count Number at last row update (Used for row-by-row computation only)

/* For Algebra Functions */				// At most, one of these may be non-NULL.
	omxAlgebra* algebra;				// If it's not an algebra, this is NULL.
	omxObjective* objective;			// If it's not an objective function, this is NULL.

/* For inclusion in(or of) other matrices */
	int numPopulateLocations;
	omxMatrix** populateFrom;
	int *populateFromRow, *populateFromCol;
	int *populateToRow, *populateToCol;

};

/* Initialize and Destroy */
	omxMatrix* omxInitMatrix(omxMatrix* om, int nrows, int ncols, unsigned short colMajor, omxState* os);	// Set up matrix
	void omxFreeMatrixData(omxMatrix* om);							// Release any held data.
	void omxFreeAllMatrixData(omxMatrix* om);						// Ditto, traversing argument trees

/* Matrix Creation Functions */
	omxMatrix* omxNewMatrixFromMxMatrix(SEXP matrix, omxState *state); 			// Create an omxMatrix from an R MxMatrix
	omxMatrix* omxNewIdentityMatrix(int nrows, omxState* state);				// Creates an Identity Matrix of a given size

/* Getters 'n Setters */
	double omxMatrixElement(omxMatrix *om, int row, int col);
	double omxVectorElement(omxMatrix *om, int index);
	void omxSetMatrixElement(omxMatrix *om, int row, int col, double value);
	void omxSetVectorElement(omxMatrix *om, int index, double value);
	double omxAliasedMatrixElement(omxMatrix *om, int row, int col);			// Element from unaliased form of the same matrix
	double* omxLocationOfMatrixElement(omxMatrix *om, int row, int col);
	void omxMarkDirty(omxMatrix *om);

/* Matrix Modification Functions */
	void omxZeroByZeroMatrix(omxMatrix *source);
	void omxResizeMatrix(omxMatrix *source, int nrows, int ncols,
							unsigned short keepMemory);									// Resize, with or without re-initialization
	omxMatrix* omxFillMatrixFromMxMatrix(omxMatrix* om, SEXP matrix, omxState *state); 	// Populate an omxMatrix from an R MxMatrix
	void omxProcessMatrixPopulationList(omxMatrix *matrix, SEXP matStruct);
	void omxCopyMatrix(omxMatrix *dest, omxMatrix *src);								// Copy across another matrix.

/* Function wrappers that switch based on inclusion of algebras */
	void omxPrint(omxMatrix *source, char* d); 											// Pretty-print a (small) matrix
	unsigned short int omxNeedsUpdate(omxMatrix *matrix);								// Does this need to be recomputed?
	void omxRecompute(omxMatrix *matrix);												// Recompute if needed
	void omxCompute(omxMatrix *matrix);													// Compute it whether needed or not

/* BLAS Wrappers */

	void omxDGEMM(unsigned short int transposeA, unsigned short int transposeB,		// result <- alpha * A %*% B + beta * C
					double* alpha, omxMatrix* a, omxMatrix *b, double* beta, omxMatrix* result);
	void omxDGEMV(unsigned short int transposeMat, double* alpha, omxMatrix* mat,	// result <- alpha * A %*% B + beta * C
					omxMatrix* vec, double* beta, omxMatrix*result, int* info);		// where B is treated as a vector

	void omxDSYMV(unsigned short int transposeMat, double* alpha, omxMatrix* mat,	// result <- alpha * A %*% B + beta * C
					omxMatrix* vec, double* beta, omxMatrix*result, int* info);		// only A is symmetric, and B is a vector

	void omxDSYMM(unsigned short int symmOnLeft, double* alpha, omxMatrix* a, 		// result <- alpha * A %*% B + beta * C
					omxMatrix *b, double* beta, omxMatrix* result);					// One of A or B is symmetric

	void omxDGETRF(omxMatrix* mat, int* ipiv, int* info);							// LUP decomposition of mat
	void omxDGETRI(omxMatrix* mat, int* ipiv, double* work, int* info);				// Invert mat from LUP decomposition

	void omxDPOTRF(omxMatrix* mat, int* info);										// Cholesky decomposition of mat
	void omxDPOTRI(omxMatrix* mat, int* info);										// Invert mat from Cholesky

/* Aliased Matrix Functions */
	void omxAliasMatrix(omxMatrix *alias, omxMatrix* const source);		// Allows aliasing for faster reset.
	void omxResetAliasedMatrix(omxMatrix *matrix);						// Reset to the original matrix
	void omxRemoveRowsAndColumns(omxMatrix* om, int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[]);

/* Matrix-Internal Helper functions */
	void omxMatrixCompute(omxMatrix *matrix);
	void omxPrintMatrix(omxMatrix *source, char* d);                    // Pretty-print a (small) matrix
	unsigned short int omxMatrixNeedsUpdate(omxMatrix *matrix);

#endif /* _OMXMATRIX_H_ */