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

#include "omxAlgebra.h"
#include "omxObjective.h"

#ifdef DEBUGMX
#define OMX_DEBUG 1
#else
#define OMX_DEBUG 0
#endif /* DEBUG */

struct omxMatrix {						// A matrix
										//TODO: Improve encapsulation

/* Actually Useful Members */
	int rows, cols;						// Matrix size  (specifically, its leading edge)
	double* data;						// Actual Data Pointer
	double* aliasedPtr;					// For now, assumes outside data if aliased.
	unsigned short colMajor;			// and column-majority.
	unsigned short isDirty;				// True if free params have been updated.

/* For Memory Administrivia */
	unsigned short localData;			// If data has been malloc'd, and must be freed.

/* For aliased matrices */	// Maybe this should be a subclass, as well.
	unsigned short originalColMajor;	// Saved for reset of aliased matrix.
	unsigned short originalRows;		// Saved for reset of aliased matrix.
	unsigned short originalCols;		// Saved for reset of aliased matrix.

/* For BLAS Multiplication Speedup */ 	// Maybe replace some of these with inlines or macros.
	const char* majority;				// Filled by compute(), included for speed
	const char* minority;				// Filled by compute(), included for speed
	int leading;						// Leading edge; depends on original majority
	int lagging;						// Non-leading edge.

/* For Algebra Functions */				// At most, one of these may be non-NULL.
	omxAlgebra* algebra;				// If it's not an algebra, this is NULL.
	omxObjective* objective;			// If it's not an objective function, this is NULL.
	
/* For inclusion in(or of) other matrices */
	int populateLocations;
	omxMatrix** populateFrom;
	double** populateTo;

};

/* Initialize and Destroy */
	omxMatrix* omxInitMatrix(omxMatrix* om, int nrows, int ncols, unsigned short colMajor);	// Set up matrix
	void omxFreeMatrixData(omxMatrix* om);							// Release any held data.
	void omxFreeAllMatrixData(omxMatrix* om);						// Ditto, traversing argument trees

/* Getters 'n Setters */
	double omxMatrixElement(omxMatrix *om, int row, int col);
	void omxSetMatrixElement(omxMatrix *om, int row, int col, double value);
	double* omxLocationOfMatrixElement(omxMatrix *om, int row, int col);
	void omxMarkDirty(omxMatrix *om);

/* Other Functions */
	void omxResizeMatrix(omxMatrix *source, int nrows, int ncols,
	 						unsigned short keepMemory);					// Resize, with or without re-initialization
	omxMatrix* omxNewMatrixFromMxMatrix(SEXP matrix); 					// Populate an omxMatrix from an R MxMatrix 
	void omxCopyMatrix(omxMatrix *dest, omxMatrix *src);				// Copy across another matrix.  
	
/* Aliased Matrix Functions */
	void omxAliasMatrix(omxMatrix *alias, omxMatrix* const source);		// Allows aliasing for faster reset.
	void omxResetAliasedMatrix(omxMatrix *matrix);						// Reset to the original majority and realias, if needed.
	void omxRemoveRowsAndColumns(omxMatrix* om, int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[]);

/* Matrix-Internal Helper functions */
	void omxComputeMatrixHelper(omxMatrix *matrix);
	void omxPrintMatrixHelper(omxMatrix *source, char* d);				// Pretty-print a (small) matrix
	unsigned short omxMatrixNeedsUpdate(omxMatrix *matrix);

/* Function wrappers that switch based on inclusion of algebras */
	void omxPrintMatrix(omxMatrix *source, char* d); 					// Pretty-print a (small) matrix
	unsigned short omxNeedsUpdate(omxMatrix *matrix);					
	void omxRecomputeMatrix(omxMatrix *matrix);
	void omxComputeMatrix(omxMatrix *matrix);

#endif /* _OMXMATRIX_H_ */
