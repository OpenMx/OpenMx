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

#ifndef _OMXDATAMATRIX_H_
#define _OMXDATAMATRIX_H_

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
#endif /* DEBUG */

class omxMatrix {						// A matrix
	
protected:								//TODO: Improve encapsulation
	bool localData;						// If data has been malloc'd, and must be freed.
	unsigned short originalColMajor;	// Saved for reset.
	unsigned short originalRows;		// Saved for memory purposes
	unsigned short originalCols;		// Saved for memory purposes
//	unsigned short isReused;			// Whether or not this data should be saved.
	bool isDirty;						// True if free params have been updated.
	bool containsDefinitions;			// True if it must be recomputed at each row.
	double* aliasedPtr;					// For now, assumes outside data for original.

public:
	static const char majorityList[3];	// For BLAS

	const char* majority;				// Filled by compute();
	const char* minority;				// Filled by compute();
	int leading;						// Leading edge; depends on original majority
	bool isAlgebra;						// <sigh> True if it's an Algebra.
	
	int rows, cols;						// Matrix size  (specifically, its leading edge)
	unsigned short colMajor;			// and column-majority.
	double* data;						// Actual Data Pointer
	
/* Initialize and Destroy */
	void init();												// Null initialization.  Just in case.
	omxMatrix();											// Null Constructor.  For when we've no idea.
	omxMatrix(const omxMatrix &in);							// Copy Constructor
	omxMatrix(int ncols, int nrows, bool colMajor=true);	// Initialize data matrix of size ncols x nrows
	void init(int ncols, int nrows, bool colMajor=true);	// Initialize data matrix of size ncols x nrows
	void freeData();										// Release any held data.
	~omxMatrix();											// Free Data

/* Getters 'n Setters */
	double element(int row, int col);
	void setElement(int row, int col, double value);
	double* locationOfElement(int row, int col);
	void markDirty() { isDirty = true; }

/* Other Functions */
	void alias(omxMatrix &dM);												// Allows aliasing for faster reset.
	void resize(int nrows, int ncols, bool keepMemory=false);				// Resize, with or without re-initialization
	void print(char* d);													// Pretty-print a (small) matrix
	void fillFromMatrix(SEXP matrix); 										// Populate a data matrix object with the values of an R matrix
	void fillFromMxMatrix(SEXP matrix) { fillFromMatrix(matrix); };							// Populate a data matrix from an Mx Matrix object
	void fillFromS3Matrix(SEXP mxMatrix); 									// Populate a data matrix to represent the $values field of an MxMatrix object
	void operator=(omxMatrix &dM);											// Copy across another element.  NOTE: Duplicates.
	
/* Will eventually be needed for evaluation optimization. */
	void recompute();
	void compute();
	bool needsUpdate();
	
	void reset();										 		// Reset to the original majority and realias, if needed.
	void inline transpose() { colMajor = !colMajor; };			// Transpose = change row/col majority.  Alters the state of the object.

	void removeRowsAndColumns(int numRowsRemoved, int numColsRemoved, int rowsRemoved[], int colsRemoved[]);
}; 

#endif /* _OMXDATAMATRIX_H_ */


