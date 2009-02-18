/***********************************************************
* 
*  omxObjective.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxObjective class
*   omxObjective objects hold necessary information to simplify
* 	objective function calculation.
*
**********************************************************/

#ifndef _OMXOBJECTIVE_H_
#define _OMXOBJECTIVE_H_

#include "R.h"
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

typedef struct omxObjective omxObjective;

#include "omxAlgebra.h"


#ifdef DEBUGMX
#define OMX_DEBUG 1
#else
#define OMX_DEBUG 0
#endif /* DEBUGMX */

struct omxObjective {					// An objective

	/* Fields unique to Objective Functions */
	void (*initFun)(omxObjective *oo, SEXP rObj, SEXP dataList);			// Wrapper for initialization function (probably not needed)
	void (*destructFun)(omxObjective* oo);									// Wrapper for the destructor object
	void (*objectiveFun)(omxObjective* oo);									// Wrapper for the objective function itself

	void* argStruct;														// Arguments to the above function
	char objType[250];														// Type of Objective Function

	omxMatrix* myMatrix;													// The (1x1) matrix populated by this objective function

};

/* Initialize and Destroy */
	void omxFillMatrixFromMxObjective(omxMatrix* om, SEXP mxobj, SEXP dataList); // Create an objective function from an R MxObjective object
	void omxFreeObjectiveArgs(omxObjective* objective);						// Frees all args

/* Algebra-specific implementations of matrix functions */
	void omxObjectiveRecompute(omxObjective *oo);
	void omxObjectiveCompute(omxObjective *oo);
	unsigned short int omxObjectiveNeedsUpdate(omxObjective *oo);

	void omxObjectivePrint(omxObjective *source, char* d);					// Pretty-print a (small) matrix

#endif /* _OMXOBJECTIVE_H_ */