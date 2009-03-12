#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_ALGEBRA_OBJECTIVE_
#define _OMX_ALGEBRA_OBJECTIVE_ TRUE

typedef struct {

	omxMatrix *algebra;

} omxAlgebraObjective;

void omxDestroyAlgebraObjective(omxObjective *oo) {

}

void omxCallAlgebraObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	omxRecomputeMatrix(((omxAlgebraObjective*)(oo->argStruct))->algebra);
	oo->myMatrix->data[0] = ((omxAlgebraObjective*)(oo->argStruct))->algebra->data[0];
	
}

unsigned short int omxNeedsUpdateAlgebraObjective(omxObjective *oo) {

	if(oo->myMatrix->data[0] != ((omxAlgebraObjective*)oo->argStruct)->algebra->data[0]) return TRUE;
	return omxNeedsUpdate(((omxAlgebraObjective*)oo->argStruct)->algebra);
}

void omxInitAlgebraObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	SEXP newptr;
	
	omxAlgebraObjective *newObj = (omxAlgebraObjective*) R_alloc(sizeof(omxAlgebraObjective), 1);
	PROTECT(newptr = GET_SLOT(rObj, install("algebra")));
	newObj->algebra = omxNewMatrixFromMxMatrixPtr(newptr);
	if(OMX_DEBUG) {Rprintf("Algebra Objective Bound to Algebra %d", newObj->algebra);}
	UNPROTECT(1);
	
	oo->needsUpdateFun = omxNeedsUpdateAlgebraObjective;
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_R_OBJECTIVE_ */