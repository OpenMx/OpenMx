#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_R_OBJECTIVE_
#define _OMX_R_OBJECTIVE_ TRUE

typedef struct {

	SEXP objfun;
	SEXP env;

} omxRObjective;

void omxDestroyRObjective(omxObjective *oo) {
	
	UNPROTECT(2); 			// objfun and env
}

void omxCallRObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	SEXP theCall, theVars, theReturn;
	PROTECT(theCall = allocList(2));
//	PROTECT(theVars = allocVector(REALSXP, *n));
//	vars = REAL(theVars);
	SET_TYPEOF(theCall, LANGSXP);
	SETCAR(theCall, ((omxRObjective*)oo->argStruct)->objfun);

//	for(k = 0; k < *n; k++) {
//		vars[k] = x[k];
//	}
//	SETCADR(theCall, theVars);
	SETCADR(theCall, NULL);

	PROTECT(theReturn = eval(theCall, R_GlobalEnv));
	oo->myMatrix->data[0] = REAL(AS_NUMERIC(theReturn))[0];
	UNPROTECT(2);
}

void omxInitRObjective(omxObjective* oo, SEXP rObj, SEXP dataList) {
	
	omxRObjective *newObj = (omxRObjective*) R_alloc(sizeof(omxRObjective), 1);
	PROTECT(newObj->objfun = GET_SLOT(rObj, install("objective")));
	PROTECT(newObj->env = GET_SLOT(rObj, install("env")));
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_R_OBJECTIVE_ */