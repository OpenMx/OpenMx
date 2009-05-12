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
 */

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
	
	omxRObjective *newObj = (omxRObjective*) R_alloc(1, sizeof(omxRObjective));
	PROTECT(newObj->objfun = GET_SLOT(rObj, install("objective")));
	PROTECT(newObj->env = GET_SLOT(rObj, install("env")));
	
	oo->argStruct = (void*) newObj;
}


#endif /* _OMX_R_OBJECTIVE_ */
