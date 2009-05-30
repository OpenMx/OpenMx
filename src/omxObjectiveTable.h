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

#ifndef _OMX_OBJECTIVE_TABLE_
#define _OMX_OBJECTIVE_TABLE_ TRUE

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

typedef struct omxObjectiveTableEntry omxObjectiveTableEntry;

#include "omxObjective.h"

struct omxObjectiveTableEntry {

	char name[250];
	void (*initFun)(omxObjective*, SEXP, SEXP) ;

};

#define omxObjectiveTableLength 5

const omxObjectiveTableEntry omxObjectiveSymbolTable[omxObjectiveTableLength];


#endif /* _OMX_OBJECTIVE_TABLE_ */
