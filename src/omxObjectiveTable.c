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

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxSymbolTable.h"
#include "omxAlgebraFunctions.h"

extern void omxInitAlgebraObjective(omxObjective *oo, SEXP rObj);
extern void omxInitFIMLObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRAMObjective(omxObjective *oo, SEXP rObj);
extern void omxInitLISRELObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRowObjective(omxObjective *oo, SEXP rObj);
extern void omxInitMLObjective(omxObjective *oo, SEXP rObj);
extern void omxInitRObjective(omxObjective *oo, SEXP rObj);	
extern void omxInitWLSObjective(omxObjective *oo, SEXP rObj);	
const omxObjectiveTableEntry omxObjectiveSymbolTable[omxObjectiveTableLength] = {
	{"MxAlgebraObjective", 			&omxInitAlgebraObjective},
	{"MxFIMLObjective",				&omxInitFIMLObjective},
	{"MxRAMObjective", 				&omxInitRAMObjective},
	{"MxWLSObjective",				&omxInitWLSObjective},
	{"MxRowObjective", 				&omxInitRowObjective},
	{"MxMLObjective", 				&omxInitMLObjective},
	{"MxRObjective",				&omxInitRObjective},
	{"MxLISRELObjective",			&omxInitLISRELObjective},
};
