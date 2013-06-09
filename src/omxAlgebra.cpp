/*
 *  Copyright 2007-2013 The OpenMx Project
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
*  omxAlgebra.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Algebras are a subclass of data matrix that evaluates
*   itself anew at each iteration, so that any changes to
*   free parameters can be incorporated into the update.
*
**********************************************************/

#include "omxMatrix.h"
#include "omxFitFunction.h"

static void
omxAlgebraAllocArgs(omxAlgebra *oa, int numArgs)
{
	if (numArgs <= 0) {
		oa->numArgs = 0;
		oa->algArgs = NULL;
		return;
	}

	if(oa->algArgs != NULL) {
		if (oa->numArgs < numArgs)
			error("omxAlgebra %p: %d args requested but %d available",
			      oa, numArgs, oa->numArgs);
		return;
	}

	oa->numArgs = numArgs;
	oa->algArgs = (omxMatrix**) R_alloc(numArgs, sizeof(omxMatrix*));
	memset(oa->algArgs, 0, sizeof(omxMatrix*) * numArgs);  //remove debug TODO
}

omxMatrix* omxInitAlgebra(omxAlgebra *oa, omxState* os) {

	omxMatrix* om = omxInitMatrix(NULL, 0, 0, TRUE, os);
	
	omxInitAlgebraWithMatrix(oa, om);

	return om;
}

void omxInitAlgebraWithMatrix(omxAlgebra *oa, omxMatrix *om) {
	
	if(oa == NULL) {
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	}
	
	if(OMX_DEBUG && om->currentState->parentState == NULL) { 
		Rprintf("Initializing algebra 0x%0x with 0x%0x.\n", oa, om);
	}
	
	omxAlgebraAllocArgs(oa, 0);
	oa->funWrapper = NULL;
	oa->matrix = om;
	om->algebra = oa;

}

void omxDuplicateAlgebra(omxMatrix* tgt, omxMatrix* src, omxState* newState) {

    if(src->algebra != NULL) {
		omxFillMatrixFromMxAlgebra(tgt, src->algebra->sexpAlgebra, src->name);
    } else if(src->fitFunction != NULL) {
        omxDuplicateFitMatrix(tgt, src, newState);
    }

}

void omxFreeAlgebraArgs(omxAlgebra *oa) {
	/* Completely destroy the algebra tree */
	
	if(OMX_DEBUG) { 
	    Rprintf("Freeing algebra at 0x%0x with %d args.\n", 
	        oa, oa->numArgs); 
	}
	
	int j;
	for(j = 0; j < oa->numArgs; j++) {
		if(OMX_DEBUG) {
			Rprintf("Freeing argument %d at 0x%0x.\n",
				j, oa->algArgs[j]);
		}
		omxFreeAllMatrixData(oa->algArgs[j]);
		oa->algArgs[j] = NULL;
	}
	omxAlgebraAllocArgs(oa, 0);
	oa->matrix = NULL;
}

void omxAlgebraRecompute(omxAlgebra *oa) {
	if(OMX_DEBUG_ALGEBRA) { 
		Rprintf("Algebra compute (%s): 0x%0x.\n", 
			oa->matrix->name, oa->matrix);
	}
		
	for(int j = 0; j < oa->numArgs; j++) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Recomputing arg %d at 0x%0x (Which %s need it).\n", j, oa->algArgs[j], (omxNeedsUpdate(oa->algArgs[j])?"does":"does not")); }
		omxRecompute(oa->algArgs[j]);
	}
   // Recompute happens in handleFreeVars, for now.
	
	if(oa->funWrapper == NULL) { 			// No-op algebra: only for algebra-is-a-matrix condition.
		if(oa->numArgs == 1) {
			if(OMX_DEBUG_ALGEBRA) { omxPrint(oa->algArgs[0], "No-op Matrix"); }
			omxCopyMatrix(oa->matrix, oa->algArgs[0]);
		} else {
			error("Internal Error: Empty algebra evaluated.\n");
		}
	} else {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Activating function with %d args.\n", oa->numArgs); }
		(*((void(*)(omxMatrix**, int, omxMatrix*))oa->funWrapper))(oa->algArgs, (oa->numArgs), oa->matrix);
	}

	omxMarkClean(oa->matrix);
	
	if(OMX_DEBUG_ALGEBRA) { omxAlgebraPrint(oa, "Result is"); }
}



void omxAlgebraCompute(omxAlgebra *oa) {
	if(OMX_DEBUG_ALGEBRA) { 
		Rprintf("Algebra compute (%s): 0x%0x.\n", 
			oa->matrix->name, oa->matrix);
	}
		
	for(int j = 0; j < oa->numArgs; j++) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Recomputing arg %d at 0x%0x (Which %s need it).\n", j, oa->algArgs[j], (omxNeedsUpdate(oa->algArgs[j])?"does":"does not")); }
		omxCompute(oa->algArgs[j]);
	}
   // Recompute happens in handleFreeVars, for now.
	
	if(oa->funWrapper == NULL) { 			// No-op algebra: only for algebra-is-a-matrix condition.
		if(oa->numArgs == 1) {
			if(OMX_DEBUG_ALGEBRA) { omxPrint(oa->algArgs[0], "No-op Matrix"); }
			omxCopyMatrix(oa->matrix, oa->algArgs[0]);
		} else {
			error("Internal Error: Empty algebra evaluated.\n");
		}
	} else {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Activating function with %d args.\n", oa->numArgs); }
		(*((void(*)(omxMatrix**, int, omxMatrix*))oa->funWrapper))(oa->algArgs, (oa->numArgs), oa->matrix);
	}

	omxMarkClean(oa->matrix);
	
	if(OMX_DEBUG_ALGEBRA) { omxAlgebraPrint(oa, "Result is"); }
}

int omxAlgebraNeedsUpdate(omxAlgebra *oa)
{
	return TRUE;
}

omxMatrix* omxNewMatrixFromMxAlgebra(SEXP alg, omxState* os, const char *name) {

	omxMatrix *om = omxInitMatrix(NULL, 0, 0, TRUE, os);

	om->hasMatrixNumber = 0;
	om->matrixNumber = 0;	

	omxFillMatrixFromMxAlgebra(om, alg, name);
	
	return om;
}

static void
omxFillAlgebraFromTableEntry(omxAlgebra *oa, const omxAlgebraTableEntry* oate, const int realNumArgs) {
	/* TODO: check for full initialization */
	if(oa == NULL) error("Internal Error: Null Algebra Detected in fillAlgebra.");

	if(OMX_DEBUG && oa->matrix->currentState->parentState == NULL) {
		Rprintf("Filling from table entry %d (%s)....", oate->number, oate->rName);
	}
	oa->funWrapper = oate->funWrapper;
	omxAlgebraAllocArgs(oa, oate->numArgs==-1? realNumArgs : oate->numArgs);
	if(OMX_DEBUG && oa->matrix->currentState->parentState == NULL) {
		Rprintf("Table Entry processed.\n");
	}
}

void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP algebra, const char *name) {

	int value;
	omxAlgebra *oa = NULL;
	SEXP algebraOperator, algebraArg, algebraElt;
	
	PROTECT(algebraOperator = AS_INTEGER(VECTOR_ELT(algebra, 0)));
	value = INTEGER(algebraOperator)[0];

	if(OMX_DEBUG && om->currentState->parentState == NULL) {Rprintf("Creating Algebra from Sexp.\n");}

	if(value > 0) { 			// This is an operator.
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
		omxInitAlgebraWithMatrix(oa, om);
		if(OMX_DEBUG && om->currentState->parentState == NULL) {Rprintf("Retrieving Table Entry %d.\n", value);}
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[value]);
		if(OMX_DEBUG && om->currentState->parentState == NULL) {Rprintf("Table Entry %d (at 0x%0x) is %s.\n", value, entry, entry->opName);}
		omxFillAlgebraFromTableEntry(oa, entry, length(algebra) - 1);
		for(int j = 0; j < oa->numArgs; j++) {
			PROTECT(algebraArg = VECTOR_ELT(algebra, j+1));
				oa->algArgs[j] = omxAlgebraParseHelper(algebraArg, om->currentState, NULL);
				if(OMX_DEBUG && om->currentState->parentState == NULL) {
					Rprintf("fillFromMxAlgebra got 0x%0x from helper, arg %d.\n", oa->algArgs[j], j);
				}
		}
	} else {		// This is an algebra pointer, and we're a No-op algebra.
		/* TODO: Optimize this by eliminating no-op algebras entirely. */
		PROTECT(algebraElt = VECTOR_ELT(algebra, 1));
		
		if(!IS_INTEGER(algebraElt)) {   			// A List: only happens if bad optimization has occurred.
			warning("Internal Error: Algebra has been passed incorrectly: detected NoOp: (Operator Arg ...)\n");
			omxFillMatrixFromMxAlgebra(om, algebraElt, NULL);		// Collapse the no-op algebra
		} else {			// Still a No-op.  Sadly, we have to keep it that way.
			
			PROTECT(algebraOperator = AS_INTEGER(algebraElt));
			value = INTEGER(algebraOperator)[0];
			
			oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
			omxInitAlgebraWithMatrix(oa, om);
			omxAlgebraAllocArgs(oa, 1);
			
			if(value < 0) {
				value = ~value;					// Bitwise reverse of number--this is a matrix index
				oa->algArgs[0] = (oa->matrix->currentState->matrixList[value]);
			} else {
				oa->algArgs[0] = (oa->matrix->currentState->algebraList[value]);
			}
		}
	}
        if (name) om->name        = name;
	oa->sexpAlgebra = algebra;
}

omxMatrix* omxAlgebraParseHelper(SEXP algebraArg, omxState* os, const char *name) {
	omxMatrix* newMat;
	if(OMX_DEBUG && os->parentState == NULL) { Rprintf("Helper: processing next arg..."); }
	
	if(!IS_INTEGER(algebraArg)) {
		if(OMX_DEBUG && os->parentState == NULL) { Rprintf("Helper detected list element.  Recursing.\n"); }
		newMat = omxNewMatrixFromMxAlgebra(algebraArg, os, name);
	} else {
		newMat = omxMatrixLookupFromState1(algebraArg, os);
	}
	
	return(newMat);
}

void omxAlgebraPrint(omxAlgebra* oa, const char* d) {
	Rprintf("(Algebra) ");
	omxPrintMatrix(oa->matrix, d);
	Rprintf("has %d args.\n", oa->numArgs);
}

static int isValidMatrixIndex(SEXP r_index)
{
	if (r_index == R_NilValue) {
		return 0;
	} else if (isString(r_index)) {
		error("String passed to omxNewMatrixFromIndex, did you forget to call imxLocateIndex?");
	}		
	return 1;
}

omxMatrix* omxMatrixLookupFromState1(SEXP matrix, omxState* os)
{
	if (!isValidMatrixIndex(matrix)) return NULL;

	// length(matrix) can be >1, but we only decode the first element

	int value = 0;
	if (IS_INTEGER(matrix)) {
		value = INTEGER(matrix)[0];
		if (value == NA_INTEGER) return NULL;
	} else if (IS_NUMERIC(matrix)) {
		value = REAL(matrix)[0];
	} else {
		error("Internal error: unknown type passed to omxMatrixLookupFromState1");
	}		
	if(OMX_DEBUG && os->parentState == NULL) {
		Rprintf("  Pointer is %d.\n", value);
	}

	omxMatrix *output = NULL;
	if (value >= 0) {
		output = os->algebraList[value];
	} else {
		output = os->matrixList[~value];		// Value invert for matrices
	}

	// omxRecompute(output); can do this here?
	
	return output;
}

void omxMatricesLookupFromState(SEXP r_index, omxState *os, int *len_out, omxMatrix **out[])
{
	if (!isValidMatrixIndex(r_index)) return;

	if (IS_INTEGER(r_index)) {
		// OK
	} else if (IS_NUMERIC(r_index)) {
		error("Matrix/algebra numbers must be given as integers (not numeric type)");
	} else {
		error("Unknown type");
	}

	*len_out = length(r_index);
	int *indices = INTEGER(r_index);

	*out = Realloc(NULL, *len_out, omxMatrix*);
	for (int mx=0; mx < *len_out; mx++) {
		int value = indices[mx];
		if (value >= 0) {
			(*out)[mx] = os->algebraList[value];
		} else {
			(*out)[mx] = os->matrixList[~value];	// Value invert for matrices
		}
		omxRecompute((*out)[mx]);
	}
}

omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix* args[], int numArgs, omxState* os) {
	
	if(OMX_DEBUG) {Rprintf("Generating new algebra from opcode %d (%s).\n", opCode, omxAlgebraSymbolTable[opCode].rName);}
	omxMatrix *om;
	omxAlgebra *oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	omxAlgebraTableEntry* entry = (omxAlgebraTableEntry*)&(omxAlgebraSymbolTable[opCode]);
	if(entry->numArgs >= 0 && entry->numArgs != numArgs) {
		error("Internal error: incorrect number of arguments passed to algebra %s.", entry->rName);
	}
	
	om = omxInitAlgebra(oa, os);
	omxFillAlgebraFromTableEntry(oa, entry, entry->numArgs);
	om->name = entry->opName;

	if(OMX_DEBUG) {Rprintf("Calculating args for %s.\n", entry->rName);}
	omxAlgebraAllocArgs(oa, numArgs);
	
	if(OMX_DEBUG) {Rprintf("Populating args for %s.\n", entry->rName);}
	
	for(int i = 0; i < numArgs;i++) {
		oa->algArgs[i] = args[i];
	}
	
	omxMarkDirty(om);
	
	return om;
	
}

