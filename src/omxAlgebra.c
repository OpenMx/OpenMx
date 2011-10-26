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

omxMatrix* omxInitAlgebra(omxAlgebra *oa, omxState* os) {

	omxMatrix* om = omxInitMatrix(NULL, 0, 0, TRUE, os);
	
	omxInitAlgebraWithMatrix(oa, om);

	return om;
}

void omxInitAlgebraWithMatrix(omxAlgebra *oa, omxMatrix *om) {
	
	if(oa == NULL) {
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	}
	
	if(OMX_DEBUG) { Rprintf("Initializing algebra 0x%0x with 0x%0x.\n", oa, om); }
	
	oa->args = NULL;
	oa->funWrapper = NULL;
	oa->numArgs = 0;
	oa->matrix = om;
	om->algebra = oa;
	oa->name = NULL;
	
}

void omxDuplicateAlgebraMatrix(omxMatrix* tgt, omxMatrix* src, omxState* newState, short fullCopy) {

    // Assumes the matrices themselves have already been duplicated.
    if(src->algebra == NULL) {
        tgt->algebra = NULL;
        return;
    }

    omxAlgebra *tgtAlg = tgt->algebra;

    if(tgtAlg == NULL) {
        omxInitAlgebraWithMatrix(NULL, tgt);
    }
}



omxAlgebra* omxDuplicateAlgebra(omxAlgebra* tgt, omxAlgebra* src, omxState* newState, short fullCopy) {

	if(src == NULL || tgt == NULL) { return NULL; }

    tgt->funWrapper = src->funWrapper;  // N.B.: This should work across processes and threads, but might have problems with more complicated structures.
    tgt->numArgs = src->numArgs;
	tgt->args = (omxMatrix**) R_alloc(tgt->numArgs, sizeof(omxMatrix*));
    for(int k = 0; k < tgt->numArgs; k++) {
	    tgt->args[k] = omxLookupDuplicateElement(newState, src->args[k]);
    }

	tgt->name = src->name;  // Name strings are constant.

	return tgt;
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
    	        j, oa->args[j]);
    	}
		omxFreeAllMatrixData(oa->args[j]);
		oa->args[j] = NULL;
	}
	oa->numArgs = 0;
	oa->matrix = NULL;
}

void omxAlgebraCompute(omxAlgebra *oa) {
	if(OMX_DEBUG_ALGEBRA) { Rprintf("Algebra compute (%s): 0x%0x (needed: %d/%d).\n", oa->name, oa->matrix, oa->matrix->lastCompute, oa->matrix->currentState->computeCount); }
		
	for(int j = 0; j < oa->numArgs; j++) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Recomputing arg %d at 0x%0x (Which %s need it).\n", j, oa->args[j], (omxNeedsUpdate(oa->args[j])?"does":"does not")); }
		omxRecompute(oa->args[j]);
	}
   // Recompute happens in handleFreeVars, for now.
	
	if(oa->funWrapper == NULL) { 			// No-op algebra: only for algebra-is-a-matrix condition.
		if(oa->numArgs == 1) {
			if(OMX_DEBUG_ALGEBRA) { omxPrint(oa->args[0], "No-op Matrix"); }
			omxCopyMatrix(oa->matrix, oa->args[0], TRUE);
		} else {
			error("Internal Error: Empty algebra evaluated.\n");
		}
	} else {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("Activating function with %d args.\n", oa->numArgs); }
		(*((void(*)(omxMatrix**, int, omxMatrix*))oa->funWrapper))(oa->args, (oa->numArgs), oa->matrix);
	}
	omxMatrixCompute(oa->matrix);
	
	if(OMX_DEBUG_ALGEBRA) { omxAlgebraPrint(oa, "Result is"); }
}

int omxAlgebraNeedsUpdate(omxAlgebra *oa)
{
	if(OMX_DEBUG_ALGEBRA) {Rprintf("AlgebraNeedsUpdate:%d?", oa->numArgs);}
	for(int j = 0; j < fabs(oa->numArgs); j++) {
		if(omxNeedsUpdate(oa->args[j])) {
			if(OMX_DEBUG_ALGEBRA) {Rprintf("Arg %d Needs Update.\n");}
			return TRUE;
		}
	}
    return FALSE;
}

omxMatrix* omxNewMatrixFromMxAlgebra(SEXP alg, omxState* os, const char *name) {

	omxMatrix *om = omxInitMatrix(NULL, 0, 0, TRUE, os);

	om->hasMatrixNumber = 0;
	om->matrixNumber = 0;	

	omxFillMatrixFromMxAlgebra(om, alg, name);
	
	return om;
}

void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP algebra, const char *name) {

	int value;
	omxAlgebra *oa = NULL;
	SEXP algebraOperator, algebraArg, algebraElt;
	
	PROTECT(algebraOperator = AS_INTEGER(VECTOR_ELT(algebra, 0)));
	value = INTEGER(algebraOperator)[0];

	if(OMX_DEBUG) {Rprintf("Creating Algebra from Sexp.\n");}

	if(value > 0) { 			// This is an operator.
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
		omxInitAlgebraWithMatrix(oa, om);
		if(OMX_DEBUG) {Rprintf("Retrieving Table Entry %d.\n", value);}
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[value]);
		if(OMX_DEBUG) {Rprintf("Table Entry %d (at 0x%0x) is %s.\n", value, entry, entry->opName);}
		omxFillAlgebraFromTableEntry(oa, entry);
		if(oa->numArgs < 0) {	// Special Case: open-ended operator.  Might want to move this section to omxFillAlgebraFromTableEntry
			oa->numArgs = length(algebra) - 1;  // Has as many arguments as there are elements after the operator
			oa->args = (omxMatrix**)R_alloc(oa->numArgs, sizeof(omxMatrix*));
		}
		for(int j = 0; j < oa->numArgs; j++) {
			PROTECT(algebraArg = VECTOR_ELT(algebra, j+1));
				oa->args[j] = omxAlgebraParseHelper(algebraArg, om->currentState, name);
				if(OMX_DEBUG) { Rprintf("fillFromMxAlgebra got 0x%0x from helper, arg %d.\n", oa->args[j-1], j); }
			UNPROTECT(1); /* algebraArg */
		}
	} else {		// This is an algebra pointer, and we're a No-op algebra.
		/* TODO: Optimize this by eliminating no-op algebras entirely. */
		PROTECT(algebraElt = VECTOR_ELT(algebra, 1));
		
		if(!IS_INTEGER(algebraElt)) {   			// A List: only happens if bad optimization has occurred.
			warning("Internal Error: Algebra has been passed incorrectly: detected NoOp: (Operator Arg ...)\n");
			omxFillMatrixFromMxAlgebra(om, algebraElt, name);		// Collapse the no-op algebra
		} else {			// Still a No-op.  Sadly, we have to keep it that way.
			
			PROTECT(algebraOperator = AS_INTEGER(algebraElt));
			value = INTEGER(algebraOperator)[0];
			
			oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
			omxInitAlgebraWithMatrix(oa, om);
			oa->args = (omxMatrix**) R_alloc(1, sizeof(omxMatrix*));
			
			if(value < 0) {
				value = ~value;					// Bitwise reverse of number--this is a matrix index
				oa->args[0] = (oa->matrix->currentState->matrixList[value]);
			} else {
				oa->args[0] = (oa->matrix->currentState->algebraList[value]);
			}
			oa->numArgs = 1;
			UNPROTECT(1); /* algebraArg */
		}
		UNPROTECT(1); /* algebraElt */
	}
	oa->name = name;

	UNPROTECT(1);	/* algebraOperator */

	// omxAlgebraCompute(oa);
	// 
	// omxMatrixCompute(om);

}

void omxFillAlgebraFromTableEntry(omxAlgebra *oa, const omxAlgebraTableEntry* oate) {
	/* TODO: check for full initialization */
	if(oa == NULL) error("Internal Error: Null Algebra Detected in fillAlgebra.");
	
	if(OMX_DEBUG) { Rprintf("Filling from table entry %d (%s)....", oate->number, oate->rName); }
	oa->funWrapper = oate->funWrapper;
	oa->numArgs = oate->numArgs;
	if(oa->numArgs >= 0) {
		oa->args = (omxMatrix**)R_alloc(oa->numArgs, sizeof(omxMatrix*));
	} else {
		oa->args = NULL;
	}
	if(OMX_DEBUG) { Rprintf("Table Entry processed.\n"); }
}

omxMatrix* omxAlgebraParseHelper(SEXP algebraArg, omxState* os, const char *name) {
	omxMatrix* newMat;
	if(OMX_DEBUG) { Rprintf("Helper: processing next arg..."); }
	
	if(!IS_INTEGER(algebraArg)) {
		if(OMX_DEBUG) { Rprintf("Helper detected list element.  Recursing.\n"); }
		newMat = omxNewMatrixFromMxAlgebra(algebraArg, os, name);
	} else {
		newMat = omxNewMatrixFromMxIndex(algebraArg, os);
	}
	
	return(newMat);
}

void omxAlgebraPrint(omxAlgebra* oa, char* d) {
	Rprintf("(Algebra) ");
	omxPrintMatrix(oa->matrix, d);
	Rprintf("has %d args.\n", oa->numArgs);
}

omxMatrix* omxNewMatrixFromMxIndex(SEXP matrix, omxState* os) {
	if(OMX_DEBUG) { Rprintf("Attaching pointer to matrix."); }

	int value = 0;
	omxMatrix* output = NULL;

	if (IS_INTEGER(matrix)) {
		SEXP intMatrix;
		PROTECT(intMatrix = AS_INTEGER(matrix));
		value = INTEGER(intMatrix)[0];
		if(value == NA_INTEGER) {
	    	if(OMX_DEBUG) {Rprintf("  Null integer matrix.  Skipping.\n");}
			UNPROTECT(1);
			return NULL;
		}
		UNPROTECT(1); // intMatrix
	} else if (IS_NUMERIC(matrix)) {
		SEXP numericMatrix;
		PROTECT(numericMatrix = AS_NUMERIC(matrix));
		value = (int) REAL(numericMatrix)[0];
		if(value == NA_INTEGER) {
	    	if(OMX_DEBUG) {Rprintf("   Null numeric matrix.  Skipping.\n");}
			UNPROTECT(1);
			return NULL;
		}
		UNPROTECT(1); // numericMatrix		
	} else {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Internal error: unknown type passed to omxNewMatrixFromMxIndex.");
		omxRaiseError(os, -1, errstr);
		free(errstr);
		return NULL;
	}		
	if(OMX_DEBUG) {Rprintf("  Pointer is %d.\n", value);}
	if (value >= 0) {										// Pre-existing algebra.  A-ok.
		output = *(os->algebraList + value);
	} else {												// Pre-existing matrix.  A-ok.
		output = os->matrixList[~value];						// Value invert for matrices.
	}
	
	return output;
}

omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix* args[], int numArgs, omxState* os) {
	
	if(OMX_DEBUG) {Rprintf("Generating new algebra from opcode %d (%s).\n", opCode, omxAlgebraSymbolTable[opCode].rName);}
	omxMatrix *om;
	omxAlgebra *oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	omxAlgebraTableEntry* entry = (omxAlgebraTableEntry*)&(omxAlgebraSymbolTable[opCode]);
	if(entry->numArgs >= 0 && entry->numArgs != numArgs) {
		char *errstr = calloc(250, sizeof(char));
		sprintf(errstr, "Internal error: incorrect number of arguments passed to algebra %s.", entry->rName);
		omxRaiseError(os, -1, errstr);
		free(errstr);
		return NULL;
	}
	
	om = omxInitAlgebra(oa, os);
	omxFillAlgebraFromTableEntry(oa, entry);
	oa->name = entry->opName;

	if(OMX_DEBUG) {Rprintf("Calculating args for %s.\n", entry->rName);}
	if(oa->args == NULL) {					// # of matrices for operator is variable
		oa->numArgs = numArgs;
		oa->args = (omxMatrix**) R_alloc(numArgs, sizeof(omxMatrix*));
	} 
	
	if(OMX_DEBUG) {Rprintf("Populating args for %s.\n", entry->rName);}
	
	for(int i = 0; i < numArgs;i++) {
		oa->args[i] = args[i];
	}
	
	omxMarkDirty(om);
	
	return om;
	
}
