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

/* Need a better way to deal with these. */
extern omxMatrix** algebraList;
extern omxMatrix** matrixList;

omxMatrix* omxInitAlgebra(omxAlgebra *oa) {

	omxMatrix* om = omxInitMatrix(NULL, 0, 0, TRUE);
	
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
	oa->myMatrix = om;
	om->algebra = oa;
	
}

void omxFreeAlgebraArgs(omxAlgebra *oa) {
	/* Completely destroy the algebra tree */
	
	int j;
	for(j = 0; j < oa->numArgs; j++) {
		omxFreeAllMatrixData(oa->args[j]);
		oa->args[j] = NULL;
	}
	oa->numArgs = 0;
	oa->myMatrix = NULL;
}

void omxAlgebraCompute(omxAlgebra *oa) {
	if(OMX_DEBUG) {Rprintf("Algebra compute: 0x%0x (needed: %d).\n", oa, oa->myMatrix->isDirty);}
		
	for(int j = 0; j < oa->numArgs; j++) {
		if(OMX_DEBUG) { Rprintf("Recomputing arg %d at 0x%0x.\n", j, oa->args[j]); }
		omxRecomputeMatrix(oa->args[j]);
	}
   // Recompute happens in handleFreeVars, for now.
	
	if(oa->funWrapper == NULL) { 			// No-op algebra: only for algebra-is-a-matrix condition.
		if(oa->numArgs == 1) {
			if(OMX_DEBUG) { omxPrintMatrix(oa->args[0], "No-op Matrix"); }
			omxCopyMatrix(oa->myMatrix, oa->args[0]);
		} else {
			error("Internal Error: Empty algebra evaluated.\n");
		}
	} else {
		if(OMX_DEBUG) { Rprintf("Activating function with %d args.\n", oa->numArgs); }
		(*((void(*)(omxMatrix**, int, omxMatrix*))oa->funWrapper))(oa->args, (oa->numArgs), oa->myMatrix);
	}
	omxComputeMatrixHelper(oa->myMatrix);
	
	if(OMX_DEBUG) { omxAlgebraPrint(oa, "Result is"); }

	oa->myMatrix->isDirty = FALSE;
}

unsigned short omxAlgebraNeedsUpdate(omxAlgebra *oa)
{
/* Permanently dirty quick fix */
    return TRUE;  

	if(OMX_DEBUG) {Rprintf("AlgebraNeedsUpdate:%d?", oa->numArgs);}
	if(oa->myMatrix->isDirty) return TRUE;  	// No need to check args if oa's dirty.
	for(int j = 0; j < fabs(oa->numArgs); j++) {
		if(omxNeedsUpdate(oa->args[j])) {
			if(OMX_DEBUG) {Rprintf("Arg Needs Update.");}
			oa->myMatrix->isDirty = TRUE;
			break;
		}
	}

	if(OMX_DEBUG) {Rprintf("Arg:%d, me:%d", oa->args[0]->isDirty, oa->myMatrix->isDirty);}
	
	return oa->myMatrix->isDirty;
}

omxMatrix* omxNewMatrixFromMxAlgebra(SEXP alg) {

	omxMatrix *om = omxInitMatrix(NULL, 0, 0, TRUE);
	
	omxFillMatrixFromMxAlgebra(om, alg);
	
	return om;
}

void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP alg) {

	int value;
	omxAlgebra *oa;
	SEXP algebraOperator, algebraArg, algebraElt;
	PROTECT(algebraOperator = AS_INTEGER(VECTOR_ELT(alg, 0)));
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
			oa->numArgs = length(alg) - 1;  // Has as many arguments as there are elements after the operator
			oa->args = (omxMatrix**)R_alloc(oa->numArgs, sizeof(omxMatrix*));
		}
		for(int j = 0; j < oa->numArgs; j++) {
			PROTECT(algebraArg = VECTOR_ELT(alg, j+1));
				oa->args[j] = omxAlgebraParseHelper(algebraArg);
				if(OMX_DEBUG) { Rprintf("fillFromMxAlgebra got 0x%0x from helper, arg %d.\n", oa->args[j-1], j); }
			UNPROTECT(1); /* algebraArg */
		}
	} else if(value == 0) {		// This is an algebra pointer, and we're a No-op algebra.
		/* TODO: Optimize this by eliminating no-op algebras entirely. */
		PROTECT(algebraElt = VECTOR_ELT(alg, 1));
		
		if(!IS_NUMERIC(algebraElt)) {   			// A List: only happens if bad optimization has occurred.
			warning("Internal Error: Algebra has been passed incorrectly: detected NoOp: (Operator Arg ...)\n");
			omxFillMatrixFromMxAlgebra(om, algebraElt);		// Collapse the no-op algebra
		} else {			// Still a No-op.  Sadly, we have to keep it that way.
			
			PROTECT(algebraOperator = AS_INTEGER(algebraElt));
			value = INTEGER(algebraOperator)[0];
			
			oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
			omxInitAlgebraWithMatrix(oa, om);
			oa->args = (omxMatrix**) R_alloc(1, sizeof(omxMatrix*));
			
			if(value < 0) {
				value = ~value;					// Bitwise reverse of number--this is a matrix index
				oa->args[0] = (matrixList[value]);
			} else {
				oa->args[0] = (algebraList[value]);
			}
			oa->numArgs = 1;
			UNPROTECT(1); /* algebraArg */
		}
		UNPROTECT(1); /* algebraElt */
	}

	UNPROTECT(1);	/* algebraOperator */

	omxAlgebraCompute(oa);
	
	omxComputeMatrixHelper(om);

}

void omxFillAlgebraFromTableEntry(omxAlgebra *oa, const omxAlgebraTableEntry* oate) {
	/* TODO: check for full initialization */
	if(oa == NULL) error("Internal Error: Null Algebra Detected in fillAlgebra.");
	
	if(OMX_DEBUG) { Rprintf("Filling from table entry %d (%s)....", oate->number, oate->rName); }
	oa->funWrapper = oate->funWrapper;
	oa->numArgs = oate->numArgs;
	if(oa->numArgs >= 0) {
		oa->args = (omxMatrix**)R_alloc(oa->numArgs, sizeof(omxMatrix*));
	}
	if(OMX_DEBUG) { Rprintf("Table Entry processed.\n"); }
}

omxMatrix* omxAlgebraParseHelper(SEXP algebraArg) {
	int value;
	omxAlgebra* newAlg;
	omxMatrix* newMat;
	SEXP argInts;
	if(OMX_DEBUG) { Rprintf("Helper: processing next arg..."); }
	
	if(!IS_NUMERIC(algebraArg)) {
		if(OMX_DEBUG) { Rprintf("Helper detected list element.  Recursing.\n"); }
		newMat = omxNewMatrixFromMxAlgebra(algebraArg);
	} else {
		newMat = omxNewMatrixFromMxIndex(algebraArg);
	}
	
	return(newMat);
}

void omxAlgebraPrint(omxAlgebra* oa, char* d) {
	Rprintf("(Algebra) ");
	omxPrintMatrixHelper(oa->myMatrix, d);
	Rprintf("has %d args.\n", oa->numArgs);
}

omxMatrix* omxNewMatrixFromMxIndex(SEXP matrix) {
	if(OMX_DEBUG) { Rprintf("Attaching pointer to matrix."); }
	SEXP intMatrix;
	int value = 0;
	omxMatrix* output = NULL;
	
	PROTECT(intMatrix = AS_INTEGER(matrix));
	value = INTEGER(intMatrix)[0];
	
	if(OMX_DEBUG) {Rprintf("  Pointer is %d.\n", value);}
	if (value >= 0) {										// Pre-existing algebra.  A-ok.
		output = *(algebraList + value);
	} else {												// Pre-existing matrix.  A-ok.
		output = matrixList[~value];						// Value invert for matrices.
	}
	
	UNPROTECT(1); // intMatrix
	return output;
}

omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix* arg1, omxMatrix* arg2) {
	/* For now, we'll be content with 2 args. */
	
	omxMatrix *om;
	omxAlgebra *oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	omxAlgebraTableEntry* entry = (omxAlgebraTableEntry*)&(omxAlgebraSymbolTable[opCode]);
	
	om = omxInitAlgebra(oa);
	omxFillAlgebraFromTableEntry(oa, entry);
	oa->args[0] = arg1;
	oa->args[1] = arg2;
	
	omxMarkDirty(om);
	
	return om;
	
}
