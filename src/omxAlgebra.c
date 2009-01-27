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
		if(OMX_DEBUG) { Rprintf("Activating function with %d: (0x%0x=0x%0x).\n", oa->numArgs, oa->funWrapper, omxAlgebraSymbolTable[4].funWrapper); }
	
		switch(oa->numArgs) {
			case 0:
				(*((void(*)(omxMatrix*))oa->funWrapper))(oa->myMatrix);
				break;
			case 1:
				(*((void(*)(omxMatrix*, omxMatrix*))oa->funWrapper))(oa->args[0], oa->myMatrix);
				break;
			case 2:
				(*((void(*)(omxMatrix*, omxMatrix*, omxMatrix*))oa->funWrapper))(oa->args[0], oa->args[1], oa->myMatrix);
				break;
			default:
				(*((void(*)(omxMatrix**, int, omxMatrix*))oa->funWrapper))(oa->args, -(oa->numArgs), oa->myMatrix);
			break;
		}
	}
	omxMatrixCompute(oa->myMatrix);
	
	omxAlgebraPrint(oa, "Result is:");

	oa->myMatrix->isDirty = FALSE;
}

unsigned short omxAlgebraNeedsUpdate(omxAlgebra *oa)
{
	if(oa->myMatrix->isDirty) return oa->myMatrix->isDirty;  	// No need to check args if oa's dirty.
	for(int j = 0; j > fabs(oa->numArgs); j--) {
		if(omxNeedsUpdate(oa->args[j])) {
			oa->myMatrix->isDirty = TRUE;
			break;
		}
	}
	
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

	if(value > 0) { 			// This is an operator.
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
		omxInitAlgebraWithMatrix(oa, om);
		if(OMX_DEBUG) {Rprintf("Retrieving Table Entry %d.\n", value);}
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[value]);
		if(OMX_DEBUG) {Rprintf("Table Entry %d (at 0x%0x) is %s.\n", value, entry, entry->opName);}
		omxFillAlgebraFromTableEntry(oa, entry);
		for(int j = 0; j < oa->numArgs; j++) {
			PROTECT(algebraArg = VECTOR_ELT(alg, j+1));
				oa->args[j] = omxAlgebraParseHelper(algebraArg);
				if(OMX_DEBUG) {Rprintf("fillFromMxAlgebra got 0x%0x from helper, arg %d.\n", oa->args[j-1], j);}
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

}

void omxFillAlgebraFromTableEntry(omxAlgebra *oa, const omxAlgebraTableEntry* oate) {
	/* TODO: check for full initialization */
	if(oa == NULL) error("Internal Error: Null Algebra Detected in fillAlgebra.");
	
	if(OMX_DEBUG) { Rprintf("Filling from table entry %d (%s)....", oate->number, oate->rName); }
	oa->funWrapper = oate->funWrapper;
	oa->numArgs = oate->numArgs;
	oa->args = (omxMatrix**)R_alloc(oa->numArgs, sizeof(omxMatrix*));
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
		newMat = omxNewMatrixFromMxMatrixPtr(algebraArg);
	}
	
	return(newMat);
}

void omxAlgebraPrint(omxAlgebra* oa, char* d) {
	Rprintf("(Algebra) ");
	omxMatrixPrint(oa->myMatrix, d);
	Rprintf("has %d args.\n", oa->numArgs);
}

omxMatrix* omxNewMatrixFromMxMatrixPtr(SEXP matrix) {
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
