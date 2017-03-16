/*
 *  Copyright 2007-2017 The OpenMx Project
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
#include "Compute.h"
#include "EnableWarnings.h"

void omxAlgebraAllocArgs(omxAlgebra *oa, int numArgs)
{
	if (numArgs <= 0) {
		oa->numArgs = 0;
		oa->algArgs = NULL;
		return;
	}

	if(oa->algArgs != NULL) {
		if (oa->numArgs < numArgs)
			Rf_error("omxAlgebra: %d args requested but %d available",
			      numArgs, oa->numArgs);
		return;
	}

	oa->numArgs = numArgs;
	oa->algArgs = (omxMatrix**) R_alloc(numArgs, sizeof(omxMatrix*));
	memset(oa->algArgs, 0, sizeof(omxMatrix*) * numArgs);  //remove debug TODO
}

omxMatrix* omxInitAlgebra(omxAlgebra *oa, omxState* os) {

	omxMatrix* om = omxInitMatrix(0, 0, TRUE, os);
	
	omxInitAlgebraWithMatrix(oa, om);

	return om;
}

void omxInitAlgebraWithMatrix(omxAlgebra *oa, omxMatrix *om) {
	
	if(oa == NULL) {
		oa = (omxAlgebra*) R_alloc(1, sizeof(omxAlgebra));
	}
	
	omxAlgebraAllocArgs(oa, 0);
	oa->funWrapper = NULL;
	oa->matrix = om;
	oa->oate = NULL;
	om->algebra = oa;
}

void omxDuplicateAlgebra(omxMatrix* tgt, omxMatrix* src, omxState* newState) {

    if(src->algebra != NULL) {
	    omxFillMatrixFromMxAlgebra(tgt, src->algebra->sexpAlgebra, src->nameStr, NULL, 0);
    } else if(src->fitFunction != NULL) {
        omxDuplicateFitMatrix(tgt, src, newState);
    }

}

void omxFreeAlgebraArgs(omxAlgebra *oa) {
	/* Completely destroy the algebra tree */
	
	int j;
	for(j = 0; j < oa->numArgs; j++) {
		omxFreeMatrix(oa->algArgs[j]);
		oa->algArgs[j] = NULL;
	}
	omxAlgebraAllocArgs(oa, 0);
	delete oa;
}

void omxAlgebraPreeval(omxMatrix *mat, FitContext *fc)
{
	if (mat->hasMatrixNumber) mat = fc->lookupDuplicate(mat);
	omxState *st = mat->currentState;
	st->setWantStage(FF_COMPUTE_PREOPTIMIZE);
	omxRecompute(mat, fc);
	auto ff = mat->fitFunction;
	if (ff) fc->fitUnits = ff->units;
	st->setWantStage(FF_COMPUTE_FIT);
}

void omxAlgebraRecompute(omxMatrix *mat, int want, FitContext *fc)
{
	omxAlgebra *oa = mat->algebra;
	if (oa->verbose >= 1) mxLog("recompute algebra '%s'", mat->name());

	if (want & FF_COMPUTE_INITIAL_FIT) {
		bool fvDep = false;
		bool dvDep = false;
		for(int j = 0; j < oa->numArgs; j++) {
			if (oa->algArgs[j]->dependsOnParameters()) {
				if ((oa->verbose + OMX_DEBUG) && !fvDep) {
					mxLog("Algebra %s depends on free parameters "
					      "because of argument[%d] %s",
					      mat->name(), j, oa->algArgs[j]->name());
				}
				fvDep = true;
			}
			if (oa->algArgs[j]->dependsOnDefinitionVariables()) {
				if ((oa->verbose + OMX_DEBUG) && !dvDep) {
					mxLog("Algebra %s depends on definition variables "
					      "because of argument[%d] %s",
					      mat->name(), j, oa->algArgs[j]->name());
				}
				dvDep = true;
			}
		}
		if (fvDep) mat->setDependsOnParameters();
		if (dvDep) mat->setDependsOnDefinitionVariables();
	}

	for(int j = 0; j < oa->numArgs; j++) {
		omxRecompute(oa->algArgs[j], fc);
	}
	if (isErrorRaised()) return;

	if(oa->funWrapper == NULL) {
		if(oa->numArgs != 1) Rf_error("Internal Error: Empty algebra evaluated");
		if(OMX_DEBUG_ALGEBRA) { omxPrint(oa->algArgs[0], "Copy no-op algebra"); }
		omxCopyMatrix(oa->matrix, oa->algArgs[0]);
	} else {
		if(OMX_DEBUG_ALGEBRA || oa->verbose >= 2) { 
			std::string buf;
			for (int ax=0; ax < oa->numArgs; ++ax) {
				if (ax) buf += ", ";
				const char *argname = oa->algArgs[ax]->name();
				buf += argname? argname : "?";
			}
			mxLog("Algebra '%s' %s(%s)", oa->matrix->name(), oa->oate->rName, buf.c_str());
		}
		(*(algebra_op_t)oa->funWrapper)(fc, oa->algArgs, (oa->numArgs), oa->matrix);
	}

	if(OMX_DEBUG_ALGEBRA || oa->verbose >= 3) {
		EigenMatrixAdaptor Emat(oa->matrix);
		int nr = std::min(10, Emat.rows());
		int nc = std::min(10, Emat.cols());
		std::string name = string_snprintf("Algebra '%s' %dx%d", oa->matrix->name(),
						   Emat.rows(), Emat.cols());
		mxPrintMat(name.c_str(), Emat.topLeftCorner(nr, nc));
	}
}

omxAlgebra::omxAlgebra()
{
	verbose = 0;
}

static omxMatrix* omxNewMatrixFromMxAlgebra(SEXP alg, omxState* os, std::string &name)
{
	omxMatrix *om = omxInitMatrix(0, 0, TRUE, os);

	om->hasMatrixNumber = 0;
	om->matrixNumber = 0;	

	omxFillMatrixFromMxAlgebra(om, alg, name, NULL, 0);
	
	return om;
}

void omxFillAlgebraFromTableEntry(omxAlgebra *oa, const omxAlgebraTableEntry* oate, const int realNumArgs)
{
	/* TODO: check for full initialization */
	if(oa == NULL) Rf_error("Internal Error: Null Algebra Detected in fillAlgebra.");

	oa->oate = oate;
	oa->funWrapper = oate->funWrapper;
	omxAlgebraAllocArgs(oa, oate->numArgs==-1? realNumArgs : oate->numArgs);
}

static omxMatrix* omxAlgebraParseHelper(SEXP algebraArg, omxState* os, std::string &name) {
	omxMatrix* newMat;
	
	if(!Rf_isInteger(algebraArg)) {
		newMat = omxNewMatrixFromMxAlgebra(algebraArg, os, name);
	} else {
		newMat = omxMatrixLookupFromState1(algebraArg, os);
	}
	
	return(newMat);
}

void omxFillMatrixFromMxAlgebra(omxMatrix* om, SEXP algebra, std::string &name, SEXP dimnames, int verbose)
{
	int value;
	omxAlgebra *oa = NULL;
	
	value = Rf_asInteger(VECTOR_ELT(algebra, 0));

	if(value > 0) { 			// This is an operator.
		oa = new omxAlgebra;
		oa->verbose = verbose;
		omxInitAlgebraWithMatrix(oa, om);
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[value]);
		if(OMX_DEBUG) {mxLog("Table Entry %d is %s.", value, entry->opName);}
		omxFillAlgebraFromTableEntry(oa, entry, Rf_length(algebra) - 1);
		for(int j = 0; j < oa->numArgs; j++) {
			SEXP algebraArg;
			{
				ScopedProtect p1(algebraArg, VECTOR_ELT(algebra, j+1));
				auto name2 = string_snprintf("%s arg %d", om->name(), j);
				oa->algArgs[j] = omxAlgebraParseHelper(algebraArg, om->currentState, name2);
			}
			if (oa->algArgs[j]->nameStr.size() == 0) {
				// A bit inefficient but invaluable for debugging
				oa->algArgs[j]->nameStr = string_snprintf("alg%03d", ++Global->anonAlgebra);
			}
		}
	} else {		// This is an algebra pointer, and we're a No-op algebra.
		/* TODO: Optimize this by eliminating no-op algebras entirely. */
		SEXP algebraElt;
		ScopedProtect p1(algebraElt, VECTOR_ELT(algebra, 1));
		
		if(!Rf_isInteger(algebraElt)) {   			// A List: only happens if bad optimization has occurred.
			Rf_error("Internal Error: Algebra has been passed incorrectly: detected NoOp: (Operator Arg ...)\n");
		} else {			// Still a No-op.  Sadly, we have to keep it that way.
			
			value = Rf_asInteger(algebraElt);
			
			oa = new omxAlgebra;
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
	om->nameStr     = name;
	oa->sexpAlgebra = algebra;
	om->loadDimnames(dimnames);
}

void omxAlgebraPrint(omxAlgebra* oa, const char* d) {
	omxPrintMatrix(oa->matrix, d);
}

omxMatrix* omxMatrixLookupFromState1(SEXP matrix, omxState* os) {
	int value = 0;

	if (Rf_length(matrix) == 0) return NULL;
	if (Rf_isInteger(matrix)) {
		value = Rf_asInteger(matrix);
		if(value == NA_INTEGER) {
			return NULL;
		}
	} else if (Rf_isReal(matrix)) {
		value = (int) Rf_asReal(matrix);  // prohibit TODO
	} else if (matrix == R_NilValue) {
		return NULL;
	} else if (Rf_isString(matrix)) {
		Rf_error("Internal Rf_error: string passed to omxMatrixLookupFromState1, did you forget to call imxLocateIndex?");
	} else {
		Rf_error("Internal Rf_error: unknown type passed to omxMatrixLookupFromState1");
	}		

	return os->getMatrixFromIndex(value);
}

omxMatrix* omxMatrixLookupFromStateByNumber(int matrix, omxState* os) {
	omxMatrix* output = NULL;
	if(matrix == NA_INTEGER){return NULL;}
	if (matrix >= 0) {
		output = os->algebraList[matrix];
	} 
	else {
		output = os->matrixList[~matrix];
	}
	return output;
}

omxMatrix* omxNewAlgebraFromOperatorAndArgs(int opCode, omxMatrix* args[], int numArgs, omxState* os) {
	
	if(OMX_DEBUG) {mxLog("Generating new algebra from opcode %d (%s).", opCode, omxAlgebraSymbolTable[opCode].rName);}
	omxMatrix *om;
	omxAlgebra *oa = new omxAlgebra;
	omxAlgebraTableEntry* entry = (omxAlgebraTableEntry*)&(omxAlgebraSymbolTable[opCode]);
	if(entry->numArgs >= 0 && entry->numArgs != numArgs) {
		Rf_error("Internal Rf_error: incorrect number of arguments passed to algebra %s.", entry->rName);
	}
	
	om = omxInitAlgebra(oa, os);
	omxFillAlgebraFromTableEntry(oa, entry, entry->numArgs);
	om->nameStr = entry->opName;

	if(OMX_DEBUG) {mxLog("Calculating args for %s.", entry->rName);}
	omxAlgebraAllocArgs(oa, numArgs);
	
	if(OMX_DEBUG) {mxLog("Populating args for %s.", entry->rName);}
	
	for(int i = 0; i < numArgs;i++) {
		oa->algArgs[i] = args[i];
	}
	
	return om;
	
}

