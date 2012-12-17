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

#include "R.h"
#include <Rinternals.h>
#include <Rdefines.h>

#include <sys/stat.h>

#include "omxDefines.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "npsolWrap.h"

void omxFinalAlgebraCalculation(omxState *currentState, SEXP matrices, SEXP algebras, SEXP expectations) {
	SEXP nextMat, algebra;
	for(int index = 0; index < currentState->numMats; index++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Matrix %d.\n", index); }
		omxMatrix* nextMatrix = currentState->matrixList[index];
		omxRecompute(nextMatrix);
		nextMat = omxExportMatrix(nextMatrix);
		SET_VECTOR_ELT(matrices, index, nextMat);
		UNPROTECT(1);	/* nextMat */
	}

	for(int index = 0; index < currentState->numAlgs; index++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation and Copy of Algebra %d.\n", index); }
		omxMatrix* nextAlgebra = currentState->algebraList[index];
		omxRecompute(nextAlgebra);
		algebra = omxExportMatrix(nextAlgebra);
		/* If an fit function, populate attributes.  Will skip if not fit function. */
		omxFitFunction* currentFit = nextAlgebra->fitFunction;
		if(currentFit != NULL) {
			if(OMX_DEBUG) { Rprintf("Algebra %d is a fit function.\n", index); }
			if(currentFit->populateAttrFun != NULL) {
				if(OMX_DEBUG) { Rprintf("Algebra %d has attribute population.\n", index); }
				currentFit->populateAttrFun(currentFit, algebra);
		    }
		}

		if(OMX_DEBUG) { Rprintf("Final Calculation of Algebra %d Complete.\n", index); }
		SET_VECTOR_ELT(algebras, index, algebra);

		UNPROTECT(1);	/* algebra */
	}
	if(OMX_DEBUG) { Rprintf("All Algebras complete.\n"); }
	
	for(int index = 0; index < currentState->numExpects; index++) {
		if(OMX_DEBUG) { Rprintf("Final Calculation of Expectation %d.\n", index); }
		omxExpectation* nextExpectation = currentState->expectationList[index];
		omxExpectationRecompute(nextExpectation);
		SEXP rExpect;
		PROTECT(rExpect = allocVector(LGLSXP, 1));
		if(nextExpectation->populateAttrFun != NULL) {
			if(OMX_DEBUG) { Rprintf("Expectation %d has attribute population.\n", index); }
			nextExpectation->populateAttrFun(nextExpectation, rExpect);
	    }
		SET_VECTOR_ELT(expectations, index, rExpect);
		UNPROTECT(1); /* rExpect */
	}
}

void omxPopulateFitFunction(omxState *currentState, int numReturns, SEXP *ans, SEXP *names) {
	omxMatrix* om = currentState->fitMatrix;
	if(om != NULL) {					// In the event of a no-fit function run.
		omxFitFunction* off = om->fitFunction;
		if(OMX_DEBUG) { Rprintf("Checking for additional fit function info.\n"); }

		if(off != NULL && off->setFinalReturns != NULL) {
			if(OMX_DEBUG) { Rprintf("Expecting fit function Info....");}
			int numEls;
			SEXP oElement;
			omxRListElement* orle = off->setFinalReturns(off, &numEls);
			PROTECT(*ans = allocVector(VECSXP, numReturns + numEls));
			PROTECT(*names = allocVector(STRSXP, numReturns + numEls));
			if(numEls != 0) {
				if(OMX_DEBUG) { Rprintf("Adding %d sets of fit function Info....", numEls);}
				for(int i = 0; i < numEls; i++) {
					PROTECT(oElement = allocVector(REALSXP, orle[i].numValues));
					for(int j = 0; j < orle[i].numValues; j++)
						REAL(oElement)[j] = orle[i].values[j];
					SET_STRING_ELT(*names, i+numReturns, mkChar(orle[i].label));
					SET_VECTOR_ELT(*ans, i+numReturns, oElement);
					UNPROTECT(1); // oElement
				}
			}
		} else {
			PROTECT(*ans = allocVector(VECSXP, numReturns));
			PROTECT(*names = allocVector(STRSXP, numReturns));
		}
		if(OMX_DEBUG) { Rprintf("Done.\n");}
	} else {
		PROTECT(*ans = allocVector(VECSXP, numReturns));
		PROTECT(*names = allocVector(STRSXP, numReturns));
	}
}

void omxPopulateHessians(int numHessians, omxMatrix* currentFit, 
		SEXP calculatedHessian, SEXP stdErrors, int calculateStdErrors, int n) {
	if(OMX_DEBUG) { Rprintf("Populating hessians for %d fit functions.\n", numHessians); }
	omxFitFunction* off = currentFit->fitFunction;
	if(off->hessian == NULL) {
		if(OMX_DEBUG) { Rprintf("Fit function 0 has no hessian. Aborting.\n");}
		return;
	}

	if(OMX_DEBUG) { Rprintf("Fit function 0 has hessian at 0x%x.\n", off->hessian);}

	double* hessian  = REAL(calculatedHessian);
	double* stdError = REAL(stdErrors);
	for(int k = 0; k < n * n; k++) {
		if(OMX_DEBUG) {Rprintf("Populating hessian at %d.\n", k);}
		hessian[k] = off->hessian[k];		// For expediency, ignore majority for symmetric matrices.
	}
	if(calculateStdErrors) {
		if(off->stdError == NULL) {
			for(int k = 0; k < n; k++) {
				if(OMX_DEBUG) {Rprintf("Populating NA standard error at %d.\n", k);}
				stdError[k] = R_NaReal;
			}
		} else {
			for(int k = 0; k < n; k++) {
				if(OMX_DEBUG) {Rprintf("Populating standard error at %d.\n", k);}
				stdError[k] = off->stdError[k];
			}
		}
	}
}

void omxPopulateConfidenceIntervals(omxState* currentState, SEXP intervals, SEXP intervalCodes) {
	int numInts = currentState->numIntervals;
	if(OMX_DEBUG) { Rprintf("Populating CIs for %d fit functions.\n", numInts); }
	double* interval = REAL(intervals);
	int* intervalCode = INTEGER(intervalCodes);
	for(int j = 0; j < numInts; j++) {
		omxConfidenceInterval *oCI = &(currentState->intervalList[j]);
		interval[j] = oCI->min;
		interval[j + numInts] = oCI->max;
		intervalCode[j] = oCI->lCode;
		intervalCode[j + numInts] = oCI->uCode;
	}
}
