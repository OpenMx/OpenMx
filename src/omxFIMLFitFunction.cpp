/*
 *  Copyright 2007-2013 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"
#include "omxFIMLSingleIteration.h"
#include "omxSadmvnWrapper.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* FIML Function body */
void omxDestroyFIMLFitFunction(omxFitFunction *off) {
	if(OMX_DEBUG) { Rprintf("Destroying FIML fit function object.\n"); }
	omxFIMLFitFunction *argStruct = (omxFIMLFitFunction*) (off->argStruct);

	if(argStruct->smallMeans != NULL) omxFreeMatrixData(argStruct->smallMeans);
	if(argStruct->ordMeans != NULL) omxFreeMatrixData(argStruct->ordMeans);
	if(argStruct->contRow != NULL) omxFreeMatrixData(argStruct->contRow);
	if(argStruct->ordRow != NULL) omxFreeMatrixData(argStruct->ordRow);
	if(argStruct->ordCov != NULL) omxFreeMatrixData(argStruct->ordCov);
	if(argStruct->ordContCov != NULL) omxFreeMatrixData(argStruct->ordContCov);
	if(argStruct->halfCov != NULL) omxFreeMatrixData(argStruct->halfCov);
	if(argStruct->reduceCov != NULL) omxFreeMatrixData(argStruct->reduceCov);

	if(argStruct->smallRow != NULL) omxFreeMatrixData(argStruct->smallRow);
	if(argStruct->smallCov != NULL) omxFreeMatrixData(argStruct->smallCov);
	if(argStruct->RCX != NULL)		omxFreeMatrixData(argStruct->RCX);
    if(argStruct->rowLikelihoods != NULL) omxFreeMatrixData(argStruct->rowLikelihoods);
    if(argStruct->rowLogLikelihoods != NULL) omxFreeMatrixData(argStruct->rowLogLikelihoods);
	if(off->expectation == NULL) {
		if(argStruct->cov != NULL) omxFreeMatrixData(argStruct->cov);
		if(argStruct->means != NULL) omxFreeMatrixData(argStruct->means);
	}
}

void omxPopulateFIMLAttributes(omxFitFunction *off, SEXP algebra) {
	omxFIMLFitFunction *argStruct = ((omxFIMLFitFunction*)off->argStruct);
	SEXP expCovExt, expMeanExt, rowLikelihoodsExt;
	omxMatrix *expCovInt, *expMeanInt, *rowLikelihoodsInt;
	expCovInt = argStruct->cov;
	expMeanInt = argStruct->means;
	rowLikelihoodsInt = argStruct->rowLikelihoods;

	PROTECT(expCovExt = allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL) {
		PROTECT(expMeanExt = allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		PROTECT(expMeanExt = allocMatrix(REALSXP, 0, 0));		
	}
	PROTECT(rowLikelihoodsExt = allocVector(REALSXP, rowLikelihoodsInt->rows));
	for(int row = 0; row < rowLikelihoodsInt->rows; row++)
		REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);

	setAttrib(algebra, install("expCov"), expCovExt);
	setAttrib(algebra, install("expMean"), expMeanExt);
	setAttrib(algebra, install("likelihoods"), rowLikelihoodsExt);

	UNPROTECT(3); // expCovExp, expCovInt, rowLikelihoodsExt
}

omxRListElement* omxSetFinalReturnsFIMLFitFunction(omxFitFunction *off, int *numReturns) {

	omxFIMLFitFunction* ofiml = (omxFIMLFitFunction *) (off->argStruct);

	omxRListElement* retVal;

	*numReturns = 1;

	if(!ofiml->returnRowLikelihoods) {
		retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));
	} else {
		retVal = (omxRListElement*) R_alloc(2, sizeof(omxRListElement));
	}

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(off->matrix, 0, 0);


	if(ofiml->returnRowLikelihoods) {
		omxData* data = ofiml->data;
		retVal[1].numValues = data->rows;
		retVal[1].values = (double*) R_alloc(data->rows, sizeof(double));
	}

	return retVal;
}

void markDefVarDependencies(omxState* os, omxDefinitionVar* defVar) {

	int numDeps = defVar->numDeps;
	int *deps = defVar->deps;

	omxMatrix** algebraList = os->algebraList;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(os->matrixList[~value]);
		} else {
			omxMarkDirty(algebraList[value]);
		}
	}

}

int handleDefinitionVarList(omxData* data, omxState *state, int row, omxDefinitionVar* defVars, double* oldDefs, int numDefs) {

	if(OMX_DEBUG_ROWS(row)) { Rprintf("Processing Definition Vars.\n"); }
	
	int numVarsFilled = 0;

	/* Fill in Definition Var Estimates */
	for(int k = 0; k < numDefs; k++) {
		if(defVars[k].source != data) {
			error("Internal error: definition variable population into incorrect data source");
		}
		double newDefVar = omxDoubleDataElement(data, row, defVars[k].column);
		if(ISNA(newDefVar)) {
			error("Error: NA value for a definition variable is Not Yet Implemented.");
		}
		if(newDefVar == oldDefs[k]) {
			continue;	// NOTE: Potential speedup vs accuracy tradeoff here using epsilon comparison
		}
		oldDefs[k] = newDefVar;
		numVarsFilled++;

		for(int l = 0; l < defVars[k].numLocations; l++) {
			if(OMX_DEBUG_ROWS(row)) {
				Rprintf("Populating column %d (value %3.2f) into matrix %d.\n", defVars[k].column, omxDoubleDataElement(defVars[k].source, row, defVars[k].column), defVars[k].matrices[l]);
			}
			int matrixNumber = defVars[k].matrices[l];
			int matrow = defVars[k].rows[l];
			int matcol = defVars[k].cols[l];
			omxMatrix *matrix = state->matrixList[matrixNumber];
			omxSetMatrixElement(matrix, matrow, matcol, newDefVar);
		}
		markDefVarDependencies(state, &(defVars[k]));
	}
	return numVarsFilled;
}

static void omxCallJointFIMLFitFunction(omxFitFunction *off, int want, double *gradient) {	
	// TODO: Figure out how to give access to other per-iteration structures.
	// TODO: Current implementation is slow: update by filtering correlations and thresholds.
	// TODO: Current implementation does not implement speedups for sorting.
	// TODO: Current implementation may fail on all-continuous-missing or all-ordinal-missing rows.
	
    if(OMX_DEBUG) { 
	    Rprintf("Beginning Joint FIML Evaluation.\n");
    }
	// Requires: Data, means, covariances, thresholds

	int numDefs;

	int returnRowLikelihoods = 0;

	omxMatrix *cov, *means, *dataColumns;

	omxThresholdColumn *thresholdCols;
	omxData* data;
	
	omxExpectation* expectation;

	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);
	omxMatrix* fitMatrix  = off->matrix;
	omxState* parentState = fitMatrix->currentState;
	int numChildren = parentState->numChildren;

	cov 		= ofiml->cov;
	means		= ofiml->means;
	data		= ofiml->data;                            //  read-only
	numDefs		= ofiml->numDefs;                         //  read-only
	dataColumns	= ofiml->dataColumns;
	thresholdCols = ofiml->thresholdCols;

	returnRowLikelihoods = ofiml->returnRowLikelihoods;   //  read-only
	expectation = off->expectation;


    if(numDefs == 0) {
        if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
		omxExpectationRecompute(expectation);
		// MCN Also do the threshold formulae!
		
		for(int j=0; j < dataColumns->cols; j++) {
			int var = omxVectorElement(dataColumns, j);
			if(omxDataColumnIsFactor(data, j) && thresholdCols[var].numThresholds > 0) { // j is an ordinal column
				omxMatrix* nextMatrix = thresholdCols[var].matrix;
				omxRecompute(nextMatrix);
				checkIncreasing(nextMatrix, thresholdCols[var].column);
				for(int index = 0; index < numChildren; index++) {
					omxMatrix *target = omxLookupDuplicateElement(parentState->childList[index], nextMatrix);
					omxCopyMatrix(target, nextMatrix);
				}
            }
        }
		for(int index = 0; index < numChildren; index++) {
			omxMatrix *childFit = omxLookupDuplicateElement(parentState->childList[index], fitMatrix);
			omxFIMLFitFunction* childOfiml = ((omxFIMLFitFunction*) childFit->fitFunction->argStruct);
			omxCopyMatrix(childOfiml->cov, cov);
			omxCopyMatrix(childOfiml->means, means);
		}
        if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
        if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
    }

	memset(ofiml->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);
    
	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				omxFIMLSingleIterationJoint(childFit, off, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIterationJoint(childFit, off, stride * i, stride);
			}
		}

		for(int i = 0; i < parallelism; i++) {
			if (parentState->childList[i]->statusCode < 0) {
				parentState->statusCode = parentState->childList[i]->statusCode;
				strncpy(parentState->statusMsg, parentState->childList[i]->statusMsg, 249);
				parentState->statusMsg[249] = '\0';
			}
		}

	} else {
		omxFIMLSingleIterationJoint(off, off, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofiml->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(off->matrix, 0, 0, sum);
	}

}

static void omxCallFIMLFitFunction(omxFitFunction *off, int want, double *gradient) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n"); }
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

	int numDefs, returnRowLikelihoods;	
	omxExpectation* expectation;
	
	omxMatrix *cov, *means;//, *oldInverse;
	omxData *data;

	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*) off->argStruct);
	omxMatrix* objMatrix  = off->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState->numChildren;

	// Locals, for readability.  Should compile out.
	cov 		= ofiml->cov;
	means		= ofiml->means;
	data		= ofiml->data;                            //  read-only
	numDefs		= ofiml->numDefs;                         //  read-only
	returnRowLikelihoods = ofiml->returnRowLikelihoods;   //  read-only
	expectation = off->expectation;

	if(numDefs == 0 && strcmp(expectation->expType, "MxExpectationStateSpace")) {
		if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
		omxExpectationCompute(expectation);
		
		for(int index = 0; index < numChildren; index++) {
			omxMatrix *childFit = omxLookupDuplicateElement(parentState->childList[index], objMatrix);
			omxFIMLFitFunction* childOfiml = ((omxFIMLFitFunction*) childFit->fitFunction->argStruct);
			omxCopyMatrix(childOfiml->cov, cov);
			omxCopyMatrix(childOfiml->means, means);
		}

		if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
		if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
	}

	memset(ofiml->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);
    
	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				omxFIMLSingleIteration(childFit, off, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIteration(childFit, off, stride * i, stride);
			}
		}

		for(int i = 0; i < parallelism; i++) {
			if (parentState->childList[i]->statusCode < 0) {
				parentState->statusCode = parentState->childList[i]->statusCode;
				strncpy(parentState->statusMsg, parentState->childList[i]->statusMsg, 249);
				parentState->statusMsg[249] = '\0';
			}
		}

	} else {
		omxFIMLSingleIteration(off, off, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofiml->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(off->matrix, 0, 0, sum);
	}
}

static void omxCallFIMLOrdinalFitFunction(omxFitFunction *off, int want, double *gradient) {	// TODO: Figure out how to give access to other per-iteration structures.
	/* TODO: Current implementation is slow: update by filtering correlations and thresholds. */
	if(OMX_DEBUG) { Rprintf("Beginning Ordinal FIML Evaluation.\n");}
	// Requires: Data, means, covariances, thresholds

	int numDefs;
	int returnRowLikelihoods = 0;

	omxMatrix *cov, *means, *dataColumns;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *corList, *weights;
	
	omxExpectation* expectation;

	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);
	omxMatrix* objMatrix  = off->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState->numChildren;

	// Locals, for readability.  Compiler should cut through this.
	cov 		= ofiml->cov;
	means		= ofiml->means;
	data		= ofiml->data;
	dataColumns	= ofiml->dataColumns;
	numDefs		= ofiml->numDefs;

	corList 	= ofiml->corList;
	weights		= ofiml->weights;
	thresholdCols = ofiml->thresholdCols;
	returnRowLikelihoods = ofiml->returnRowLikelihoods;
	
	expectation = off->expectation;
	
	if(numDefs == 0) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("No Definition Vars: precalculating."); }
		omxExpectationCompute(expectation);
		for(int j = 0; j < dataColumns->cols; j++) {
			if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
				omxMatrix* nextMatrix = thresholdCols[j].matrix;
				omxRecompute(nextMatrix);
				checkIncreasing(nextMatrix, thresholdCols[j].column);
				for(int index = 0; index < numChildren; index++) {
					omxMatrix *target = omxLookupDuplicateElement(parentState->childList[index], nextMatrix);
					omxCopyMatrix(target, nextMatrix);
				}
			}
		}
		omxStandardizeCovMatrix(cov, corList, weights);	// Calculate correlation and covariance
		for(int index = 0; index < numChildren; index++) {
			omxMatrix *childFit = omxLookupDuplicateElement(parentState->childList[index], objMatrix);
			omxFIMLFitFunction* childOfiml = ((omxFIMLFitFunction*) childFit->fitFunction->argStruct);
			omxCopyMatrix(childOfiml->cov, cov);
			omxCopyMatrix(childOfiml->means, means);
			memcpy(childOfiml->weights, weights, sizeof(double) * cov->rows);
			memcpy(childOfiml->corList, corList, sizeof(double) * (cov->rows * (cov->rows - 1)) / 2);
		}
	}

	memset(ofiml->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);

	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				omxFIMLSingleIterationOrdinal(childFit, off, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIterationOrdinal(childFit, off, stride * i, stride);
			}
		}

		for(int i = 0; i < parallelism; i++) {
			if (parentState->childList[i]->statusCode < 0) {
				parentState->statusCode = parentState->childList[i]->statusCode;
				strncpy(parentState->statusMsg, parentState->childList[i]->statusMsg, 249);
				parentState->statusMsg[249] = '\0';
			}
		}

	} else {
		omxFIMLSingleIterationOrdinal(off, off, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofiml->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(off->matrix, 0, 0, sum);
	}
}

void omxInitFIMLFitFunction(omxFitFunction* off, SEXP rObj) {

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Initializing FIML fit function function.\n");
	}

	SEXP nextMatrix;
	int numOrdinal = 0, numContinuous = 0;
	omxMatrix *cov, *means;

	omxFIMLFitFunction *newObj = (omxFIMLFitFunction*) R_alloc(1, sizeof(omxFIMLFitFunction));
	omxExpectation* expectation = off->expectation;
	if(expectation == NULL) {
		omxRaiseError(off->matrix->currentState, -1, "FIML cannot fit without model expectations.");
		return;
	}

	cov = omxGetExpectationComponent(expectation, off, "cov");
	if(cov == NULL) { 
		omxRaiseError(off->matrix->currentState, -1, "No covariance expectation in FIML evaluation.");
		return;
	}

	means = omxGetExpectationComponent(expectation, off, "means");
	
	if(means == NULL) { 
		omxRaiseError(off->matrix->currentState, -1, "No means model in FIML evaluation.");
		return;
	}

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("FIML Initialization Completed.");
	}
	
    newObj->cov = cov;
    newObj->means = means;
    newObj->smallMeans = NULL;
    newObj->ordMeans   = NULL;
    newObj->contRow    = NULL;
    newObj->ordRow     = NULL;
    newObj->ordCov     = NULL;
    newObj->ordContCov = NULL;
    newObj->halfCov    = NULL;
    newObj->reduceCov  = NULL;
    
    /* Set default FitFunction calls to FIML FitFunction Calls */
	off->computeFun = omxCallFIMLFitFunction;
	off->fitType = "omxFIMLFitFunction";
	off->setFinalReturns = omxSetFinalReturnsFIMLFitFunction;
	off->destructFun = omxDestroyFIMLFitFunction;
	off->populateAttrFun = omxPopulateFIMLAttributes;
	off->repopulateFun = NULL;

	off->usesChildModels = TRUE;

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing data source.\n");
	}
	newObj->data = off->expectation->data;

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing row likelihood option.\n");
	}
	PROTECT(nextMatrix = AS_INTEGER(GET_SLOT(rObj, install("vector")))); // preparing the object by using the vector to populate and the flag
	newObj->returnRowLikelihoods = INTEGER(nextMatrix)[0];
	if(newObj->returnRowLikelihoods) {
	   omxResizeMatrix(off->matrix, newObj->data->rows, 1, FALSE); // 1=column matrix, FALSE=discards memory as this is a one time resize
    }
    newObj->rowLikelihoods = omxInitMatrix(NULL, newObj->data->rows, 1, TRUE, off->matrix->currentState);
    newObj->rowLogLikelihoods = omxInitMatrix(NULL, newObj->data->rows, 1, TRUE, off->matrix->currentState);
	UNPROTECT(1); // nextMatrix

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing variable mapping structure.\n");
	}
	newObj->dataColumns = off->expectation->dataColumns;

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing Threshold matrix.\n");
	}
	newObj->thresholdCols = off->expectation->thresholds;
	numOrdinal = off->expectation->numOrdinal;
	numContinuous = newObj->dataColumns->cols - off->expectation->numOrdinal;

	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);
	
	newObj->numDefs = off->expectation->numDefs;
	newObj->defVars = off->expectation->defVars;
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars

	if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing definition variables structure.\n");
	}
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars
	memset(newObj->oldDefs, NA_REAL, sizeof(double) * newObj->numDefs);			// Does this work?
	// for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
	// 	newObj->oldDefs[nextDef] = NA_REAL;					// Def Vars default to NA
	// }

    /* Temporary storage for calculation */
    int covCols = newObj->cov->cols;
	if(OMX_DEBUG){Rprintf("Number of columns found is %d", covCols);}
    // int ordCols = omxDataNumFactor(newObj->data);        // Unneeded, since we don't use it.
    // int contCols = omxDataNumNumeric(newObj->data);
    newObj->smallRow = omxInitMatrix(NULL, 1, covCols, TRUE, off->matrix->currentState);
    newObj->smallCov = omxInitMatrix(NULL, covCols, covCols, TRUE, off->matrix->currentState);
    newObj->RCX = omxInitMatrix(NULL, 1, covCols, TRUE, off->matrix->currentState);
//  newObj->zeros = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, off->matrix->currentState);

    omxAliasMatrix(newObj->smallCov, newObj->cov);          // Will keep its aliased state from here on.
    off->argStruct = (void*)newObj;

    if(numOrdinal > 0 && numContinuous <= 0) {
        if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
            Rprintf("Ordinal Data detected.  Using Ordinal FIML.");
        }
        newObj->weights = (double*) R_alloc(covCols, sizeof(double));
        newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, off->matrix->currentState);
        omxAliasMatrix(newObj->smallMeans, newObj->means);
        newObj->corList = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
        newObj->smallCor = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
        newObj->lThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->uThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->Infin = (int*) R_alloc(covCols, sizeof(int));

        off->computeFun = omxCallFIMLOrdinalFitFunction;
    } else if(numOrdinal > 0) {
        if(OMX_DEBUG && off->matrix->currentState->parentState == NULL) {
            Rprintf("Ordinal and Continuous Data detected.  Using Joint Ordinal/Continuous FIML.");
        }

        newObj->weights = (double*) R_alloc(covCols, sizeof(double));
        newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, off->matrix->currentState);
        newObj->ordMeans = omxInitMatrix(NULL, covCols, 1, TRUE, off->matrix->currentState);
        newObj->contRow = omxInitMatrix(NULL, covCols, 1, TRUE, off->matrix->currentState);
        newObj->ordRow = omxInitMatrix(NULL, covCols, 1, TRUE, off->matrix->currentState);
        newObj->ordCov = omxInitMatrix(NULL, covCols, covCols, TRUE, off->matrix->currentState);
        newObj->ordContCov = omxInitMatrix(NULL, covCols, covCols, TRUE, off->matrix->currentState);
        newObj->halfCov = omxInitMatrix(NULL, covCols, covCols, TRUE, off->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(NULL, covCols, covCols, TRUE, off->matrix->currentState);
        omxAliasMatrix(newObj->smallMeans, newObj->means);
        omxAliasMatrix(newObj->ordMeans, newObj->means);
        omxAliasMatrix(newObj->contRow, newObj->smallRow );
        omxAliasMatrix(newObj->ordRow, newObj->smallRow );
        omxAliasMatrix(newObj->ordCov, newObj->cov);
        omxAliasMatrix(newObj->ordContCov, newObj->cov);
        omxAliasMatrix(newObj->smallMeans, newObj->means);
        newObj->corList = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
        newObj->lThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->uThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->Infin = (int*) R_alloc(covCols, sizeof(int));

        off->computeFun = omxCallJointFIMLFitFunction;
    }
}
