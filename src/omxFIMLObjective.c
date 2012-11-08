/*
 *  Copyright 2007-2012 The OpenMx Project
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
#include "omxFIMLObjective.h"
#include "omxFIMLSingleIteration.h"
#include "omxSadmvnWrapper.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/* FIML Function body */
void omxDestroyFIMLObjective(omxObjective *oo) {
	if(OMX_DEBUG) { Rprintf("Destroying FIML objective object.\n"); }
	omxFIMLObjective *argStruct = (omxFIMLObjective*) (oo->argStruct);

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
	if(oo->subObjective == NULL) {
		if(argStruct->cov != NULL) omxFreeMatrixData(argStruct->cov);
		if(argStruct->means != NULL) omxFreeMatrixData(argStruct->means);
	}
}

void omxPopulateFIMLAttributes(omxObjective *oo, SEXP algebra) {
	omxFIMLObjective *argStruct = ((omxFIMLObjective*)oo->argStruct);
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

	UNPROTECT(3);
}

omxRListElement* omxSetFinalReturnsFIMLObjective(omxObjective *oo, int *numReturns) {

	omxFIMLObjective* ofo = (omxFIMLObjective *) (oo->argStruct);

	omxRListElement* retVal;

	*numReturns = 1;

	if(!ofo->returnRowLikelihoods) {
		retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));
	} else {
		retVal = (omxRListElement*) R_alloc(2, sizeof(omxRListElement));
	}

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);


	if(ofo->returnRowLikelihoods) {
		omxData* data = ofo->data;
		retVal[1].numValues = data->rows;
		retVal[1].values = (double*) R_alloc(data->rows, sizeof(double));
	}

	return retVal;
}

void markDefVarDependencies(omxState* os, omxDefinitionVar* defVar) {

	int numDeps = defVar->numDeps;
	int *deps = defVar->deps;

	omxMatrix** matrixList = os->matrixList;
	omxMatrix** algebraList = os->algebraList;

	for (int i = 0; i < numDeps; i++) {
		int value = deps[i];

		if(value < 0) {
			omxMarkDirty(matrixList[~value]);
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
			omxRaiseError(data->currentState, -1, 
					"Internal error: definition variable population into incorrect data source");
			error("Internal error: definition variable population into incorrect data source"); // Kept for historical reasons
			continue; //Do not populate this variable.
		}
		double newDefVar = omxDoubleDataElement(data, row, defVars[k].column);
		if(ISNA(newDefVar)) {
			omxRaiseError(data->currentState, -1, "Error: NA value for a definition variable is Not Yet Implemented.");
			error("Error: NA value for a definition variable is Not Yet Implemented."); // Kept for historical reasons
			return -1;
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

void omxCallJointFIMLObjective(omxObjective *oo) {	
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
	
	omxObjective* subObjective;

    omxFIMLObjective* ofo = ((omxFIMLObjective*)oo->argStruct);
	omxMatrix* objMatrix  = oo->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState->numChildren;

    cov 		= ofo->cov;
	means		= ofo->means;
	data		= ofo->data;                            //  read-only
	numDefs		= ofo->numDefs;                         //  read-only
	dataColumns	= ofo->dataColumns;
	thresholdCols = ofo->thresholdCols;

	returnRowLikelihoods = ofo->returnRowLikelihoods;   //  read-only
	subObjective = oo->subObjective;


    if(numDefs == 0) {
        if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
        if(subObjective != NULL) {
            omxObjectiveCompute(subObjective);
        } else {
            omxRecompute(cov);			// Only recompute this here if there are no definition vars
            omxRecompute(means);
        }
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
			omxMatrix *childObjective = omxLookupDuplicateElement(parentState->childList[index], objMatrix);
			omxFIMLObjective* childOfo = ((omxFIMLObjective*) childObjective->objective->argStruct);
			omxCopyMatrix(childOfo->cov, cov);
			omxCopyMatrix(childOfo->means, means);
		}
        if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
        if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
    }

	memset(ofo->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);
    
	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxObjective *childObjective = childMatrix->objective;
			if (i == parallelism - 1) {
				omxFIMLSingleIterationJoint(childObjective, oo, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIterationJoint(childObjective, oo, stride * i, stride);
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
		omxFIMLSingleIterationJoint(oo, oo, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofo->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(oo->matrix, 0, 0, sum);
	}

}

void omxCallFIMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n"); }
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

	int numDefs, returnRowLikelihoods;	
	omxObjective* subObjective;
	
	omxMatrix *cov, *means;//, *oldInverse;
	omxData *data;

	omxFIMLObjective* ofo = ((omxFIMLObjective*) oo->argStruct);
	omxMatrix* objMatrix  = oo->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState->numChildren;

	// Locals, for readability.  Should compile out.
	cov 		= ofo->cov;
	means		= ofo->means;
	data		= ofo->data;                            //  read-only
	numDefs		= ofo->numDefs;                         //  read-only
	returnRowLikelihoods = ofo->returnRowLikelihoods;   //  read-only
	subObjective = oo->subObjective;

	if(numDefs == 0) {
		if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
		if(subObjective != NULL) {
			omxObjectiveCompute(subObjective);
		} else {
			omxRecompute(cov);			// Only recompute this here if there are no definition vars
			omxRecompute(means);
		}
		for(int index = 0; index < numChildren; index++) {
			omxMatrix *childObjective = omxLookupDuplicateElement(parentState->childList[index], objMatrix);
			omxFIMLObjective* childOfo = ((omxFIMLObjective*) childObjective->objective->argStruct);
			omxCopyMatrix(childOfo->cov, cov);
			omxCopyMatrix(childOfo->means, means);
		}

		if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
		if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
	}

	memset(ofo->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);
    
	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxObjective *childObjective = childMatrix->objective;
			if (i == parallelism - 1) {
				omxFIMLSingleIteration(childObjective, oo, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIteration(childObjective, oo, stride * i, stride);
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
		omxFIMLSingleIteration(oo, oo, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofo->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(oo->matrix, 0, 0, sum);
	}
}

void omxCallFIMLOrdinalObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
	/* TODO: Current implementation is slow: update by filtering correlations and thresholds. */
	if(OMX_DEBUG) { Rprintf("Beginning Ordinal FIML Evaluation.\n");}
	// Requires: Data, means, covariances, thresholds

	int numDefs;
	int returnRowLikelihoods = 0;

	omxMatrix *cov, *means, *dataColumns;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *corList, *weights;
	
	omxObjective* subObjective;	

	omxFIMLObjective* ofo = ((omxFIMLObjective*)oo->argStruct);
	omxMatrix* objMatrix  = oo->matrix;
	omxState* parentState = objMatrix->currentState;
	int numChildren = parentState->numChildren;

	// Locals, for readability.  Compiler should cut through this.
	cov 		= ofo->cov;
	means		= ofo->means;
	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	numDefs		= ofo->numDefs;

	corList 	= ofo->corList;
	weights		= ofo->weights;
	thresholdCols = ofo->thresholdCols;
	returnRowLikelihoods = ofo->returnRowLikelihoods;
	
	subObjective = oo->subObjective;
	
	if(numDefs == 0) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("No Definition Vars: precalculating."); }
		if(subObjective != NULL) {
			omxObjectiveCompute(subObjective);
		} else {
			omxRecompute(cov);			// Only recompute this here if there are no definition vars
			omxRecompute(means);
		}
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
			omxMatrix *childObjective = omxLookupDuplicateElement(parentState->childList[index], objMatrix);
			omxFIMLObjective* childOfo = ((omxFIMLObjective*) childObjective->objective->argStruct);
			omxCopyMatrix(childOfo->cov, cov);
			omxCopyMatrix(childOfo->means, means);
			memcpy(childOfo->weights, weights, sizeof(double) * cov->rows);
			memcpy(childOfo->corList, corList, sizeof(double) * (cov->rows * (cov->rows - 1)) / 2);
		}
	}

	memset(ofo->rowLogLikelihoods->data, 0, sizeof(double) * data->rows);

	int parallelism = (numChildren == 0) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		#pragma omp parallel for num_threads(parallelism) 
		for(int i = 0; i < parallelism; i++) {
			omxMatrix *childMatrix = omxLookupDuplicateElement(parentState->childList[i], objMatrix);
			omxObjective *childObjective = childMatrix->objective;
			if (i == parallelism - 1) {
				omxFIMLSingleIterationOrdinal(childObjective, oo, stride * i, data->rows - stride * i);
			} else {
				omxFIMLSingleIterationOrdinal(childObjective, oo, stride * i, stride);
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
		omxFIMLSingleIterationOrdinal(oo, oo, 0, data->rows);
	}

	if(!returnRowLikelihoods) {
		double val, sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			val = omxVectorElement(ofo->rowLogLikelihoods, i);
//			Rprintf("%d , %f, %llx\n", i, val, *((unsigned long long*) &val));
			sum += val;
		}	
		if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
		omxSetMatrixElement(oo->matrix, 0, 0, sum);
	}
}

void omxInitFIMLObjective(omxObjective* oo, SEXP rObj) {

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Initializing FIML objective function.\n");
	}

	SEXP nextMatrix;
    omxMatrix *cov, *means;
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	means = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(means == NULL) { 
		omxRaiseError(oo->matrix->currentState, -1, "No means model in FIML evaluation.");
	}
	UNPROTECT(1);	// UNPROTECT(means)

	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	cov = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);	// UNPROTECT(covariance)
	
	omxCreateFIMLObjective(oo, rObj, cov, means);

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("FIML Initialization Completed.");
	}
}

void omxCreateFIMLObjective(omxObjective* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means) {

	SEXP nextMatrix, itemList, nextItem, threshMatrix;
	int nextDef, index, numOrdinal = 0, numContinuous = 0, numCols;

    omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(1, sizeof(omxFIMLObjective));

	numCols = cov->rows;
	
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
    
    /* Set default Objective calls to FIML Objective Calls */
	oo->objectiveFun = omxCallFIMLObjective;
	strncpy(oo->objType, "omxFIMLObjective", 17);
	oo->setFinalReturns = omxSetFinalReturnsFIMLObjective;
	oo->destructFun = omxDestroyFIMLObjective;
	oo->populateAttrFun = omxPopulateFIMLAttributes;
	oo->repopulateFun = NULL;
	

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing data source.\n");
	}
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data"))); // TODO: Need better way to process data elements.
	newObj->data = omxNewDataFromMxDataPtr(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing row likelihood option.\n");
	}
	PROTECT(nextMatrix = AS_INTEGER(GET_SLOT(rObj, install("vector")))); // preparing the object by using the vector to populate and the flag
	newObj->returnRowLikelihoods = INTEGER(nextMatrix)[0];
	if(newObj->returnRowLikelihoods) {
	   omxResizeMatrix(oo->matrix, newObj->data->rows, 1, FALSE); // 1=column matrix, FALSE=discards memory as this is a one time resize
    }
    newObj->rowLikelihoods = omxInitMatrix(NULL, newObj->data->rows, 1, TRUE, oo->matrix->currentState);
    newObj->rowLogLikelihoods = omxInitMatrix(NULL, newObj->data->rows, 1, TRUE, oo->matrix->currentState);
	UNPROTECT(1);

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing variable mapping structure.\n");
	}
	PROTECT(nextMatrix = GET_SLOT(rObj, install("dataColumns")));
	newObj->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, oo->matrix->currentState, 0, 0);
	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		omxPrint(newObj->dataColumns, "Variable mapping");
	}
	UNPROTECT(1);

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing Threshold matrix.\n");
	}
	PROTECT(threshMatrix = GET_SLOT(rObj, install("thresholds")));
    
    if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
        if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
			Rprintf("Accessing Threshold Mappings.\n");
		}
        
        /* Process the data and threshold mapping structures */
    	/* if (threshMatrix == NA_INTEGER), then we could ignore the slot "thresholdColumns"
         * and fill all the thresholdCols with {NULL, 0, 0}.
    	 * However the current path does not have a lot of overhead. */
    	PROTECT(nextMatrix = GET_SLOT(rObj, install("thresholdColumns")));
    	PROTECT(itemList = GET_SLOT(rObj, install("thresholdLevels")));
        int* thresholdColumn, *thresholdNumber;
        thresholdColumn = INTEGER(nextMatrix);
        thresholdNumber = INTEGER(itemList);
    	newObj->thresholdCols = (omxThresholdColumn *) R_alloc(numCols, sizeof(omxThresholdColumn));
    	for(index = 0; index < numCols; index++) {
    		if(thresholdColumn[index] == NA_INTEGER) {	// Continuous variable
    			if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
					Rprintf("Column %d is continuous.\n", index);
				}
    			newObj->thresholdCols[index].matrix = NULL;
    			newObj->thresholdCols[index].column = 0;
    			newObj->thresholdCols[index].numThresholds = 0;
                numContinuous++;
    		} else {
    			newObj->thresholdCols[index].matrix = omxNewMatrixFromMxIndex(threshMatrix, 
    				oo->matrix->currentState);
    			newObj->thresholdCols[index].column = thresholdColumn[index];
    			newObj->thresholdCols[index].numThresholds = thresholdNumber[index];
    			if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
    				Rprintf("Column %d is ordinal with %d thresholds in threshold column %d.\n", 
    				    index, thresholdColumn[index], thresholdNumber[index]);
    			}
    			numOrdinal++;
    		}
    	}
    	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
			Rprintf("%d threshold columns processed.\n", numOrdinal);
		}
    	UNPROTECT(2); /* nextMatrix and itemList ("thresholds" and "thresholdColumns") */
    } else {
        if (OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
			Rprintf("No thresholds matrix; not processing thresholds.");
		}
        numContinuous = numCols;
        newObj->thresholdCols = NULL;
        numOrdinal = 0;
    }
    UNPROTECT(1); /* threshMatrix */

	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);

	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Accessing definition variables structure.\n");
	}
	PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
	newObj->numDefs = length(nextMatrix);
	if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
		Rprintf("Number of definition variables is %d.\n", newObj->numDefs);
	}
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars
	for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
		SEXP dataSource, columnSource, depsSource; 
		int nextDataSource, numDeps;

		PROTECT(itemList = VECTOR_ELT(nextMatrix, nextDef));
		PROTECT(dataSource = VECTOR_ELT(itemList, 0));
		nextDataSource = INTEGER(dataSource)[0];
		if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
			Rprintf("Data source number is %d.\n", nextDataSource);
		}
		newObj->defVars[nextDef].data = nextDataSource;
		newObj->defVars[nextDef].source = oo->matrix->currentState->dataList[nextDataSource];
		PROTECT(columnSource = VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
			Rprintf("Data column number is %d.\n", INTEGER(columnSource)[0]);
		}
		newObj->defVars[nextDef].column = INTEGER(columnSource)[0];
		PROTECT(depsSource = VECTOR_ELT(itemList, 2));
		numDeps = LENGTH(depsSource);
		newObj->defVars[nextDef].numDeps = numDeps;
		newObj->defVars[nextDef].deps = (int*) R_alloc(numDeps, sizeof(int));
		for(int i = 0; i < numDeps; i++) {
			newObj->defVars[nextDef].deps[i] = INTEGER(depsSource)[i];
		}
		UNPROTECT(3); // unprotect dataSource, columnSource, and depsSource



		newObj->defVars[nextDef].numLocations = length(itemList) - 3;
		newObj->defVars[nextDef].matrices = (int *) R_alloc(length(itemList) - 3, sizeof(int));
		newObj->defVars[nextDef].rows = (int *) R_alloc(length(itemList) - 3, sizeof(int));
		newObj->defVars[nextDef].cols = (int *) R_alloc(length(itemList) - 3, sizeof(int));
		newObj->oldDefs[nextDef] = NA_REAL;					// Def Vars default to NA
		for(index = 3; index < length(itemList); index++) {
			PROTECT(nextItem = VECTOR_ELT(itemList, index));
			newObj->defVars[nextDef].matrices[index-3] = INTEGER(nextItem)[0];
			newObj->defVars[nextDef].rows[index-3] = INTEGER(nextItem)[1];
			newObj->defVars[nextDef].cols[index-3] = INTEGER(nextItem)[2];
			UNPROTECT(1); // unprotect nextItem
		}
		UNPROTECT(1); // unprotect itemList
	}
	UNPROTECT(1); // unprotect nextMatrix

    /* Temporary storage for calculation */
    int covCols = newObj->cov->cols;
    // int ordCols = omxDataNumFactor(newObj->data);        // Unneeded, since we don't use it.
    // int contCols = omxDataNumNumeric(newObj->data);
    newObj->smallRow = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
    newObj->smallCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
    newObj->RCX = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
//  newObj->zeros = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);

    omxAliasMatrix(newObj->smallCov, newObj->cov);          // Will keep its aliased state from here on.
    oo->argStruct = (void*)newObj;

    if(numOrdinal > 0 && numContinuous <= 0) {
        if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
            Rprintf("Ordinal Data detected.  Using Ordinal FIML.");
        }
        newObj->weights = (double*) R_alloc(covCols, sizeof(double));
        newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        omxAliasMatrix(newObj->smallMeans, newObj->means);
        newObj->corList = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
        newObj->smallCor = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
        newObj->lThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->uThresh = (double*) R_alloc(covCols, sizeof(double));
        newObj->Infin = (int*) R_alloc(covCols, sizeof(int));

        oo->objectiveFun = omxCallFIMLOrdinalObjective;
    } else if(numOrdinal > 0) {
        if(OMX_DEBUG && oo->matrix->currentState->parentState == NULL) {
            Rprintf("Ordinal and Continuous Data detected.  Using Joint Ordinal/Continuous FIML.");
        }

        newObj->weights = (double*) R_alloc(covCols, sizeof(double));
        newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->ordMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->contRow = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->ordRow = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->ordCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
        newObj->ordContCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
        newObj->halfCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
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

        oo->objectiveFun = omxCallJointFIMLObjective;
    }
}
