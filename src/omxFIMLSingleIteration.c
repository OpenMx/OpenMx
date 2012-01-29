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
#include "omxSadmvnWrapper.h"

/**
 * The localobj reference is used to access read-only variables,
 * or variables that can be modified but whose state cannot be
 * accessed from other threads.
 *
 * The sharedobj reference is used to access write-only variables,
 * where the memory writes of any two threads are non-overlapping.
 * No synchronization mechanisms are employed to maintain consistency
 * of sharedobj references.
 *
 *
 * Because (1) these functions may be invoked with arbitrary 
 * rowbegin and rowcount values, and (2) the log-likelihood
 * values for all data rows must be calculated (even in cases
 * of errors), this function is forbidden from return()-ing early.
 *
 * As another consequence of (1) and (2), if "rowbegin" is in
 * the middle of a sequence of identical rows, then defer
 * move "rowbegin" to after the sequence of identical rows.
 * Grep for "[[Comment 4]]" in source code.
 */
void omxFIMLSingleIterationOrdinal(omxObjective *localobj, omxObjective *sharedobj, int rowbegin, int rowcount) {

    omxFIMLObjective* ofo = ((omxFIMLObjective*) localobj->argStruct);
    omxFIMLObjective* shared_ofo = ((omxFIMLObjective*) sharedobj->argStruct);

	double Q = 0.0;
	double* oldDefs;
	int numDefs;
	int numRemoves;
	int returnRowLikelihoods;
	int keepCov = 0, keepInverse = 0;

	omxObjective* subObjective;
	
	omxMatrix *cov, *means, *smallCov, *dataColumns;//, *oldInverse;
    omxMatrix *rowLikelihoods, *rowLogLikelihoods;;
    omxThresholdColumn *thresholdCols;
    double *lThresh, *uThresh, *corList, *weights;
	int *Infin;
	omxDefinitionVar* defVars;
	omxData *data;

	// Locals, for readability.  Should compile out.
	cov 		     = ofo->cov;
	means		     = ofo->means;
	smallCov 	     = ofo->smallCov;
	oldDefs		     = ofo->oldDefs;
	data		     = ofo->data;                       //  read-only
	dataColumns	     = ofo->dataColumns;                //  read-only
	defVars		     = ofo->defVars;                    //  read-only
	numDefs		     = ofo->numDefs;                    //  read-only
	returnRowLikelihoods = ofo->returnRowLikelihoods;   //  read-only
	rowLikelihoods    = shared_ofo->rowLikelihoods;     // write-only
	rowLogLikelihoods = shared_ofo->rowLogLikelihoods;  // write-only

    corList          = ofo->corList;
    weights          = ofo->weights;
    lThresh          = ofo->lThresh;
    uThresh          = ofo->uThresh;
    thresholdCols    = ofo->thresholdCols;

    Infin            = ofo->Infin;

	subObjective = localobj->subObjective;

	int firstRow = 1;
    int row = rowbegin;

    resetDefinitionVariables(oldDefs, numDefs);

	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}

	while(row < data->rows && (row - rowbegin) < rowcount) {
        if(OMX_DEBUG_ROWS(row)) {Rprintf("Row %d.\n", row);}
        localobj->matrix->currentState->currentRow = row;		// Set to a new row.
		int numIdentical = omxDataNumIdenticalRows(data, row);
		if(numIdentical == 0) numIdentical = 1; 
		// N.B.: numIdentical == 0 means an error occurred and was not properly handled;
		// it should never be the case.

        Q = 0.0;

        // Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
        numRemoves = 0;

        // Handle Definition Variables.
        if(numDefs != 0) {
			if(keepCov <= 0) {  // If we're keeping covariance from the previous row, do not populate 
				if(OMX_DEBUG_ROWS(row)) { Rprintf("Handling Definition Vars.\n"); }
				if(handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs) || firstRow) {
					// Use firstrow instead of rows == 0 for the case where the first row is all NAs
					// N.B. handling of definition var lists always happens, regardless of firstRow.
					if(!(subObjective == NULL)) {
						omxObjectiveCompute(subObjective);
					} else {
						omxRecompute(cov);
						omxRecompute(means);
					}
					for(int j=0; j < dataColumns->cols; j++) {
						if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
							omxRecompute(thresholdCols[j].matrix);
							checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column);
						}
					}
					// Calculate correlation matrix from covariance
					omxStandardizeCovMatrix(cov, corList, weights);
				}
			}
		}

		// Filter down correlation matrix and calculate thresholds

		for(int j = 0; j < dataColumns->cols; j++) {
			int var = omxVectorElement(dataColumns, j);
			int value = omxIntDataElement(data, row, var); // Indexing correction means this is the index of the upper bound +1.
			if(ISNA(value) || value == NA_INTEGER) {  // Value is NA, therefore filter.
				numRemoves++;
				// toRemove[j] = 1;
				Infin[j] = -1;
				continue;
			} else {			// For joint, check here for continuousness
				value--;		// Correct for C indexing: value is now the index of the upper bound.
				// Note : Tested subsampling of the corList and thresholds for speed. 
				//			Doesn't look like that's much of a speedup.
				double mean = (means == NULL) ? 0 : omxVectorElement(means, j);
				double weight = weights[j];
				if(OMX_DEBUG_ROWS(row)) { 
					Rprintf("Row %d, column %d. Mean is %f and weight is %f\n", row, j, mean, weight);
				}
				if(value == 0) { 									// Lowest threshold = -Inf
					lThresh[j] = (omxMatrixElement(thresholdCols[j].matrix, 0, thresholdCols[j].column) - mean) / weight;
					uThresh[j] = lThresh[j];
					Infin[j] = 0;
				} else {
					lThresh[j] = (omxMatrixElement(thresholdCols[j].matrix, value-1, thresholdCols[j].column) - mean) / weight;
					if(thresholdCols[j].numThresholds > value) {	// Highest threshold = Inf
						double tmp = (omxMatrixElement(thresholdCols[j].matrix, value, thresholdCols[j].column) - mean) / weight;
						uThresh[j] = tmp;
						Infin[j] = 2;
					} else {
						uThresh[j] = NA_INTEGER; // NA is a special to indicate +Inf
						Infin[j] = 1;
					}
				}
				
				if(uThresh[j] == NA_INTEGER || isnan(uThresh[j])) { // for matrix-style specification.
					uThresh[j] = lThresh[j];
					Infin[j] = 1;
				}

				if(OMX_DEBUG_ROWS(row)) { 
					Rprintf("Row %d, column %d.  Thresholds for data column %d and row %d are %f -> %f. (Infin=%d)\n", 
						row, j, var, value, lThresh[j], uThresh[j], Infin[j]);
				}
			}
		}

		if(numRemoves >= smallCov->rows) {
			for(int nid = 0; nid < numIdentical; nid++) {
				if(returnRowLikelihoods) {
					omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 1.0);
				}
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 1.0);
			}
    		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
    		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
            row += numIdentical;
    		keepCov -= numIdentical;
    		keepInverse -= numIdentical;
			continue;
		}

   		double likelihood;
		int inform;

		omxSadmvnWrapper(localobj, cov, smallCov, corList, lThresh, uThresh, Infin, &likelihood, &inform);

		if(inform == 2) {
			if(!returnRowLikelihoods) {
				char helperstr[200];
				char *errstr = calloc(250, sizeof(char));
				sprintf(helperstr, "Improper value detected by integration routine in data row %d: Most likely the expected covariance matrix is not positive-definite", omxDataIndex(data, row));
				if(localobj->matrix->currentState->computeCount <= 0) {
					sprintf(errstr, "%s at starting values.\n", helperstr);
				} else {
					sprintf(errstr, "%s at major iteration %d.\n", helperstr, localobj->matrix->currentState->majorIteration);
				}
				omxRaiseError(localobj->matrix->currentState, -1, errstr);
				free(errstr);
			}
			for(int nid = 0; nid < numIdentical; nid++) {
				if (returnRowLikelihoods)
					omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
			}
       		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
       		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
			if(OMX_DEBUG) {Rprintf("Improper input to sadmvn in row likelihood.  Skipping Row.");}
			row += numIdentical;
			keepCov -= numIdentical;
			keepInverse -= numIdentical;
			continue;
		}

		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS(row)) { 
				Rprintf("Row %d likelihood is %3.3f.\n", row, likelihood);
			} 
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, j), 0, likelihood);
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, likelihood);
			}
		} else {
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, likelihood);
			}
			double logDet = -2 * log(likelihood);
			logDet *= numIdentical;

			omxSetMatrixElement(rowLogLikelihoods, row, 0, logDet);

			if(OMX_DEBUG_ROWS(row)) { 
				Rprintf("-2 Log Likelihood this row is %3.3f, total change %3.3f\n",
				    logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
			}
		}
		
		if(firstRow) firstRow = 0;
		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);

		row += numIdentical;		// Step forward by the number of identical rows
		keepCov -= numIdentical;
		keepInverse -= numIdentical;
	}
}



/**
 * The localobj reference is used to access read-only variables,
 * or variables that can be modified but whose state cannot be
 * accessed from other threads.
 *
 * The sharedobj reference is used to access write-only variables,
 * where the memory writes of any two threads are non-overlapping.
 * No synchronization mechanisms are employed to maintain consistency
 * of sharedobj references.
 *
 * Because (1) these functions may be invoked with arbitrary 
 * rowbegin and rowcount values, and (2) the log-likelihood
 * values for all data rows must be calculated (even in cases
 * of errors), this function is forbidden from return()-ing early.
 *
 * As another consequence of (1) and (2), if "rowbegin" is in
 * the middle of a sequence of identical rows, then defer
 * move "rowbegin" to after the sequence of identical rows.
 * Grep for "[[Comment 4]]" in source code.
 * 
 */
void omxFIMLSingleIteration(omxObjective *localobj, omxObjective *sharedobj, int rowbegin, int rowcount) {
    
    omxFIMLObjective* ofo = ((omxFIMLObjective*) localobj->argStruct);
    omxFIMLObjective* shared_ofo = ((omxFIMLObjective*) sharedobj->argStruct);

	char u = 'U';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	int onei = 1;
	double determinant = 0.0;
	double Q = 0.0;
	double* oldDefs;
	int numDefs;
	int isContiguous, contiguousStart, contiguousLength;
	int numRemoves;
	int returnRowLikelihoods;
	int keepCov = 0, keepInverse = 0;

	omxObjective* subObjective;
	
	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *dataColumns;//, *oldInverse;
	omxMatrix *rowLikelihoods, *rowLogLikelihoods;
	omxDefinitionVar* defVars;
	omxData *data;

	// Locals, for readability.  Should compile out.
	cov 		     = ofo->cov;
	means		     = ofo->means;
	smallRow 	     = ofo->smallRow;
	smallCov 	     = ofo->smallCov;
	oldDefs		     = ofo->oldDefs;
	RCX 		     = ofo->RCX;
	data		     = ofo->data;                       //  read-only
	dataColumns	     = ofo->dataColumns;                //  read-only
	defVars		     = ofo->defVars;                    //  read-only
	numDefs		     = ofo->numDefs;                    //  read-only
	returnRowLikelihoods = ofo->returnRowLikelihoods;   //  read-only
	rowLikelihoods   = shared_ofo->rowLikelihoods;      // write-only
	rowLogLikelihoods = shared_ofo->rowLogLikelihoods;  // write-only
	isContiguous     = ofo->contiguous.isContiguous;    //  read-only
	contiguousStart  = ofo->contiguous.start;           //  read-only
	contiguousLength = ofo->contiguous.length;          //  read-only

	subObjective = localobj->subObjective;

	int toRemove[cov->cols];
	int zeros[cov->cols];

	int firstRow = 1;
	int row = rowbegin;

	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}

	resetDefinitionVariables(oldDefs, numDefs);

	while(row < data->rows && (row - rowbegin) < rowcount) {
        if (OMX_DEBUG_ROWS(row)) {Rprintf("Row %d.\n", row);} //:::DEBUG:::
		localobj->matrix->currentState->currentRow = row;		// Set to a new row.

		int numIdentical = omxDataNumIdenticalRows(data, row);
		// N.B.: numIdentical == 0 means an error occurred and was not properly handled;
		// it should never be the case.
		if (numIdentical == 0) numIdentical = 1; 
		
		Q = 0.0;

		numRemoves = 0;
		omxResetAliasedMatrix(smallRow); 			// Resize smallrow
		if (isContiguous) {
			omxContiguousDataRow(data, row, contiguousStart, contiguousLength, smallRow);
		} else {
			omxDataRow(data, row, dataColumns, smallRow);	// Populate data row
		}

		// Handle Definition Variables.
		if(numDefs != 0) {
			if(keepCov <= 0) {  // If we're keeping covariance from the previous row, do not populate
				if(OMX_DEBUG_ROWS(row)) { Rprintf("Handling Definition Vars.\n"); }
				if(handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs) || firstRow) {
				// Use firstrow instead of rows == 0 for the case where the first row is all NAs
				// N.B. handling of definition var lists always happens, regardless of firstRow.
					if(!(subObjective == NULL)) {
						omxObjectiveCompute(subObjective);
					} else {
						omxRecompute(cov);
						omxRecompute(means);
					}
				}
			} else if(OMX_DEBUG_ROWS(row)){ Rprintf("Identical def vars: Not repopulating"); }
		}

		if(OMX_DEBUG_ROWS(row)) { omxPrint(means, "Local Means"); }
		if(OMX_DEBUG_ROWS(row)) {
			char note[50];
			sprintf(note, "Local Data Row %d", row);
			omxPrint(smallRow, note); 
		}
		
		/* Censor row and censor and invert cov. matrix. */
		// Determine how many rows/cols to remove.
		memset(zeros, 0, sizeof(int) * dataColumns->cols);
		memset(toRemove, 0, sizeof(int) * dataColumns->cols);
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxVectorElement(smallRow, j);
			if(isnan(dataValue) || dataValue == NA_INTEGER || !R_FINITE(dataValue)) {
				numRemoves++;
				toRemove[j] = 1;
			} else if(means != NULL) {
				omxSetVectorElement(smallRow, j, (dataValue -  omxVectorElement(means, j)));
			}
		}

		if(cov->cols <= numRemoves) {
			for(int nid = 0; nid < numIdentical; nid++) {
				if(returnRowLikelihoods) {
					omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 1);
				}
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 1);
			}
			if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
    		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
            // Rprintf("Incrementing Row."); //:::DEBUG:::
    		row += numIdentical;
    		keepCov -= numIdentical;
    		keepInverse -= numIdentical;
			continue;
		}
		
		omxRemoveRowsAndColumns(smallRow, 0, numRemoves, zeros, toRemove); 	// Reduce it.
		
		if(OMX_DEBUG_ROWS(row)) { Rprintf("Keeper codes: inverse: %d, cov:%d, identical:%d\n", keepInverse, keepCov, omxDataNumIdenticalRows(data, row)); }

		if(keepInverse <= 0 || keepCov <= 0 || firstRow) { // If defs and missingness don't change, skip.
			omxResetAliasedMatrix(smallCov);				// Re-sample covariance matrix
			omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);

			if(OMX_DEBUG_ROWS(row)) { omxPrint(smallCov, "Local Covariance Matrix"); }

			/* Calculate derminant and inverse of Censored Cov matrix */
			// TODO : Speed this up.
			F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
			if(info != 0) {
				int skip;
				if(!returnRowLikelihoods) {
					char helperstr[200];
					char *errstr = calloc(250, sizeof(char));
					sprintf(helperstr, "Expected covariance matrix is not positive-definite in data row %d", omxDataIndex(data, row));
					if(localobj->matrix->currentState->computeCount <= 0) {
						sprintf(errstr, "%s at starting values.\n", helperstr);
					} else {
						sprintf(errstr, "%s at major iteration %d (minor iteration %d).\n", helperstr, 
							localobj->matrix->currentState->majorIteration, 
							localobj->matrix->currentState->minorIteration);
					}
					omxRaiseError(localobj->matrix->currentState, -1, errstr);
					free(errstr);
				}
				if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
				if(keepInverse <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
				
				skip = (keepCov < keepInverse) ? keepCov : keepInverse;

				for(int nid = 0; nid < skip; nid++) {
					if (returnRowLikelihoods) 
						omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
					omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
				}
				row += skip;
				keepCov -= skip;
				keepInverse -= skip;
				continue;
			}
			
			// Calculate determinant: squared product of the diagonal of the decomposition
			// For speed, we'll take the sum of the logs, rather than the log of the product
			determinant = 0.0;
			for(int diag = 0; diag < (smallCov->rows); diag++) {
				determinant += log(fabs(omxMatrixElement(smallCov, diag, diag)));
                // if(OMX_DEBUG_ROWS) { Rprintf("Next det is: %3.3d\n", determinant);} //:::DEBUG:::
			}
            // determinant = determinant * determinant; // Delayed for now.
			
			F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
			if(info != 0) {
				if(!returnRowLikelihoods) {
					char *errstr = calloc(250, sizeof(char));
					for(int nid = 0; nid < numIdentical; nid++) {
						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
					}
					sprintf(errstr, "Cannot invert expected covariance matrix. Error %d.", info);
					omxRaiseError(localobj->matrix->currentState, -1, errstr);
					free(errstr);
				} else {
					for(int nid = 0; nid < numIdentical; nid++) {
						omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
					}
					if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
					if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
					// Rprintf("Incrementing Row."); //:::DEBUG:::
					row += numIdentical;
					keepCov -= numIdentical;
					keepInverse -= numIdentical;
					continue;
				}
			}
		}

		/* Calculate Row Likelihood */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);

		double likelihood = pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/exp(determinant)) * exp(-.5 * Q);
		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS(row)) {Rprintf("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f\n", pow(2 * M_PI, -.5 * smallRow->cols), (1.0/exp(determinant)), exp(-.5 * Q), pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/exp(determinant)) * exp(-.5 * Q));}

			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, j), 0, likelihood);
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, likelihood);
			}
		} else {
			double logLikelihood = ((2.0*determinant) + Q + (log(2 * M_PI) * smallRow->cols)) * numIdentical;
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, likelihood);
			}
			omxSetMatrixElement(rowLogLikelihoods, row, 0, logLikelihood);

			if(OMX_DEBUG_ROWS(row)) {
				Rprintf("Change in Total Likelihood for row %d is %3.3f + %3.3f + %3.3f = %3.3f", localobj->matrix->currentState->currentRow, (2.0*determinant), Q, (log(2 * M_PI) * smallRow->cols), (2.0*determinant) + Q + (log(2 * M_PI) * smallRow->cols));
			}
		}
		if(firstRow) firstRow = 0;
		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);

		row += numIdentical;
		keepCov -= numIdentical;
		keepInverse -= numIdentical;
	}
}
