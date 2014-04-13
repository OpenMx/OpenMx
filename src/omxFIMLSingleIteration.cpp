/*
 *  Copyright 2007-2014 The OpenMx Project
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


#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"
#include "omxSadmvnWrapper.h"


void omxFIMLAdvanceRow(int *keepCov, int *keepInverse, int *row,
	omxData *data, int numIdentical) {
	int rowVal = *row;
	if(*keepCov <= 0) *keepCov = omxDataNumIdenticalDefs(data, rowVal);
	if(*keepInverse  <= 0) *keepInverse = omxDataNumIdenticalMissingness(data, rowVal);
	// mxLog("Incrementing Row."); //:::DEBUG:::
	*row += numIdentical;
	*keepCov -= numIdentical;
	*keepInverse -= numIdentical;
}

void omxFIMLAdvanceJointRow(int *row, int *numIdenticalDefs, 
	int *numIdenticalContinuousMissingness,
	int *numIdenticalOrdinalMissingness, 
	int *numIdenticalContinuousRows,
	int *numIdenticalOrdinalRows,
	omxData *data, int numDefs, int numIdentical) {

	int rowVal = *row;

    if(numDefs != 0 && *numIdenticalDefs <= 0) *numIdenticalDefs = 
		omxDataNumIdenticalDefs(data, rowVal);
	if(*numIdenticalContinuousMissingness <= 0) *numIdenticalContinuousMissingness =
		omxDataNumIdenticalContinuousMissingness(data, rowVal);
	if(*numIdenticalOrdinalMissingness <= 0) *numIdenticalOrdinalMissingness = 
		omxDataNumIdenticalOrdinalMissingness(data, rowVal);
	if(*numIdenticalContinuousRows <= 0) *numIdenticalContinuousRows = 
		omxDataNumIdenticalContinuousRows(data, rowVal);
	if(*numIdenticalOrdinalRows <= 0) *numIdenticalOrdinalRows = 
		omxDataNumIdenticalOrdinalRows(data, rowVal);

	*row += numIdentical;
	*numIdenticalDefs -= numIdentical;
	*numIdenticalContinuousMissingness -= numIdentical;
	*numIdenticalContinuousRows -= numIdentical;
	*numIdenticalOrdinalMissingness -= numIdentical;
	*numIdenticalOrdinalRows -= numIdentical;
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
 *
 * Because (1) these functions may be invoked with arbitrary 
 * rowbegin and rowcount values, and (2) the log-likelihood
 * values for all data rows must be calculated (even in cases
 * of Rf_errors), this function is forbidden from return()-ing early.
 *
 * As another consequence of (1) and (2), if "rowbegin" is in
 * the middle of a sequence of identical rows, then defer
 * move "rowbegin" to after the sequence of identical rows.
 * Grep for "[[Comment 4]]" in source code.
 */
bool omxFIMLSingleIterationJoint(FitContext *fc, omxFitFunction *localobj, omxFitFunction *sharedobj, int rowbegin, int rowcount) {

    omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) localobj->argStruct);
    omxFIMLFitFunction* shared_ofo = ((omxFIMLFitFunction*) sharedobj->argStruct);

	double Q = 0.0;
	int numDefs;
	int numOrdRemoves = 0, numContRemoves=0;
	int returnRowLikelihoods = 0;
    int numIdenticalDefs = 0, numIdenticalOrdinalMissingness = 0, numIdenticalOrdinalRows = 0, numOrdinal = 1,
	                    numIdenticalContinuousMissingness = 0, numIdenticalContinuousRows = 0, numContinuous = 1;

	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns;
	omxMatrix *rowLikelihoods, *rowLogLikelihoods;
    omxMatrix *ordMeans, *ordCov, *ordRow, *contRow;
    omxMatrix *halfCov, *reduceCov, *ordContCov;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *lThresh, *uThresh, *corList, *weights, *oldDefs;
	int *Infin;
	omxDefinitionVar* defVars;
	
	omxExpectation* expectation;
	

	// Locals, for readability.  Compiler should cut through this.
	cov 		= ofo->cov;
	means		= ofo->means;
	smallRow 	= ofo->smallRow;
	smallCov 	= ofo->smallCov;
	smallMeans	= ofo->smallMeans;
    ordMeans    = ofo->ordMeans;
    ordCov      = ofo->ordCov;
    ordRow      = ofo->ordRow;
    contRow     = ofo->contRow;
    halfCov     = ofo->halfCov;
    reduceCov   = ofo->reduceCov;
    ordContCov  = ofo->ordContCov;
	RCX 		= ofo->RCX;

	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	defVars		= ofo->defVars;
	oldDefs		= ofo->oldDefs;
	numDefs		= ofo->numDefs;

	corList 	= ofo->corList;
	weights		= ofo->weights;
	lThresh		= ofo->lThresh;
	uThresh		= ofo->uThresh;
	thresholdCols = ofo->thresholdCols;
	returnRowLikelihoods = ofo->returnRowLikelihoods;
	rowLikelihoods = shared_ofo->rowLikelihoods;		// write-only
	rowLogLikelihoods = shared_ofo->rowLogLikelihoods;  // write-only

	Infin			= ofo->Infin;
	expectation 	= localobj->expectation;

    int ordRemove[cov->cols], contRemove[cov->cols];
    char u = 'U', l = 'L';
    int info;
    double determinant = 0.0;
    double oned = 1.0, zerod = 0.0, minusoned = -1.0;
    int onei = 1;
    double likelihood;
	int inform;

	int firstRow = 1;
    int row = rowbegin;

    resetDefinitionVariables(oldDefs, numDefs);

	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}

	while(row < data->rows && (row - rowbegin) < rowcount) {
        localobj->matrix->currentState->currentRow = row;		// Set to a new row.
        int numIdentical = omxDataNumIdenticalRows(data, row);
        if(numIdentical == 0) numIdentical = 1; 
        // N.B.: numIdentical == 0 means an Rf_error occurred and was not properly handled;
        // it should never be the case.
        
        omxDataRow(data, row, dataColumns, smallRow);                               // Populate data row
        
        if(OMX_DEBUG_ROWS(row)) {
            mxLog("Identicality check. Is %sfirst. Total: %d rows identical, %d defs, %d missingness: Continuous: %d rows, %d missingness; Ordinal: %d rows, %d missingness.", 
                    (firstRow?"":"not "), numIdentical, numIdenticalDefs, omxDataNumIdenticalRows(data, row), 
                    numIdenticalContinuousRows, numIdenticalContinuousMissingness, 
                    numIdenticalOrdinalRows, numIdenticalOrdinalMissingness);
        }

        if(numIdenticalDefs <= 0 || numIdenticalContinuousMissingness <= 0 || numIdenticalOrdinalMissingness <= 0 || firstRow ) {  // If we're keeping covariance from the previous row, do not populate 
            // Handle Definition Variables.
            if((numDefs != 0 && numIdenticalDefs > 0) || firstRow) {
				int numVarsFilled = 0;
				if(OMX_DEBUG_ROWS(row)) { mxLog("Handling Definition Vars."); }
				numVarsFilled = handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs);
				if (numVarsFilled < 0) {
					for(int nid = 0; nid < numIdentical; nid++) {
						if(returnRowLikelihoods) omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
					}
					omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
						&numIdenticalContinuousMissingness,
						&numIdenticalOrdinalMissingness, 
						&numIdenticalContinuousRows,
						&numIdenticalOrdinalRows,
						data, numDefs, numIdentical);
					continue;
				} else if (numVarsFilled || firstRow) { 
					// Use firstrow instead of rows == 0 for the case where the first row is all NAs
					// N.B. handling of definition var lists always happens, regardless of firstRow.
					// Recalculate means and covariances.
					omxExpectationCompute(expectation, NULL);
				}
			}
            // Filter down correlation matrix and calculate thresholds.
            // TODO: If identical ordinal or continuous missingness, ignore only the appropriate columns.
            numOrdRemoves = 0;
            numContRemoves = 0;
			for(int j = 0; j < dataColumns->cols; j++) {
    			int var = omxVectorElement(dataColumns, j);
    			int value = omxIntDataElement(data, row, var);// Indexing correction means this is the index of the upper bound +1.
    			// TODO: Might save time by preseparating ordinal from continuous.
    			if(std::isnan(value) || value == NA_INTEGER) {  // Value is NA, therefore filter.
    				numOrdRemoves++;
                    numContRemoves++;
                    ordRemove[j] = 1;
                    contRemove[j] = 1;
    				Infin[j] = -1;
                    if(OMX_DEBUG_ROWS(row)) { 
    				    mxLog("Row %d, column %d, value %d.  NA.", row, j, value);
                    }
    				continue;
    			} else if(omxDataColumnIsFactor(data, var)) {             // Ordinal column.
                    numContRemoves++;
                    ordRemove[j] = 0;
                    contRemove[j] = 1;
                    if(OMX_DEBUG_ROWS(row)) { 
    			        mxLog("Row %d, column %d, value %d.  Ordinal.", row, j, value);
                    }
    			} else {
    			    numOrdRemoves++;
                    ordRemove[j] = 1;
                    contRemove[j] = 0;
                    if(OMX_DEBUG_ROWS(row)) { 
    			        mxLog("Row %d, column %d, value %d.  Continuous.", row, j, value);
                    }
    			}
    		}
    		
            if(OMX_DEBUG_ROWS(row)) {
                mxLog("\n\nRemovals: %d ordinal, %d continuous out of %d total.", numOrdRemoves, numContRemoves, dataColumns->cols);
            }
    		
			for(int j=0; j < dataColumns->cols; j++) {
				int var = omxVectorElement(dataColumns, j);
				if(omxDataColumnIsFactor(data, j) && thresholdCols[var].numThresholds > 0) { // j is an ordinal column
					omxRecompute(thresholdCols[var].matrix); // Only one of these--save time by only doing this once
					checkIncreasing(thresholdCols[var].matrix, thresholdCols[var].column);
				}
			}
            numContinuous = dataColumns->cols - numContRemoves;
            numOrdinal = dataColumns->cols - numOrdRemoves;

		}

        // TODO: Possible solution here: Manually record threshold column and index from data 
        //   during this initial reduction step.  Since all the rest is algebras, it'll filter 
        //   naturally.  Calculate offsets from continuous data, then dereference actual 
        //   threshold values from the threshold matrix in its original state.  
        //   Alternately, rearrange the thresholds matrix (and maybe data matrix) to split
        //    ordinal and continuous variables.
        //   Requirement: colNum integer vector

		if(numContinuous <= 0 && numOrdinal <= 0) {
		    // All elements missing.  Skip row.
			for(int nid = 0; nid < numIdentical; nid++) {	
				if(returnRowLikelihoods) {
					omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 1.0);
				}
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 1.0);
			}
			omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
				&numIdenticalContinuousMissingness,
				&numIdenticalOrdinalMissingness, 
				&numIdenticalContinuousRows,
				&numIdenticalOrdinalRows,
				data, numDefs, numIdentical);
            continue;
		}

		//  smallCov <- cov[!contRemove, !contRemove] : covariance of continuous elements
		//  smallMeans <- means[ALL, !contRemove] : continuous means
		//  smallRow <- data[ALL, !contRemove]  : continuous data
		//              ordCov <- cov[!ordRemove, !ordRemove]
		//              ordMeans <- means[NULL, !ordRemove]
		//              ordData <- data[NULL, !ordRemove]
		//              ordContCov <- cov[!contRemove, !ordRemove]

        // TODO: Data handling is confusing.  Maybe set two self-aliased row-reduction "datacolumns" elements?
        
        // SEPARATION: 
        // Catch here: If continuous columns are all missing, skip everything except the ordCov calculations
        //              in this case, log likelihood of the continuous is 1 (likelihood is 0)
        // Do not recompute ordcov if missingness is identical and no def vars

        // SEPARATION: 
        //  Unprojected covariances only need to reset and re-filter if there are def vars or the appropriate missingness pattern changes
        //  Also, if each one is not all-missing.

		if(numContinuous <= 0) {
		    // All continuous missingness.  Populate some stuff.
            Q = 0.0;
            determinant = 0.0;
            if(numIdenticalDefs <= 0 || numIdenticalOrdinalRows <= 0 || firstRow) {
                // Recalculate Ordinal covariance matrix
                omxResetAliasedMatrix(ordCov);				// Re-sample covariance and means matrices for ordinal
                omxRemoveRowsAndColumns(ordCov, numOrdRemoves, numOrdRemoves, ordRemove, ordRemove);
                
                // Recalculate ordinal fs
                omxResetAliasedMatrix(ordMeans);
                omxRemoveElements(ordMeans, numOrdRemoves, ordRemove); 	    // Reduce the row to just ordinal.
                
                // These values pass through directly without modification by continuous variables
                
                // Calculate correlation matrix, correlation list, and weights from covariance
    		    omxStandardizeCovMatrix(ordCov, corList, weights);
            }
        } else if(numIdenticalDefs <= 0 || numIdenticalContinuousRows <= 0 || firstRow) {

            /* Reset and Resample rows if necessary. */
            // First Cov and Means (if they've changed)
            if(numIdenticalDefs <= 0 || numIdenticalContinuousMissingness <= 0 || firstRow) {
                omxResetAliasedMatrix(smallMeans);
                omxRemoveElements(smallMeans, numContRemoves, contRemove);
                omxResetAliasedMatrix(smallCov);				// Re-sample covariance matrix
                omxRemoveRowsAndColumns(smallCov, numContRemoves, numContRemoves, contRemove, contRemove);

                /* Calculate derminant and inverse of Censored continuousCov matrix */
                if(OMX_DEBUG_ROWS(row)) { 
                    omxPrint(smallCov, "Cont Cov to Invert"); 
                }
                
			    F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);

                if(info != 0) {
                    if(!returnRowLikelihoods) {
                        for(int nid = 0; nid < numIdentical; nid++) {
                            omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
                        }
			if (fc) fc->recordIterationError("Expected covariance matrix for continuous variables "
							 "is not positive-definite in data row %d", omxDataIndex(data, row));
			return TRUE;
                    }
                    for(int nid = 0; nid < numIdentical; nid++) {
                        if (returnRowLikelihoods)
					        omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
                        omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
                    }
                    if(OMX_DEBUG) {mxLog("Non-positive-definite covariance matrix in row likelihood.  Skipping Row.");}
   					omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
						&numIdenticalContinuousMissingness,
						&numIdenticalOrdinalMissingness, 
						&numIdenticalContinuousRows,
						&numIdenticalOrdinalRows,
						data, numDefs, numIdentical);
                    continue;
                }
                // Calculate determinant: squared product of the diagonal of the decomposition
    			// For speed, use sum of logs rather than log of product.
    			
                determinant = 0.0;
    			for(int diag = 0; diag < (smallCov->rows); diag++) {
    				determinant += log(fabs(omxMatrixElement(smallCov, diag, diag)));
    			}
                // determinant = determinant * determinant;  // Delayed.
    			F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
    			if(info != 0) {
    				if(!returnRowLikelihoods) {
    					char *errstr = (char*) calloc(250, sizeof(char));
    					sprintf(errstr, "Cannot invert expected continuous covariance matrix. Error %d.", info);
    					omxRaiseError(localobj->matrix->currentState, -1, errstr);
    					free(errstr);
                    }
   					for(int nid = 0; nid < numIdentical; nid++) {
                        if (returnRowLikelihoods)
   						    omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
   						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
   					}
 					omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
						&numIdenticalContinuousMissingness,
						&numIdenticalOrdinalMissingness, 
						&numIdenticalContinuousRows,
						&numIdenticalOrdinalRows,
						data, numDefs, numIdentical);
                    continue;
    			}
            }

            // Reset continuous data row (always needed)
            omxResetAliasedMatrix(contRow);                                            // Reset smallRow
            omxRemoveElements(contRow, numContRemoves, contRemove); 	// Reduce the row to just continuous.
            F77_CALL(daxpy)(&(contRow->cols), &minusoned, smallMeans->data, &onei, contRow->data, &onei);

		    /* Calculate Row Likelihood */
		    /* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
    		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), contRow->data, &onei, &zerod, RCX->data, &onei);       // RCX is the continuous-column mahalanobis distance.
    		Q = F77_CALL(ddot)(&(contRow->cols), contRow->data, &onei, RCX->data, &onei); //Q is the total mahalanobis distance
    		
    		if(numOrdinal > 0) {

                // Precalculate Ordinal things that change with continuous changes
                // Reserve: 1) Inverse continuous covariance (smallCov)
                //          2) Columnwise Mahalanobis distance (contCov^-1)%*%(Data - Means) (RCX)
                //          3) Overall Mahalanobis distance (FIML likelihood of data) (Q)
                //Calculate:4) Cont/ord covariance %*% Mahalanobis distance  (halfCov)
                //          5) ordCov <- ordCov - Cont/ord covariance %*% Inverse continuous cov

                if(numIdenticalContinuousMissingness <= 0 || firstRow) {
                    // Re-sample covariance between ordinal and continuous only if the continuous missingness changes.
                    omxResetAliasedMatrix(ordContCov);
                    omxRemoveRowsAndColumns(ordContCov, numContRemoves, numOrdRemoves, contRemove, ordRemove);

                    // TODO: Make this less of a hack.
                    halfCov->rows = smallCov->rows;
                    halfCov->cols = ordContCov->cols;
                    omxMatrixLeadingLagging(halfCov);
                    reduceCov->rows = ordContCov->cols;
                    reduceCov->cols = ordContCov->cols;
                    omxMatrixLeadingLagging(reduceCov);

                    F77_CALL(dsymm)(&l, &u, &(smallCov->rows), &(ordContCov->cols), &oned, smallCov->data, &(smallCov->leading), ordContCov->data, &(ordContCov->leading), &zerod, halfCov->data, &(halfCov->leading));          // halfCov is inverse continuous %*% cont/ord covariance
                    F77_CALL(dgemm)((ordContCov->minority), (halfCov->majority), &(ordContCov->cols), &(halfCov->cols), &(ordContCov->rows), &oned, ordContCov->data, &(ordContCov->leading), halfCov->data, &(halfCov->leading), &zerod, reduceCov->data, &(reduceCov->leading));      // reduceCov is cont/ord^T %*% (contCov^-1 %*% cont/ord)
                }

                if(numIdenticalOrdinalMissingness <= 0 || firstRow) {
                    // Means, projected covariance, and Columnwise mahalanobis distance must be recalculated
                    //   unless there are no ordinal variables or the continuous variables are identical
                    
                    // Recalculate Ordinal and Ordinal/Continuous covariance matrices.
                    if(OMX_DEBUG_ROWS(row)) {
                        mxLog("Resetting Ordinal Covariance Matrix.");
                        omxPrint(ordCov, "Was:");
                    }

                    omxResetAliasedMatrix(ordCov);				// Re-sample covariance and means matrices for ordinal
                    if(OMX_DEBUG_ROWS(row)) {
                        mxLog("Resetting/Filtering Ordinal Covariance Matrix.");
                        omxPrint(ordCov, "Reset to:");
                    }

                    omxRemoveRowsAndColumns(ordCov, numOrdRemoves, numOrdRemoves, ordRemove, ordRemove);
                    if(OMX_DEBUG_ROWS(row)) {
                        mxLog("Resetting/Filtering Ordinal Covariance Matrix.");
                        omxPrint(ordCov, "Filtered to:");
                    }
                
                    // FIXME: This assumes that ordCov and reducCov have the same row/column majority.
                    int vlen = reduceCov->rows * reduceCov->cols;
                    F77_CALL(daxpy)(&vlen, &minusoned, reduceCov->data, &onei, ordCov->data, &onei); // ordCov <- (ordCov - reduceCov) %*% cont/ord
            
                }
                
                // Projected means must be recalculated if the continuous variables change at all.
                omxResetAliasedMatrix(ordMeans);
                omxRemoveElements(ordMeans, numOrdRemoves, ordRemove); 	    // Reduce the row to just ordinal.
                F77_CALL(dgemv)((smallCov->minority), &(halfCov->rows), &(halfCov->cols), &oned, halfCov->data, &(halfCov->leading), contRow->data, &onei, &oned, ordMeans->data, &onei);                      // ordMeans += halfCov %*% contRow
            }
            
        } // End of continuous likelihood values calculation
        
        if(numOrdinal <= 0) {       // No Ordinal Vars at all.
            likelihood = 1;
        } else {  
            // There are ordinal vars, and not everything is identical, so we're recalculating
            // Calculate correlation matrix, correlation list, and weights from covariance
            if(numIdenticalDefs <=0 || numIdenticalContinuousMissingness <= 0 || numIdenticalOrdinalMissingness <= 0 || firstRow) {
                // if(OMX_DEBUG_ROWS(row)) {omxPrint(ordCov, "Ordinal cov matrix for standardization."); } //:::DEBUG:::
                omxStandardizeCovMatrix(ordCov, corList, weights);
            }
            
            omxResetAliasedMatrix(ordRow);                                              // Propagate to ordinal row
            omxRemoveElements(ordRow, numOrdRemoves, ordRemove); 	    // Reduce the row to just ordinal.

            // omxPrint(ordMeans, "Ordinal Projected Means"); //:::DEBUG:::
            // omxPrint(ordRow, "Filtered Ordinal Row"); //:::DEBUG:::


            // Inspect elements, reweight, and set to 
            int count = 0;
    		for(int j = 0; j < dataColumns->cols; j++) {
                if(ordRemove[j]) continue;         // NA or non-ordinal
                int var = omxVectorElement(dataColumns, j);
    			int value = omxIntDataElement(data, row, var); //  TODO: Compare with extraction from dataRow.
                // mxLog("Row %d, Column %d, value %d+1\n", row, j, value); // :::DEBUG:::
    	        value--;		// Correct for C indexing: value is now the index of the upper bound.
                // mxLog("Row %d, Column %d, value %d+1\n", row, j, value); // :::DEBUG:::
    			double offset;
    			if(means == NULL) offset = 0;
    			else offset = omxVectorElement(ordMeans, count);
    			double weight = weights[count];
    			if(value == 0) { 									// Lowest threshold = -Inf
    				lThresh[count] = (omxMatrixElement(thresholdCols[j].matrix, 0, thresholdCols[j].column) - offset) / weight;
    				uThresh[count] = lThresh[count];
    				Infin[count] = 0;
    			} else {
    				lThresh[count] = (omxMatrixElement(thresholdCols[j].matrix, value-1, thresholdCols[j].column) - offset) / weight;
    				if(thresholdCols[j].numThresholds > value) {	// Highest threshold = Inf
    					double tmp = (omxMatrixElement(thresholdCols[j].matrix, value, thresholdCols[j].column) - offset) / weight;
    					uThresh[count] = tmp;
    					Infin[count] = 2;
    				} else {
    					uThresh[count] = NA_INTEGER; // NA is a special to indicate +Inf
    					Infin[count] = 1;
    				}
    			}

    			if(uThresh[count] == NA_INTEGER || std::isnan(uThresh[count])) { // for matrix-style specification.
    				uThresh[count] = lThresh[count];
    				Infin[count] = 1;
    			}
    			if(OMX_DEBUG) { 
    			    mxLog("Row %d, column %d.  Thresholds for data column %d and threshold column %d are %f -> %f. (Infin=%d).  Offset is %f and weight is %f",
    			            row, count, j, value, lThresh[count], uThresh[count], Infin[count], offset, weight);
    			    mxLog("       Thresholds were %f -> %f, scaled by weight %f and shifted by mean %f and total offset %f.",
    			            omxMatrixElement(thresholdCols[j].matrix, (Infin[count]==0?0:value-1), thresholdCols[j].column), 
    			            omxMatrixElement(thresholdCols[j].matrix, (Infin[count]==1?value-1:value), thresholdCols[j].column), 
                            weight, (means==NULL?0:omxVectorElement(ordMeans, count)), offset);
    			}
                count++;
    		}

			omxSadmvnWrapper(localobj, cov, ordCov, corList, lThresh, uThresh, Infin, &likelihood, &inform);

    		if(inform == 2) {
    			if(!returnRowLikelihoods) {
				if (fc) fc->recordIterationError("Improper value detected by integration routine in data row %d: Most likely the maximum number of ordinal variables (20) has been exceeded.  \n Also check that expected covariance matrix is not positive-definite", omxDataIndex(data, row));
				return TRUE;
			}
  				for(int nid = 0; nid < numIdentical; nid++) {
                    if (returnRowLikelihoods)
						omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
  					omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
   				}
                if(OMX_DEBUG) {mxLog("Improper input to sadmvn in row likelihood.  Skipping Row.");}
				omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
					&numIdenticalContinuousMissingness,
					&numIdenticalOrdinalMissingness, 
					&numIdenticalContinuousRows,
					&numIdenticalOrdinalRows,
					data, numDefs, numIdentical);
                continue;
    		}
		}

		double rowLikelihood = pow(2 * M_PI, -.5 * numContinuous) * (1.0/exp(determinant)) * exp(-.5 * Q) * likelihood;

		if(returnRowLikelihoods) {
		    if(OMX_DEBUG_ROWS(row)) {mxLog("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f", 
		        pow(2 * M_PI, -.5 * numContinuous), (1.0/exp(determinant)), exp(-.5 * Q), 
		        pow(2 * M_PI, -.5 * numContinuous) * (1.0/exp(determinant)) * exp(-.5 * Q));}
            
			if(OMX_DEBUG_ROWS(row)) {mxLog("Row %d likelihood is %3.3f.", row, rowLikelihood);}
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, j), 0, rowLikelihood);
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, rowLikelihood);
			}
		} else {
			double logLikelihood = -2 * log(likelihood);       // -2 Log of ordinal likelihood
            logLikelihood += ((2 * determinant) + Q + (log(2 * M_PI) * numContinuous));    // -2 Log of continuous likelihood
            logLikelihood *= numIdentical;

			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, rowLikelihood);
			}
			omxSetMatrixElement(rowLogLikelihoods, row, 0, logLikelihood);
			
			if(OMX_DEBUG_ROWS(row)) { 
				mxLog("Change in Total log Likelihood for row %ld is %3.3f + %3.3f + %3.3f + %3.3f= %3.3f, ", 
				    localobj->matrix->currentState->currentRow, (2.0*determinant), Q, (log(2 * M_PI) * numContinuous), 
				    -2  * log(rowLikelihood), (2.0 *determinant) + Q + (log(2 * M_PI) * numContinuous));
			} 

        }
        if(firstRow) firstRow = 0;
		omxFIMLAdvanceJointRow(&row, &numIdenticalDefs, 
			&numIdenticalContinuousMissingness,
			&numIdenticalOrdinalMissingness, 
			&numIdenticalContinuousRows,
			&numIdenticalOrdinalRows,
			data, numDefs, numIdentical);
        continue;

	}
	return FALSE;
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
 *
 * Because (1) these functions may be invoked with arbitrary 
 * rowbegin and rowcount values, and (2) the log-likelihood
 * values for all data rows must be calculated (even in cases
 * of Rf_errors), this function is forbidden from return()-ing early.
 *
 * As another consequence of (1) and (2), if "rowbegin" is in
 * the middle of a sequence of identical rows, then defer
 * move "rowbegin" to after the sequence of identical rows.
 * Grep for "[[Comment 4]]" in source code.
 */
bool omxFIMLSingleIterationOrdinal(FitContext *fc, omxFitFunction *localobj, omxFitFunction *sharedobj, int rowbegin, int rowcount) {

    omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) localobj->argStruct);
    omxFIMLFitFunction* shared_ofo = ((omxFIMLFitFunction*) sharedobj->argStruct);

	double Q = 0.0;
	double* oldDefs;
	int numDefs;
	int numRemoves;
	int returnRowLikelihoods;
	int keepCov = 0, keepInverse = 0;

	omxExpectation* expectation;
	
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

	expectation 	 = localobj->expectation;

	int firstRow = 1;
    int row = rowbegin;

    resetDefinitionVariables(oldDefs, numDefs);

	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}

	while(row < data->rows && (row - rowbegin) < rowcount) {
        if(OMX_DEBUG_ROWS(row)) {mxLog("Row %d.", row);}
        localobj->matrix->currentState->currentRow = row;		// Set to a new row.
		int numIdentical = omxDataNumIdenticalRows(data, row);
		if(numIdentical == 0) numIdentical = 1; 
		// N.B.: numIdentical == 0 means an Rf_error occurred and was not properly handled;
		// it should never be the case.

        Q = 0.0;

        // Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
        numRemoves = 0;

        // Handle Definition Variables.
        if(numDefs != 0) {
			if(keepCov <= 0) {  // If we're keeping covariance from the previous row, do not populate 
				int numVarsFilled = 0;
				if(OMX_DEBUG_ROWS(row)) { mxLog("Handling Definition Vars."); }
				numVarsFilled = handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs);
				if (numVarsFilled < 0) {
					for(int nid = 0; nid < numIdentical; nid++) {
						if(returnRowLikelihoods) omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
					}
					omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
					continue;
				} else if (numVarsFilled || firstRow) {
					// Use firstrow instead of rows == 0 for the case where the first row is all NAs
					// N.B. handling of definition var lists always happens, regardless of firstRow.
					omxExpectationCompute(expectation, NULL);
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
					mxLog("Row %d, column %d. Mean is %f and weight is %f", row, j, mean, weight);
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
				
				if(uThresh[j] == NA_INTEGER || std::isnan(uThresh[j])) { // for matrix-style specification.
					uThresh[j] = lThresh[j];
					Infin[j] = 1;
				}

				if(OMX_DEBUG_ROWS(row)) { 
					mxLog("Row %d, column %d.  Thresholds for data column %d and row %d are %f -> %f. (Infin=%d)", 
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
			omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
			continue;
		}

   		double likelihood;
		int inform;

		omxSadmvnWrapper(localobj, cov, smallCov, corList, lThresh, uThresh, Infin, &likelihood, &inform);

		if(inform == 2) {
			if(!returnRowLikelihoods) {
				if (fc) fc->recordIterationError("Improper value detected by integration "
							 "routine in data row %d: \n Most likely the maximum "
							 "number of ordinal variables (20) has been exceeded.\n "
							 "Also check that the expected covariance matrix is positive-definite",
							 omxDataIndex(data, row));
				return TRUE;
			}
			for(int nid = 0; nid < numIdentical; nid++) {
				if (returnRowLikelihoods)
					omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
			}
			omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
			continue;
		}

		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS(row)) { 
				mxLog("Row %d likelihood is %3.3f.", row, likelihood);
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
				mxLog("-2 Log Likelihood this row is %3.3f, total change %3.3f",
				    logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
			}
		}
		
		if(firstRow) firstRow = 0;
		omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
	}
	return FALSE;
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
 * of Rf_errors), this function is forbidden from return()-ing early.
 *
 * As another consequence of (1) and (2), if "rowbegin" is in
 * the middle of a sequence of identical rows, then defer
 * move "rowbegin" to after the sequence of identical rows.
 * Grep for "[[Comment 4]]" in source code.
 * 
 */
bool omxFIMLSingleIteration(FitContext *fc, omxFitFunction *localobj, omxFitFunction *sharedobj, int rowbegin, int rowcount) {
    if(OMX_DEBUG_ALGEBRA) {mxLog("Entering FIML Single Iteration"); }
    omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) localobj->argStruct);
    omxFIMLFitFunction* shared_ofo = ((omxFIMLFitFunction*) sharedobj->argStruct);

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

	omxExpectation* expectation;
	
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

	expectation = localobj->expectation;

	int toRemove[cov->cols];
	int dataColumnCols = dataColumns->cols;

	int firstRow = 1;
	int row = rowbegin;
	

	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}
	
	if(row == 0 && !strcmp(expectation->expType, "MxExpectationStateSpace") ) {
		if(OMX_DEBUG){ mxLog("Resetting State Space state (x) and Rf_error cov (P)."); }
		omxSetExpectationComponent(expectation, localobj, "Reset", NULL);
	}

	resetDefinitionVariables(oldDefs, numDefs);

	while(row < data->rows && (row - rowbegin) < rowcount) {
        if (OMX_DEBUG_ROWS(row)) {mxLog("Row %d.", row);} //:::DEBUG:::
		localobj->matrix->currentState->currentRow = row;		// Set to a new row.

		int numIdentical = omxDataNumIdenticalRows(data, row);
		// N.B.: numIdentical == 0 means an Rf_error occurred and was not properly handled;
		// it should never be the case.
		if (numIdentical == 0) numIdentical = 1; 
		
		Q = 0.0;
		
		numRemoves = 0;
		omxResizeMatrix(smallRow, 1, cov->cols, FALSE);
		if (isContiguous) {
			omxContiguousDataRow(data, row, contiguousStart, contiguousLength, smallRow);
		} else {
			omxDataRow(data, row, dataColumns, smallRow);	// Populate data row
		}
		if(OMX_DEBUG_ROWS(row)){omxPrint(smallRow, "...smallRow"); }
		if(!strcmp(expectation->expType, "MxExpectationStateSpace")) {
			omxSetExpectationComponent(expectation, localobj, "y", smallRow);
		}
		//If the expectation is a state space model then
		// set the y attribute of the state space expectation to smallRow.
		
		// Handle Definition Variables.
		if(numDefs != 0 || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
			if(keepCov <= 0 || !strcmp(expectation->expType, "MxExpectationStateSpace")) {  // If we're keeping covariance from the previous row, do not populate
				int numVarsFilled = 0;
				if(OMX_DEBUG_ROWS(row)) { mxLog("Handling Definition Vars."); }
				numVarsFilled = handleDefinitionVarList(data, localobj->matrix->currentState, row, defVars, oldDefs, numDefs);
				if (numVarsFilled < 0) {
					for(int nid = 0; nid < numIdentical; nid++) {
						if(returnRowLikelihoods) omxSetMatrixElement(sharedobj->matrix, omxDataIndex(data, row+nid), 0, 0.0);
						omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, row+nid), 0, 0.0);
					}
					if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
					if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
					// mxLog("Incrementing Row."); //:::DEBUG:::
					row += numIdentical;
					keepCov -= numIdentical;
					keepInverse -= numIdentical;
					continue;
				} else if (numVarsFilled || firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
				// Use firstrow instead of rows == 0 for the case where the first row is all NAs
				// N.B. handling of definition var lists always happens, regardless of firstRow.
					omxExpectationCompute(expectation, NULL);
					if(OMX_DEBUG_ROWS(row)){ mxLog("Successfully finished expectation compute."); }
				}
			} else if(OMX_DEBUG_ROWS(row)){ mxLog("Identical def vars: Not repopulating"); }
		}

		if(OMX_DEBUG_ROWS(row)) { omxPrint(means, "Local Means"); }
		if(OMX_DEBUG_ROWS(row)) {
			char note[50];
			sprintf(note, "Local Data Row %d", row);
			omxPrint(smallRow, note); 
		}
		
		/* Censor row and censor and invert cov. matrix. */
		// Determine how many rows/cols to remove.
		memset(toRemove, 0, sizeof(int) * dataColumnCols);
		for(int j = 0; j < dataColumnCols; j++) {
			double dataValue = omxVectorElement(smallRow, j);
			int dataValuefpclass = std::fpclassify(dataValue);
			if(dataValuefpclass == FP_NAN || dataValuefpclass == FP_INFINITE) {
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
			omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
			continue;
		}
		
		omxRemoveElements(smallRow, numRemoves, toRemove); 	// Reduce it.
		
		if(OMX_DEBUG_ROWS(row)) { mxLog("Keeper codes: inverse: %d, cov:%d, identical:%d", keepInverse, keepCov, omxDataNumIdenticalRows(data, row)); }

		if((keepInverse <= 0 || keepCov <= 0 || firstRow) && strcmp(expectation->expType, "MxExpectationStateSpace")) { // If defs and missingness don't change, skip.
			// also skip if this is a state space expectation
			if(OMX_DEBUG_ROWS(row)) { mxLog("Beginning to recompute inverse cov for standard models"); }
			omxResetAliasedMatrix(smallCov);				// Re-sample covariance matrix
			omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);

			if(OMX_DEBUG_ROWS(row)) { omxPrint(smallCov, "Local Covariance Matrix"); }

			/* Calculate derminant and inverse of Censored Cov matrix */
			// TODO : Speed this up.
			F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
			if(info != 0) {
				int skip;
				if(!returnRowLikelihoods) {
					if (fc) fc->recordIterationError("Expected covariance matrix is not "
								 "positive-definite in data row %d", omxDataIndex(data, row));
					return TRUE;
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
                // if(OMX_DEBUG_ROWS) { mxLog("Next det is: %3.3d", determinant);} //:::DEBUG:::
			}
            // determinant = determinant * determinant; // Delayed for now.
			
			F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
			if(info != 0) {
				if(!returnRowLikelihoods) {
					char *errstr = (char*) calloc(250, sizeof(char));
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
					omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
					continue;
				}
			}
		}
		
		/* If it's a state space expectation, extract the inverse rather than recompute it */
		if(!strcmp(expectation->expType, "MxExpectationStateSpace")) {
			if(OMX_DEBUG_ROWS(row)) { mxLog("Beginning to extract inverse cov for state space models"); }
			smallCov = omxGetExpectationComponent(expectation, localobj, "inverse");
			if(OMX_DEBUG_ROWS(row)) { omxPrint(smallCov, "Inverse of Local Covariance Matrix in state space model"); }
			//Get covInfo from state space expectation
			info = (int) omxGetExpectationComponent(expectation, localobj, "covInfo")->data[0];
			if(info!=0) {
				if (fc) fc->recordIterationError("Expected covariance matrix is not positive-definite in data row %d", omxDataIndex(data, row));
				return TRUE;
			}
			
			determinant = *omxGetExpectationComponent(expectation, localobj, "determinant")->data;
			if(OMX_DEBUG_ROWS(row)) { mxLog("0.5*log(det(Cov)) is: %3.3f", determinant);}
			
		}

		/* Calculate Row Likelihood */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);

		double likelihood = pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/exp(determinant)) * exp(-.5 * Q);
		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS(row)) {mxLog("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f", pow(2 * M_PI, -.5 * smallRow->cols), (1.0/exp(determinant)), exp(-.5 * Q), pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/exp(determinant)) * exp(-.5 * Q));}

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
				mxLog("Change in Total Likelihood for row %ld is %3.3f + %3.3f + %3.3f = %3.3f", localobj->matrix->currentState->currentRow, (2.0*determinant), Q, (log(2 * M_PI) * smallRow->cols), (2.0*determinant) + Q + (log(2 * M_PI) * smallRow->cols));
			}
		}
		if(firstRow) firstRow = 0;
		omxFIMLAdvanceRow(&keepCov, &keepInverse, &row, data, numIdentical);
	}
	return FALSE;
}
