/*
 *  Copyright 2007-2015 The OpenMx Project
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
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"
#include "omxSadmvnWrapper.h"


void omxFIMLAdvanceJointRow(int *row, int *numIdenticalDefs, 
	int *numIdenticalContinuousMissingness,
	int *numIdenticalOrdinalMissingness, 
	int *numIdenticalContinuousRows,
	int *numIdenticalOrdinalRows,
	omxData *data, int numDefs, int numIdentical) {

	int rowVal = *row;

	if(*numIdenticalDefs <= 0) *numIdenticalDefs = 
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
 * of errors), this function is forbidden from return()-ing early.
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
	int returnRowLikelihoods = 0;
	int numIdenticalDefs = 0, numIdenticalOrdinalMissingness = 0, numIdenticalOrdinalRows = 0,
		numIdenticalContinuousMissingness = 0, numIdenticalContinuousRows = 0;
	
	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns;
	omxMatrix *rowLikelihoods, *rowLogLikelihoods;
	omxMatrix *ordMeans, *ordCov, *ordRow, *contRow;
	omxMatrix *halfCov, *reduceCov, *ordContCov;
	omxData* data;
	double *lThresh, *uThresh, *corList, *weights;
	int *Infin;
	
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
	int numDefs = data->defVars.size();
	
	corList 	= ofo->corList;
	weights		= ofo->weights;
	lThresh		= ofo->lThresh;
	uThresh		= ofo->uThresh;
	returnRowLikelihoods = ofo->returnRowLikelihoods;
	rowLikelihoods = shared_ofo->rowLikelihoods;		// write-only
	rowLogLikelihoods = shared_ofo->rowLogLikelihoods;  // write-only
	
	Infin			= ofo->Infin;
	omxExpectation* expectation = localobj->expectation;
	std::vector< omxThresholdColumn > &thresholdCols = expectation->thresholds;
	
	Eigen::VectorXi ordRemove(cov->cols);
	Eigen::VectorXi contRemove(cov->cols);
	char u = 'U', l = 'L';
	int info;
	double determinant = 0.0;
	double oned = 1.0, zerod = 0.0, minusoned = -1.0;
	int onei = 1;
	double likelihood;
	int inform;
	
	bool firstRow = true;
	int row = rowbegin;
	
	Eigen::VectorXd oldDefs;
	oldDefs.resize(data->defVars.size());
	oldDefs.setConstant(NA_REAL);
	
	// [[Comment 4]] moving row starting position
	if (row > 0) {
		int prevIdentical = omxDataNumIdenticalRows(data, row - 1);
		row += (prevIdentical - 1);
	}
	
	int numContinuous = 0;
	int numOrdinal = 0;
	
	while(row < data->rows && (row - rowbegin) < rowcount) {
		mxLogSetCurrentRow(row);
		int numIdentical = omxDataNumIdenticalRows(data, row);
		if(numIdentical == 0) numIdentical = 1; 
		// N.B.: numIdentical == 0 means an error occurred and was not properly handled;
		// it should never be the case.
		
		omxDataRow(data, row, dataColumns, smallRow);                               // Populate data row
		
		if(OMX_DEBUG_ROWS(row)) {
			mxLog("Identicality check. Is %sfirst row of data. Total: %d rows identical, %d identical definition vars, %d identical missingness patterns. Continuous: %d rows, %d missingness patterns; Ordinal: %d rows, %d missingness patterns.", 
			(firstRow?"":"not "), numIdentical, numIdenticalDefs, omxDataNumIdenticalRows(data, row), 
			numIdenticalContinuousRows, numIdenticalContinuousMissingness, 
			numIdenticalOrdinalRows, numIdenticalOrdinalMissingness);
		}
		if(!strcmp(expectation->expType, "MxExpectationStateSpace")) {
			omxSetExpectationComponent(expectation, localobj, "y", smallRow);
		}
		//If the expectation is a state space model then
		// set the y attribute of the state space expectation to smallRow.
		
		if(numIdenticalDefs <= 0 || numIdenticalContinuousMissingness <= 0 || numIdenticalOrdinalMissingness <= 0 || 
			firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {  // If we're keeping covariance from the previous row, do not populate 
			// Handle Definition Variables.
			if((numDefs && numIdenticalDefs <= 0) || firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
				int numVarsFilled = 0;
				if(OMX_DEBUG_ROWS(row)) { mxLog("Handling Definition Vars."); }
				numVarsFilled = data->handleDefinitionVarList(localobj->matrix->currentState, row, oldDefs.data());
				if (numVarsFilled || firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
					if(row == 0 && !strcmp(expectation->expType, "MxExpectationStateSpace") ) {
						if(OMX_DEBUG){ mxLog("Resetting State Space state (x) and error cov (P)."); }
						omxSetExpectationComponent(expectation, localobj, "Reset", NULL);
					}
					omxExpectationCompute(expectation, NULL);
				}
			}
			// Filter down correlation matrix and calculate thresholds.
			// TODO: If identical ordinal or continuous missingness, ignore only the appropriate columns.
			numContinuous = 0;
			numOrdinal = 0;
			for(int j = 0; j < dataColumns->cols; j++) {
				int var = omxVectorElement(dataColumns, j);
				int value = omxIntDataElement(data, row, var);// Indexing correction means this is the index of the upper bound +1.
				// TODO: Might save time by preseparating ordinal from continuous.
				if(std::isnan(value) || value == NA_INTEGER) {  // Value is NA, therefore filter.
					ordRemove[j] = 1;
					contRemove[j] = 1;
					Infin[j] = -1;
					if(OMX_DEBUG_ROWS(row)) { 
						mxLog("Row %d, column %d, value %d.  NA.", row, j, value);
					}
					continue;
				}
				else if(omxDataColumnIsFactor(data, var)) {             // Ordinal column.
					++numOrdinal;
					ordRemove[j] = 0;
					contRemove[j] = 1;
					if(OMX_DEBUG_ROWS(row)) { 
						mxLog("Row %d, column %d, value %d.  Ordinal.", row, j, value);
					}
				} 
				else {
					++numContinuous;
					ordRemove[j] = 1;
					contRemove[j] = 0;
					if(OMX_DEBUG_ROWS(row)) { 
						mxLog("Row %d, column %d, value %d.  Continuous.", row, j, value);
					}
				}
			}
			
			if(OMX_DEBUG_ROWS(row)) {
				mxLog("Removals: %d ordinal, %d continuous out of %d total.",
				      dataColumns->cols - numOrdinal, dataColumns->cols - numContinuous, dataColumns->cols);
			}

			for(int j=0; j < dataColumns->cols; j++) {
				int var = omxVectorElement(dataColumns, j);
				if(!omxDataColumnIsFactor(data, var) || j >= int(thresholdCols.size()) || thresholdCols[j].numThresholds == 0) continue;
				omxRecompute(thresholdCols[j].matrix, fc); // Only one of these--save time by only doing this once
				checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column, thresholdCols[j].numThresholds, fc);
			}
			
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
					omxCopyMatrix(ordCov, cov);
					omxRemoveRowsAndColumns(ordCov, ordRemove.data(), ordRemove.data());
					
					// Recalculate ordinal fs
					omxCopyMatrix(ordMeans, means);
					omxRemoveElements(ordMeans, ordRemove.data()); 	    // Reduce the row to just ordinal.
					
					// These values pass through directly without modification by continuous variables
					
					// Calculate correlation matrix, correlation list, and weights from covariance
					omxStandardizeCovMatrix(ordCov, corList, weights);
				}
			} 
			else if( numIdenticalDefs <= 0 || numIdenticalContinuousRows <= 0 || firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
				
				/* Reset and Resample rows if necessary. */
				// First Cov and Means (if they've changed)
				if( numIdenticalDefs <= 0 || numIdenticalContinuousMissingness <= 0 || firstRow || !strcmp(expectation->expType, "MxExpectationStateSpace")) {
					if(OMX_DEBUG_ROWS(row)) { mxLog("Beginning to recompute inverse cov for standard models"); }
					
					/* If it's a state space expectation, extract the inverse rather than recompute it */
					if(!strcmp(expectation->expType, "MxExpectationStateSpace")) {
						smallMeans = omxGetExpectationComponent(expectation, localobj, "means");
						omxRemoveElements(smallMeans, contRemove.data());
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
					} //If it's a GREML expectation, extract the inverse rather than recompute it:
					else if(!strcmp(expectation->expType, "MxExpectationGREML")){
						smallMeans = omxGetExpectationComponent(expectation, localobj, "means");
						smallCov = omxGetExpectationComponent(expectation, localobj, "invcov");
						info = (int) omxGetExpectationComponent(expectation, localobj, "cholV_fail_om")->data[0];
						if(info!=0) {
							if (fc) fc->recordIterationError("expected covariance matrix is not positive-definite in data row %d", omxDataIndex(data, row));
							return TRUE;
						}
						determinant = 0.5 * omxGetExpectationComponent(expectation, localobj, "logdetV_om")->data[0];
						if(OMX_DEBUG_ROWS(row)) { mxLog("0.5*log(det(Cov)) is: %3.3f", determinant);}
					}
					else {
						/* Calculate derminant and inverse of Censored continuousCov matrix */
						omxCopyMatrix(smallMeans, means);
						omxRemoveElements(smallMeans, contRemove.data());
						omxCopyMatrix(smallCov, cov);
						omxRemoveRowsAndColumns(smallCov, contRemove.data(), contRemove.data());
						
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
					}
					
					if(info != 0) {
						if(!returnRowLikelihoods) {
							char *errstr = (char*) calloc(250, sizeof(char));
							sprintf(errstr, "Cannot invert expected continuous covariance matrix. Error %d.", info);
							omxRaiseError(errstr);
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
				omxCopyMatrix(contRow, smallRow);
				omxRemoveElements(contRow, contRemove.data()); 	// Reduce the row to just continuous.
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
						omxCopyMatrix(ordContCov, cov);
						omxRemoveRowsAndColumns(ordContCov, contRemove.data(), ordRemove.data());
						
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
						
						omxCopyMatrix(ordCov, cov);
						if(OMX_DEBUG_ROWS(row)) {
							mxLog("Resetting/Filtering Ordinal Covariance Matrix.");
							omxPrint(ordCov, "Reset to:");
						}
						
						omxRemoveRowsAndColumns(ordCov, ordRemove.data(), ordRemove.data());
						if(OMX_DEBUG_ROWS(row)) {
							mxLog("Resetting/Filtering Ordinal Covariance Matrix.");
							omxPrint(ordCov, "Filtered to:");
						}
						
						// FIXME: This assumes that ordCov and reducCov have the same row/column majority.
						int vlen = reduceCov->rows * reduceCov->cols;
						F77_CALL(daxpy)(&vlen, &minusoned, reduceCov->data, &onei, ordCov->data, &onei); // ordCov <- (ordCov - reduceCov) %*% cont/ord
						
					}
					
					// Projected means must be recalculated if the continuous variables change at all.
					omxCopyMatrix(ordMeans, means);
					omxRemoveElements(ordMeans, ordRemove.data()); 	    // Reduce the row to just ordinal.
					F77_CALL(dgemv)((smallCov->minority), &(halfCov->rows), &(halfCov->cols), &oned, halfCov->data, &(halfCov->leading), contRow->data, &onei, &oned, ordMeans->data, &onei);                      // ordMeans += halfCov %*% contRow
				}
				
			} // End of continuous likelihood values calculation
			
			if(numOrdinal <= 0) {       // No Ordinal Vars at all.
			likelihood = 1;
			} 
			else {  
				// There are ordinal vars, and not everything is identical, so we're recalculating
				// Calculate correlation matrix, correlation list, and weights from covariance
				if(numIdenticalDefs <=0 || numIdenticalContinuousMissingness <= 0 || numIdenticalOrdinalMissingness <= 0 || firstRow) {
					// if(OMX_DEBUG_ROWS(row)) {omxPrint(ordCov, "Ordinal cov matrix for standardization."); } //:::DEBUG:::
					omxStandardizeCovMatrix(ordCov, corList, weights);
				}
				
				omxCopyMatrix(ordRow, smallRow);
				omxRemoveElements(ordRow, ordRemove.data()); 	    // Reduce the row to just ordinal.
				
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
					} 
					else {
						lThresh[count] = (omxMatrixElement(thresholdCols[j].matrix, value-1, thresholdCols[j].column) - offset) / weight;
						if(thresholdCols[j].numThresholds > value) {	// Highest threshold = Inf
						double tmp = (omxMatrixElement(thresholdCols[j].matrix, value, thresholdCols[j].column) - offset) / weight;
						uThresh[count] = tmp;
						Infin[count] = 2;
						} 
						else {
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
						if (0) {
							// This diagnostic triggers an omxMatrixElement by models/passing/JointFIMLTest.R
							mxLog("       Thresholds were %f -> %f, scaled by weight %f and shifted by mean %f and total offset %f.",
							omxMatrixElement(thresholdCols[j].matrix, (Infin[count]==0?0:value-1), thresholdCols[j].column), 
							omxMatrixElement(thresholdCols[j].matrix, (Infin[count]==1?value-1:value), thresholdCols[j].column), 
							weight, (means==NULL?0:omxVectorElement(ordMeans, count)), offset);
						}
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
			} 
			else {
				double logLikelihood = -2 * log(likelihood);       // -2 Log of ordinal likelihood
				logLikelihood += ((2 * determinant) + Q + (log(2 * M_PI) * numContinuous));    // -2 Log of continuous likelihood
				logLikelihood *= numIdentical;
				
				for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(rowLikelihoods, omxDataIndex(data, j), 0, rowLikelihood);
				}
				omxSetMatrixElement(rowLogLikelihoods, row, 0, logLikelihood);
				
				if(OMX_DEBUG_ROWS(row)) { 
					mxLog("row log Likelihood %3.3f + %3.3f + %3.3f + %3.3f= %3.3f, ", 
					      (2.0*determinant), Q, (log(2 * M_PI) * numContinuous), 
					-2  * log(rowLikelihood), (2.0 *determinant) + Q + (log(2 * M_PI) * numContinuous));
				} 
				
			}
			if(firstRow) firstRow = false;
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
