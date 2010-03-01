/*
 *  Copyright 2007-2009 The OpenMx Project
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

#ifndef _OMX_FIML_OBJECTIVE_
#define _OMX_FIML_OBJECTIVE_ TRUE

extern omxMatrix** matrixList;
extern void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*, double*, double*, double*, double*, int*);

typedef struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
	omxData* source;		// Data source
	int numLocations;		// Num locations
	double** location;		// And where it goes
	omxMatrix** matrices;	// Matrix numbers for dirtying

} omxDefinitionVar;

typedef struct omxThresholdColumn {		 	// Definition Var
	omxMatrix* matrix;		// Which Matrix/Algebra it comes from
	int column;				// Which column has the thresholds
	int numThresholds;		// And how many thresholds
} omxThresholdColumn;

typedef struct omxFIMLRowOutput {  // Output object for each row of estimation.  Mirrors the Mx1 output vector
	double Minus2LL;		// Minus 2 Log Likelihood
	double Mahalanobis;		// Square root of Mahalanobis distance Q = (row - means) %*% solve(sigma) %*% (row - means)
	double Z;				// Z-score of Mahalnobis.  This is ((Q/n )^(1/3) - 1 + 2/(9n)) (9n/2)^(.5)
	double Nrows;			// Number of rows in the data set--Unclear why this is returned
	double Ncols;			// Number of columns in the (censored) data set
	int finalMissed;		// Whether or not likelihood was calculable in the final case
	int modelNumber;		// Not used
} omxFIMLRowOutput;

typedef struct omxFIMLObjective {

	/* Parts of the R  MxFIMLObjective Object */
	omxMatrix* cov;				// Covariance Matrix
	omxMatrix* means;			// Vector of means
	omxData* data;				// The data
	omxMatrix* dataColumns;		// The order of columns in the data matrix
	omxMatrix* dataRow;			// One row of data
	int returnRowLikelihoods;   // Whether or not to return row-by-row likelihoods

//	double* zeros;

	/* Structures determined from info in the MxFIMLObjective Object*/
	omxDefinitionVar* defVars;	// A list of definition variables
	int numDefs;				// The length of the defVars list

	/* Reserved memory for faster calculation */
	omxMatrix* smallRow;		// Memory reserved for operations on each data row
	omxMatrix* smallCov;		// Memory reserved for operations on covariance matrix
	omxMatrix* smallMeans;		// Memory reserved for operations on the means matrix
	omxMatrix* RCX;				// Memory reserved for computations

	/* Structures for FIMLOrdinalObjective Objects */
	omxMatrix* cor;				// To calculate correlation matrix from covariance
	double* weights;			// Covariance weights to shift parameter estimates
	omxMatrix* smallThresh;		// Memory reserved for reduced threshold matrix
	omxThresholdColumn* thresholdCols;		// List of column thresholds

	/* Argument space for SADMVN function */
	double* lThresh;			// Specific list of lower thresholds
	double* uThresh;			// Specific list of upper thresholds
	double* corList;			// SADMVN-specific list of correlations
	double* smallCor;			// Reduced SADMVN-specific list of correlations
	int* Infin;					// Which thresholds to use
	int maxPts;					// From MxOptions (?)
	double absEps;				// From MxOptions
	double relEps;				// From MxOptions

} omxFIMLObjective;

void omxDestroyFIMLObjective(omxObjective *oo) {

}

omxRListElement* omxSetFinalReturnsFIMLObjective(omxObjective *oo, int *numReturns) {

	omxFIMLObjective* ofo = (omxFIMLObjective *)(oo->argStruct);

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

void handleDefinitionVarList(omxData* data, int row, omxDefinitionVar* defVars, int numDefs) {

	if(OMX_DEBUG_ROWS) { Rprintf("Processing Definition Vars.\n"); }

	/* Fill in Definition Var Estimates */
	for(int k = 0; k < numDefs; k++) {
		if(defVars[k].source != data) {
			error("Internal error: definition variable population into incorrect data source");
			continue; // don't populate from the wrong data frame
		}
		for(int l = 0; l < defVars[k].numLocations; l++) {
			if(OMX_DEBUG_ROWS) {
				Rprintf("Populating column %d (value %3.2f) into matrix %d.\n", defVars[k].column, omxDoubleDataElement(defVars[k].source, row, defVars[k].column), defVars[k].matrices[l]);
			}
			*(defVars[k].location[l]) = omxDoubleDataElement(data, row, defVars[k].column);
			omxMarkDirty(defVars[k].matrices[l]);
			if(ISNA(omxDoubleDataElement(data, row, defVars[k].column))) {
				error("Error NYI: Missing Definition Vars Not Yet Implemented.");
			}
		}
	}
}

void omxCallFIMLOrdinalObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
	/* TODO: Current implementation is slow: update by filtering correlations and thresholds. */
	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n");}
	// Requires: Data, means, covariances, thresholds

	double sum;
	double Q = 0.0;
	double logDet = 0;
	int numDefs;
	int numCols, numRemoves = 0;
	int returnRowLikelihoods = 0;

	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns;
	omxMatrix *cor, *smallThresh;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *lThresh, *uThresh, maxPts, absEps, relEps, *corList, *smallCor, *weights;
	int *Infin;
	omxDefinitionVar* defVars;

	// Locals, for readability.  Compiler should cut through this.
	cov 		= ((omxFIMLObjective*)oo->argStruct)->cov;
	means		= ((omxFIMLObjective*)oo->argStruct)->means;
	smallRow 	= ((omxFIMLObjective*)oo->argStruct)->smallRow;
	smallCov 	= ((omxFIMLObjective*)oo->argStruct)->smallCov;
	smallMeans	= ((omxFIMLObjective*)oo->argStruct)->smallMeans;
	RCX 		= ((omxFIMLObjective*)oo->argStruct)->RCX;
	data		= ((omxFIMLObjective*)oo->argStruct)->data;
	dataColumns	= ((omxFIMLObjective*)oo->argStruct)->dataColumns;
	defVars		= ((omxFIMLObjective*)oo->argStruct)->defVars;
	numDefs		= ((omxFIMLObjective*)oo->argStruct)->numDefs;

	cor 		= ((omxFIMLObjective*)oo->argStruct)->cor;
	corList 	= ((omxFIMLObjective*)oo->argStruct)->corList;
	smallCor	= ((omxFIMLObjective*)oo->argStruct)->smallCor;
	weights		= ((omxFIMLObjective*)oo->argStruct)->weights;
	smallThresh	= ((omxFIMLObjective*)oo->argStruct)->smallThresh;
	lThresh		= ((omxFIMLObjective*)oo->argStruct)->lThresh;
	uThresh		= ((omxFIMLObjective*)oo->argStruct)->uThresh;
	thresholdCols = ((omxFIMLObjective*)oo->argStruct)->thresholdCols;
	returnRowLikelihoods = ((omxFIMLObjective*)oo->argStruct)->returnRowLikelihoods;

	Infin		= ((omxFIMLObjective*)oo->argStruct)->Infin;
	maxPts		= ((omxFIMLObjective*)oo->argStruct)->maxPts;
	absEps		= ((omxFIMLObjective*)oo->argStruct)->absEps;
	relEps		= ((omxFIMLObjective*)oo->argStruct)->relEps;

	if(numDefs == 0) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("No Definition Vars: precalculating."); }
		omxRecompute(cov);			// Only recompute this here if there are no definition vars
		omxRecompute(means);
		for(int j = 0; j < data->cols; j++) {
			if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
				omxRecompute(thresholdCols[j].matrix);
				checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column);
			}
		}
		omxStandardizeCovMatrix(cov, corList, weights);	// Calculate correlation and covariance
	}

	sum = 0.0;

	for(int row = 0; row < data->rows; row++) {
		oo->matrix->currentState->currentRow = row;		// Set to a new row.
		logDet = 0.0;
		Q = 0.0;

		// Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
		numCols = 0;
		numRemoves = 0;

		// Handle Definition Variables.
		if(numDefs != 0) {
			if(OMX_DEBUG_ROWS) { Rprintf("Handling Definition Vars.\n"); }
			handleDefinitionVarList(data, row, defVars, numDefs);
			omxRecompute(cov);
			omxRecompute(means);
			for(int j=0; j < data->cols; j++) {
				if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
					omxRecompute(thresholdCols[j].matrix);
					checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column);
				}
			}
			omxStandardizeCovMatrix(cov, corList, weights);	// Calculate correlation and covariance
		}

		int nextRow = 0;

		for(int j = 0; j < dataColumns->cols; j++) {
			int var = omxVectorElement(dataColumns, j);
			int value = omxIntDataElement(data, row, var);// Indexing correction means this is the index of the upper bound +1.
			if(ISNA(value) || value == NA_INTEGER) {  // Value is NA, therefore filter.
				numRemoves++;
				// toRemove[j] = 1;
				Infin[j] = -1;
				continue;
			} else {									// For joint, check here for continuousness
				value--;								// Correct for C indexing: value is now the index of the upper bound.
				// TODO: Subsample the corList and thresholds for speed. Doesn't look like that's much of a speedup.
				// TODO: Check for high and low thresholds of NA
				double mean;
				if(means == NULL) mean = 0;
				else mean = omxVectorElement(means, var);
				double weight = weights[var];
				// toRemove[j] = 0;
				// Rprintf(":::Data column %d, matrix column %d, values %d & %d, matrix size (%d, %d), with %d thresholds:::\n", j, thresholdCols[j].column, value, value-1, thresholdCols[j].matrix->rows, thresholdCols[j].matrix->cols, thresholdCols[j].numThresholds);
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

			if(OMX_DEBUG_ROWS) { Rprintf("Row %d, column %d.  Thresholds for data column %d and row %d are %f -> %f. (Infin=%d)\n", row, j, var, value-1, lThresh[j], uThresh[j], Infin[j]);}
			nextRow++;
			}
		}

		if(numRemoves >= smallCov->rows) {
			if(returnRowLikelihoods) {
			    omxSetMatrixElement(oo->matrix, row, 0, 1.0);
			}
			continue;
		}

		// SADMVN calls Alan Genz's sadmvn.f--see appropriate file for licensing info.
		// TODO: Check with Genz: should we be using sadmvn or sadmvn?
		// Parameters are:
		// 	N 		int			# of vars
		//	Lower	double*		Array of lower bounds
		//	Upper	double*		Array of upper bounds
		//	Infin	int*		Array of flags: 0 = (-Inf, upper] 1 = [lower, Inf), 2 = [lower, upper]
		//	Correl	double*		Array of correlation coeffs: in row-major lower triangular order
		//	MaxPts	int			Maximum # of function values (use 1000*N or 1000*N*N)
		//	Abseps	double		Absolute error tolerance.  Yick.
		//	Releps	double		Relative error tolerance.  Use EPSILON.
		//	Error	&double		On return: absolute real error, 99% confidence
		//	Value	&double		On return: evaluated value
		//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
		// TODO: Separate block diagonal covariance matrices into pieces for integration separately
		double Error;
		double absEps = 1e-3;
		double relEps = 0;
		int MaxPts = 100000*cov->rows;
		double likelihood;
		int inform;
		int numVars = smallCov->rows;
		/* FOR DEBUGGING PURPOSES */
/*		numVars = 2;
		lThresh[0] = -2;
		uThresh[0] = -1.636364;
		Infin[0] = 2;
		lThresh[1] = 0;
		uThresh[1] = 0;
		Infin[1] = 0;
		smallCor[0] = 1.0; smallCor[1] = 0; smallCor[2] = 1.0; */
		F77_CALL(sadmvn)(&numVars, lThresh, uThresh, Infin, corList, &MaxPts, &absEps, &relEps, &Error, &likelihood, &inform);

		if(!oo->matrix->currentState->currentRow && OMX_DEBUG) {
			char infinCodes[3][20];
			strcpy(infinCodes[0], "(-INF, upper]");
			strcpy(infinCodes[1], "[lower, INF)");
			strcpy(infinCodes[2], "[lower, upper]");
			Rprintf("Input to sadmvn is (%d rows):\n", numVars);

			omxPrint(omxDataRow(data, row, dataColumns, smallRow), "Data Row");

			for(int i = 0; i < numVars; i++) {
				Rprintf("Row %d: %f, %f, %d(%s)\n", i, lThresh[i], uThresh[i], Infin[i], infinCodes[Infin[i]]);
			}

			Rprintf("Cor: (Lower %d x %d):", cov->rows, cov->cols);
			for(int i = 0; i < cov->rows*(cov->rows-1)/2; i++) {
				// Rprintf("Row %d of Cor: ", i);
				// for(int j = 0; j < i; j++)
				Rprintf(" %f", corList[i]); // (i*(i-1)/2) + j]);
				// Rprintf("\n");
			}
			Rprintf("\n");
		}

		if(OMX_DEBUG_ROWS) { Rprintf("Output of sadmvn is %f, %f, %d.\n", Error, likelihood, inform); }

		if(inform == 2) {
			error("Improper input to sadmvn.");
		}

		if(returnRowLikelihoods) {
            if(OMX_DEBUG_ROWS) {Rprintf("Row %d likelihood is %3.3f.\n", row, likelihood);}
            omxSetMatrixElement(oo->matrix, row, 0, likelihood);
        } else {
            logDet = -2 * log(likelihood);

            sum += logDet;// + (log(2 * M_PI) * (cov->cols - numRemoves));

            if(OMX_DEBUG_ROWS) {
                Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n",
				    sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
            }
        }
	}

    if(!returnRowLikelihoods) {
        if(OMX_DEBUG) {
            Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n",
                sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
        }

        oo->matrix->data[0] = sum;
    }
}

void omxCallFIMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n");}
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

	double sum;
	char u = 'U';
	int info = 0;
	double oned = 1.0;
	double zerod = 0.0;
	int onei = 1;
	double Q = 0.0;
	double determinant = 1.0;
	int numDefs;
	int numCols, numRemoves;
	int returnRowLikelihoods;

	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *dataColumns;
	omxDefinitionVar* defVars;
	omxData *data;

    omxFIMLObjective* ofo = ((omxFIMLObjective*)oo->argStruct);

	cov 		= ofo->cov;		// Locals, for readability.  Should compile out.
	means		= ofo->means;
	smallRow 	= ofo->smallRow;
	smallCov 	= ofo->smallCov;
	RCX 		= ofo->RCX;
	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	defVars		= ofo->defVars;
	numDefs		= ofo->numDefs;
    returnRowLikelihoods = ofo->returnRowLikelihoods;

	if(numDefs == 0) {
		if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
		omxRecompute(cov);			// Only recompute this here if there are no definition vars
		omxRecompute(means);
		if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
		if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
	}
	int toRemove[cov->cols];
	int zeros[cov->cols];

	sum = 0.0;

	for(int row = 0; row < data->rows; row++) {
		oo->matrix->currentState->currentRow = row;		// Set to a new row.
		determinant = 1.0;
		Q = 0.0;

		numCols = 0;
		numRemoves = 0;
		omxResetAliasedMatrix(smallRow); 									// Resize smallrow
		omxDataRow(data, row, dataColumns, smallRow);						// Populate data row

		// Handle Definition Variables.
		if(numDefs != 0) {
			// omxStateNextRow(oo->matrix->currentState);						// Advance Row
			if(OMX_DEBUG_ROWS) {Rprintf("Handling definition vars and calculating cov and means for row %d.\n", oo->matrix->currentState->currentRow);}
			handleDefinitionVarList(data, row, defVars, numDefs);
			omxRecompute(cov);
			omxRecompute(means);
		}

		if(OMX_DEBUG_ROWS) { omxPrint(means, "Local Means"); }
		if(OMX_DEBUG_ROWS) {
			char note[50];
			sprintf(note, "Local Data Row %d", row);
			omxPrint(smallRow, note); }

		// Determine how many rows/cols to remove.
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxMatrixElement(smallRow, 0, j);
			if(isnan(dataValue) || dataValue == NA_REAL) {
				numRemoves++;
				toRemove[j] = 1;
			} else {
				if(means != NULL) {
					omxSetMatrixElement(smallRow, 0, j, (dataValue -  omxVectorElement(means, j)));
				}
				numCols++;
				toRemove[j] = 0;
			}
			zeros[j] = 0;
		}

		if(cov->cols <= numRemoves) {
            if(returnRowLikelihoods) {
                omxSetMatrixElement(oo->matrix, row, 0, 1);
            }
            continue;
        }

		omxRemoveRowsAndColumns(smallRow, 0, numRemoves, zeros, toRemove); 	// Reduce it.

		omxResetAliasedMatrix(smallCov);						// Subsample covariance matrix
		omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);

		if(OMX_DEBUG_ROWS) { omxPrint(smallCov, "Local Covariance Matrix"); }

		/* The Calculation */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (t(dataRow) %*% (solve(ExpectedCov)) %*% dataRow)^(1/2) */
		F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			if(!returnRowLikelihoods) {
				char *errstr = calloc(250, sizeof(char));
				sprintf(errstr, "Expected covariance matrix is not positive-definite in row %d.\n", row);
				omxRaiseError(oo->matrix->currentState, -1, errstr);
				free(errstr);
				return;
			} else {
				omxSetMatrixElement(oo->matrix, row, 0, 0.0);
				continue;
			}
		}
		for(int diag = 0; diag < (smallCov->rows); diag++) {
			determinant *= omxMatrixElement(smallCov, diag, diag);
		}
		determinant = determinant * determinant;

		F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			if(!returnRowLikelihoods) {
				char *errstr = calloc(250, sizeof(char));
				sprintf(errstr, "Cannot invert expected covariance matrix. Error %d.", info);
				omxRaiseError(oo->matrix->currentState, -1, errstr);
				free(errstr);
				return;
			} else {
				omxSetMatrixElement(oo->matrix, row, 0, 0.0);
				continue;
			}
		}
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);

        if(returnRowLikelihoods) {
            if(OMX_DEBUG_ROWS) {Rprintf("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f\n", pow(2 * M_PI, -.5 * smallRow->cols), (1.0/sqrt(determinant)), exp(-.5 * Q), pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q));}
            sum = pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q);
            omxSetMatrixElement(oo->matrix, row, 0, sum);
        } else {
            sum += log(determinant) + Q + (log(2 * M_PI) * smallRow->cols);
            if(OMX_DEBUG_ROWS) {Rprintf("Change in Total Likelihood for row %d is %3.3f + %3.3f + %3.3f = %3.3f, total Likelihood is %3.3f\n", oo->matrix->currentState->currentRow, log(determinant), Q, (log(2 * M_PI) * smallRow->cols), log(determinant) + Q + (log(2 * M_PI) * smallRow->cols), sum);}
        }
	}

    if(!returnRowLikelihoods) {
	   if(OMX_VERBOSE || OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
	   omxSetMatrixElement(oo->matrix, 0, 0, sum);
    }

}

unsigned short int omxNeedsUpdateFIMLObjective(omxObjective* oo) {
	return omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->cov)
		|| omxMatrixNeedsUpdate(((omxFIMLObjective*)oo->argStruct)->means);
}

void omxInitFIMLObjective(omxObjective* oo, SEXP rObj) {

	if(OMX_DEBUG) { Rprintf("Initializing FIML objective function.\n"); }

	SEXP nextMatrix, itemList, nextItem, dataSource, columnSource;
	int nextDef, index, hasOrdinal = FALSE;
	int *nextInt;
	omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(1, sizeof(omxFIMLObjective));

	PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
	newObj->means = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	if(newObj->means == NULL) { error("No means model in FIML evaluation.");}
	UNPROTECT(1);

	PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
	newObj->cov = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);

	PROTECT(nextMatrix = GET_SLOT(rObj, install("thresholdColumns")));
	if(OMX_DEBUG) {Rprintf("Processing Thresholds.\n");}
	newObj->thresholdCols = (omxThresholdColumn *) R_alloc(length(nextMatrix), sizeof(omxThresholdColumn));
	for(index = 0; index < length(nextMatrix); index++) {
		PROTECT(itemList = VECTOR_ELT(nextMatrix, index));
		if(INTEGER(itemList)[0] == NA_INTEGER) {	// Continuous variable
			if(OMX_DEBUG) {Rprintf("Threshold %d fail.\n", index);}
			newObj->thresholdCols[index].matrix = NULL;
			newObj->thresholdCols[index].column = 0;
			newObj->thresholdCols[index].numThresholds = 0;
		} else {
			// PROTECT(nextItem = VECTOR_ELT(itemList, 0));
			nextInt = INTEGER(itemList);
			if(OMX_DEBUG) {Rprintf("Threshold %d processed.\n", index);}
			newObj->thresholdCols[index].matrix = omxNewMatrixFromMxIndex(itemList, oo->matrix->currentState);
			newObj->thresholdCols[index].column = nextInt[1];
			newObj->thresholdCols[index].numThresholds = newObj->thresholdCols[index].matrix->rows;
			hasOrdinal++;
			// UNPROTECT(1); /* nextItem */
		}
		UNPROTECT(1); /* itemList */
	}
	if(OMX_DEBUG) {Rprintf("%d thresholds processed.\n", index);}
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing data source.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("data"))); // TODO: Need better way to process data elements.
	newObj->data = omxNewDataFromMxDataPtr(nextMatrix, oo->matrix->currentState);
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing row likelihood option.\n"); }
	PROTECT(nextMatrix = AS_INTEGER(GET_SLOT(rObj, install("vector")))); // preparing the object by using the vector to populate and the flag
	newObj->returnRowLikelihoods = INTEGER(nextMatrix)[0];
	if(newObj->returnRowLikelihoods) {
	   omxResizeMatrix(oo->matrix, newObj->data->rows, 1, FALSE); // 1=column matrix, FALSE=discards memory as this is a one time resize
    }
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing variable mapping structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("dataColumns")));
	newObj->dataColumns = omxNewMatrixFromMxMatrix(nextMatrix, oo->matrix->currentState);
	if(OMX_DEBUG) {omxPrint(newObj->dataColumns, "Variable mapping"); }
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing definition variables structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
	newObj->numDefs = length(nextMatrix);
	if(OMX_DEBUG) {Rprintf("Number of definition variables is %d.\n", newObj->numDefs); }
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	for(nextDef = 0; nextDef < newObj->numDefs; nextDef++) {
		PROTECT(itemList = VECTOR_ELT(nextMatrix, nextDef));
		PROTECT(dataSource = VECTOR_ELT(itemList, 0));
		if(OMX_DEBUG) {Rprintf("Data source number is %d.\n", INTEGER(dataSource)[0]); }
		newObj->defVars[nextDef].data = INTEGER(dataSource)[0];
		newObj->defVars[nextDef].source = oo->matrix->currentState->dataList[INTEGER(dataSource)[0]];
		PROTECT(columnSource = VECTOR_ELT(itemList, 1));
		if(OMX_DEBUG) {Rprintf("Data column number is %d.\n", INTEGER(columnSource)[0]); }
		newObj->defVars[nextDef].column = INTEGER(columnSource)[0];
		UNPROTECT(2); // unprotect dataSource and columnSource
		newObj->defVars[nextDef].numLocations = length(itemList) - 2;
		newObj->defVars[nextDef].location = (double **) R_alloc(length(itemList) - 2, sizeof(double*));
		newObj->defVars[nextDef].matrices = (omxMatrix **) R_alloc(length(itemList) - 2, sizeof(omxMatrix*));
		for(index = 2; index < length(itemList); index++) {
			PROTECT(nextItem = VECTOR_ELT(itemList, index));
			newObj->defVars[nextDef].location[index-2] = omxLocationOfMatrixElement(
				oo->matrix->currentState->matrixList[INTEGER(nextItem)[0]],
				INTEGER(nextItem)[1], INTEGER(nextItem)[2]);
			newObj->defVars[nextDef].matrices[index-2] = oo->matrix->currentState->matrixList[INTEGER(nextItem)[0]];
			UNPROTECT(1); // unprotect nextItem
		}
		UNPROTECT(1); // unprotect itemList
	}
	UNPROTECT(1); // unprotect nextMatrix

	/* Temporary storage for calculation */
	int covCols = newObj->cov->cols;
	newObj->smallRow = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
	newObj->smallCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
	newObj->RCX = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
//	newObj->zeros = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);

	omxAliasMatrix(newObj->smallCov, newObj->cov);					// Will keep its aliased state from here on.

	oo->objectiveFun = omxCallFIMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateFIMLObjective;
	oo->setFinalReturns = omxSetFinalReturnsFIMLObjective;
	oo->destructFun = omxDestroyFIMLObjective;
	oo->repopulateFun = NULL;

//	Rprintf("Checking 0x%d.", newObj->data);
	if(OMX_DEBUG) { Rprintf("Checking %d.", omxDataColumnIsFactor(newObj->data, 0)); }

	if(hasOrdinal) {
		if(OMX_DEBUG) { Rprintf("Ordinal Data detected.  Using Ordinal FIML."); }
		newObj->weights = (double*) R_alloc(covCols, sizeof(double));
		newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		newObj->corList = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
		newObj->smallCor = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
		newObj->lThresh = (double*) R_alloc(covCols, sizeof(double));
		newObj->uThresh = (double*) R_alloc(covCols, sizeof(double));
		newObj->Infin = (int*) R_alloc(covCols, sizeof(int));

		oo->objectiveFun = omxCallFIMLOrdinalObjective;
	}

	oo->argStruct = (void*) newObj;
}

#endif /* _OMX_FIML_OBJECTIVE_ */
