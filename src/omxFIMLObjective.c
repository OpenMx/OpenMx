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
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"

#ifndef _OMX_FIML_OBJECTIVE_
#define _OMX_FIML_OBJECTIVE_ TRUE
#define OMX_DEBUG_ROWS FALSE

extern omxMatrix** matrixList;
extern void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*, double*, double*, double*, double*, int*);

typedef struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
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
	int misses;				// How many times likelihood broke on this row
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
	omxMatrix* weights;			// Covariance weights to shift parameter estimates
	omxMatrix* smallWeights;	// Memory reserved for reduced weight matrix	
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

void handleDefinitionVarList(omxData* data, int row, omxDefinitionVar* defVars, int numDefs) {

	if(OMX_DEBUG_ROWS) { Rprintf("Processing Definition Vars.\n"); }

	/* Fill in Definition Var Estimates */
	for(int k = 0; k < numDefs; k++) {
		for(int l = 0; l < defVars[k].numLocations; l++) {
			if(OMX_DEBUG_ROWS) {
				Rprintf("Populating column %d (value %3.2f) into matrix %d.\n", defVars[k].column, omxDoubleDataElement(data, row, defVars[k].column), defVars[k].matrices[l]);
			}
			*(defVars[k].location[l]) = omxDoubleDataElement(data, row, defVars[k].column);
			omxMarkDirty(defVars[k].matrices[l]);
			if(ISNA(omxDoubleDataElement(data, row, defVars[k].column))) {
				error("Error NYI: Missing Definition Vars Not Yet Implemented.");
			}
		}
	}
}

void omxDestroyFIMLObjective(omxObjective *oo) {

}

void omxStandardizeCovMatrix(omxMatrix* cov, double* corList, omxMatrix* weights) {
	// Maybe coerce this into an algebra or sequence of algebras?
	
	if(OMX_DEBUG) { Rprintf("Standardizing matrix."); }
	
	int rows = cov->rows;

	for(int i = 0; i < rows; i++) {
		omxSetMatrixElement(weights, 0, i, sqrt(omxMatrixElement(cov, i, i)));
	}

	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < i; j++) {
			corList[((i*(i-1))/2) + j] = omxMatrixElement(cov, i, j) / (omxVectorElement(weights, i) * omxVectorElement(weights, j));
		}
	}
}

void checkIncreasing(omxMatrix* om, int column) {
	double previous = - INFINITY;
	double current;
	for(int j = 0; j < om->rows; j++ ) {
		current = omxMatrixElement(om, j, column);
		if(current == NA_REAL || current == NA_INTEGER) {
			continue;
		}
		if(current <= previous) {
			char errstr[250];
			sprintf(errstr, "Thresholds are not strictly increasing.");
			omxRaiseError(om->currentState, -1, errstr);
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

	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns, *smallWeights;
	omxMatrix *cor, *weights, *smallThresh;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *lThresh, *uThresh, maxPts, absEps, relEps, *corList, *smallCor;
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
	smallWeights = ((omxFIMLObjective*)oo->argStruct)->smallWeights;
	thresholdCols = ((omxFIMLObjective*)oo->argStruct)->thresholdCols;

	Infin		= ((omxFIMLObjective*)oo->argStruct)->Infin;
	maxPts		= ((omxFIMLObjective*)oo->argStruct)->maxPts;
	absEps		= ((omxFIMLObjective*)oo->argStruct)->absEps;
	relEps		= ((omxFIMLObjective*)oo->argStruct)->relEps;
	
	if(numDefs == 0) {
		if(OMX_DEBUG) { Rprintf("No Definition Vars: precalculating."); }
		omxRecompute(cov);			// Only recompute this here if there are no definition vars
		for(int j = 0; j < data->cols; j++) {
			if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
				omxRecompute(thresholdCols[j].matrix);
				checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column);
			}
		}
		omxResetAliasedMatrix(weights);
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
				double weight = omxVectorElement(weights, var);
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
		
		if(numRemoves >= smallCov->rows) continue;

		// SADMVN calls Alan Gentz's sadmvn.f--see appropriate file for licensing info.
		// TODO: Check with Gentz: should we be using sadmvn or sadmvn?
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
			Rprintf("Input to sadmvn is (%d rows):\n", numVars);
			
			omxPrint(omxDataRow(data, row, dataColumns, smallRow), "Data Row");
			
			for(int i = 0; i < numVars; i++) {
				Rprintf("Row %d: %f, %f, %d\n", i, lThresh[i], uThresh[i], Infin[i]);
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
		
		if(likelihood < 10e-10) {
			char errstr[250];
			sprintf(errstr, "Likelihood 0 for row %d.", row);
			if(OMX_DEBUG) { Rprintf(errstr); }
			omxRaiseError(smallCov->currentState, -1, errstr);
			return;			
		}

		logDet = -2 * log(likelihood);

		sum += logDet;// + (log(2 * M_PI) * (cov->cols - numRemoves));

		if(!oo->matrix->currentState->currentRow && OMX_DEBUG) {
			Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n", 
				sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
		}
	}

	if(OMX_DEBUG) {
		Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n", 
			sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
	}

	oo->matrix->data[0] = sum;

}

omxRListElement* omxSetFinalReturnsFIMLObjective(omxObjective *oo, int *numReturns) {
	*numReturns = 1;
	omxRListElement* retVal = (omxRListElement*) R_alloc(1, sizeof(omxRListElement));

	retVal[0].numValues = 1;
	retVal[0].values = (double*) R_alloc(1, sizeof(double));
	strncpy(retVal[0].label, "Minus2LogLikelihood", 20);
	retVal[0].values[0] = omxMatrixElement(oo->matrix, 0, 0);
	
	return retVal;
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
	double logDet = 0;
	int numDefs;
	int numCols, numRemoves;

	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *dataColumns;
	omxDefinitionVar* defVars;
	omxData *data;

	cov 		= ((omxFIMLObjective*)oo->argStruct)->cov;		// Locals, for readability.  Should compile out.
	means		= ((omxFIMLObjective*)oo->argStruct)->means;
	smallRow 	= ((omxFIMLObjective*)oo->argStruct)->smallRow;
	smallCov 	= ((omxFIMLObjective*)oo->argStruct)->smallCov;
	RCX 		= ((omxFIMLObjective*)oo->argStruct)->RCX;
	data		= ((omxFIMLObjective*)oo->argStruct)->data;
	dataColumns	= ((omxFIMLObjective*)oo->argStruct)->dataColumns;
	defVars		= ((omxFIMLObjective*)oo->argStruct)->defVars;
	numDefs		= ((omxFIMLObjective*)oo->argStruct)->numDefs;

	if(numDefs == 0) {
		omxRecompute(cov);			// Only recompute this here if there are no definition vars
		omxRecompute(means);
	}
	
	if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }

	int toRemove[cov->cols];
	int zeros[cov->cols];

	sum = 0.0;

	for(int row = 0; row < data->rows; row++) {
		oo->matrix->currentState->currentRow = row;		// Set to a new row.
		logDet = 1.0;
		Q = 0.0;

		numCols = 0;
		numRemoves = 0;
		omxResetAliasedMatrix(smallRow); 									// Resize smallrow
		omxDataRow(data, row, dataColumns, smallRow);						// Populate data row

		// Handle Definition Variables.
		if(numDefs != 0) {
			handleDefinitionVarList(data, row, defVars, numDefs);
			omxStateNextRow(oo->matrix->currentState);							// Advance Row
			omxRecompute(cov);
			omxRecompute(means);
		}

		if(OMX_DEBUG) { omxPrint(means, "Local Means"); }
		if(OMX_DEBUG) { omxPrint(smallRow, "Local Data Row"); }
		
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

		if(cov->cols <= numRemoves) continue;
		omxRemoveRowsAndColumns(smallRow, 0, numRemoves, zeros, toRemove); 	// Reduce it.
		
		omxResetAliasedMatrix(smallCov);						// Subsample covariance matrix
		omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);
		
		if(OMX_DEBUG) { omxPrint(smallCov, "Local Covariance Matrix"); }

		/* The Calculation */
		F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			char errStr[250];
			sprintf(errStr, "Covariance matrix is not positive-definite.");
			omxRaiseError(oo->matrix->currentState, -1, errStr);
			return;
		}
		for(int diag = 0; diag < (smallCov->rows); diag++) {
			logDet *= omxMatrixElement(smallCov, diag, diag);
		}
		logDet = log(logDet * logDet);
		
		F77_CALL(dpotri)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
		if(info != 0) {
			char errstr[250];
			sprintf(errstr, "Cannot invert covariance matrix. Error %d.", info);
			omxRaiseError(oo->matrix->currentState, -1, errstr);
			return;
		}
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);

		sum += logDet + Q + (log(2 * M_PI) * smallRow->cols);
		if(OMX_DEBUG_ROWS) {Rprintf("Change in Total Likelihood is %3.3f + %3.3f + %3.3f = %3.3f, total Likelihood is %3.3f\n", logDet, Q, (log(2 * M_PI) * smallRow->cols), logDet + Q + (log(2 * M_PI) * smallRow->cols), sum);}
	}

	if(OMX_DEBUG) {Rprintf("Total Likelihood is %3.3f\n", sum);}
	oo->matrix->data[0] = sum;

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
	newObj->smallRow = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);
	newObj->smallCov = omxInitMatrix(NULL, newObj->cov->rows, newObj->cov->cols, TRUE, oo->matrix->currentState);
	newObj->RCX = omxInitMatrix(NULL, 1, newObj->data->cols, TRUE, oo->matrix->currentState);
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
		newObj->weights = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);
		newObj->smallWeights = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);
		newObj->smallMeans = omxInitMatrix(NULL, newObj->cov->cols, 1, TRUE, oo->matrix->currentState);
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		omxAliasMatrix(newObj->smallWeights, newObj->weights);
		newObj->corList = (double*) R_alloc(newObj->cov->cols * (newObj->cov->cols + 1) / 2, sizeof(double));
		newObj->smallCor = (double*) R_alloc(newObj->cov->cols * (newObj->cov->cols + 1) / 2, sizeof(double));
		newObj->lThresh = (double*) R_alloc(newObj->cov->cols, sizeof(double));
		newObj->uThresh = (double*) R_alloc(newObj->cov->cols, sizeof(double));
		newObj->Infin = (int*) R_alloc(newObj->cov->cols, sizeof(int));
		
		oo->objectiveFun = omxCallFIMLOrdinalObjective;
	}

	oo->argStruct = (void*) newObj;
}

#endif /* _OMX_FIML_OBJECTIVE_ */
