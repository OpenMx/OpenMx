/*
 *  Copyright 2007-2010 The OpenMx Project
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
#include "omxObjectiveMetadata.h"

#ifndef _OMX_FIML_OBJECTIVE_
#define _OMX_FIML_OBJECTIVE_ TRUE

extern omxMatrix** matrixList;
extern void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*, double*, double*, double*, double*, int*);

/* FIML Computation Structures */
typedef struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
	omxData* source;		// Data source
	int numLocations;		// Num locations
	double** location;		// And where it goes
	omxMatrix** matrices;	// Matrix numbers for dirtying

} omxDefinitionVar;

typedef struct omxThresholdColumn {		 	// Threshold
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
	omxContiguousData contiguous;		// Are the dataColumns contiguous within the data set

//	double* zeros;

	/* Structures determined from info in the MxFIMLObjective Object*/
	omxDefinitionVar* defVars;	// A list of definition variables
	double* oldDefs;			// Stores definition variables between rows
	int numDefs;				// The length of the defVars list

	/* Reserved memory for faster calculation */
	omxMatrix* smallRow;		// Memory reserved for operations on each data row
	omxMatrix* smallCov;		// Memory reserved for operations on covariance matrix
	omxMatrix* smallMeans;		// Memory reserved for operations on the means matrix

	omxMatrix* RCX;				// Memory reserved for computationxs
		
	/* Structures for FIMLOrdinalObjective Objects */
	omxMatrix* cor;				// To calculate correlation matrix from covariance
	double* weights;			// Covariance weights to shift parameter estimates
	omxMatrix* smallThresh;		// Memory reserved for reduced threshold matrix
	omxThresholdColumn* thresholdCols;		// List of column thresholds
	
	/* Structures for JointFIMLObjective */
	omxMatrix* ordRow;		    // Memory reserved for ordinal data row
	omxMatrix* ordCov;	    	// Memory reserved for ordinal covariance matrix
	omxMatrix* ordMeans;		// Memory reserved for ordinal column means    
    omxMatrix* ordContCov;      // Memory reserved for ordinal/continuous covariance
	omxMatrix* halfCov;         // Memory reserved for computations
    omxMatrix* reduceCov;       // Memory reserved for computations
    
    

	/* Argument space for SADMVN function */
	double* lThresh;			// Specific list of lower thresholds
	double* uThresh;			// Specific list of upper thresholds
	double* corList;			// SADMVN-specific list of correlations
	double* smallCor;			// Reduced SADMVN-specific list of correlations
	int* Infin;					// Which thresholds to use
	int maxPts;					// From MxOptions (?)
	double absEps;				// From MxOptions
	double relEps;				// From MxOptions

	/* Space for inner sub-objective */
	void* subObjective;			// Inner Objective Object
	void (*covarianceMeansFunction)(void* subObjective, omxMatrix* cov, omxMatrix* means);
								// Inner Objective Function
	void (*destroySubObjective)(void* subObjective, omxObjectiveMetadataContainer* oomc);

} omxFIMLObjective;


/* FIML Function body */
void omxDestroyFIMLObjective(omxObjective *oo) {
	if(OMX_DEBUG) { Rprintf("Destroying FIML objective object.\n"); }
	omxFIMLObjective *argStruct = (omxFIMLObjective*) (oo->argStruct);

	if(argStruct->smallRow != NULL) omxFreeMatrixData(argStruct->smallRow);
	if(argStruct->smallCov != NULL) omxFreeMatrixData(argStruct->smallCov);
	if(argStruct->RCX != NULL)		omxFreeMatrixData(argStruct->RCX);
	if(argStruct->subObjective != NULL) {
		omxObjectiveMetadataContainer oomc = {argStruct->cov, argStruct->means,
			 argStruct->subObjective, argStruct->covarianceMeansFunction,
			 argStruct->destroySubObjective};
		argStruct->destroySubObjective(argStruct->subObjective, &oomc);
		argStruct->cov = oomc.cov;
		argStruct->means = oomc.means;
		argStruct->subObjective = oomc.subObjective;
		argStruct->covarianceMeansFunction = oomc.covarianceMeansFunction;
		argStruct->destroySubObjective = oomc.destroySubObjective;
	}
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

int handleDefinitionVarList(omxData* data, int row, omxDefinitionVar* defVars, double* oldDefs, int numDefs) {

	if(OMX_DEBUG_ROWS) { Rprintf("Processing Definition Vars.\n"); }
	
	int numVarsFilled = 0;

	/* Fill in Definition Var Estimates */
	for(int k = 0; k < numDefs; k++) {
		if(defVars[k].source != data) {
			error("Internal error: definition variable population into incorrect data source");
			continue; // don't populate from the wrong data frame
		}
		double newDefVar = omxDoubleDataElement(data, row, defVars[k].column);
		if(newDefVar == oldDefs[k]) continue;	// NOTE: Potential speedup vs accuracy tradeoff here using epsilon comparison
		oldDefs[k] = newDefVar;
		numVarsFilled++;

		for(int l = 0; l < defVars[k].numLocations; l++) {
			if(OMX_DEBUG_ROWS) {
				Rprintf("Populating column %d (value %3.2f) into matrix %d.\n", defVars[k].column, omxDoubleDataElement(defVars[k].source, row, defVars[k].column), defVars[k].matrices[l]);
			}
			*(defVars[k].location[l]) = newDefVar;
			omxMarkDirty(defVars[k].matrices[l]);
			if(ISNA(omxDoubleDataElement(data, row, defVars[k].column))) {
				error("Error: NA value for a definition variable is Not Yet Implemented.");
			}
		}
	}
	return numVarsFilled;
}

void omxCallFIMLOrdinalObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.
	/* TODO: Current implementation is slow: update by filtering correlations and thresholds. */
	if(OMX_DEBUG) { Rprintf("Beginning Ordinal FIML Evaluation.\n");}
	// Requires: Data, means, covariances, thresholds

	double sum;
	double Q = 0.0;
	double logDet = 0;
	int numDefs;
	int numCols, numRemoves = 0;
	int returnRowLikelihoods = 0;
	int keepCov = 0, keepInverse = 0;

	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns;
	omxMatrix *cor, *smallThresh;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *lThresh, *uThresh, maxPts, absEps, relEps, *corList, *smallCor, *weights, *oldDefs;
	int *Infin;
	omxDefinitionVar* defVars;
	int firstRow = 1;
	
	void* subObjective;
	void (*covMeans)(void* subObjective, omxMatrix* cov, omxMatrix* means);
	

	// Locals, for readability.  Compiler should cut through this.
	omxFIMLObjective* ofo = (omxFIMLObjective*)oo->argStruct;
	cov 		= ofo->cov;
	means		= ofo->means;
	smallRow 	= ofo->smallRow;
	smallCov 	= ofo->smallCov;
	smallMeans	= ofo->smallMeans;
	RCX 		= ofo->RCX;
	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	defVars		= ofo->defVars;
	oldDefs		= ofo->oldDefs;
	numDefs		= ofo->numDefs;

	cor 		= ofo->cor;
	corList 	= ofo->corList;
	smallCor	= ofo->smallCor;
	weights		= ofo->weights;
	smallThresh	= ofo->smallThresh;
	lThresh		= ofo->lThresh;
	uThresh		= ofo->uThresh;
	thresholdCols = ofo->thresholdCols;
	returnRowLikelihoods = ofo->returnRowLikelihoods;

	Infin		= ofo->Infin;
	maxPts		= ofo->maxPts;
	absEps		= ofo->absEps;
	relEps		= ofo->relEps;
	
	subObjective = ofo->subObjective;
	
	covMeans	= ofo->covarianceMeansFunction;

	if(numDefs == 0) {
		if(OMX_DEBUG_ALGEBRA) { Rprintf("No Definition Vars: precalculating."); }
		if(!(covMeans == NULL)) {
			covMeans(subObjective, cov, means);
		} else {
			omxRecompute(cov);			// Only recompute this here if there are no definition vars
			omxRecompute(means);
		}
		for(int j = 0; j < dataColumns->cols; j++) {
			if(thresholdCols[j].numThresholds > 0) { // Actually an ordinal column
				omxRecompute(thresholdCols[j].matrix);
				checkIncreasing(thresholdCols[j].matrix, thresholdCols[j].column);
			}
		}
		omxStandardizeCovMatrix(cov, corList, weights);	// Calculate correlation and covariance
	}

	sum = 0.0;
	int row = 0;

    while(row < data->rows) {
        if(OMX_DEBUG_ROWS) {Rprintf("Row %d.\n", row);}
        oo->matrix->currentState->currentRow = row;		// Set to a new row.
        logDet = 0.0;
        Q = 0.0;

        // Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
        numCols = 0;
        numRemoves = 0;

        // Handle Definition Variables.
        if(numDefs != 0) {
			if(keepCov <= 0) {  // If we're keeping covariance from the previous row, do not populate 
				if(OMX_DEBUG_ROWS) { Rprintf("Handling Definition Vars.\n"); }
				if(handleDefinitionVarList(data, row, defVars, oldDefs, numDefs) || firstRow) {
					// Use firstrow instead of rows == 0 for the case where the first row is all NAs
					// N.B. handling of definition var lists always happens, regardless of firstRow.
					if(!(covMeans == NULL)) {
						covMeans(subObjective, cov, means);
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
				double mean;
				if(means == NULL) mean = 0;
				else mean = omxVectorElement(means, j);
				double weight = weights[j];
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
			}
		}
		
		if(numRemoves >= smallCov->rows) {
			if(returnRowLikelihoods) {
			    omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 1.0);
			}
    		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
    		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
            row += 1;
    		keepCov -= 1;
    		keepInverse -= 1;
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

		if(OMX_DEBUG && !oo->matrix->currentState->currentRow) {
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

		int numIdentical = omxDataNumIdenticalRows(data, row);

		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS) { 
				Rprintf("Row %d likelihood is %3.3f.\n", row, likelihood);
			} 
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(oo->matrix, omxDataIndex(data, j), 0, likelihood);
			}
		} else {
			logDet = -2 * log(likelihood);
			logDet *= numIdentical;

			sum += logDet;// + (log(2 * M_PI) * (cov->cols - numRemoves));

			if(OMX_DEBUG_ROWS) { 
				Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n",
				    sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
            } 
        }
		if(firstRow) firstRow = 0;
		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);

		row += numIdentical;		// Step forward by the number of identical rows
		keepCov -= numIdentical;
		keepInverse -= numIdentical;
	}

    if(!returnRowLikelihoods) {
        if(OMX_DEBUG) {
            Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n",
                sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols - numRemoves)));
        }

        oo->matrix->data[0] = sum;
    }
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

	double sum = 0.0;
	double Q = 0.0;
	double logDet = 0;
	int numDefs;
	int numCols, numOrdRemoves = 0, numContRemoves=0;
	int returnRowLikelihoods = 0;
	int keepCov = 0, keepInverse = 0;

	omxMatrix *cov, *means, *smallRow, *smallCov, *smallMeans, *RCX, *dataColumns, *ordColumns;
	omxMatrix *cor, *smallThresh;
    omxMatrix *ordMeans, *ordCov, *ordRow;
    omxMatrix *halfCov, *reduceCov, *ordContCov;
	omxThresholdColumn *thresholdCols;
	omxData* data;
	double *lThresh, *uThresh, maxPts, absEps, relEps, *corList, *smallCor, *weights, *oldDefs;
	int *Infin;
	omxDefinitionVar* defVars;
	int firstRow = 1;
	
	void* subObjective;
	void (*covMeans)(void* subObjective, omxMatrix* cov, omxMatrix* means);
	

	// Locals, for readability.  Compiler should cut through this.
	omxFIMLObjective* ofo = (omxFIMLObjective*)oo->argStruct;
	cov 		= ofo->cov;
	means		= ofo->means;
	smallRow 	= ofo->smallRow;
	smallCov 	= ofo->smallCov;
	smallMeans	= ofo->smallMeans;
    ordMeans    = ofo->ordMeans;
    ordCov      = ofo->ordCov;
    ordRow      = ofo->ordRow;
    halfCov     = ofo->halfCov;
    reduceCov   = ofo->reduceCov;
    ordContCov  = ofo->ordContCov;
	RCX 		= ofo->RCX;
	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	defVars		= ofo->defVars;
	oldDefs		= ofo->oldDefs;
	numDefs		= ofo->numDefs;

	cor 		= ofo->cor;
	corList 	= ofo->corList;
	smallCor	= ofo->smallCor;
	weights		= ofo->weights;
	smallThresh	= ofo->smallThresh;
	lThresh		= ofo->lThresh;
	uThresh		= ofo->uThresh;
	thresholdCols = ofo->thresholdCols;
	returnRowLikelihoods = ofo->returnRowLikelihoods;

	Infin		= ofo->Infin;
	maxPts		= ofo->maxPts;
	absEps		= ofo->absEps;
	relEps		= ofo->relEps;
	
	subObjective = ofo->subObjective;
	
	covMeans	= ofo->covarianceMeansFunction;

	// if(numDefs == 0) {
	//         if(OMX_DEBUG_ALGEBRA) { Rprintf("No Definition Vars: precalculating."); }
	//         if(!(covMeans == NULL)) {
	//             covMeans(subObjective, cov, means);
	//         } else {
	//             omxRecompute(cov);          // Only recompute this here if there are no definition vars
	//             omxRecompute(means);
	//         }
	//         for(int j = 0; j < dataColumns->cols; j++) {
	//             int var = omxVectorElement(dataColumns, j);
	//             if(thresholdCols[var].numThresholds > 0) { // Actually an ordinal column
	//                 omxRecompute(thresholdCols[var].matrix);
	//                 checkIncreasing(thresholdCols[var].matrix, thresholdCols[var].column);
	//             }
	//         }
	//     }

    int row = 0;
    int ordRemove[cov->cols], contRemove[cov->cols];
    int zeros[cov->cols];
    char u = 'U';
    int info;
    double determinant;
    double oned = 1.0, zerod = 0.0, minusoned = -1.0;
    int onei = 1, zeroi = 0;

    while(row < data->rows) {
        if(OMX_DEBUG_ROWS) { Rprintf("Row %d.\n", row); } //:::DEBUG:::
        oo->matrix->currentState->currentRow = row;		// Set to a new row.
        logDet = 0.0;
        Q = 0.0;

        // Note:  This next bit really aught to be done using a matrix multiply.  Why isn't it?
        numCols = 0;
        numOrdRemoves = 0;
        numContRemoves = 0;

        // Handle Definition Variables.
        if(numDefs != 0) {
			if(keepCov <= 0) {  // If we're keeping covariance from the previous row, do not populate 
				if(OMX_DEBUG_ROWS) { Rprintf("Handling Definition Vars.\n"); }
				if(handleDefinitionVarList(data, row, defVars, oldDefs, numDefs) || firstRow || 1) { // TODO: Implement sorting-based speedups here.
					// Use firstrow instead of rows == 0 for the case where the first row is all NAs
					// N.B. handling of definition var lists always happens, regardless of firstRow.
					if(!(covMeans == NULL)) {
						covMeans(subObjective, cov, means);
					} else {
						omxRecompute(cov);
						omxRecompute(means);
					}
					for(int j=0; j < dataColumns->cols; j++) {
						int var = omxVectorElement(dataColumns, j);
						if(omxDataColumnIsFactor(data, j) && thresholdCols[var].numThresholds > 0) { // j is an ordinal column
							omxRecompute(thresholdCols[var].matrix); // Only one of these--save time by only doing this once
							checkIncreasing(thresholdCols[var].matrix, thresholdCols[var].column);
						}
					}
				}
			}
		}

        // TODO: Possible solution here: Manually record threshold column and index from data during this initial reduction step.  Since all the rest is algebras, it'll filter naturally.  Calculate offsets from continuous data, then dereference actual threshold values from the threshold matrix in its original state.  Requirement: colNum integer vector

		// Filter down correlation matrix and calculate thresholds
        // if(keepInverse <= 0 || keepCov <= 0 || firstRow) { // If defs and missingness don't change, skip. // TODO: Add in handling of sort speedup data
    		for(int j = 0; j < dataColumns->cols; j++) {
                zeros[j] = 0;
    			int var = omxVectorElement(dataColumns, j);
    			int value = omxIntDataElement(data, row, var);// Indexing correction means this is the index of the upper bound +1.
    			// TODO: Might save time by preseparating ordinal from continuous.
    			if(ISNA(value) || value == NA_INTEGER || value == NA_REAL) {  // Value is NA, therefore filter.
    				numOrdRemoves++;
                    numContRemoves++;
                    ordRemove[j] = 1;
                    contRemove[j] = 1;
    				// toRemove[j] = 1;
    				Infin[j] = -1;
    				if(OMX_DEBUG_ROWS) { Rprintf("Row %d, column %d.  Not a factor.\n", row, j);}
    				continue;
    			} else if(omxDataColumnIsFactor(data, j)) {             // Ordinal column.
                    numContRemoves++;
                    ordRemove[j] = 0;
                    contRemove[j] = 1;
    			    if(OMX_DEBUG_ROWS) { Rprintf("Row %d, column %d.  Thresholds for data column %d and row %d are %f -> %f. (Infin=%d)\n", row, j, var, value-1, lThresh[j], uThresh[j], Infin[j]);}
    			} else {
    			    numOrdRemoves++;
                    ordRemove[j] = 1;
                    contRemove[j] = 0;
    			}
    		}

    		if(numOrdRemoves >= dataColumns->cols && numContRemoves >=  dataColumns->cols) {
    		    // All elements missing.  Skip row.
    			if(returnRowLikelihoods) {
    			    omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 1.0);
    			}
        		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
        		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
                if(OMX_DEBUG) { Rprintf("All elements missing.  Skipping row."); } // WAS: OMX_DEBUG_ROWS
                row += 1;
        		keepCov -= 1;
        		keepInverse -= 1;
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
            omxResetAliasedMatrix(smallRow);				// Reset smallRow
            omxDataRow(data, row, dataColumns, smallRow);						        // Populate data row
            if(OMX_DEBUG) { // :::DEBUG:::
				omxPrint(smallRow, "Original Data elements");
			} // :::DEBUG:::
            omxResetAliasedMatrix(ordRow);                                              // Propagate to ordinal row
            omxRemoveRowsAndColumns(ordRow, 0, numOrdRemoves, zeros, ordRemove); 	    // Reduce the row to just ordinal.
    		omxRemoveRowsAndColumns(smallRow, 0, numContRemoves, zeros, contRemove); 	// Reduce the row to just continuous.
    		omxResetAliasedMatrix(smallMeans);
    		omxResetAliasedMatrix(ordMeans);
            omxRemoveRowsAndColumns(smallMeans, 0, numContRemoves, zeros, contRemove);
            omxRemoveRowsAndColumns(ordMeans, 0, numOrdRemoves, zeros, ordRemove); 	    // Reduce the row to just ordinal.
            

    		if(OMX_DEBUG_ROWS) { Rprintf("Keeper codes: inverse: %d, cov:%d, identical:%d\n", keepInverse, keepCov, omxDataNumIdenticalRows(data, row)); }

			omxResetAliasedMatrix(smallCov);				// Re-sample covariance matrix
			omxRemoveRowsAndColumns(smallCov, numContRemoves, numContRemoves, contRemove, contRemove);
			omxResetAliasedMatrix(ordCov);				// Re-sample covariance matrix for ordinal
			omxRemoveRowsAndColumns(ordCov, numOrdRemoves, numOrdRemoves, ordRemove, ordRemove);
			omxResetAliasedMatrix(ordContCov);				// Re-sample covariance between ordinal and continuous
			omxRemoveRowsAndColumns(ordContCov, numContRemoves, numOrdRemoves, contRemove, ordRemove);

            /* :::DEBUG::: */
            if(OMX_DEBUG) { Rprintf("Removed %d continuous and %d ordinal cols from length %d(%d) data row.\n", numContRemoves, numOrdRemoves, dataColumns->cols, cov->cols);}
			if(OMX_DEBUG) { omxPrint(cov, "Original Covariance Matrix"); }
			if(OMX_DEBUG) { 
				omxPrint(smallCov, "Continuous Covariance Matrix"); 
				}
            if(OMX_DEBUG) { 
				omxPrint(smallRow, "Continuous elements");
				}
			if(OMX_DEBUG) { omxPrint(ordCov, "Ordinal Covariance Matrix"); }
			if(OMX_DEBUG) { 
				omxPrint(ordRow, "Ordinal elements");
				}
			if(OMX_DEBUG) { 
				omxPrint(ordContCov, "Ordinal/Continuous Covariance Matrix");
				 }
            /* :::DEBUG::: */
            
            // if(smallCov->cols < 1) {     // TODO: Implement catch for all-continuous-missing case

			/* Calculate derminant and inverse of Censored continuousCov matrix */
			// TODO : Speed this up.
			F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);

			if(info != 0) {
				if(!returnRowLikelihoods) {
					char helperstr[200];
					char *errstr = calloc(250, sizeof(char));
					sprintf(helperstr, "Expected covariance matrix is not positive-definite in data row %d", omxDataIndex(data, row));
					if(oo->matrix->currentState->computeCount <= 0) {
						sprintf(errstr, "%s at starting values.\n", helperstr);
					} else {
						sprintf(errstr, "%s at iteration %d.\n", helperstr, oo->matrix->currentState->majorIteration);
					}
					omxRaiseError(oo->matrix->currentState, -1, errstr);
					free(errstr);
					return;
				} else {
					omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 0.0);
            		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
            		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
                    if(OMX_DEBUG) {Rprintf("Non-positive-definite covariance matrix in row likelihood.  Skipping Row.");}
                    row += 1;
            		keepCov -= 1;
            		keepInverse -= 1;
					continue;
				}
			}

			// Calculate determinant: squared product of the diagonal of the decomposition
			determinant = 1.0;
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
					omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 0.0);
            		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
            		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
                    // Rprintf("Incrementing Row."); //:::DEBUG:::
                    row += 1;
            		keepCov -= 1;
            		keepInverse -= 1;
					continue;
				}
			}
        // }
        
		/* Calculate Row Likelihood */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);                          // RCX is the continuous-column mahalanobis distance.
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei); //Q is the total mahalanobis distance
	
		if(OMX_DEBUG) {Rprintf("Continuous, row %d: %3.3f %3.3f (%d cols)", row, determinant, Q, smallRow->cols);} // :::DEBUG:::
		
		// Reserve: 1) Inverse continuous covariance (smallCov)
		//          2) Columnwise Mahalanobis distance (contCov^-1)%*%(Data - Means) (RCX)
		//          3) Overall Mahalanobis distance (FIML likelihood of data) (Q)
		//Calculate:4) Cont/ord covariance %*% Mahalanobis distance  (halfCov)
		//          5) ordCov <- ordCov - Cont/ord covariance %*% Inverse continuous cov
		
		// TODO: Make this less of a hack.
        halfCov->rows = smallCov->rows;
        halfCov->cols = ordContCov->cols;
        omxMatrixCompute(halfCov);
        reduceCov->rows = ordContCov->rows;
        reduceCov->cols = ordContCov->rows;
        omxMatrixCompute(reduceCov);

        F77_CALL(dsymv)(&u, &(ordContCov->rows), &oned, ordContCov->data, (&ordContCov->cols), RCX->data, &onei, &zerod, RCX->data, &onei);                             // RCX is the influence of the continuous on the thresholds
        F77_CALL(dgemm)((smallCov->majority), (ordContCov->majority), &(smallCov->rows), &(smallCov->cols), &(ordContCov->rows), &oned, smallCov->data, &(smallCov->leading), ordContCov->data, &(ordContCov->leading), &zerod, halfCov->data, &(halfCov->leading));          // halfCov is inverse continuous %*% cont/ord covariance
        F77_CALL(dgemm)((ordContCov->minority), (halfCov->majority), &(ordContCov->rows), &(ordContCov->cols), &(halfCov->rows), &oned, ordContCov->data, &(ordContCov->leading), halfCov->data, &(halfCov->leading), &zerod, reduceCov->data, &(reduceCov->leading));      // reduceCov is  cont/ord %&% contCov^-1
        int vlen = reduceCov->rows * reduceCov->cols;
        // FIXME: This assumes that ordCov and reducCov have the same row/column majority.
        F77_CALL(daxpy)(&vlen, &minusoned, reduceCov->data, &onei, ordCov->data, &onei); // ordCov <- (ordCov - reduceCov)
        F77_CALL(dgemv)((smallCov->minority), &(halfCov->rows), &(halfCov->cols), &oned, halfCov->data, &(halfCov->leading), smallRow->data, &onei, &oned, ordMeans->data, &onei);                      // ordMeans += halfCov %*% contRow

		// TODO: Implement all-ordinal-missing case
        // if(numOrdRemoves < dataColumns->cols) {    // Ordinal all missing.
        //             likelihood = 1;
        // 
        //         } else {

		    // Calculate correlation matrix from covariance
		    omxStandardizeCovMatrix(ordCov, corList, weights);

            int count = 0;
    		for(int j = 0; j < dataColumns->cols; j++) {
                if(ordRemove[j]) continue;         // NA or non-ordinal
    		    // FIXME: Once we connect ordinal columns to thresholds, use var to correct indexing.
                // int var = omxVectorElement(ordColumns, j);
    			int value = omxIntDataElement(data, row, j); //  TODO: Compare with extraction from dataRow.
                // Rprintf("Row %d, Column %d, value %d+1\n", row, j, value); // :::DEBUG:::
    	        value--;		// Correct for C indexing: value is now the index of the upper bound.
                // Rprintf("Row %d, Column %d, value %d+1\n", row, j, value); // :::DEBUG:::
    			double offset;
    			if(means == NULL) offset = 0;
    			else offset = omxVectorElement(ordMeans, count);
    			double weight = weights[count];
                offset += omxVectorElement(RCX, j);          // Offset adjustment now covers mahalnobis adjustment as well
    			if(value == 0) { 									// Lowest threshold = -Inf
    			    // Rprintf("0 %d = %d, %x, %d, %3.3f, %3.3f.\n", j, count, thresholdCols[j].matrix, thresholdCols[j].column, offset, weight); //:::DEBUG::::
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

    			if(uThresh[count] == NA_INTEGER || isnan(uThresh[count])) { // for matrix-style specification.
    				uThresh[count] = lThresh[count];
    				Infin[count] = 1;
    			}
    			if(OMX_DEBUG) { Rprintf("Row %d, column %d.  Thresholds for data column %d and threshold column %d are %f -> %f. (Infin=%d)\n", row, count, j, value, lThresh[count], uThresh[count], Infin[count]);} // :::DEBUG:::
                count++;
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
    		int numVars = ordCov->rows;
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

    		if(OMX_DEBUG && !oo->matrix->currentState->currentRow) {
    			char infinCodes[3][20];
    			strcpy(infinCodes[0], "(-INF, upper]");
    			strcpy(infinCodes[1], "[lower, INF)");
    			strcpy(infinCodes[2], "[lower, upper]");
    			Rprintf("Input to sadmvn is (%d rows):\n", numVars); //:::DEBUG:::

    			omxPrint(ordCov, "Ordinal Covariance Matrix"); //:::DEBUG:::

    			for(int i = 0; i < numVars; i++) {
    				Rprintf("Row %d: %f, %f, %d(%s)\n", i, lThresh[i], uThresh[i], Infin[i], infinCodes[Infin[i]]);
    			}

    			Rprintf("Cor: (Lower %d x %d):", cov->rows, cov->cols); //:::DEBUG:::
    			for(int i = 0; i < cov->rows*(cov->rows-1)/2; i++) {
    				// Rprintf("Row %d of Cor: ", i);
    				// for(int j = 0; j < i; j++)
    				Rprintf(" %f", corList[i]); // (i*(i-1)/2) + j]);
    				// Rprintf("\n");
    			}
    			Rprintf("\n");
    		}

    		if(OMX_DEBUG) {
				Rprintf("Output of sadmvn is %f, %f, %d.\n", Error, likelihood, inform); 
			} 

    		if(inform == 2) {
    			error("Improper input to sadmvn.");
    		}
        // }
        
		int numIdentical = omxDataNumIdenticalRows(data, row);

		if(returnRowLikelihoods) {
		    if(OMX_DEBUG_ROWS) {Rprintf("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f\n", pow(2 * M_PI, -.5 * smallRow->cols), (1.0/sqrt(determinant)), exp(-.5 * Q), pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q));}
			sum = pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q) * likelihood;
            
			if(OMX_DEBUG_ROWS) {Rprintf("Row %d likelihood is %3.3f.\n", row, likelihood);}
			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(oo->matrix, omxDataIndex(data, j), 0, sum);
			}
		} else {
			logDet = -2 * log(likelihood);       // -2 Log of ordinal likelihood
            logDet += (log(determinant) + Q + (log(2 * M_PI) * smallRow->cols));    // -2 Log of continuous likelihood
            // logDet *= numIdentical;

            sum += logDet;
			
			if(OMX_DEBUG_ROWS) { 
				Rprintf("Change in Total log Likelihood for row %d is %3.3f + %3.3f + %3.3f + %3.3f= %3.3f, total Likelihood is %3.3f\n", oo->matrix->currentState->currentRow, log(determinant), Q, (log(2 * M_PI) * smallRow->cols), -2  * log(likelihood), log(determinant) + Q + (log(2 * M_PI) * smallRow->cols), sum);
			} 

			if(OMX_DEBUG_ROWS) {
				Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change \n",
				    sum, logDet);
            }
        }
		if(firstRow) firstRow = 0;
		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
        // Rprintf("Incrementing Row."); //:::DEBUG:::
		row += numIdentical;		// Step forward by the number of identical rows
		keepCov -= numIdentical;
		keepInverse -= numIdentical;
	}

    if(!returnRowLikelihoods) {
        if(OMX_DEBUG) {
            Rprintf("Total over all rows is %3.3f. -2 Log Likelihood this row is %3.3f, total change %3.3f\n",
                sum, logDet, logDet + Q + (log(2 * M_PI) * (cov->cols)));
        }

        oo->matrix->data[0] = sum;
    }
}

void omxCallFIMLObjective(omxObjective *oo) {	// TODO: Figure out how to give access to other per-iteration structures.

	if(OMX_DEBUG) { Rprintf("Beginning FIML Evaluation.\n"); }
	// Requires: Data, means, covariances.
	// Potential Problem: Definition variables currently are assumed to be at the end of the data matrix.

	double sum;
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
	int numCols, numRemoves;
	int returnRowLikelihoods;
	int keepCov = 0, keepInverse = 0;
	
	void* subObjective;
	
	void (*covMeans)(void* subObjective, omxMatrix* cov, omxMatrix* means);
	
	omxMatrix *cov, *means, *smallRow, *smallCov, *RCX, *dataColumns;//, *oldInverse;
	omxDefinitionVar* defVars;
	omxData *data;

    omxFIMLObjective* ofo = ((omxFIMLObjective*)oo->argStruct);

	cov 		= ofo->cov;		// Locals, for readability.  Should compile out.
	means		= ofo->means;
	smallRow 	= ofo->smallRow;
	smallCov 	= ofo->smallCov;
	oldDefs		= ofo->oldDefs;
	RCX 		= ofo->RCX;
	data		= ofo->data;
	dataColumns	= ofo->dataColumns;
	defVars		= ofo->defVars;
	numDefs		= ofo->numDefs;
	returnRowLikelihoods = ofo->returnRowLikelihoods;
	isContiguous    = ofo->contiguous.isContiguous;
	contiguousStart = ofo->contiguous.start;
	contiguousLength = ofo->contiguous.length;

	subObjective = ofo->subObjective;

	covMeans	= ofo->covarianceMeansFunction;

	if(numDefs == 0) {
		if(OMX_DEBUG) {Rprintf("Precalculating cov and means for all rows.\n");}
		if(!(covMeans == NULL)) {
			covMeans(subObjective, cov, means);
		} else {
			omxRecompute(cov);			// Only recompute this here if there are no definition vars
			omxRecompute(means);
		}
		if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
		if(OMX_DEBUG) { omxPrintMatrix(means, "Means"); }
	}
	int toRemove[cov->cols];
	int zeros[cov->cols];

	sum = 0.0;
	int firstRow = 1;
	int row = 0;

	while(row < data->rows) {
        // Rprintf("Row %d.\n", row); //:::DEBUG:::
		oo->matrix->currentState->currentRow = row;		// Set to a new row.
		Q = 0.0;

		numCols = 0;
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
				if(OMX_DEBUG_ROWS) { Rprintf("Handling Definition Vars.\n"); }
				if(handleDefinitionVarList(data, row, defVars, oldDefs, numDefs) || firstRow) {
				// Use firstrow instead of rows == 0 for the case where the first row is all NAs
				// N.B. handling of definition var lists always happens, regardless of firstRow.
					if(!(covMeans == NULL)) {
						covMeans(subObjective, cov, means);
					} else {
						omxRecompute(cov);
						omxRecompute(means);
					}
				}
			} else if(OMX_DEBUG_ROWS){ Rprintf("Identical def vars: Not repopulating"); }
		}

		if(OMX_DEBUG_ROWS) { omxPrint(means, "Local Means"); }
		if(OMX_DEBUG_ROWS) {
			char note[50];
			sprintf(note, "Local Data Row %d", row);
			omxPrint(smallRow, note); 
		}
		
		/* Censor row and censor and invert cov. matrix. */
		// Determine how many rows/cols to remove.
		memset(zeros, 0, sizeof(int) * dataColumns->cols);
		memset(toRemove, 0, sizeof(int) * dataColumns->cols);
		for(int j = 0; j < dataColumns->cols; j++) {
			double dataValue = omxMatrixElement(smallRow, 0, j);
			if(isnan(dataValue) || dataValue == NA_REAL) {
				numRemoves++;
				toRemove[j] = 1;
			} else if(means != NULL) {
				omxSetMatrixElement(smallRow, 0, j, (dataValue -  omxVectorElement(means, j)));
			}
		}
		numCols = dataColumns->cols - numRemoves;

		if(cov->cols <= numRemoves) {
			if(returnRowLikelihoods) {
				omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 1);
			}
			if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
    		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
            // Rprintf("Incrementing Row."); //:::DEBUG:::
    		row += 1;
    		keepCov -= 1;
    		keepInverse -= 1;
			continue;
		}
		
		omxRemoveRowsAndColumns(smallRow, 0, numRemoves, zeros, toRemove); 	// Reduce it.
		
		if(OMX_DEBUG_ROWS) { Rprintf("Keeper codes: inverse: %d, cov:%d, identical:%d\n", keepInverse, keepCov, omxDataNumIdenticalRows(data, row)); }

		if(keepInverse <= 0 || keepCov <= 0 || firstRow) { // If defs and missingness don't change, skip.
			omxResetAliasedMatrix(smallCov);				// Re-sample covariance matrix
			omxRemoveRowsAndColumns(smallCov, numRemoves, numRemoves, toRemove, toRemove);

			if(OMX_DEBUG_ROWS) { omxPrint(smallCov, "Local Covariance Matrix"); }

			/* Calculate derminant and inverse of Censored Cov matrix */
			// TODO : Speed this up.
			F77_CALL(dpotrf)(&u, &(smallCov->rows), smallCov->data, &(smallCov->cols), &info);
			if(info != 0) {
				if(!returnRowLikelihoods) {
					char helperstr[200];
					char *errstr = calloc(250, sizeof(char));
					sprintf(helperstr, "Expected covariance matrix is not positive-definite in data row %d", omxDataIndex(data, row));
					if(oo->matrix->currentState->computeCount <= 0) {
						sprintf(errstr, "%s at starting values.\n", helperstr);
					} else {
						sprintf(errstr, "%s at iteration %d.\n", helperstr, oo->matrix->currentState->majorIteration);
					}
					omxRaiseError(oo->matrix->currentState, -1, errstr);
					free(errstr);
					return;
				} else {
					omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 0.0);
            		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
            		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
                    // Rprintf("Incrementing Row."); //:::DEBUG:::
                    row += 1;
            		keepCov -= 1;
            		keepInverse -= 1;
					continue;
				}
			}
			
			// Calculate determinant: squared product of the diagonal of the decomposition
			determinant = 1.0;
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
					omxSetMatrixElement(oo->matrix, omxDataIndex(data, row), 0, 0.0);
            		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
            		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);
                    // Rprintf("Incrementing Row."); //:::DEBUG:::
                    row += 1;
            		keepCov -= 1;
            		keepInverse -= 1;
					continue;
				}
			}
		}

		/* Calculate Row Likelihood */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), smallRow->data, &onei, &zerod, RCX->data, &onei);
		Q = F77_CALL(ddot)(&(smallRow->cols), smallRow->data, &onei, RCX->data, &onei);
		int numIdentical = omxDataNumIdenticalRows(data, row);

		if(returnRowLikelihoods) {
			if(OMX_DEBUG_ROWS) {Rprintf("Change in Total Likelihood is %3.3f * %3.3f * %3.3f = %3.3f\n", pow(2 * M_PI, -.5 * smallRow->cols), (1.0/sqrt(determinant)), exp(-.5 * Q), pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q));}
			sum = pow(2 * M_PI, -.5 * smallRow->cols) * (1.0/sqrt(determinant)) * exp(-.5 * Q);

			for(int j = numIdentical + row - 1; j >= row; j--) {  // Populate each successive identical row
				omxSetMatrixElement(oo->matrix, omxDataIndex(data, j), 0, sum);
			}
		} else {
			sum += (log(determinant) + Q + (log(2 * M_PI) * smallRow->cols)) * numIdentical;
			if(OMX_DEBUG_ROWS) {
				Rprintf("Change in Total Likelihood for row %d is %3.3f + %3.3f + %3.3f = %3.3f, total Likelihood is %3.3f\n", oo->matrix->currentState->currentRow, log(determinant), Q, (log(2 * M_PI) * smallRow->cols), log(determinant) + Q + (log(2 * M_PI) * smallRow->cols), sum);
			}
		}
		if(firstRow) firstRow = 0;
		if(keepCov <= 0) keepCov = omxDataNumIdenticalDefs(data, row);
		if(keepInverse  <= 0) keepInverse = omxDataNumIdenticalMissingness(data, row);

		row += numIdentical;
		keepCov -= numIdentical;
		keepInverse -= numIdentical;
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

	SEXP nextMatrix, itemList, nextItem, dataSource, columnSource, threshMatrix, objectiveClass;
    SEXP levelList;
	int nextDef, index, numOrdinal = 0, numContinuous = 0, numCols;
	int *nextInt;
	omxFIMLObjective *newObj = (omxFIMLObjective*) R_alloc(1, sizeof(omxFIMLObjective));
	
	newObj->subObjective = NULL;
	
	/* Set default Objective calls to FIML Objective Calls */
	oo->objectiveFun = omxCallFIMLObjective;
	oo->needsUpdateFun = omxNeedsUpdateFIMLObjective;
	oo->setFinalReturns = omxSetFinalReturnsFIMLObjective;
	oo->destructFun = omxDestroyFIMLObjective;
	oo->repopulateFun = NULL;
	
	PROTECT(nextMatrix = GET_SLOT(rObj, install("metadata")));
	if(IS_S4_OBJECT(nextMatrix)) {
		if(OMX_DEBUG) { Rprintf("Initializing subobjective metadata.\n"); }
		// Metadata exists and must be processed.
		PROTECT(objectiveClass = STRING_ELT(getAttrib(nextMatrix, install("class")), 0));
		char const * objType;
		objType = CHAR(objectiveClass);
		
		if(OMX_DEBUG) { Rprintf("Subobjective metadata is class %s.\n", objType); }
		
		omxObjectiveMetadataContainer oomc = {NULL, NULL, NULL, NULL, NULL};
		omxObjectiveMetadataContainer *poomc = &oomc;
		
		for(int i = 0; i < numObjectiveMetadatas; i++) {
			if(strncmp(objType, omxObjectiveMetadataTable[i].name, 50) == 0) {
				omxObjectiveMetadataTable[i].initFunction(nextMatrix, poomc, oo->matrix->currentState);
				break;
			}
		}
		
		if(oomc.cov != NULL) {
			newObj->cov = oomc.cov;
			newObj->means = oomc.means;
			newObj->subObjective = oomc.subObjective;
			newObj->covarianceMeansFunction = oomc.covarianceMeansFunction;
			newObj->destroySubObjective = oomc.destroySubObjective;
		}
		
		UNPROTECT(1); // objectiveClass

		if(newObj->subObjective == NULL) {
			error("Unrecognized subObjective Metadata passed to back-end.");
		}
	} else {
		// No metadata.  Continue as usual.
		newObj->subObjective = NULL;
		newObj->covarianceMeansFunction = NULL;

		PROTECT(nextMatrix = GET_SLOT(rObj, install("means")));
		newObj->means = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
		if(newObj->means == NULL) { error("No means model in FIML evaluation.");}
		UNPROTECT(1);	// UNPROTECT(means)

		PROTECT(nextMatrix = GET_SLOT(rObj, install("covariance")));
		newObj->cov = omxNewMatrixFromMxIndex(nextMatrix, oo->matrix->currentState);
		UNPROTECT(1);	// UNPROTECT(covariance)
	}
	UNPROTECT(1);	// UNPROTECT(metadata)
	
	numCols = newObj->cov->rows;

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
	newObj->dataColumns = omxNewMatrixFromRPrimitive(nextMatrix, oo->matrix->currentState);
	if(OMX_DEBUG) { omxPrint(newObj->dataColumns, "Variable mapping"); }
	UNPROTECT(1);

	if(OMX_DEBUG) {Rprintf("Accessing Threshold matrix.\n"); }
	PROTECT(threshMatrix = GET_SLOT(rObj, install("thresholds")));
    
    if(INTEGER(threshMatrix)[0] != NA_INTEGER) {
        if(OMX_DEBUG) {Rprintf("Accessing Threshold Mappings.\n"); }
        
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
    			if(OMX_DEBUG) {Rprintf("Column %d is continuous.\n", index);}
    			newObj->thresholdCols[index].matrix = NULL;
    			newObj->thresholdCols[index].column = 0;
    			newObj->thresholdCols[index].numThresholds = 0;
                numContinuous++;
    		} else {
    			newObj->thresholdCols[index].matrix = omxNewMatrixFromMxIndex(threshMatrix, 
    				oo->matrix->currentState);
    			newObj->thresholdCols[index].column = thresholdColumn[index];
    			newObj->thresholdCols[index].numThresholds = thresholdNumber[index];
    			if(OMX_DEBUG) {
    				Rprintf("Column %d is ordinal with %d thresholds in threshold column %d.\n", 
    				    thresholdColumn[index], thresholdNumber[index]);
    			}
    			numOrdinal++;
    		}
    	}
    	if(OMX_DEBUG) {Rprintf("%d threshold columns processed.\n", numOrdinal);}
    	UNPROTECT(2); /* nextMatrix and itemList ("thresholds" and "thresholdColumns") */
    } else {
        if (OMX_DEBUG) Rprintf("No thresholds matrix; not processing thresholds.");
        numContinuous = newObj->dataColumns->rows;
        newObj->thresholdCols = NULL;
        numOrdinal = 0;
    }
    UNPROTECT(1); /* threshMatrix */

	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, newObj->dataColumns);

	if(OMX_DEBUG) {Rprintf("Accessing definition variables structure.\n"); }
	PROTECT(nextMatrix = GET_SLOT(rObj, install("definitionVars")));
	newObj->numDefs = length(nextMatrix);
	if(OMX_DEBUG) {Rprintf("Number of definition variables is %d.\n", newObj->numDefs); }
	newObj->defVars = (omxDefinitionVar *) R_alloc(newObj->numDefs, sizeof(omxDefinitionVar));
	newObj->oldDefs = (double *) R_alloc(newObj->numDefs, sizeof(double));		// Storage for Def Vars
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
		newObj->oldDefs[nextDef] = NA_REAL;					// Def Vars default to NA
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
    // int ordCols = omxDataNumFactor(newObj->data);        // Unneeded, since we don't use it.
    // int contCols = omxDataNumNumeric(newObj->data);
	newObj->smallRow = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
	newObj->smallCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
	newObj->RCX = omxInitMatrix(NULL, 1, covCols, TRUE, oo->matrix->currentState);
//	newObj->zeros = omxInitMatrix(NULL, 1, newObj->cov->cols, TRUE, oo->matrix->currentState);

	omxAliasMatrix(newObj->smallCov, newObj->cov);					// Will keep its aliased state from here on.

	if(numOrdinal > 0 && numContinuous <= 0) {
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
	} else if(numOrdinal > 0) {
		if(OMX_DEBUG) { Rprintf("Ordinal and Continuous Data detected.  Using Joint Ordinal/Continuous FIML."); }
		newObj->weights = (double*) R_alloc(covCols, sizeof(double));
		newObj->smallMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		newObj->ordCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		newObj->ordMeans = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->ordRow = omxInitMatrix(NULL, covCols, 1, TRUE, oo->matrix->currentState);
        newObj->ordContCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
        newObj->halfCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(NULL, covCols, covCols, TRUE, oo->matrix->currentState);
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		omxAliasMatrix(newObj->ordMeans, newObj->means);
		omxAliasMatrix(newObj->ordCov, newObj->cov);
		omxAliasMatrix(newObj->ordContCov, newObj->cov);
		omxAliasMatrix(newObj->ordRow, newObj->smallRow );
		omxAliasMatrix(newObj->smallMeans, newObj->means);
		omxAliasMatrix(newObj->ordMeans, newObj->means);
		newObj->corList = (double*) R_alloc(covCols * (covCols + 1) / 2, sizeof(double));
		newObj->lThresh = (double*) R_alloc(covCols, sizeof(double));
		newObj->uThresh = (double*) R_alloc(covCols, sizeof(double));
		newObj->Infin = (int*) R_alloc(covCols, sizeof(int));

		oo->objectiveFun = omxCallJointFIMLObjective;
	}

	oo->argStruct = (void*) newObj;
	if(OMX_DEBUG) { Rprintf("FIML Initialization Completed."); }
}

#endif /* _OMX_FIML_OBJECTIVE_ */
