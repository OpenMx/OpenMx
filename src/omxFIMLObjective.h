/*
 *  Copyright 2007-2011 The OpenMx Project
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
 *
 */
 
#ifndef _OMXFIMLOBJECTIVE_H_
#define _OMXFIMLOBJECTIVE_H_

/* FIML Computation Structures */
typedef struct omxDefinitionVar {		 	// Definition Var

	int data, column;		// Where it comes from
	omxData* source;		// Data source
	int numLocations;		// Num locations
	int* rows;				// row positions
	int* cols;				// column positions
	int* matrices;			// matrix numbers

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
    omxMatrix* rowLikelihoods;  // The row-by-row likelihoods
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

} omxFIMLObjective;

int handleDefinitionVarList(omxData* data, omxState *state, int row, 
	omxDefinitionVar* defVars, double* oldDefs, int numDefs);

void omxCallFIMLObjective(omxObjective *oo);
void omxCallFIMLOrdinalObjective(omxObjective *oo);
unsigned short int omxNeedsUpdateFIMLObjective(omxObjective* oo);
omxRListElement* omxSetFinalReturnsFIMLObjective(omxObjective *oo, int *numReturns);
void omxDestroyFIMLObjective(omxObjective *oo);
void omxPopulateFIMLObjective(omxObjective *oo, SEXP algebra);
void omxInitFIMLObjective(omxObjective* oo, SEXP rObj);
void omxCreateFIMLObjective(omxObjective* oo, SEXP rObj, omxMatrix* cov, omxMatrix* means);


#endif /* _OMXFIMLOBJECTIVE_H_ */
