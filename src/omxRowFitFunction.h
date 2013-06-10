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

#ifndef _OMX_ROW_FITFUNCTION_
#define _OMX_ROW_FITFUNCTION_ TRUE

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

typedef struct omxRowFitFunction {

	/* Parts of the R  MxRowFitFunction Object */
	omxMatrix* rowAlgebra;		// Row-by-row algebra
	omxMatrix* rowResults;		// Aggregation of row algebra results
	omxMatrix* reduceAlgebra;	// Algebra performed after row-by-row computation
    omxMatrix* filteredDataRow; // Data row minus NAs
    omxMatrix* existenceVector; // Set of NAs
    omxMatrix* dataColumns;		// The order of columns in the data matrix

    /* Contiguous data note for contiguity speedup */
	omxContiguousData contiguous;		// Are the dataColumns contiguous within the data set

	/* Structures determined from info in the MxRowFitFunction Object*/
    double* oldDefs;            // The previous defVar vector.  To avoid recalculations where possible.         // NYI: This element currently unused.
	omxDefinitionVar* defVars;	// A list of definition variables
	int numDefs;				// The length of the defVars list
	omxMatrix* dataRow;         // One row of data, kept for aliasing only
	omxData*   data;			// The data

	int numDataRowDeps;         // number of algebra/matrix dependencies
	int *dataRowDeps;           // indices of algebra/matrix dependencies

} omxRowFitFunction;



void omxDestroyRowFitFunction(omxFitFunction *oo);

omxRListElement* omxSetFinalReturnsRowFitFunction(omxFitFunction *oo, int *numReturns);


void omxCopyMatrixToRow(omxMatrix* source, int row, omxMatrix* target);

void omxInitRowFitFunction(omxFitFunction* oo);





#endif /* _OMX_ROW_FITFUNCTION_ */
