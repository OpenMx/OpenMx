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

/***********************************************************
*
*  omxAlgebraFunctions.c
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Includes the functions required for omxAlgebra statements.
*   These functions should take a number of values that
*   evenly matches the number of args requested by the
*   omxSymbolTable.
*
**********************************************************/

#include "omxMatrix.h"
#include "merge.h"
#include "omxBLAS.h"
#include "omxOpenmpWrap.h"
#include "omxSadmvnWrapper.h"
#include "matrix.h"

// TODO: Implement wrappers for BLAS functions used here.

/* omxAlgebraFunction Wrappers */

static void omxMatrixTranspose(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	omxCopyMatrix(result, inMat);
	result->colMajor = !result->colMajor;
	int rowtemp = result->rows;
	result->rows = result->cols;
	result->cols = rowtemp;
	result->transposePopulate();
	omxMatrixLeadingLagging(result);
}

static void omxMatrixInvert(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxCopyMatrix(result, inMat);

	Matrix resultMat(result);
	int info = MatrixInvert1(result);
	if (info) {
		result->data[0] = nan("singular");
		// recordIterationError not available here
	}
}

static int BroadcastIndex = 0;

static void nameBroadcastAlg(omxMatrix *bc)
{
	std::string str = string_snprintf("broadcast%03d", ++BroadcastIndex);
	SEXP name;
	Rf_protect(name = Rf_mkChar(str.c_str()));
	bc->name = CHAR(name);
}

static void ensureElemConform(const char *op, omxMatrix **matList, omxMatrix *result)
{
	omxMatrix *mat0 = matList[0];
	omxMatrix *mat1 = matList[1];

	// This is for backward compatibility with opcodes that
	// do not do anything special for FF_COMPUTE_DIMS.
	if (mat0->cols == 0 || mat0->rows == 0) {
		omxRecompute(mat0, FF_COMPUTE_INITIAL_FIT, NULL);
	}
	if (mat1->cols == 0 || mat1->rows == 0) {
		omxRecompute(mat1, FF_COMPUTE_INITIAL_FIT, NULL);
	}

	if (mat0->cols == mat1->cols && mat0->rows == mat1->rows) {
		if(OMX_DEBUG_ALGEBRA) { 
			mxLog("Resize %s to %dx%d", result->name, mat0->rows, mat0->cols);
		}
		omxResizeMatrix(result, mat0->rows, mat0->cols);
		return;
	}

	if (mat0->cols == 1 && mat0->rows == 1) {
		omxResizeMatrix(result, mat1->rows, mat1->cols);
		omxMatrix* om = omxInitMatrix(mat1->rows, mat1->cols, TRUE, result->currentState);
		nameBroadcastAlg(om);
		omxAlgebra *oa = new omxAlgebra;
		omxInitAlgebraWithMatrix(oa, om);
		omxAlgebraAllocArgs(oa, 1);
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[62]); // broadcast
		omxFillAlgebraFromTableEntry(oa, entry, 1);
		oa->algArgs[0] = mat0;
		matList[0] = om;
		return;
	}
	if (mat1->cols == 1 && mat1->rows == 1) {
		omxResizeMatrix(result, mat0->rows, mat0->cols);
		omxMatrix* om = omxInitMatrix(mat0->rows, mat0->cols, TRUE, result->currentState);
		nameBroadcastAlg(om);
		omxAlgebra *oa = new omxAlgebra;
		omxInitAlgebraWithMatrix(oa, om);
		omxAlgebraAllocArgs(oa, 1);
		const omxAlgebraTableEntry* entry = &(omxAlgebraSymbolTable[62]); // broadcast
		omxFillAlgebraFromTableEntry(oa, entry, 1);
		oa->algArgs[0] = mat1;
		matList[1] = om;
		return;
	}

	Rf_error("Matrices %s and %s are non-conformable in %s; rows %d != %d or cols %d != %d",
		 mat0->name, mat1->name, op, mat0->rows, mat1->rows, mat0->cols, mat1->cols);
}

static void omxBroadcast(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) return;

	omxMatrix *src = matList[0];

	if (src->rows != 1 || src->cols != 1) {
		Rf_error("Don't know how to broadcast from a non 1x1 source matrix");
	}

	int size = result->rows * result->cols;
	for (int dx=0; dx < size; ++dx) {
		result->data[dx] = src->data[0];
	}
}

static void omxMatrixMult(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];

	if(preMul == NULL || postMul == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Null matrix pointer detected.\n");
		free(errstr);
		return;
	}

	/* Conformability Check! */
	if(preMul->cols != postMul->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-conformable matrices [(%d x %d) and (%d x %d)] in Matrix Multiply.", preMul->rows, preMul->cols, postMul->rows, postMul->cols);
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if(result->rows != preMul->rows || result->cols != postMul->cols)
		omxResizeMatrix(result, preMul->rows, postMul->cols);

	omxDGEMM(FALSE, FALSE, 1.0, preMul, postMul, 0.0, result);

	result->colMajor = TRUE;

	omxMatrixLeadingLagging(result);
}

static void omxElementPower(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("elementwise power", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
		for(int i = 0; i < size; i++) {
			omxSetVectorElement(result, i,
				pow(omxVectorElement(first, i),
					omxVectorElement(second, i)));
		}
		result->colMajor = first->colMajor;
		omxMatrixLeadingLagging(result);
	} else {
		for(int i = 0; i < rows; i++) {
			for(int j = 0; j < cols; j++) {
				omxSetMatrixElement(result, i, j,
					pow(omxMatrixElement(first, i, j),
						omxMatrixElement(second, i, j)));
			}
		}
	}
}

static void omxMatrixElementMult(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("elementwise multiplication", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;
	
	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
		for(int i = 0; i < size; i++) {
			omxSetVectorElement(result, i,
				omxVectorElement(first, i) *
				omxVectorElement(second, i));
		}
		result->colMajor = first->colMajor;
		omxMatrixLeadingLagging(result);
	} else {
		for(int i = 0; i < rows; i++) {
			for(int j = 0; j < cols; j++) {
				omxSetMatrixElement(result, i, j,
					omxMatrixElement(first, i, j) *
					omxMatrixElement(second, i, j));
			}
		}
	}
}


static void omxKroneckerProd(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];

	int preMulRows = preMul->rows;
	int preMulCols = preMul->cols;
	int postMulRows = postMul->rows;
	int postMulCols = postMul->cols;
	int rows = preMulRows * postMulRows;
	int cols = preMulCols * postMulCols;

	if(result->rows != rows || result->cols != cols)
		omxResizeMatrix(result, rows, cols);

	for(int preRow = 0; preRow < preMulRows; preRow++)
		for(int postRow = 0; postRow < postMulRows; postRow++)
			for(int preCol = 0; preCol < preMulCols; preCol++)
				for(int postCol = 0; postCol < postMulCols; postCol++)
					omxSetMatrixElement(result, preRow * postMulRows + postRow,
						preCol * postMulCols + postCol,
						omxMatrixElement(preMul, preRow, preCol) * omxMatrixElement(postMul, postRow, postCol));
}

static void omxKroneckerPower(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];

	int rows = preMul->rows * postMul->rows;
	int cols = preMul->cols * postMul->cols;

	if(result->rows != rows || result->cols != cols)
		omxResizeMatrix(result, rows, cols);

	for(int preRow = 0; preRow < preMul->rows; preRow++)
		for(int postRow = 0; postRow < postMul->rows; postRow++)
			for(int preCol = 0; preCol < preMul->cols; preCol++)
				for(int postCol = 0; postCol < postMul->cols; postCol++)
					omxSetMatrixElement(result, preRow*postMul->rows + postRow,
						preCol*postMul->cols + postCol,
						pow(omxMatrixElement(preMul, preRow, preCol), omxMatrixElement(postMul, postRow, postCol)));
}

static void omxQuadraticProd(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];
	/* A %&% B = ABA' */

	static double zero = 0.0;
	static double one = 1.0;

	/* Conformability Check! */
	if(preMul->cols != postMul->rows || postMul->rows != postMul->cols) {
		omxRaiseError("Non-conformable matrices in Matrix Quadratic Product.");
		return;
	}

	omxMatrix* intermediate = NULL;
	intermediate = omxInitMatrix(preMul->rows, postMul->cols, TRUE, preMul->currentState);

	if(result->rows != preMul->rows || result->cols != preMul->rows)
		omxResizeMatrix(result, preMul->rows, preMul->rows);

	/* The call itself */
	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic: premul.");}
	F77_CALL(omxunsafedgemm)((preMul->majority), (postMul->majority), &(preMul->rows), &(postMul->cols), &(preMul->cols), &one, preMul->data, &(preMul->leading), postMul->data, &(postMul->leading), &zero, intermediate->data, &(intermediate->leading));

	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic: postmul.");}
//	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic postmul: result is (%d x %d), %d leading, inter is (%d x %d), prem is (%d x %d), post is (%d x %d).", result->rows, result->cols, result->leading, intermediate->rows, intermediate->cols, preMul->rows, preMul->cols, postMul->rows, postMul->cols);}
	F77_CALL(omxunsafedgemm)((intermediate->majority), (preMul->minority), &(intermediate->rows), &(preMul->rows), &(intermediate->cols), &one, intermediate->data, &(intermediate->leading), preMul->data, &(preMul->leading), &zero, result->data, &(result->leading));
	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic: clear.");}

	omxFreeMatrix(intermediate);

}

static void omxElementDivide(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("elementwise divide", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
		for(int i = 0; i < size; i++) {
			omxSetVectorElement(result, i,
				omxVectorElement(first, i) /
				omxVectorElement(second, i));
		}
		result->colMajor = first->colMajor;
		omxMatrixLeadingLagging(result);
	} else {
		for(int i = 0; i < rows; i++) {
			for(int j = 0; j < cols; j++) {
				omxSetMatrixElement(result, i, j,
					omxMatrixElement(first, i, j) /
					omxMatrixElement(second, i, j));
			}
		}
	}
}

static void omxUnaryNegation(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int rows = inMat->rows;
	int cols = inMat->cols;

	if((rows != result->rows) || (cols != result->cols)){
		omxResizeMatrix(result, rows, cols);
	}

	int vec_Rf_length = rows * cols;
	for (int i=0; i < vec_Rf_length; i++){
		double ith_value = omxVectorElement(inMat, i);
		if (ith_value == 0.0){
			omxSetVectorElement(result, i, 1.0);
		}
		else {
			omxSetVectorElement(result, i, 0.0);
		}
	}
	result->colMajor = inMat->colMajor;
	omxMatrixLeadingLagging(result);
}

static void omxBinaryOr(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("binary or", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
		int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

		if (first->colMajor == second->colMajor) {
	        	for(int i = 0; i < size; i++) {
					double ith_first  = omxVectorElement(first, i);
					double ith_second =omxVectorElement(second, i);
					if ((ith_first == 0.0) && (ith_second == 0.0)){
						omxSetVectorElement(result, i, 0.0);
					}
					else {
						omxSetVectorElement(result, i, 1.0);
					}
	        	}
	        result->colMajor = first->colMajor;
	        omxMatrixLeadingLagging(result);
		} else {
	        	for(int i = 0; i < rows; i++) {
	                	for(int j = 0; j < cols; j++) {
							double ith_first  = omxMatrixElement(first, i, j);
							double ith_second = omxMatrixElement(second, i, j);
	                        	if ((ith_first == 0.0) && (ith_second == 0.0)){
	                                	omxSetMatrixElement(result, i, j, 0.0);
	                        	}
	                        	else {
	                                	omxSetMatrixElement(result, i, j, 1.0);
	                        	}
	                	}
	        	}
		}
}

static void omxBinaryAnd(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("binary and", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
		int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

		if (first->colMajor == second->colMajor) {
	        	for(int i = 0; i < size; i++) {
					double ith_first  = omxVectorElement(first, i);
					double ith_second =omxVectorElement(second, i);
					if ((ith_first == 0.0) || (ith_second == 0.0)){
						omxSetVectorElement(result, i, 0.0);
					}
					else {
						omxSetVectorElement(result, i, 1.0);
					}
	        	}
	        result->colMajor = first->colMajor;
	        omxMatrixLeadingLagging(result);
		} else {
	        	for(int i = 0; i < rows; i++) {
	                	for(int j = 0; j < cols; j++) {
							double ith_first  = omxMatrixElement(first, i, j);
							double ith_second = omxMatrixElement(second, i, j);
	                        	if ((ith_first == 0.0) || (ith_second == 0.0)){
	                                	omxSetMatrixElement(result, i, j, 0.0);
	                        	}
	                        	else {
	                                	omxSetMatrixElement(result, i, j, 1.0);
	                        	}
	                	}
	        	}
		}
}

static void omxBinaryLessThan(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("binary less than", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
		int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

		if (first->colMajor == second->colMajor) {
	        	for(int i = 0; i < size; i++) {
	                	double ith_value = omxVectorElement(first, i) -
	                        		   omxVectorElement(second, i);
						if (ith_value < 0.0){
							omxSetVectorElement(result, i, 1.0);
						}
						else {
							omxSetVectorElement(result, i, 0.0);
						}
	        	}
	        result->colMajor = first->colMajor;
	        omxMatrixLeadingLagging(result);
		} else {
	        	for(int i = 0; i < rows; i++) {
	                	for(int j = 0; j < cols; j++) {
	 				double ith_value = omxMatrixElement(first, i, j) -
	                                   omxMatrixElement(second, i, j);

	                        	if (ith_value < 0.0){
	                                	omxSetMatrixElement(result, i, j, 1.0);
	                        	}
	                        	else {
	                                	omxSetMatrixElement(result, i, j, 0.0);
	                        	}
	                	}
	        	}
		}
}

static void omxBinaryGreaterThan(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("binary greater than", matList, result);
		return;
	}
	
        omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
        	for(int i = 0; i < size; i++) {
                	double ith_value = omxVectorElement(first, i) -
                        		   omxVectorElement(second, i);
			if (ith_value > 0.0){
				omxSetVectorElement(result, i, 1.0);
			}
			else {
				omxSetVectorElement(result, i, 0.0);
			}
        	}
        result->colMajor = first->colMajor;
        omxMatrixLeadingLagging(result);
	} else {
        	for(int i = 0; i < rows; i++) {
                	for(int j = 0; j < cols; j++) {
 				double ith_value = omxMatrixElement(first, i, j) -
                                           	   omxMatrixElement(second, i, j);

                        	if (ith_value > 0.0){
                                	omxSetMatrixElement(result, i, j, 1.0);
                        	}
                        	else {
                                	omxSetMatrixElement(result, i, j, 0.0);
                        	}
                	}
        	}
	}
}

static void omxBinaryApproxEquals(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
        omxMatrix* first  = matList[0];
	omxMatrix* second = matList[1];
	omxMatrix* epsilon = matList[2]; 
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;
	double negativeOne = -1.0;

	omxResizeMatrix(result, first->rows, first->cols);

	if (first->colMajor == second->colMajor && second->colMajor == epsilon->colMajor) {
        	for(int i = 0; i < size; i++) {
                double ith_value = omxVectorElement(first, i) -
                        		   omxVectorElement(second, i);
				double epsilon_value = omxVectorElement(epsilon, i);
				
				if (ith_value < 0.0){
					ith_value = ith_value * negativeOne;
				}
				if (ith_value < epsilon_value){
					omxSetVectorElement(result, i, 1.0);
				}
				else {
					omxSetVectorElement(result, i, 0.0);
				}
        	}
        result->colMajor = first->colMajor;
        omxMatrixLeadingLagging(result);
	} else {
        	for(int i = 0; i < rows; i++) {
                	for(int j = 0; j < cols; j++) {
 						    double ith_value = omxMatrixElement(first, i, j) -
                                           	   omxMatrixElement(second, i, j);

							double epsilon_value = omxMatrixElement(epsilon, i, j);
							if (ith_value < 0.0){
								ith_value = ith_value * negativeOne;
							}
                        	if (ith_value < epsilon_value){
                                	omxSetMatrixElement(result, i, j, 1.0);
                        	}
                        	else {
                                	omxSetMatrixElement(result, i, j, 0.0);
                        	}
                	}
        	}
	}

}

static void omxMatrixAdd(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("matrix add", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
		for(int i = 0; i < size; i++) {
			omxSetVectorElement(result, i,
				omxVectorElement(first, i) +
				omxVectorElement(second, i));
		}
		result->colMajor = first->colMajor;
		omxMatrixLeadingLagging(result);
	} else {
		for(int i = 0; i < rows; i++) {
			for(int j = 0; j < cols; j++) {
				omxSetMatrixElement(result, i, j,
					omxMatrixElement(first, i, j) +
					omxMatrixElement(second, i, j));
			}
		}
	}
}

static int matrixExtractIndices(omxMatrix *source, int dimLength, int **indices, omxMatrix *result) {

	int *retval;
	/* Case 1: the source vector contains no elements */
	if (source->rows == 0 || source->cols == 0) {
		retval = (int*) calloc(dimLength, sizeof(int));
		for(int i = 0; i < dimLength; i++) {
			retval[i] = i;
		}
		*indices = retval;
		return(dimLength);
	}
	int zero = 0, positive = 0, negative = 0;
	/* Count the number of zero, positive, and negative elements */
	for(int i = 0; i < source->rows * source->cols; i++) {
		double delement = omxVectorElement(source, i);
		if (!R_finite(delement)) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "non-finite value in '[' operator.\n");
			omxRaiseError(errstr);
			free(errstr);
			return(0);
		}
		int element = (int) delement;
		if (element < 0) {
			/* bounds checking */
			if (element < - dimLength) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "index %d is out of bounds in '[' operator.", element);
				omxRaiseError(errstr);
				free(errstr);
				return(0);
			}
			negative++;
		} else if (element == 0) {
			zero++;
		} else {
			/* bounds checking */
			if (element > dimLength) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "index %d is out of bounds in '[' operator.", element);
				omxRaiseError(errstr);
				free(errstr);
				return(0);
			}
			positive++;
		}
	}
	/* It is illegal to mix positive and negative elements */
	if (positive > 0 && negative > 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Positive and negative indices together in '[' operator.");
		omxRaiseError(errstr);
		free(errstr);
		return(0);
	}
	/* convert negative indices into a list of positive indices */
	if (negative > 0) {
		int *track = (int*) calloc(dimLength, sizeof(int));
		int Rf_length = dimLength;
		for(int i = 0; i < source->rows * source->cols; i++) {
			int element = (int) omxVectorElement(source, i);
			if (element < 0) {
				if (!track[-element - 1]) Rf_length--;
				track[-element - 1]++;
			}
		}
		if (Rf_length == 0) {
			free(track);
			return(0);
		}
		retval = (int*) calloc(Rf_length, sizeof(int));
		int j = 0;
		for(int i = 0; i < dimLength; i++) {
			if(!track[i]) {
				retval[j++] = i;
			}
		}
		free(track);
		*indices = retval;
		return(Rf_length);
	}
	/* convert positive indices with offset of zero instead of one */
	if (positive > 0) {
		int Rf_length = positive - zero;
		retval = (int*) calloc(Rf_length, sizeof(int));
		int j = 0;
		for(int i = 0; i < source->rows * source->cols; i++) {
			int element = (int) omxVectorElement(source, i);
			if (element > 0) {
				retval[j++] = element - 1;
			}
		}
		*indices = retval;
		return(Rf_length);
	}
	/* return zero Rf_length if no positive or negative elements */
	return(0);
}

static void omxMatrixExtract(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* rowMatrix = matList[1];
	omxMatrix* colMatrix = matList[2];

	if(OMX_DEBUG_ALGEBRA) { omxPrint(rowMatrix, "Row matrix: "); }
	if(OMX_DEBUG_ALGEBRA) { omxPrint(colMatrix, "Col matrix: "); }

	int *rowIndices, *colIndices;
	int rowIndexLength, colIndexLength;

	rowIndexLength = matrixExtractIndices(rowMatrix, inMat->rows, &rowIndices, result);
	colIndexLength = matrixExtractIndices(colMatrix, inMat->cols, &colIndices, result);

	if (result->rows != rowIndexLength || result->cols != colIndexLength) {
		omxResizeMatrix(result, rowIndexLength, colIndexLength);
	}

	for(int row = 0; row < rowIndexLength; row++) {
		for(int col = 0; col < colIndexLength; col++) {
			if(OMX_DEBUG_ALGEBRA) { mxLog("ALGEBRA: Matrix Extract: (%d, %d)[%d, %d] <- (%d, %d)[%d,%d].", result->rows, result->cols, row, col, rowIndexLength, colIndexLength, rowIndices[row], colIndices[col]);}
			double element = omxMatrixElement(inMat, rowIndices[row], colIndices[col]);
			omxSetMatrixElement(result, row, col, element);
		}
	}

	if (rowIndexLength > 0) free(rowIndices);
	if (colIndexLength > 0) free(colIndices);

}

static void omxMatrixSubtract(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (want == FF_COMPUTE_DIMS) {
		ensureElemConform("matrix subtract", matList, result);
		return;
	}

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

	omxResizeMatrix(result, rows, cols);

	if (first->colMajor == second->colMajor) {
		for(int i = 0; i < size; i++) {
			omxSetVectorElement(result, i,
				omxVectorElement(first, i) -
				omxVectorElement(second, i));
		}
		result->colMajor = first->colMajor;
		omxMatrixLeadingLagging(result);
	} else {
		for(int i = 0; i < rows; i++) {
			for(int j = 0; j < cols; j++) {
				omxSetMatrixElement(result, i, j,
					omxMatrixElement(first, i, j) -
					omxMatrixElement(second, i, j));
			}
		}
	}
}

static void omxUnaryMinus(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int rows = inMat->rows;
	int cols = inMat->cols;
	int size = rows * cols;

	if((rows != result->rows) || (cols != result->cols)) {
		omxResizeMatrix(result, rows, cols);
	}

	for(int i = 0; i < size; i++) {
		omxSetVectorElement(result, i,
			- omxVectorElement(inMat, i));
	}
	result->colMajor = inMat->colMajor;
	omxMatrixLeadingLagging(result);

}

static void omxMatrixHorizCatOp(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrixHorizCat(matList, numArgs, result);
}

static void omxMatrixVertCatOp(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrixVertCat(matList, numArgs, result);
}

static void omxMatrixDeterminant(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* calcMat;					// This should be preallocated.

	int rows = inMat->rows;
	int cols = inMat->cols;
	double det = 1;
	int info;

	if(rows != cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Determinant of non-square matrix cannot be found.\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	calcMat = omxInitMatrix(rows, cols, TRUE, inMat->currentState);
	omxCopyMatrix(calcMat, inMat);

	int* ipiv = (int*) calloc(inMat->rows, sizeof(int));

	F77_CALL(dgetrf)(&(calcMat->rows), &(calcMat->cols), calcMat->data, &(calcMat->cols), ipiv, &info);

	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Determinant Calculation: Nonsingular matrix (at row %d) on LUP decomposition.", info);
		omxRaiseError(errstr);
		free(errstr);
		free(ipiv);
		omxFreeMatrix(calcMat);
		return;
	}

	if(OMX_DEBUG_ALGEBRA) {
		omxPrint(calcMat, "LU Decomp");
		mxLog("info is %d.", info);
	}

	for(int i = 0; i < rows; i++) {
		det *= omxMatrixElement(calcMat, i, i);
		if(ipiv[i] != (i+1)) det *= -1;
	}

	if(OMX_DEBUG_ALGEBRA) {
		mxLog("det is %f.", det);
	}

	omxFreeMatrix(calcMat);

	omxSetMatrixElement(result, 0, 0, det);

	free(ipiv);
}

static void omxMatrixTraceOp(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrixTrace(matList, numArgs, result);
}

static void omxMatrixTotalSum(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	double sum = 0.0;

	/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matRf_length = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matRf_length; k++) {
			sum += data[k];
		}
	}

	omxSetMatrixElement(result, 0, 0, sum);
}

static void omxMatrixTotalProduct(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	double product = 1.0;

	/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matRf_length = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matRf_length; k++) {
			product *= data[k];
		}
	}

	omxSetMatrixElement(result, 0, 0, product);
}

static void omxMatrixArithmeticMean(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	omxMatrix *input = matList[0];
	int matLength = input->rows * input->cols;
	if (matLength == 0) return;
	double mean = omxVectorElement(input, 0);
	for(int i = 1; i < matLength; i++) {
		double val = omxVectorElement(input, i);
		mean += (val - mean) / (i + 1);	
	}

	omxSetMatrixElement(result, 0, 0, mean);
}

static void omxMatrixMinimum(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	double min = DBL_MAX; // DBL_MAX is the maximum possible DOUBLE value, usually 10e37.
						  // We could change this to use NPSOL's INFINITY, but why bother?

	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matRf_length = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matRf_length; k++) {
			if(data[k] < min) min = data[k];
		}
	}

	omxSetMatrixElement(result, 0, 0, min);
}

static void omxMatrixMaximum(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1);
	}

	double max = -DBL_MAX;

	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matRf_length = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matRf_length; k++) {
			if(data[k] > max) max = data[k];
		}
	}

	omxSetMatrixElement(result, 0, 0, max);
}

static void omxMatrixAbsolute(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = fabs(data[j]);
	}

}

static void omxMatrixDiagonal(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	int diags = inMat->cols;
	if(inMat->cols > inMat->rows) {
		diags = inMat->rows;
	}

	if (result->cols != 1 || result->rows != diags) {
		omxResizeMatrix(result, diags, 1);
	}

	for(int j = 0; j < diags; j++) {
		omxSetMatrixElement(result, j, 0, omxMatrixElement(inMat, j, j));
	}

}

static void omxMatrixFromDiagonal(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	int diags = inMat->cols;

	if(inMat->cols < inMat->rows) {
		diags = inMat->rows;
	}

	if(inMat->cols != 1 && inMat->rows != 1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "To generate a matrix from a diagonal that is not 1xN or Nx1.");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if (result->cols != diags || result->rows != diags) {
			omxResizeMatrix(result, diags, diags);
	}

	for(int j = 0; j < diags; j++) {
		for(int k = 0; k < diags; k++) {
			if(j == k) {
				omxSetMatrixElement(result, j, k, omxVectorElement(inMat, j));
			} else {
				omxSetMatrixElement(result, j, k, 0);
			}
		}
	}
}

static void omxElementCosine(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cos(data[j]);
	}

}

static void omxElementCosh(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cosh(data[j]);
	}

}

static void omxElementSine(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sin(data[j]);
	}

}

static void omxElementSinh(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sinh(data[j]);
	}

}

static void omxElementTangent(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = tan(data[j]);
	}

}

static void omxElementTanh(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = tanh(data[j]);
	}

}

static void omxElementExponent(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = exp(data[j]);
	}

}

static void omxElementNaturalLog(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = log(data[j]);
	}

}

static void omxElementSquareRoot(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sqrt(data[j]);
	}
}

static void omxMatrixVech(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {
	omxMatrix *inMat = matList[0];

	int size;
	if (inMat->rows > inMat->cols) {
		size = inMat->cols * (2 * inMat->rows - inMat->cols + 1) / 2;
	} else {
		size = inMat->rows * (inMat->rows + 1) / 2;
	}

	/* Consistency check: */
	if(result->rows != size || result->cols != 1) {
		omxResizeMatrix(result, size, 1);
	}

	int counter = 0;
	for(int i = 0; i < inMat->cols; i++) {
		for(int j = i; j < inMat->rows; j++) {
			omxSetMatrixElement(result, counter, 0, omxMatrixElement(inMat, j, i));
			counter++;
		}
	}

	if(counter != size) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Internal Rf_error in vech().\n");
		omxRaiseError(errstr);
		free(errstr);
	}

}

static void omxMatrixVechs(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {
	omxMatrix *inMat = matList[0];

	int size;
	if (inMat->rows > inMat->cols) {
		size = inMat->cols * (2 * inMat->rows - inMat->cols + 1) / 2 - inMat->cols;
	} else {
		size = inMat->rows * (inMat->rows + 1) / 2 - inMat->rows;
	}

	/* Consistency check: */
	if(result->rows != size || result->cols != 1) {
		omxResizeMatrix(result, size, 1);
	}

	int counter = 0;
	for(int i = 0; i < inMat->cols; i++) {
		for(int j = i + 1; j < inMat->rows; j++) {
			omxSetMatrixElement(result, counter, 0, omxMatrixElement(inMat, j, i));
			counter++;
		}
	}

	if(counter != size) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Internal Rf_error in vechs().\n");
		omxRaiseError(errstr);
		free(errstr);
	}

}

static void omxRowVectorize(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int size = (inMat->rows * inMat->cols);

	/* Consistency Check */
	if(result->rows != size || result->cols != 1)
		omxResizeMatrix(result, size, 1);

	if(!inMat->colMajor) {		// Special case: we can just memcpy.
		memcpy(result->data, inMat->data, size*sizeof(double));
	} else {
		int next = 0;
		for(int i = 0; i < inMat->rows; i++) {
			for(int j = 0; j < inMat->cols; j++) {
				omxSetMatrixElement(result, next++, 0, omxMatrixElement(inMat, i, j));
			}
		}
	}
}

static void omxColVectorize(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int size = (inMat->rows * inMat->cols);

	/* Consistency Check */
	if(result->rows != size || result->cols != 1)
		omxResizeMatrix(result, size, 1);
	if(inMat->colMajor) {		// Special case: we can just memcpy.
		memcpy(result->data, inMat->data, size * sizeof(double));
	} else {
		int next = 0;
		for(int i = 0; i < inMat->cols; i++) {
			for(int j = 0; j < inMat->rows; j++) {
				omxSetMatrixElement(result, next++, 0, omxMatrixElement(inMat, j, i));
			}
		}
	}
}


static void omxSequenceGenerator(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {

	double start = omxVectorElement(matList[0], 0);
	double stop = omxVectorElement(matList[1], 0);

	if (!R_finite(start)) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-finite start value in ':' operator.\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if (!R_finite(stop)) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-finite stop value in ':' operator.\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	double difference = stop - start;
	if (difference < 0) difference = - difference;

	int size = ((int) difference) + 1;

	/* Consistency check: */
	if(result->rows != size || result->cols != 1) {
		omxResizeMatrix(result, size, 1);
	}

	/* Sanity-checking.  This loop can be eliminated */
	for(int i = 0; i < size; i++) {
		omxSetVectorElement(result, i, 0);
	}

	int count = 0;
	if ((stop - start) >= 0) {
		while (start <= stop) {
			omxSetVectorElement(result, count, start);
			start = start + 1.0;
			count++;
		}
	} else {
		while (start >= stop) {
			omxSetVectorElement(result, count, start);
			start = start - 1.0;
			count++;
		}
	}
}

static void omxMultivariateNormalIntegration(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {

	omxMatrix* cov = matList[0];
	omxMatrix* means = matList[1];
	omxMatrix* lBoundMat = matList[2];
	omxMatrix* uBoundMat = matList[3];

	/* Conformance checks: */
	if (result->rows != 1 || result->cols != 1) omxResizeMatrix(result, 1, 1);

	if (cov->rows != cov->cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "covariance is not a square matrix");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if (means->rows > 1 && means->cols > 1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "means is neither row nor column vector");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if (lBoundMat->rows > 1 && lBoundMat->cols > 1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "lbound is neither row nor column vector");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if (uBoundMat->rows > 1 && uBoundMat->cols > 1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "ubound is neither row nor column vector");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	int nElements = (cov->cols > 1) ? cov->cols : cov->rows;
	double *lBounds, *uBounds;
	double *weights;
	double *corList;
	lBounds = (double*) malloc(nElements * sizeof(double));
	uBounds = (double*) malloc(nElements * sizeof(double));
	weights = (double*) malloc(nElements * sizeof(double));
	corList = (double*) malloc((nElements * (nElements + 1) / 2) * sizeof(double));

	omxStandardizeCovMatrix(cov, corList, weights);

	// SADMVN calls Alan Genz's sadmvn.f--see appropriate file for licensing info.
	// TODO: Check with Genz: should we be using sadmvn or sadmvn?
	// Parameters are:
	// 	N 		int			# of vars
	//	Lower	double*		Array of lower bounds
	//	Upper	double*		Array of upper bounds
	//	Infin	int*		Array of flags: <0 = (-Inf, Inf) 0 = (-Inf, upper] 1 = [lower, Inf), 2 = [lower, upper]
	//	Correl	double*		Array of correlation coeffs: in row-major lower triangular order
	//	MaxPts	int			Maximum # of function values (use 1000*N or 1000*N*N)
	//	Abseps	double		Absolute Rf_error tolerance.  Yick.
	//	Releps	double		Relative Rf_error tolerance.  Use EPSILON.
	//	Error	&double		On return: absolute real Rf_error, 99% confidence
	//	Value	&double		On return: evaluated value
	//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
	// TODO: Separate block diagonal covariance matrices into pieces for integration separately
	double Error;
	double likelihood;
	int inform;
	int numVars = cov->rows;
	Eigen::VectorXi Infin(cov->rows);
	int fortranThreadId = omx_absolute_thread_num() + 1;

	for(int i = 0; i < nElements; i++) {
		lBounds[i] = (omxVectorElement(lBoundMat, i) - omxVectorElement(means, i))/weights[i];
		uBounds[i] = (omxVectorElement(uBoundMat, i) - omxVectorElement(means, i))/weights[i];
		Infin[i] = 2; // Default to both thresholds
		if(uBounds[i] <= lBounds[i]) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Thresholds are not strictly increasing: %3.3f >= %3.3f.", lBounds[i], uBounds[i]);
			omxRaiseError(errstr);
			free(errstr);
			free(corList);
			free(weights);
			free(uBounds);
			free(lBounds);
			return;
		}
		if(!R_finite(lBounds[i]) ) {
			Infin[i] -= 2;	// NA or INF or -INF means no lower threshold.
		} else {

		}
		if(!R_finite(uBounds[i]) ) {
			Infin[i] -= 1; // NA or INF or -INF means no upper threshold.
		}

	}


	double absEps = Global->absEps;
	double relEps = Global->relEps;
	int MaxPts = Global->maxptsa + Global->maxptsb * cov->rows + Global->maxptsc * cov->rows * cov->rows;
	F77_CALL(sadmvn)(&numVars, &(lBounds[0]), &(*uBounds), Infin.data(), corList, 
		&MaxPts, &absEps, &relEps, &Error, &likelihood, &inform, &fortranThreadId);

	if(OMX_DEBUG_ALGEBRA) { mxLog("Output of sadmvn is %f, %f, %d.", Error, likelihood, inform); }

	if(inform == 2) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Improper input to sadmvn.");
		omxRaiseError(errstr);
		free(errstr);
		free(corList);
		free(weights);
		free(uBounds);
		free(lBounds);
		return;
	}

	free(corList);
	free(weights);
	free(uBounds);
	free(lBounds);

	omxSetMatrixElement(result, 0, 0, likelihood);

}

static void omxAllIntegrationNorms(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* cov = matList[0];
	omxMatrix* means = matList[1];
	int nCols = cov->cols;
	int i,j,k;

	int totalLevels = 1;
	omxMatrix **thresholdMats = (omxMatrix **) malloc(nCols * sizeof(omxMatrix*));
	int *numThresholds = (int*) malloc(nCols * sizeof(int));
	int *matNums = (int*) malloc(nCols * sizeof(int));
	int *thresholdCols = (int*) malloc(nCols * sizeof(int));
	int *currentThresholds = (int*) malloc(nCols * sizeof(int));

	int currentMat = 0;

	for(i = currentMat; i < nCols;) {							// Map out the structure of levels.
	if(OMX_DEBUG_ALGEBRA) {
		mxLog("All-part multivariate normal integration: Examining threshold column %d.", i);
	}
		thresholdMats[currentMat] = matList[currentMat+2];		// Get the thresholds for this covariance column

		for(j = 0; j < thresholdMats[currentMat]->cols; j++) {	// We walk along the columns of this threshold matrix
			double ubound, lbound = omxMatrixElement(thresholdMats[currentMat], 0, j);
			if(ISNA(lbound)) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "Invalid lowest threshold for dimension %d of Allint.", j);
				omxRaiseError(errstr);
				free(errstr);
				return;
			}

			thresholdCols[i] = j;

			for(k = 1; k < thresholdMats[currentMat]->rows; k++) {
				ubound = omxMatrixElement(thresholdMats[currentMat], k, j);
				if(ISNA(ubound)) {
					numThresholds[i] = k-1;
					totalLevels *= numThresholds[i];
					break;
				}

				if(!(ubound > lbound)) {
					char *errstr = (char*) calloc(250, sizeof(char));
					sprintf(errstr, "Thresholds (%f and %f) are not strictly increasing for dimension %d of Allint.", lbound, ubound, j+1);
					omxRaiseError(errstr);
					free(errstr);
					return;
				}

				if(!R_finite(ubound)) {					// Infinite bounds must be last.
					numThresholds[i] = k;
					totalLevels *= numThresholds[i];
					break;
				}

				if(k == (thresholdMats[currentMat]->rows -1)) { // In case the highest threshold isn't Infinity
					numThresholds[i] = k;
					totalLevels *= numThresholds[i];
				}
			}
			currentThresholds[i] = 1;
			matNums[i] = currentMat;
			if(++i >= nCols) {							// We have all we need
				break;
			}
		}
		currentMat++;
	}

	/* Conformance checks: */
	if(result->rows != totalLevels || result->cols != 1) omxResizeMatrix(result, totalLevels, 1);

	double *weights = (double*) malloc(nCols * sizeof(double));
	double *corList = (double*) malloc((nCols * (nCols + 1) / 2) * sizeof(double));

	omxStandardizeCovMatrix(cov, &(*corList), &(*weights));

	// SADMVN calls Alan Genz's sadmvn.f--see appropriate file for licensing info.
	// TODO: Check with Genz: should we be using sadmvn or sadmvn?
	// Parameters are:
	// 	N 		int			# of vars
	//	Lower	double*		Array of lower bounds
	//	Upper	double*		Array of upper bounds
	//	Infin	int*		Array of flags: <0 = (-Inf, Inf) 0 = (-Inf, upper] 1 = [lower, Inf), 2 = [lower, upper]
	//	Correl	double*		Array of correlation coeffs: in row-major lower triangular order
	//	MaxPts	int			Maximum # of function values (use 1000*N or 1000*N*N)
	//	Abseps	double		Absolute Rf_error tolerance.  Yick.
	//	Releps	double		Relative Rf_error tolerance.  Use EPSILON.
	//	Error	&double		On return: absolute real Rf_error, 99% confidence
	//	Value	&double		On return: evaluated value
	//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
	// TODO: Separate block diagonal covariance matrices into pieces for integration separately
	double Error;
	double likelihood;
	int inform;
	int numVars = nCols;
	int* Infin = (int*) malloc(nCols * sizeof(int));
	double* lBounds = (double*) malloc(nCols * sizeof(double));
	double* uBounds = (double*) malloc(nCols * sizeof(double));
	int fortranThreadId = omx_absolute_thread_num() + 1;

	/* Set up first row */
	for(j = (nCols-1); j >= 0; j--) {					// For each threshold set, starting from the fastest

		Infin[j] = 2; 									// Default to using both thresholds
		lBounds[j] = (omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j]-1, thresholdCols[j]) - omxVectorElement(means, j))/weights[j];
		if(!R_finite(lBounds[j])) { 					// Inifinite lower bounds = -Inf to ?
				Infin[j] -= 2;
		}

		uBounds[j] = (omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j], thresholdCols[j]) - omxVectorElement(means, j))/weights[j];

		if(!R_finite(uBounds[j])) { 					// Inifinite lower bounds = -Inf to ?
				Infin[j] -= 1;
		}

		if(Infin[j] < 0) { Infin[j] = 3; }			// Both bounds infinite.
	}

	double absEps = Global->absEps;
	double relEps = Global->relEps;
	int MaxPts = Global->maxptsa + Global->maxptsb * cov->rows + Global->maxptsc * cov->rows * cov->rows;
	F77_CALL(sadmvn)(&numVars, &(lBounds[0]), &(*uBounds), Infin, corList, 
		&MaxPts, &absEps, &relEps, &Error, &likelihood, &inform, &fortranThreadId);

	if(OMX_DEBUG_ALGEBRA) { mxLog("Output of sadmvn is %f, %f, %d.", Error, likelihood, inform); }

	if(inform == 2) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Improper input to sadmvn.");
		omxRaiseError(errstr);
		free(errstr);
		goto AllIntCleanup;
	}

	omxSetMatrixElement(result, 0, 0, likelihood);


	/* And repeat with increments for all other rows. */
	for(i = 1; i < totalLevels; i++) {
		for(j = (nCols-1); j >= 0; j--) {							// For each threshold set, starting from the fastest
			currentThresholds[j]++;									// Move to the next threshold set.
			if(currentThresholds[j] > numThresholds[j]) {			// Hit the end; cycle to the next.
				currentThresholds[j] = 1;
			}

			/* Update only the rows that need it. */
			Infin[j] = 2; // Default to both thresholds
			lBounds[j] = (omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j]-1, thresholdCols[j]) - omxVectorElement(means, j))/weights[j];
			if(!R_finite(lBounds[j])) {								// Inifinite lower bounds = -Inf to ?
				Infin[j] -= 2;
			}
			uBounds[j] = (omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j], thresholdCols[j]) - omxVectorElement(means, j))/weights[j];

			if(!R_finite(uBounds[j])) { 							// Inifinite lower bounds = -Inf to ?
				Infin[j] -= 1;
			}

			if(Infin[j] < 0) { Infin[j] = 3; }						// Both bounds infinite.

			if(currentThresholds[j] != 1) {							// If we just cycled, we need to see the next set.
				break;
			}

		}

		F77_CALL(sadmvn)(&numVars, &(lBounds[0]), &(*uBounds), Infin, corList,
			&MaxPts, &absEps, &relEps, &Error, &likelihood, &inform, &fortranThreadId);

		if(OMX_DEBUG_ALGEBRA) { mxLog("Output of sadmvn is %f, %f, %d.", Error, likelihood, inform); }

		if(inform == 2) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Improper input to sadmvn.");
			omxRaiseError(errstr);
			free(errstr);
			goto AllIntCleanup;
		}

		omxSetMatrixElement(result, i, 0, likelihood);
	}

AllIntCleanup:
	free(Infin);
	free(lBounds);
	free(uBounds);
	free(weights);
	free(corList);
	free(thresholdMats);
	free(numThresholds);
	free(matNums);
	free(thresholdCols);
	free(currentThresholds);
}

static int omxComparePointerContentsHelper(const void* one, const void* two, void *ign) {
	double diff = (*(*(double**) two)) - (*(*(double**) one));
	if(diff > EPSILON) {
		return 1;
	} else if(diff < -EPSILON) {
		return -1;
	} else return 0;
}

static void omxSortHelper(double* sortOrder, omxMatrix* original, omxMatrix* result) {
	/* Sorts the columns of a matrix or the rows of a column vector
					in decreasing order of the elements of sortOrder. */

	if(OMX_DEBUG) {mxLog("SortHelper:Original is (%d x %d), result is (%d x %d).", original->rows, original->cols, result->rows, result->cols);}

	if(!result->colMajor || !original->colMajor
		|| result->cols != original->cols || result->rows != original->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Incorrect input to omxRowSortHelper: %d %d %d %d", result->cols, original->cols, result->rows, original->rows);
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	std::vector<double*> sortArray(original->rows);
	int numElements = original->cols;
	int numRows = original->rows;

	if(numElements == 1)  numElements = numRows;		// Special case for column vectors

	for(int i = 0; i < numElements; i++) {
		sortArray[i] = sortOrder + i;
	}

	freebsd_mergesort(sortArray.data(), numElements, sizeof(double*), omxComparePointerContentsHelper, NULL);

	if(OMX_DEBUG) {mxLog("Original is (%d x %d), result is (%d x %d).", original->rows, original->cols, result->rows, result->cols);}


	for(int i = 0; i < numElements; i++) {
		if(original->cols == 1) {
			omxSetMatrixElement(result, i, 0, omxMatrixElement(original, (sortArray[i] - sortOrder), 0));
		} else {
			memcpy(omxLocationOfMatrixElement(result, 0, i), omxLocationOfMatrixElement(original, 0, sortArray[i]-sortOrder), numRows * sizeof(double));
		}
	}

	return;
}

static void omxRealEigenvalues(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* A = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxMatrix* B = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxCopyMatrix(B, matList[0]);
	omxResizeMatrix(A, B->rows, 1);

	/* Conformability Check! */
	if(B->cols != B->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-square matrix in eigenvalue decomposition.\n");
		omxRaiseError(errstr);
		free(errstr);
		omxFreeMatrix(A);
		omxFreeMatrix(B);
		return;
	}

	if(result->rows != B->rows || result->cols != 1)
		omxResizeMatrix(result, B->rows, 1);

	char N = 'N';						// Indicators for BLAS
	// char V = 'V';						// Indicators for BLAS

	int One = 1;
	int lwork = 10*B->rows;

	int info;

	double* work = (double*) malloc(lwork * sizeof(double));
	double* WI = (double*) malloc(B->cols * sizeof(double));

	F77_CALL(dgeev)(&N, &N, &(B->rows), B->data, &(B->leading), A->data, WI, NULL, &One, NULL, &One, work, &lwork, &info);
	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "DGEEV returned %d in (real) eigenvalue decomposition:", info);
		if(info > 0)
			sprintf(errstr, "%s argument %d had an illegal value.  Post this to the OpenMx wiki.\n", errstr, info);
		else
			sprintf(errstr, "%s Unable to decompose matrix: Not of full rank.\n", errstr);
		omxRaiseError(errstr);
		free(errstr);
		goto RealEigenValCleanup;
	}

	result->colMajor = TRUE;

	// Calculate Eigenvalue modulus.
	for(int i = 0; i < A->rows; i++) {
		double value = omxMatrixElement(A, i, 0);
		if(WI[i] != 0) {				// FIXME: Might need to be abs(WI[i] > EPSILON)
			value = sqrt(WI[i]*WI[i] + value*value);				// Sort by eigenvalue modulus
			WI[i] = value;
			WI[++i] = value;										// Conjugate pair.
		} else {
			WI[i] = fabs(value); 									// Modulus of a real is its absolute value
		}
	}

	omxSortHelper(WI, A, result);

RealEigenValCleanup:
	omxFreeMatrix(A);				// FIXME: State-keeping for algebras would save significant time in memory allocation/deallocation
	omxFreeMatrix(B);
	omxMatrixLeadingLagging(result);

	free(work);
	free(WI);
}

static void omxRealEigenvectors(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* A = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxCopyMatrix(result, matList[0]);
	omxResizeMatrix(A, result->rows, result->cols);


	if(A == NULL) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Null matrix pointer detected.\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	/* Conformability Check! */
	if(A->cols != A->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-square matrix in (real) eigenvalue decomposition.\n");
		omxRaiseError(errstr);
		free(errstr);
		omxFreeMatrix(A);
		return;
	}

	char N = 'N';						// Indicators for BLAS
	char V = 'V';						// Indicators for BLAS

	int One = 1;
	int lwork = 10*A->rows;

	int info;

	double *WR = (double*) malloc(A->cols * sizeof(double));
	double *WI = (double*) malloc(A->cols * sizeof(double));
	double *work = (double*) malloc(lwork * sizeof(double));

	F77_CALL(dgeev)(&N, &V, &(result->rows), result->data, &(result->leading), WR, WI, NULL, &One, A->data, &(A->leading), work, &lwork, &info);
	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "DGEEV returned %d in eigenvalue decomposition:", info);
		if(info > 0)
			sprintf(errstr, "%s argument %d had an illegal value.  Post this to the OpenMx wiki.\n", errstr, info);
		else
			sprintf(errstr, "%s Unable to decompose matrix: Not of full rank.\n", errstr);
		omxRaiseError(errstr);
		free(errstr);
		goto RealEigenVecCleanup;
	}

	// Filter real and imaginary eigenvectors.  Real ones have no WI.
	for(int i = 0; i < A->cols; i++) {
		if(fabs(WI[i]) > EPSILON) {									// If this is part of a conjugate pair
			memcpy(omxLocationOfMatrixElement(A, 0, i+1), omxLocationOfMatrixElement(A, 0, i), A->rows * sizeof(double));
				// ^^ This is column-major, so we can clobber columns over one another.
			WR[i] = sqrt(WR[i] *WR[i] + WI[i]*WI[i]);				// Sort by eigenvalue modulus
			WR[i+1] = WR[i];										// Identical--conjugate pair
			i++; 	// Skip the next one; we know it's the conjugate pair.
		} else {
			WR[i] = fabs(WR[i]); 									// Modulus of a real is its absolute value
		}
	}

	result->colMajor = TRUE;

	// Sort results
	omxSortHelper(WR, A, result);

RealEigenVecCleanup:
	omxFreeMatrix(A);		// FIXME: State-keeping for algebras would save significant time in memory allocation/deallocation
	omxMatrixLeadingLagging(result);

	free(WR);
	free(WI);
	free(work);	
}

static void omxImaginaryEigenvalues(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* A = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxMatrix* B = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxCopyMatrix(B, matList[0]);
	omxResizeMatrix(A, B->rows, 1);

	/* Conformability Check! */
	if(B->cols != B->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-square matrix in eigenvalue decomposition.\n");
		omxRaiseError(errstr);
		free(errstr);
		omxFreeMatrix(A);
		omxFreeMatrix(B);
		return;
	}

	if(result->cols != 1 || result->rows != A->rows)
		omxResizeMatrix(result, B->rows, 1);

	char N = 'N';						// Indicators for BLAS

	int One = 1;
	int lwork = 10*B->rows;

	int info;

	double *WR = (double*) malloc(B->cols * sizeof(double));
	double *VR = (double*) malloc(B->rows * B->cols * sizeof(double));
	double *work = (double*) malloc(lwork * sizeof(double));

	F77_CALL(dgeev)(&N, &N, &(B->rows), B->data, &(B->leading), WR, A->data, NULL, &One, NULL, &One, work, &lwork, &info);
	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "DGEEV returned %d in (real) eigenvalue decomposition:", info);
		if(info > 0)
			sprintf(errstr, "%s argument %d had an illegal value.  Post this to the OpenMx wiki.\n", errstr, info);
		else
			sprintf(errstr, "%s Unable to decompose matrix: Not of full rank.\n", errstr);
		omxRaiseError(errstr);
		free(errstr);
		goto ImagEigenValCleanup;
	}

	// Calculate Eigenvalue modulus.
	for(int i = 0; i < result->rows; i++) {
		double value = omxMatrixElement(A, i, 0);					// A[i] is the ith imaginary eigenvalue
		value *= value;												// Squared imaginary part
		if(value > EPSILON) {
			value = sqrt(WR[i] *WR[i] + value);				// Sort by eigenvalue modulus
			WR[i] = value;
			WR[++i] = value;										// Conjugate pair.
		} else {
			WR[i] = fabs(WR[i]);
		}
	}

	result->colMajor = TRUE;

	// Sort results
	omxSortHelper(WR, A, result);

ImagEigenValCleanup:
	omxFreeMatrix(A);		// FIXME: State-keeping for algebras would save significant time in memory allocation/deallocation
	omxFreeMatrix(B);
	omxMatrixLeadingLagging(result);

	free(WR);
	free(VR);
	free(work);
}

static void omxImaginaryEigenvectors(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* A = omxInitMatrix(0, 0, TRUE, result->currentState);
	omxCopyMatrix(result, matList[0]);
	omxResizeMatrix(A, result->rows, result->cols);

	/* Conformability Check! */
	if(A->cols != A->rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-square matrix in (imaginary) eigenvalue decomposition.\n");
		omxRaiseError(errstr);
		free(errstr);
		omxFreeMatrix(A);
		return;
	}

	char N = 'N';						// Indicators for BLAS
	char V = 'V';						// Indicators for BLAS

	int One = 1;
	int lwork = 10*A->rows;

	int info;

	double *WR = (double*) malloc(A->cols * sizeof(double));
	double *WI = (double*) malloc(A->cols * sizeof(double));
	double *work = (double*) malloc(lwork * sizeof(double));

	if(result->rows != A->rows || result->cols != A->cols)
		omxResizeMatrix(result, A->rows, A->cols);

	F77_CALL(dgeev)(&N, &V, &(result->rows), result->data, &(result->leading), WR, WI, NULL, &One, A->data, &(A->leading), work, &lwork, &info);
	if(info != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "DGEEV returned %d in eigenvalue decomposition:", info);
		if(info > 0)
			sprintf(errstr, "%s argument %d had an illegal value.  Post this to the OpenMx wiki.\n", errstr, info);
		else
			sprintf(errstr, "%s Unable to decompose matrix: Not of full rank.\n", errstr);
		omxRaiseError(errstr);
		free(errstr);
		goto ImagEigenVecCleanup;
	}

	// Filter real and imaginary eigenvectors.  Imaginary ones have a WI.
	for(int i = 0; i < result->cols; i++) {
		if(WI[i] != 0) {				// FIXME: Might need to be abs(WI[i] > EPSILON)
			// memcpy(omxLocationOfMatrixElement(A, 0, i), omxLocationOfMatrixElement(A, 0, i+1), A->rows * sizeof(double));
			for(int j = 0; j < result->rows; j++) {
				double value = omxMatrixElement(A, j, i+1);			// Conjugate pair
				omxSetMatrixElement(A, j, i, value);				// Positive first,
				omxSetMatrixElement(A, j, i+1, -value);				// Negative second
			}
			WR[i] = sqrt(WR[i] *WR[i] + WI[i]*WI[i]);				// Sort by eigenvalue modulus
			WR[i+1] = WR[i];										// Identical--conjugate pair
			i++; 	// Skip the next one; we know it's the conjugate pair.
		} else {						// If it's not imaginary, it's zero.
			for(int j = 0; j < A->rows; j++) {
				omxSetMatrixElement(A, j, i, 0.0);
			}
			WR[i] = fabs(WR[i]); 									// Modulus of a real is its absolute value

		}
	}

	result->colMajor = TRUE;

	omxSortHelper(WR, A, result);

ImagEigenVecCleanup:
	omxFreeMatrix(A);			// FIXME: State-keeping for algebras would save significant time in memory allocation/deallocation
	omxMatrixLeadingLagging(result);

	free(WR);
	free(WI);
	free(work);

}

static void omxSelectRows(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int rows = inMat->rows;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(rows);
    int numRemoves = 0;
    
    if((selector->cols != 1) && selector->rows !=1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Selector must have a single row or a single column.\n");
        omxRaiseError(errstr);
		free(errstr);
		return;
    }

	if(selectLength != rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-conformable matrices for row selection.\n");
        omxRaiseError(errstr);
		free(errstr);
		return;
	}
	
	omxCopyMatrix(result, inMat);
	
    for(int index = 0; index < selectLength; index++) {
        if(omxVectorElement(selector, index) == 0) {
            numRemoves++;
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    if(numRemoves >= rows) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Attempted to select zero columns.\n");
        omxRaiseError(errstr);
		free(errstr);        
		return;
    }
    
    std::vector<int> zeros(inMat->cols);
    omxRemoveRowsAndColumns(result, numRemoves, 0, toRemove.data(), zeros.data());

}

static void omxSelectCols(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int cols = inMat->cols;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(cols);
    int numRemoves = 0;

    if((selector->cols != 1) && selector->rows !=1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Selector must have a single row or a single column.\n");
        omxRaiseError(errstr);
		free(errstr);        
		return;
    }

	if(selectLength != cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-conformable matrices for row selection.\n");
        omxRaiseError(errstr);
		free(errstr);
		return;
	}
	
	omxCopyMatrix(result, inMat);
	
    for(int index = 0; index < selectLength; index++) {
        if(omxVectorElement(selector, index) == 0) {
            numRemoves++;
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    if(numRemoves >= cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Attempted to select zero columns.\n");
        omxRaiseError(errstr);
		free(errstr);        
		return;
    }
    
    std::vector<int> zeros(inMat->rows);
    omxRemoveRowsAndColumns(result, 0, numRemoves, zeros.data(), toRemove.data());
}

static void omxSelectRowsAndCols(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int rows = inMat->rows;
	int cols = inMat->cols;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(cols);
    int numRemoves = 0;

    if((selector->cols != 1) && selector->rows !=1) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Selector must have a single row or a single column.\n");
        omxRaiseError(errstr);
		free(errstr);        
		return;
    }

	if(rows != cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Can only select rows and columns from square matrices.\n");
        omxRaiseError(errstr);
		free(errstr);
		return;
	}

	if(selectLength != cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-conformable matrices for row selection.\n");
        omxRaiseError(errstr);
		free(errstr);
		return;
	}
	
	omxCopyMatrix(result, inMat);
	
    for(int index = 0; index < selectLength; index++) {
        if(omxVectorElement(selector, index) == 0) {
            numRemoves++;
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    if(numRemoves >= cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Attempted to select zero columns.\n");
        omxRaiseError(errstr);
		free(errstr);        
		return;
    }
    
    omxRemoveRowsAndColumns(result, numRemoves, numRemoves, toRemove.data(), toRemove.data());
}

static void omxCovToCor(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{

    omxMatrix* inMat = matList[0];
    int rows = inMat->rows;

	omxMatrix* intermediate;

    if(inMat->rows != inMat->cols) {
        char *errstr = (char*) calloc(250, sizeof(char));
        sprintf(errstr, "cov2cor of non-square matrix cannot even be attempted\n");
        omxRaiseError(errstr);
        free(errstr);
		return;
	}

	if(result->rows != rows || result->cols != rows) {
        if(OMX_DEBUG_ALGEBRA) { mxLog("ALGEBRA: cov2cor resizing result.");}
        omxResizeMatrix(result, rows, rows);
	}

    intermediate = omxInitMatrix(1, rows, TRUE, inMat->currentState);

    for(int i = 0; i < rows; i++) {
        intermediate->data[i] = sqrt(1.0 / omxMatrixElement(inMat, i, i));
    }

    if (inMat->colMajor) {
        for(int col = 0; col < rows; col++) {
            for(int row = 0; row < rows; row++) {
                result->data[col * rows + row] = inMat->data[col * rows + row] * 
                    intermediate->data[row] * intermediate->data[col];
            }
        }
    } else {
        for(int col = 0; col < rows; col++) {
            for(int row = 0; row < rows; row++) {
                result->data[col * rows + row] = inMat->data[row * rows + col] * 
                    intermediate->data[row] * intermediate->data[col];
            }
        }
    }

    for(int i = 0; i < rows; i++) {
        result->data[i * rows + i] = 1.0;
    }

    omxFreeMatrix(intermediate);
}

static void omxCholesky(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

    int l = 0; char u = 'U';
	omxCopyMatrix(result, inMat);
	if(result->rows != result->cols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Cholesky decomposition of non-square matrix cannot even be attempted\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

    F77_CALL(dpotrf)(&u, &(result->rows), result->data, &(result->cols), &l);
	if(l != 0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Attempted to Cholesky decompose non-invertable matrix.\n");
		omxRaiseError(errstr);
		free(errstr);
		return;
	} else {
	    for(int i = 0; i < result->rows; i++) {
			for(int j = 0; j < i; j++) {
                omxSetMatrixElement(result, i, j, 0.0);
			}
		}
	}
}

static void omxVechToMatrix(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {

	omxMatrix *inMat = matList[0];
	
	int dim = (inMat->cols > inMat->rows) ? inMat->cols : inMat->rows;

	int size = sqrt(2.0 * dim + 0.25) - 0.5;

	int counter = 0;

    if(inMat->cols > 1 && inMat->rows > 1) {
        char *errstr = (char*) calloc(250, sizeof(char));
        sprintf(errstr, "vech2full input has %d rows and %d columns\n", inMat->rows, inMat->cols);
        omxRaiseError(errstr);
        free(errstr);
		return;
	}

	/* Consistency check: */
	if(result->rows != size || result->cols != size) {
		omxResizeMatrix(result, size, size);
	}

	for(int i = 0; i < size; i++) {
		for(int j = i; j < size; j++) {

			double next = omxVectorElement(inMat, counter);

			omxSetMatrixElement(result, i, j, next);
			omxSetMatrixElement(result, j, i, next);

			counter++;
		}
	}

}


static void omxVechsToMatrix(FitContext *fc, int want, omxMatrix** matList, int numArgs, omxMatrix* result) {
	omxMatrix *inMat = matList[0];
	
	int dim = (inMat->cols > inMat->rows) ? inMat->cols : inMat->rows;
	
	int size = sqrt(2.0 * dim + 0.25) + 0.5; //note the plus 0.5

	int counter = 0;

    if(inMat->cols > 1 && inMat->rows > 1) {
        char *errstr = (char*) calloc(250, sizeof(char));
        sprintf(errstr, "vechs2full input has %d rows and %d columns\n", inMat->rows, inMat->cols);
        omxRaiseError(errstr);
        free(errstr);
		return;
	}

	/* Consistency check: */
	if(result->rows != size || result->cols != size) {
		omxResizeMatrix(result, size, size);
	}

	for(int i = 0; i < size; i++) {

		omxSetMatrixElement(result, i, i, 0.0);

		for(int j = i + 1; j < size; j++) {

			double next = omxVectorElement(inMat, counter);

			omxSetMatrixElement(result, i, j, next);
			omxSetMatrixElement(result, j, i, next);

			counter++;
		}
	}

}


