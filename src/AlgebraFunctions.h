/*
 *  Copyright 2007-2016 The OpenMx Project
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

#include <limits>
#include "omxMatrix.h"
#include "merge.h"
#include "omxSadmvnWrapper.h"
#include "matrix.h"
#include "omxState.h"

// TODO: Implement wrappers for BLAS functions used here.

/* omxAlgebraFunction Wrappers */

static void omxMatrixTranspose(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxMatrixInvert(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxCopyMatrix(result, inMat);

	Matrix resultMat(result);
	int info = MatrixInvert1(result);
	if (info) {
		result->data[0] = nan("singular");
		// recordIterationError TODO
	}
}

static int BroadcastIndex = 0;

static void nameBroadcastAlg(omxMatrix *bc)
{
	bc->nameStr = string_snprintf("broadcast%03d", ++BroadcastIndex);
}

static void ensureElemConform(FitContext *fc, omxMatrix **matList, omxMatrix *result)
{
	omxMatrix *mat0 = matList[0];
	omxMatrix *mat1 = matList[1];

	if (mat0->cols == mat1->cols && mat0->rows == mat1->rows) {
		omxResizeMatrix(result, mat0->rows, mat0->cols);
		return;
	}

	if (mat0->cols == 1 && mat0->rows == 1 && mat1->rows != 0 && mat1->cols != 0) {
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
		omxAlgebraRecompute(om, FF_COMPUTE_INITIAL_FIT, fc);
		return;
	}
	if (mat1->cols == 1 && mat1->rows == 1 && mat0->rows != 0 && mat0->cols != 0) {
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
		omxAlgebraRecompute(om, FF_COMPUTE_INITIAL_FIT, fc);
		return;
	}
}

static void omxBroadcast(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *src = matList[0];

	if (src->rows == result->rows && src->cols == result->cols) {
		omxCopyMatrix(result, src);
		return;
	}

	if (src->rows != 1 || src->cols != 1) {
		Rf_error("Don't know how to broadcast from %dx%d source "
			 "matrix '%s' to %dx%d result matrix '%s'",
			 src->rows, src->cols, src->name(),
			 result->rows, result->cols, result->name());
	}

	int size = result->rows * result->cols;
	for (int dx=0; dx < size; ++dx) {
		result->data[dx] = src->data[0];
	}
}

static void omxMatrixMult(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxElementPower(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxMatrixElementMult(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;
	
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


static void omxKroneckerProd(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxKroneckerPower(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxQuadraticProd(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];
	/* A %&% B = ABA' */

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
	omxDGEMM(0, 0, 1, preMul, postMul, 0, intermediate);
	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic: postmul.");}
	omxDGEMM(0, 1, 1, intermediate, preMul, 0, result);
	if(OMX_DEBUG_ALGEBRA) { mxLog("Quadratic: clear.");}

	omxFreeMatrix(intermediate);

}

static void omxElementDivide(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxUnaryNegation(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxBinaryOr(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxBinaryAnd(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxBinaryLessThan(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxBinaryGreaterThan(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

        omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxBinaryApproxEquals(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxMatrixAdd(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

template <typename T> 
static void matrixExtractIndices(omxMatrix *source, int dimLength, Eigen::ArrayBase<T> &out)
{
	/* Case 1: the source vector contains no elements */
	if (source->rows == 0 || source->cols == 0) {
		out.derived().resize(dimLength);
		for(int i = 0; i < dimLength; i++) {
			out[i] = i;
		}
		return;
	}
	int positive = 0, negative = 0;
	/* Count the number of zero, positive, and negative elements */
	for(int i = 0; i < source->rows * source->cols; i++) {
		double delement = omxVectorElement(source, i);
		if (!R_finite(delement)) return;
		int element = (int) delement;
		if (element < 0) {
			/* bounds checking */
			if (element < - dimLength) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "index %d is out of bounds in '[' operator.", element);
				omxRaiseError(errstr);
				free(errstr);
				return;
			}
			negative++;
		} else if (element == 0) {
			// OK
		} else {
			/* bounds checking */
			if (element > dimLength) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "index %d is out of bounds in '[' operator.", element);
				omxRaiseError(errstr);
				free(errstr);
				return;
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
		return;
	}
	/* convert negative indices into a list of positive indices */
	if (negative > 0) {
		Eigen::ArrayXi track(dimLength);
		track.setZero();
		for(int i = 0; i < source->rows * source->cols; i++) {
			int element = (int) omxVectorElement(source, i);
			if (element < 0) {
				track[-element - 1]++;
			}
		}
		Eigen::DenseIndex numSelected = (track == 0).count();
		out.derived().resize(numSelected);
		int j = 0;
		for(int i = 0; i < dimLength; i++) {
			if(!track[i]) {
				out[j++] = i;
			}
		}
		return;
	}
	/* convert positive indices with offset of zero instead of one */
	if (positive > 0) {
		Eigen::Map< Eigen::ArrayXd > mask(source->data, source->rows * source->cols);
		Eigen::DenseIndex numSelected = (mask > 0).count();
		out.derived().resize(numSelected);
		int j = 0;
		for(int i = 0; i < numSelected; i++) {
			int element = (int) omxVectorElement(source, i);
			if (element > 0) {
				out[j++] = element - 1;
			}
		}
	}
}

static void omxMatrixExtract(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* rowMatrix = matList[1];
	omxMatrix* colMatrix = matList[2];

	if(OMX_DEBUG_ALGEBRA) { omxPrint(rowMatrix, "Row matrix: "); }
	if(OMX_DEBUG_ALGEBRA) { omxPrint(colMatrix, "Col matrix: "); }

	Eigen::ArrayXi rowIndices;
	matrixExtractIndices(rowMatrix, inMat->rows, rowIndices);

	Eigen::ArrayXi colIndices;
	matrixExtractIndices(colMatrix, inMat->cols, colIndices);

	omxResizeMatrix(result, rowIndices.size(), colIndices.size());

	for(int row = 0; row < (int) rowIndices.size(); row++) {
		for(int col = 0; col < (int) colIndices.size(); col++) {
			if(OMX_DEBUG_ALGEBRA) {
				mxLog("ALGEBRA: Matrix Extract: (%d, %d)[%d, %d] <- (%d, %d)[%d,%d].",
				      result->rows, result->cols, row, col,
				      (int) rowIndices.size(), (int) colIndices.size(), rowIndices[row], colIndices[col]);
			}
			double element = omxMatrixElement(inMat, rowIndices[row], colIndices[col]);
			omxSetMatrixElement(result, row, col, element);
		}
	}
}

static void omxMatrixSubtract(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	ensureElemConform(fc, matList, result);

	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	int rows = first->rows;
	int cols = first->cols;
	int size = rows * cols;

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

static void omxUnaryMinus(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxMatrixHorizCatOp(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrixHorizCat(matList, numArgs, result);
}

static void omxMatrixVertCatOp(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrixVertCat(matList, numArgs, result);
}

static void omxMatrixDeterminant(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);

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

static void omxMatrixTraceOp(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);
	omxMatrixTrace(matList, numArgs, result);
}

static void omxMatrixTotalSum(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);

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

static void omxMatrixTotalProduct(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);

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

static void omxMatrixArithmeticMean(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);
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

static void omxMatrixMinimum(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);

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

static void omxMatrixMaximum(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxResizeMatrix(result, 1, 1);

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

static void omxMatrixAbsolute(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = fabs(data[j]);
	}

}

static void omxMatrixDiagonal(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxMatrixFromDiagonal(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxElementCosine(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cos(data[j]);
	}

}

static void omxElementCosh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cosh(data[j]);
	}

}

static void omxElementSine(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sin(data[j]);
	}

}

static void omxElementSinh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sinh(data[j]);
	}

}

static void omxElementTangent(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = tan(data[j]);
	}

}

static void omxElementTanh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = tanh(data[j]);
	}

}

static void omxElementArcCosine(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = acos(data[j]);
	}
	
}

static void omxElementArcSine(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = asin(data[j]);
	}
	
}

static void omxElementArcTangent(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = atan(data[j]);
	}
	
}

static void omxElementAsinh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = asinh(data[j]);
	}
	
}

static void omxElementAcosh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = acosh(data[j]);
	}
	
}

static void omxElementAtanh(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = atanh(data[j]);
	}
	
}

static void omxElementExponent(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = exp(data[j]);
	}

}

static void omxElementNaturalLog(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = log(data[j]);
	}

}

static void omxElementSquareRoot(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sqrt(data[j]);
	}
}

static void omxElementPtoZ(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = Rf_qnorm5(data[j], 0, 1, 1, 0);
	}
}

static void omxElementLogPtoZ(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = Rf_qnorm5(data[j], 0, 1, 1, 1);
	}
}

static void omxElementLgamma(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int max = inMat->cols * inMat->rows;

	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = Rf_lgammafn(data[j]);
	}
}

static void omxElementLgamma1p(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = Rf_lgamma1p(data[j]);
	}
}

static void omxElementDbeta(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *a = matList[1];
	omxMatrix *b = matList[2];
	omxMatrix *ncp = matList[3];
	omxMatrix *give_log = matList[4];
	
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(a);
	omxEnsureColumnMajor(b);
	omxEnsureColumnMajor(ncp);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int aDataSize = a->rows * a->cols;
	int bDataSize = b->rows * b->cols;
	int ncpDataSize = ncp->rows * ncp->cols;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		if( Rf_sign(ncp->data[j%ncpDataSize]) == -1 ){
			data[j] = Rf_dbeta(data[j],a->data[j%aDataSize],b->data[j%bDataSize],give_log_arg);
		}
		else{
			data[j] = Rf_dnbeta(data[j],a->data[j%aDataSize],b->data[j%bDataSize],ncp->data[j%ncpDataSize],give_log_arg);
		}
	}
}

static void omxElementPbeta(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *a = matList[1];
	omxMatrix *b = matList[2];
	omxMatrix *ncp = matList[3];
	omxMatrix *lower_tail = matList[4];
	omxMatrix *give_log = matList[5];
	
	int lower_tail_arg = (int)(lower_tail->data[0] != 0);
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(a);
	omxEnsureColumnMajor(b);
	omxEnsureColumnMajor(ncp);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int aDataSize = a->rows * a->cols;
	int bDataSize = b->rows * b->cols;
	int ncpDataSize = ncp->rows * ncp->cols;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		if( Rf_sign(ncp->data[j%ncpDataSize]) == -1 ){
			data[j] = Rf_pbeta(data[j],a->data[j%aDataSize],b->data[j%bDataSize],lower_tail_arg,give_log_arg);
		}
		else{
			data[j] = Rf_pnbeta(data[j],a->data[j%aDataSize],b->data[j%bDataSize],ncp->data[j%ncpDataSize],lower_tail_arg,give_log_arg);
		}
	}
}

static void omxElementDpois(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *lambda = matList[1];
	omxMatrix *give_log = matList[2];
	
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(lambda);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int lambdaDataSize = lambda->rows * lambda->cols;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		data[j] = Rf_dpois(data[j],lambda->data[j%lambdaDataSize],give_log_arg);
	}
}

static void omxElementPpois(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *lambda = matList[1];
	omxMatrix *lower_tail = matList[2];
	omxMatrix *give_log = matList[3];
	
	int lower_tail_arg = (int)(lower_tail->data[0] != 0);
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(lambda);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int lambdaDataSize = lambda->rows * lambda->cols;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		data[j] = Rf_ppois(data[j],lambda->data[j%lambdaDataSize],lower_tail_arg,give_log_arg);
	}
}

static void omxElementDnbinom(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *size = matList[1];
	omxMatrix *prob = matList[2];
	omxMatrix *mu = matList[3];
	omxMatrix *give_log = matList[4];
	
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(size);
	omxEnsureColumnMajor(prob);
	omxEnsureColumnMajor(mu);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int sizeDataSize = size->rows * size->cols;
	int probDataSize = prob->rows * prob->cols;
	int muDataSize = mu->rows * mu->cols;
	int isSizeNeg=0, isProbNeg=0, isMuNeg=0, jumpVal=0;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	double sizecurr=0, probcurr=0, mucurr=0;
	
	for(int j = 0; j < inMatDataSize; j++) {
		sizecurr = size->data[j%sizeDataSize];
		probcurr = prob->data[j%probDataSize];
		mucurr = mu->data[j%muDataSize];
		isSizeNeg = (Rf_sign(sizecurr) == -1) ? 1L : 0L ;
		isProbNeg = (Rf_sign(probcurr) == -1) ? 3L : 0L ;
		isMuNeg = (Rf_sign(mucurr) == -1) ? 5L : 0L ;
		jumpVal = isSizeNeg + isProbNeg + isMuNeg;
		switch(jumpVal){
		case 1:
			data[j] = Rf_dnbinom(data[j],mucurr*probcurr/(1-probcurr),probcurr,give_log_arg);
			break;
		case 3:
			data[j] = Rf_dnbinom_mu(data[j],sizecurr,mucurr,give_log_arg);
			break;
		case 5:
			data[j] = Rf_dnbinom(data[j],sizecurr,probcurr,give_log_arg);
			break;
		default:
			Rf_warning("exactly one of arguments 'size', 'prob', and 'mu' must be negative (and therefore ignored)\n");
			data[j] = Rf_dnbinom(data[j],sizecurr,probcurr,give_log_arg); //let the R API handle bad inputs
		}
	}
}

static void omxElementPnbinom(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *size = matList[1];
	omxMatrix *prob = matList[2];
	omxMatrix *mu = matList[3];
	omxMatrix *lower_tail = matList[4];
	omxMatrix *give_log = matList[5];
	
	int lower_tail_arg = (int)(lower_tail->data[0] != 0);
	int give_log_arg = (int)(give_log->data[0] != 0);
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(size);
	omxEnsureColumnMajor(prob);
	omxEnsureColumnMajor(mu);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int sizeDataSize = size->rows * size->cols;
	int probDataSize = prob->rows * prob->cols;
	int muDataSize = mu->rows * mu->cols;
	int isSizeNeg=0, isProbNeg=0, isMuNeg=0, jumpVal=0;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	double sizecurr=0, probcurr=0, mucurr=0;
	
	for(int j = 0; j < inMatDataSize; j++) {
		sizecurr = size->data[j%sizeDataSize];
		probcurr = prob->data[j%probDataSize];
		mucurr = mu->data[j%muDataSize];
		isSizeNeg = (Rf_sign(sizecurr) == -1) ? 1L : 0L ;
		isProbNeg = (Rf_sign(probcurr) == -1) ? 3L : 0L ;
		isMuNeg = (Rf_sign(mucurr) == -1) ? 5L : 0L ;
		jumpVal = isSizeNeg + isProbNeg + isMuNeg;
		switch(jumpVal){
		case 1:
			data[j] = Rf_pnbinom(data[j],mucurr*probcurr/(1-probcurr),probcurr,lower_tail_arg,give_log_arg);
			break;
		case 3:
			data[j] = Rf_pnbinom_mu(data[j],sizecurr,mucurr,lower_tail_arg,give_log_arg);
			break;
		case 5:
			data[j] = Rf_pnbinom(data[j],sizecurr,probcurr,lower_tail_arg,give_log_arg);
			break;
		default:
			Rf_warning("exactly one of arguments 'size', 'prob', and 'mu' must be negative (and therefore ignored)\n");
			data[j] = Rf_pnbinom(data[j],sizecurr,probcurr,lower_tail_arg,give_log_arg); //let the R API handle bad inputs
		}
	}
}

static void omxElementBesselI(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *nu = matList[1];
	omxMatrix *expo = matList[2];
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(nu);
	omxEnsureColumnMajor(expo);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int nuDataSize = nu->rows * nu->cols;
	int expoDataSize = expo->rows * expo->cols;
	double expocurr;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		expocurr = (expo->data[j%expoDataSize] != 0) ? 2 : 1;
		data[j] = Rf_bessel_i(data[j],nu->data[j%nuDataSize],expocurr);
	}
}

static void omxElementBesselJ(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *nu = matList[1];
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(nu);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int nuDataSize = nu->rows * nu->cols;

	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		data[j] = Rf_bessel_j(data[j],nu->data[j%nuDataSize]);
	}
}

static void omxElementBesselK(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *nu = matList[1];
	omxMatrix *expo = matList[2];
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(nu);
	omxEnsureColumnMajor(expo);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int nuDataSize = nu->rows * nu->cols;
	int expoDataSize = expo->rows * expo->cols;
	double expocurr;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		expocurr = (expo->data[j%expoDataSize] != 0) ? 2 : 1;
		data[j] = Rf_bessel_k(data[j],nu->data[j%nuDataSize],expocurr);
	}
}

static void omxElementBesselY(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];
	omxMatrix *nu = matList[1];
	
	omxEnsureColumnMajor(inMat);
	omxEnsureColumnMajor(nu);
	
	int inMatDataSize = inMat->rows * inMat->cols;
	int nuDataSize = nu->rows * nu->cols;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < inMatDataSize; j++) {
		data[j] = Rf_bessel_y(data[j],nu->data[j%nuDataSize]);
	}
}

static void omxMatrixVech(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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

static void omxMatrixVechs(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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

static void omxRowVectorize(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxColVectorize(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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


static void omxSequenceGenerator(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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

static void omxMultivariateNormalIntegration(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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

	EigenVectorAdaptor Emean(means);

	if (lBoundMat->rows > 1 && lBoundMat->cols > 1) {
		omxRaiseErrorf("lbound must be a vector of length %d (not %dx%d)",
			       Emean.size(), lBoundMat->rows, lBoundMat->cols);
		return;
	}

	if (uBoundMat->rows > 1 && uBoundMat->cols > 1) {
		omxRaiseErrorf("ubound must be a vector of length %d (not %dx%d)",
			       Emean.size(), uBoundMat->rows, uBoundMat->cols);
		return;
	}

	OrdinalLikelihood ol;
	EigenMatrixAdaptor Ecov(cov);
	ol.setCovariance(Ecov, fc);
	ol.setMean(Emean);

	if(!R_finite(omxMatrixElement(cov, 0, 0))) {
		omxSetMatrixElement(result, 0, 0, NA_REAL);
		return;
	}

	EigenVectorAdaptor elb(lBoundMat);
	if (elb.size() != Emean.size()) {
		omxRaiseErrorf("lBound vector is length %d, not matching mean vector length %d",
			       elb.size(), Emean.size());
		return;
	}

	EigenVectorAdaptor eub(uBoundMat);
	if (eub.size() != Emean.size()) {
		omxRaiseErrorf("uBound vector is length %d, not matching mean vector length %d",
			       eub.size(), Emean.size());
		return;
	}

	omxSetMatrixElement(result, 0, 0, ol.likelihood(elb, eub));
}

static void omxAllIntegrationNorms(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* cov = matList[0];
	omxMatrix* means = matList[1];
	int nCols = cov->cols;
	int totalLevels = 1;
	std::vector<omxMatrix *> thresholdMats;
	thresholdMats.reserve(nCols);
	Eigen::ArrayXi numThresholds(nCols);
	Eigen::ArrayXi matNums(nCols);
	Eigen::ArrayXi thresholdCols(nCols);
	Eigen::ArrayXi currentThresholds(nCols);

	int currentMat = 0;

	for(int i = 0; i < nCols;) {
		if(OMX_DEBUG_ALGEBRA) {
			mxLog("All-part multivariate normal integration: Examining threshold column %d.", i);
		}
		thresholdMats.push_back(matList[currentMat+2]);		// Get the thresholds for this covariance column

		for(int j = 0; j < thresholdMats[currentMat]->cols; j++) {
			double ubound, lbound = omxMatrixElement(thresholdMats[currentMat], 0, j);
			if(ISNA(lbound)) {
				char *errstr = (char*) calloc(250, sizeof(char));
				sprintf(errstr, "Invalid lowest threshold for dimension %d of Allint.", j);
				omxRaiseError(errstr);
				free(errstr);
				return;
			}

			thresholdCols[i] = j;

			for(int k = 1; k < thresholdMats[currentMat]->rows; k++) {
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

	OrdinalLikelihood ol;
	EigenVectorAdaptor Emean(means);
	EigenMatrixAdaptor Ecov(cov);
	ol.setCovariance(Ecov, fc);
	ol.setMean(Emean);

	Eigen::VectorXd lBounds(nCols);
	Eigen::VectorXd uBounds(nCols);

	/* Set up first row */
	for(int j = (nCols-1); j >= 0; j--) {
		lBounds[j] = omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j]-1, thresholdCols[j]);
		uBounds[j] = omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j], thresholdCols[j]);
	}

	double likelihood = ol.likelihood(lBounds, uBounds);

	if (likelihood == 0.0) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Improper input to sadmvn.");
		omxRaiseError(errstr);
		free(errstr);
		return;
	}

	omxSetMatrixElement(result, 0, 0, likelihood);

	/* And repeat with increments for all other rows. */
	for(int i = 1; i < totalLevels; i++) {
		for(int j = (nCols-1); j >= 0; j--) {							// For each threshold set, starting from the fastest
			currentThresholds[j]++;									// Move to the next threshold set.
			if(currentThresholds[j] > numThresholds[j]) {			// Hit the end; cycle to the next.
				currentThresholds[j] = 1;
			}

			/* Update only the rows that need it. */
			lBounds[j] = omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j]-1, thresholdCols[j]);
			uBounds[j] = omxMatrixElement(thresholdMats[matNums[j]], currentThresholds[j], thresholdCols[j]);

			if (currentThresholds[j] != 1) break;
		}

		likelihood = ol.likelihood(lBounds, uBounds);

		if(likelihood == 0.0) {
			char *errstr = (char*) calloc(250, sizeof(char));
			sprintf(errstr, "Improper input to sadmvn.");
			omxRaiseError(errstr);
			free(errstr);
			return;
		}

		omxSetMatrixElement(result, i, 0, likelihood);
	}
}

static int omxComparePointerContentsHelper(const void* one, const void* two, void *ign) {
	double diff = (*(*(double**) two)) - (*(*(double**) one));
	if(diff > std::numeric_limits<double>::epsilon()) {
		return 1;
	} else if(diff < -std::numeric_limits<double>::epsilon()) {
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

static void omxRealEigenvalues(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxRealEigenvectors(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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
		if(fabs(WI[i]) > std::numeric_limits<double>::epsilon()) {									// If this is part of a conjugate pair
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

static void omxImaginaryEigenvalues(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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
		if(value > std::numeric_limits<double>::epsilon()) {
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

static void omxImaginaryEigenvectors(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxSelectRows(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int rows = inMat->rows;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(rows);
    
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
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    std::vector<int> zeros(inMat->cols);
    omxRemoveRowsAndColumns(result, toRemove.data(), zeros.data());

}

static void omxSelectCols(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int cols = inMat->cols;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(cols);

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
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    std::vector<int> zeros(inMat->rows);
    omxRemoveRowsAndColumns(result, zeros.data(), toRemove.data());
}

static void omxSelectRowsAndCols(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* selector = matList[1];

	int rows = inMat->rows;
	int cols = inMat->cols;
    int selectLength = selector->rows * selector->cols;
    Eigen::VectorXi toRemove(cols);

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
            toRemove[index] = 1;
        } else {
            toRemove[index] = 0;
        }
    }
    
    omxRemoveRowsAndColumns(result, toRemove.data(), toRemove.data());
}

static void omxCovToCor(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxCholesky(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
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

static void omxVechToMatrix(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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


static void omxVechsToMatrix(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
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

static void omxExponential(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	if (result->currentState->getWantStage() == FF_COMPUTE_INITIAL_FIT) {
		if (numArgs == 2) {
			Rf_warning("The second argument to omxExponential is ignored");
		}
	}

	omxMatrix* inMat = matList[0];
	if (inMat->rows != inMat->cols) Rf_error("omxExponential requires a symmetric matrix");
	omxEnsureColumnMajor(inMat);
	omxResizeMatrix(result, inMat->rows, inMat->cols);
	result->colMajor = true;

	expm_eigen(inMat->rows, inMat->data, result->data);
}

static void mxMatrixLog(FitContext *fc, omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	if (inMat->rows != inMat->cols) Rf_error("logm requires a symmetric matrix");
	omxEnsureColumnMajor(inMat);
	omxResizeMatrix(result, inMat->rows, inMat->cols);
	result->colMajor = true;

	logm_eigen(inMat->rows, inMat->data, result->data);
}
