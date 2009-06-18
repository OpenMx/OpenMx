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

#include "omxAlgebraFunctions.h"
#include "omxMatrix.h"

void omxMatrixTranspose(omxMatrix** matList, int numArgs, omxMatrix* result) {
	
	omxMatrix* inMat = matList[0];
	
	omxCopyMatrix(result, inMat);
	result->colMajor = !result->colMajor;
	int rowtemp = result->rows;
	result->rows = result->cols;
	result->cols = rowtemp;
	omxComputeMatrixHelper(result);
}

void omxMatrixInvert(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int ipiv[inMat->rows], lwork;
	lwork = 4 * inMat->rows * inMat->cols;
	double work[lwork];
	int l = 0;

	omxCopyMatrix(result, inMat);
	F77_CALL(dgetrf)(&(result->cols), &(result->rows), result->data, &(result->leading), ipiv, &l);
	if(l != 0) {
		error("Attempted to invert non-invertable matrix.");
	}
	F77_CALL(dgetri)(&(result->cols), result->data, &(result->leading), ipiv, work, &lwork, &l);
	
}

void omxElementPower(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* power = matList[1];
	
	/* Conformability Check! */
	if(inMat->cols != power->cols || inMat->rows != power->rows) 
		error("Non-conformable matrices in Element Division.");	

	int rows = inMat->rows;
	int cols = inMat->cols;
	
	if(result->rows != rows || result->cols != cols)
		omxResizeMatrix(result, rows, cols, FALSE);	

	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			omxSetMatrixElement(result, i, j, 
				pow(omxMatrixElement(inMat, i, j), 
					omxMatrixElement(power, i, j)));
		}
	}
		
}

void omxMatrixMult(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];

	if(OMX_DEBUG) { Rprintf("Multiplying two matrices.\n");}
	
	if(preMul == NULL || postMul == NULL) {
		error("Null matrix pointer detected.\n");
	}
	
	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul->cols != postMul->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	if(result->rows != preMul->rows || result->cols != preMul->cols)
		omxResizeMatrix(result, preMul->rows, postMul->cols, FALSE);
	
	
	/* For debugging */
	if(OMX_DEBUG) {
		omxPrintMatrix(result, "NewMatrix");
		Rprintf("DGEMM: %c, %c, %d, %d, %d, %f, %0x %d %0x %d %f %0x %d\n", *(preMul->majority), *(postMul->majority), (preMul->rows), (postMul->cols), (preMul->cols), one, (preMul->data), (preMul->leading), (postMul->data), (postMul->leading), zero, (result->data), (result->leading)); 
	}
	
	/* The call itself */
	F77_CALL(dgemm)((preMul->majority), (postMul->majority), &(preMul->rows), &(postMul->cols), &(preMul->cols), &one, preMul->data, &(preMul->leading), postMul->data, &(postMul->leading), &zero, result->data, &(result->leading));
	result->colMajor = TRUE;
	omxComputeMatrixHelper(result);

}

void omxMatrixElementMult(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* first = matList[0];
	omxMatrix* second = matList[1];
	
	/* Conformability Check! */
	if(first->cols != second->cols || first->rows != second->rows) 
		error("Non-conformable matrices in Element Multiplication.");
	
	omxCopyMatrix(result, first);
	
	int rows = first->rows;
	int cols = first->cols;
	
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			omxSetMatrixElement(result, i, j, 
				omxMatrixElement(result, i, j) * 
				omxMatrixElement(second, i, j));
		}
	}


}
void omxKroneckerProd(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];
	
	int rows = preMul->rows * postMul->rows;
	int cols = preMul->cols * postMul->cols;

	if(result->rows != rows || result->cols != cols)
		omxResizeMatrix(result, rows, cols, FALSE);

	for(int preRow = 0; preRow < preMul->rows; preRow++)
		for(int postRow = 0; postRow < postMul->rows; postRow++)
			for(int preCol = 0; preCol < preMul->cols; preCol++)
				for(int postCol = 0; postCol < postMul->cols; postCol++)
					omxSetMatrixElement(result, preRow*postMul->rows + postRow, 
						preCol*postMul->cols + postCol, 
						omxMatrixElement(preMul, preRow, preCol) * omxMatrixElement(postMul, postRow, postCol));
	
}

void omxQuadraticProd(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* preMul = matList[0];
	omxMatrix* postMul = matList[1];
	/* A %&% B = ABA' */

	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul->cols != postMul->rows) 
		error("Non-conformable matrices in Matrix Quadratic Product.");
		
	omxMatrix* intermediate = NULL;
	omxInitMatrix(intermediate, preMul->rows, postMul->cols, TRUE, preMul->currentState);
		
	if(result->rows != preMul->rows || result->cols != preMul->rows)
		omxResizeMatrix(result, preMul->cols, preMul->cols, FALSE);
	
	
	/* The call itself */
	F77_CALL(dgemm)((preMul->majority), (postMul->majority), &(preMul->rows), &(postMul->cols), &(preMul->cols), &one, preMul->data, &(preMul->leading), postMul->data, &(postMul->leading), &zero, intermediate->data, &(intermediate->leading));
	F77_CALL(dgemm)((intermediate->majority), (preMul->minority), &(intermediate->rows), &(preMul->cols), &(intermediate->cols), &one, intermediate->data, &(intermediate->leading), preMul->data, &(preMul->leading), &zero, result->data, &(result->leading));


}

void omxElementDivide(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* divisor = matList[1];
	
	/* Conformability Check! */
	if(inMat->cols != divisor->cols || inMat->rows != divisor->rows) 
		error("Non-conformable matrices in Element Division.");
		
	omxCopyMatrix(result, inMat);
	
	int rows = inMat->rows;
	int cols = inMat->cols;
	
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			omxSetMatrixElement(result, i, j, 
				omxMatrixElement(result, i, j) / 
				omxMatrixElement(divisor, i, j));
		}
	}
}

void omxMatrixAdd(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* addend = matList[1];
	
	/* Conformability Check! */
	if(inMat->cols != addend->cols || inMat->rows != addend->rows) 
		error("Non-conformable matrices in Matrix Addition.");
		
	omxCopyMatrix(result, inMat);
	
	int rows = inMat->rows;
	int cols = inMat->cols;
	
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			omxSetMatrixElement(result, i, j, 
				omxMatrixElement(result, i, j) +  
				omxMatrixElement(addend, i, j));
		}
	}	
}

void omxMatrixSubtract(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	omxMatrix* subtrahend = matList[1];
	
	if(inMat->cols != subtrahend->cols || inMat->rows != subtrahend->rows) 
		error("Non-conformable matrices in Matrix Subtract.");
		
	omxCopyMatrix(result, inMat);
	
	int max = inMat->cols * inMat->rows;
	
	//F77_CALL(dgemm)((inMat->majority), (result->majority), &(inMat->rows), &(result->cols), &(inMat->cols), &zero, inMat->data, &(inMat->leading), result->data, &(result->leading), &one, result->data, &(result->leading));   // TODO: Compare performance on BLAS vs for.
	if(OMX_DEBUG) {
		omxPrintMatrix(subtrahend, "Subtracting");
		omxPrintMatrix(inMat, "From");
	}

	int rows = inMat->rows;
	int cols = inMat->cols;
	
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			omxSetMatrixElement(result, i, j, 
				omxMatrixElement(result, i, j) -  
				omxMatrixElement(subtrahend, i, j));
		}
	}

	if(OMX_DEBUG) {
		omxPrintMatrix(inMat, "And Got");
	}
}

void omxUnaryMinus(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = - data[j];
	}
	
}

void omxMatrixHorizCat(omxMatrix** matList, int numArgs, omxMatrix* result) {

	int totalRows = 0, totalCols = 0, currentCol=0;
	
	if(numArgs == 0) return;
	
	totalRows = matList[0]->rows;			// Assumed constant.  Assert this below.

	for(int j = 0; j < numArgs; j++) {
		if(totalRows != matList[j]->rows) {
			error("Non-conformable matrices in horizontal concatenation.");
		}
		totalCols += matList[j]->cols;
	}
	
	if(result->rows != totalRows || result->cols != totalCols) {
		omxResizeMatrix(result, totalRows, totalCols, FALSE);
	}
	
	for(int j = 0; j < numArgs; j++) {
		for(int k = 0; k < matList[j]->cols; j++) {
			for(int l = 0; l < totalRows; l++) {		// Gotta be a faster way to do this.
				omxSetMatrixElement(result, l, currentCol, omxMatrixElement(matList[j], l, k));
			}
			currentCol++;
		}
	}

}

void omxMatrixVertCat(omxMatrix** matList, int numArgs, omxMatrix* result) {
	
	int totalRows = 0, totalCols = 0, currentRow=0;
	
	if(numArgs == 0) return;
	
	totalCols = matList[0]->cols;			// Assumed constant.  Assert this below.

	for(int j = 0; j < numArgs; j++) {
		if(totalCols != matList[j]->cols) {
			error("Non-conformable matrices in horizontal concatenation.");
		}
		totalRows += matList[j]->rows;
	}
	
	if(result->rows != totalRows || result->cols != totalCols) {
		omxResizeMatrix(result, totalRows, totalCols, FALSE);
	}
	
	for(int j = 0; j < numArgs; j++) {
		for(int k = 0; k < matList[j]->rows; j++) {
			for(int l = 0; l < totalCols; l++) {		// Gotta be a faster way to do this.
				omxSetMatrixElement(result, l, currentRow, omxMatrixElement(matList[j], k, l));
			}
			currentRow++;
		}
	}

}

void omxMatrixDeterminant(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
		error("NYI: Not yet Implemented.\n");
	}

void omxMatrixTrace(omxMatrix** matList, int numArgs, omxMatrix* result)
	{
		omxMatrix* inMat = matList[0];

		/* Consistency check: */
		if(result->rows != 1 || result->cols != 1) {
			omxResizeMatrix(result, 1, 1, FALSE);
		}

		if(inMat->rows != inMat->cols) {
			error("Non-square matrix in Trace().\n");
		}
		
		double trace = 0.0;

		/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
		for(int j = 0; j < inMat->rows; j++) {
			trace += omxMatrixElement(inMat, j, j);
		}

		omxSetMatrixElement(result, 0, 0, trace);
	}

void omxMatrixTotalSum(omxMatrix** matList, int numArgs, omxMatrix* result) {
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1, FALSE);
	}
	
	double sum = 0.0;
	
	/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matlength = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matlength; k++) {
			sum += data[k];
		}
	}
	
	omxSetMatrixElement(result, 0, 0, sum);	
}

void omxMatrixTotalProduct(omxMatrix** matList, int numArgs, omxMatrix* result) {
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1, FALSE);
	}
	
	double product = 1.0;
	
	/* Note: This algorithm is numerically unstable.  Sorry, dudes. */
	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matlength = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matlength; k++) {
			product *= data[k];
		}
	}
	
	omxSetMatrixElement(result, 0, 0, product);	
}

void omxMatrixMinimum(omxMatrix** matList, int numArgs, omxMatrix* result){
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1, FALSE);
	}
	
	double min = DBL_MAX; // DBL_MAX is the maximum possible DOUBLE value, usually 10e37.  
						  // We could change this to use NPSOL's INFINITY, but why bother?
	
	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;
		int matlength = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matlength; k++) {
			if(data[k] < min) min = data[k];
		}
	}
	
	omxSetMatrixElement(result, 0, 0, min);	
}

void omxMatrixMaximum(omxMatrix** matList, int numArgs, omxMatrix* result){
	/* Consistency check: */
	if(result->rows != 1 || result->cols != 1) {
		omxResizeMatrix(result, 1, 1, FALSE);
	}
	
	double max = -DBL_MAX;

	for(int j = 0; j < numArgs; j++) {
		double* data = matList[j]->data;		
		int matlength = matList[j]->rows * matList[j]->cols;
		for(int k = 0; k < matlength; k++) {
			if(data[k] > max) max = data[k];
		}
	}
	
	omxSetMatrixElement(result, 0, 0, max);	
}

void omxMatrixAbsolute(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = fabs(data[j]);
	}
	
}

void omxElementCosine(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cos(data[j]);
	}
	
}

void omxElementCosh(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = cosh(data[j]);
	}
	
}

void omxElementSine(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sin(data[j]);
	}
	
}

void omxElementSinh(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];

	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sinh(data[j]);
	}
	
}

void omxElementTangent(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);

	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = tan(data[j]);
	}
	
}

void omxElementTanh(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;	
	for(int j = 0; j < max; j++) {
		data[j] = tanh(data[j]);
	}
	
}

void omxElementExponent(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;	
	for(int j = 0; j < max; j++) {
		data[j] = exp(data[j]);
	}
	
}

void omxElementNaturalLog(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix* inMat = matList[0];
	
	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);

	double* data = result->data;	
	for(int j = 0; j < max; j++) {
		data[j] = log(data[j]);
	}
	
}

void omxElementSquareRoot(omxMatrix** matList, int numArgs, omxMatrix* result)
{
	omxMatrix *inMat = matList[0];

	int max = inMat->cols * inMat->rows;
	
	omxCopyMatrix(result, inMat);
	
	double* data = result->data;
	for(int j = 0; j < max; j++) {
		data[j] = sqrt(data[j]);
	}

}
