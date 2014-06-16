/* Copyright (c) 2001, 2002 Enthought, Inc.
 *  All rights reserved.
 *
 *  Copyright (c) 2003-2009 SciPy Developers.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *    a. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *    b. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    c. Neither the name of the Enthought nor the names of its contributors
 *       may be used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
 *  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 *  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 *  DAMAGE.
 */

// Translated from the function linalg.expm in the SciPy library

#include "omxAlgebraFunctions.h"
#include "omxMatrix.h"
#include "omxBLAS.h"

void matrixExponential(omxMatrix* inMat, int order, omxMatrix* result) {
	int nrow = inMat->rows;
	int ncol = inMat->cols;
	if (nrow != ncol) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Non-square matrix for matrix exponential.\n");
		omxRaiseError(errstr);
		free(errstr);
	}

	if((nrow != result->rows) || (ncol != result->cols)) {
		omxResizeMatrix(result, nrow, ncol);
	}

	double maxRowSum = 0.0;
	for(int col = 0; col < ncol; col++) {
		maxRowSum += fabs(omxMatrixElement(inMat, 0, col));
	}
	for(int row = 1; row < nrow; row++) {
		double rowSum = 0.0;
		for(int col = 0; col < ncol; col++) {
			rowSum += fabs(omxMatrixElement(inMat, row, col));
		}
		if (rowSum > maxRowSum) {
			maxRowSum = rowSum;
		}
	}
	if (maxRowSum == 0.0) {
		for(int i = 0; i < nrow; i++) {
			for(int j = 0; j < ncol; j++) {
				if (i == j) {
					omxSetMatrixElement(result, i, j, 1.0);
				} else {
					omxSetMatrixElement(result, i, j, 0.0);
				}
			}
		}
		return;
	}
	static double one = 1.0;
	static double zero = 0.0;
	double e = floor(log2(maxRowSum));
	int j = (int) (e + 1);
	if (j < 0) j = 0;
	int Rf_length = nrow * ncol;

	omxMatrix *normInMat = omxInitMatrix(nrow, ncol, 1, inMat->currentState);
	omxCopyMatrix(normInMat, inMat);

	double multipleOfTwo = pow(2.0, j);
	for(int i = 0; i < Rf_length; i++) {
		omxSetVectorElement(normInMat, i, omxVectorElement(normInMat, i) / multipleOfTwo);
	}

	omxMatrix *tempA = omxInitMatrix(nrow, ncol, 1, inMat->currentState);
	omxMatrix *tempResult = omxInitMatrix(nrow, ncol, 1, inMat->currentState);
	omxMatrix *N = omxNewIdentityMatrix(nrow, inMat->currentState);
	omxMatrix *D = omxNewIdentityMatrix(nrow, inMat->currentState);
	for(int row = 0; row < nrow; row++) {
		for(int col = 0; col < ncol; col++) {
			double element = omxMatrixElement(normInMat, row, col);
			omxSetMatrixElement(N, row, col, omxMatrixElement(N, row, col) + 0.5 * element);
			omxSetMatrixElement(D, row, col, omxMatrixElement(D, row, col) - 0.5 * element);
		}
	}
	for(int i = 0; i < Rf_length; i++) {
		omxSetVectorElement(normInMat, i, omxVectorElement(inMat, i) / multipleOfTwo);
	}

	omxCopyMatrix(tempA, normInMat);
	double constant = 0.5;


	for(int k = 2; k <= order; k++) {
		constant *= (order - k + 1.0) / (k * (2.0 * order - k + 1.0));
		F77_CALL(omxunsafedgemm)((normInMat->majority), (tempA->majority), &(normInMat->rows), 
			&(tempA->cols), &(normInMat->cols), &one, normInMat->data, &(normInMat->leading), 
			tempA->data, &(tempA->leading), &zero, tempResult->data, &(tempResult->leading));
		omxCopyMatrix(tempA, tempResult);
		for(int i = 0; i < Rf_length; i++) {
			omxSetVectorElement(tempResult, i, omxVectorElement(tempResult, i) * constant);
		}
		for(int i = 0; i < nrow; i++) {
			for(int j = 0; j < ncol; j++) {
				omxSetMatrixElement(N, i, j, omxMatrixElement(N, i, j) + omxMatrixElement(tempResult, i, j));
				if (k % 2 == 0) {
					omxSetMatrixElement(D, i, j, omxMatrixElement(D, i, j) + omxMatrixElement(tempResult, i, j));
				} else {
					omxSetMatrixElement(D, i, j, omxMatrixElement(D, i, j) - omxMatrixElement(tempResult, i, j));
				}
			}
		}	
	}


	Eigen::VectorXi ipiv(nrow);
	int l = 0;
	const char trans = 'N';

	F77_CALL(dgetrf)(&nrow, &nrow, D->data, &(D->leading), ipiv.data(), &l);

	F77_CALL(dgetrs)(&trans, &nrow, &nrow, D->data, &(D->leading), ipiv.data(), 
		N->data, &(N->leading), &l);

	if (j > 0) {  
		for(int k = 1; k <= j; k++) {
			F77_CALL(omxunsafedgemm)((N->majority), (N->majority), &nrow, 
				&nrow, &nrow, &one, N->data, &(N->leading), 
				N->data, &(N->leading), &zero, result->data, &(result->leading));
			omxCopyMatrix(N, result);		
		}
	} else {
		omxCopyMatrix(result, N);
	}

	omxFreeMatrix(normInMat);
	omxFreeMatrix(tempA);
	omxFreeMatrix(tempResult);
	omxFreeMatrix(N);
	omxFreeMatrix(D);
}


void omxExponential(omxMatrix** matList, int numArgs, omxMatrix* result) {
	omxMatrix* inMat = matList[0];
	matrixExponential(inMat, 7, result);
}

void omxExponentialOrder(omxMatrix** matList, int numArgs, omxMatrix* result) {
	omxMatrix* inMat = matList[0];
	omxMatrix* orderMat = matList[1];
	int order = (int) omxVectorElement(orderMat, 0);
	matrixExponential(inMat, order, result);
}



