/***********************************************************
* 
*  omxAlgebraFunctions.h
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

void omxMatrixTranspose(omxMatrix *inMat, omxMatrix* result) {
	*result = *inMat;
	result->transpose();
}

void omxMatrixInvert(omxMatrix *inMat, omxMatrix* result)
{
	int ipiv[inMat->rows];
	int l = 0;
/* TODO: Does this interact with transposition? */
	*result = *inMat;
	F77_CALL(dgetrf)(&(result->cols), &(result->rows), result->data, &(result->leading), ipiv, &l);
}

void omxElementPower(omxMatrix *inMat, omxMatrix *power, omxMatrix* result)
{
	*result = *inMat;
	
	if(power->rows == 1 && power->cols ==1) {
		for(int j = 0; j < result->cols * result->rows; j++) {
			result->data[j] = pow(result->data[j], power->data[0]);
		}
	} else if(power->rows == inMat->rows && power->cols == inMat->cols) {
		for(int j = 0; j < result->rows; j++) 
			for(int k = 0; k < result->cols; k++)
				result->setElement(j, k, pow(result->element(j,k), power->element(j,k)));
	} else {
		error("Non-conformable matrices in elementPower.");
	}
};

void omxMatrixMult(omxMatrix *preMul, omxMatrix *postMul, omxMatrix* result) {

	if(OMX_DEBUG) { Rprintf("Multiplying two matrices.\n"); preMul->print("This one"); postMul->print("by this");}
	
	if(preMul == NULL || postMul == NULL) {
		error("Null matrix pointer detected.\n");
	}
	
	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul->cols != postMul->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	if(result->rows != preMul->rows || result->cols != preMul->cols)
		result->resize(preMul->rows, postMul->cols);
	
	
	/* For debugging */
	if(OMX_DEBUG) {
		result->print("NewMatrix");
		Rprintf("DGEMM: %c, %c, %d, %d, %d, %d\n", *(preMul->majority), *(postMul->majority), (preMul->rows), (postMul->cols), preMul->leading, postMul->leading);
		Rprintf("Note: majorityList is [%c, %c]\n", *(preMul->majorityList), *(preMul->majorityList + 1));
	}
	
	/* The call itself */
	F77_CALL(dgemm)((preMul->majority), (postMul->majority), &(preMul->rows), &(postMul->cols), &(preMul->cols), &one, preMul->data, &(preMul->leading), postMul->data, &(postMul->leading), &zero, result->data, &(result->leading));

};

void omxMatrixDot(omxMatrix* preDot, omxMatrix* postDot, omxMatrix* result) {
// Actually, this is element-by-element multiplication.
	
	/* Conformability Check! */
	if(preDot->cols != postDot->cols || preDot->rows != postDot->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	result = postDot;
	
	int max = preDot->cols * preDot->rows;
	
	/* The call itself */
	//F77_CALL(dgemm)((preDot->majority), (result->majority), &(preDot->rows), &(result->cols), &(preDot->cols), &zero, preDot->data, &(preDot->leading), result->data, &(result->leading), &one, result->data, &(result->leading));
	for(int j = 0; j < max; j++) {
		result->data[j] = preDot->data[j] * postDot->data[j];
	}


};
void omxKroneckerProd(omxMatrix* preMul, omxMatrix* postMul, omxMatrix* result) {

	int rows = preMul->rows * postMul->rows;
	int cols = preMul->cols * postMul->cols;

	if(result->rows != rows || result->cols != cols)
		result->resize(rows, cols);

	for(int preRow = 0; preRow < preMul->rows; preRow++)
		for(int postRow = 0; postRow < postMul->rows; postRow++)
			for(int preCol = 0; preCol < preMul->cols; preCol++)
				for(int postCol = 0; postCol < postMul->cols; postCol++)
					result->setElement(preRow*postMul->rows + postRow, preCol*postMul->cols + postCol, preMul->element(preRow, preCol)*postMul->element(postRow, postCol));
					
};

void omxQuadraticProd(omxMatrix* preMul, omxMatrix* postMul, omxMatrix* result) {
	/* A %&% B = ABA' */

	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul->cols != postMul->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	omxMatrix* intermediate = new omxMatrix(preMul->rows, postMul->cols);
		
	if(result->rows != preMul->rows || result->cols != preMul->rows)
		result->resize(preMul->cols, preMul->cols);
	
	
	/* The call itself */
	F77_CALL(dgemm)((preMul->majority), (postMul->majority), &(preMul->rows), &(postMul->cols), &(preMul->cols), &one, preMul->data, &(preMul->leading), postMul->data, &(postMul->leading), &zero, intermediate->data, &(intermediate->leading));
	F77_CALL(dgemm)((intermediate->majority), (preMul->minority), &(intermediate->rows), &(preMul->cols), &(intermediate->cols), &one, intermediate->data, &(intermediate->leading), preMul->data, &(preMul->leading), &zero, result->data, &(result->leading));
	
};

void omxElementDivide(omxMatrix* inMat, omxMatrix* divisor, omxMatrix* result) {
	
	
	/* Conformability Check! */
	if(inMat->cols != divisor->cols || inMat->rows != divisor->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	result = divisor;
	
	int max = inMat->cols * inMat->rows;
	
	for(int j = 0; j < max; j++) {
		result->data[j] = inMat->data[j] / divisor->data[j];
	}
	
};

void omxMatrixAdd(omxMatrix* inMat, omxMatrix* addend, omxMatrix* result) {
	
	/* Conformability Check! */
	if(inMat->cols != addend->cols || inMat->rows != addend->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	result = addend;
	
	int max = inMat->cols * inMat->rows;
	
	for(int j = 0; j < max; j++) {
		result->data[j] = inMat->data[j] + addend->data[j];
	}

};

void omxMatrixSubtract(omxMatrix* inMat, omxMatrix* subtrahend, omxMatrix* result) {	/* Conformability Check! */
	if(inMat->cols != subtrahend->cols || inMat->rows != subtrahend->rows) 
		error("Non-conformable matrices in Matrix Multiply.");
		
	result = subtrahend;
	
	int max = inMat->cols * inMat->rows;
	
	//F77_CALL(dgemm)((inMat->majority), (result->majority), &(inMat->rows), &(result->cols), &(inMat->cols), &zero, inMat->data, &(inMat->leading), result->data, &(result->leading), &one, result->data, &(result->leading));   // TODO: Compare performance on BLAS vs for.
	
	for(int j = 0; j < max; j++) {
		result->data[j] = inMat->data[j] - subtrahend->data[j];
	}

};

void omxUnaryMinus(omxMatrix* inMat, omxMatrix* result) {
	
	int max = inMat->cols * inMat->rows;
	
	result = inMat;
	
	for(int j = 0; j < max; j++) {
		result->data[j] = -result->data[j];
	}
	
};

void omxMatrixHorizCat(omxMatrix** matList, double numArgs, omxMatrix* result) {

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
		result->resize(totalRows, totalCols);
	}
	
	for(int j = 0; j < numArgs; j++) {
		for(int k = 0; k < matList[j]->cols; j++) {
			for(int l = 0; l < totalRows; l++) {		// Gotta be a faster way to do this.
				result->setElement(l, currentCol, matList[j]->element(l, k));
			}
			currentCol++;
		}
	}

};

void omxMatrixVertCat(omxMatrix** matList, double numArgs, omxMatrix* result) {
	
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
		result->resize(totalRows, totalCols);
	}
	
	for(int j = 0; j < numArgs; j++) {
		for(int k = 0; k < matList[j]->rows; j++) {
			for(int l = 0; l < totalCols; l++) {		// Gotta be a faster way to do this.
				result->setElement(currentRow, l, matList[j]->element(k, l));
			}
			currentRow++;
		}
	}

};
