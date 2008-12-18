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

omxMatrix* omxMatrixTranspose(omxMatrix *inMat) {
	omxMatrix* result = new omxMatrix(*inMat);
	result->transpose();
	return result;
}

omxMatrix* omxMatrixInvert(omxMatrix *inMat)
{
	int ipiv[inMat->rows];
	int l = 0;
/* TODO: Does this interact with transposition? */
	omxMatrix* result = new omxMatrix(*inMat);
	F77_CALL(dgetrf)(&(result->cols), &(result->rows), result->data, &(result->cols), ipiv, &l);
	return result;
}

omxMatrix* omxElementPower(omxMatrix *inMat, omxMatrix *power)
{
	omxMatrix* result = new omxMatrix(*inMat);
	
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

omxMatrix* omxMatrixMult(omxMatrix *preMul, omxMatrix *postMul) {
	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul->lagging != postMul->leading) {
		error("Non-conformable matrices in Matrix Multiply.");
	}
	
	omxMatrix* newMat;
	newMat = new omxMatrix(preMul->rows, postMul->cols, TRUE);
	
	F77_CALL(dgemm)(&(preMul->majority), &(postMul->majority), &(preMul->leading), &(postMul->lagging), &(preMul->lagging), &one, preMul->data, &(postMul->leading), postMul->data, &(postMul->cols), &zero, newMat->data, &(newMat->leading));
	
	return newMat;
};

omxMatrix* omxMatrixDot(omxMatrix* preDot, omxMatrix* postDot) {};
omxMatrix* omxKroneckerProd(omxMatrix* preMul, omxMatrix* postMul) {};
omxMatrix* omxQuadraticProd(omxMatrix* preMul, omxMatrix* postMul) {};
omxMatrix* omxElementDivide(omxMatrix* inMat, omxMatrix* divisor) {};
omxMatrix* omxMatrixAdd(omxMatrix* inMat, omxMatrix* addend) {};
omxMatrix* omxMatrixSubtract(omxMatrix* inMat, omxMatrix* subtrahend) {};
omxMatrix* omxUnaryMinus(omxMatrix* inMat) {};
omxMatrix* omxMatrixHorizCat(omxMatrix* matList, double numArgs) {};
omxMatrix* omxMatrixVertCat(omxMatrix* matList, double numArgs) {};
