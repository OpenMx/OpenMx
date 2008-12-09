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

omxDataMatrix* omxMatrixTranspose(omxDataMatrix &inMat) {
	omxDataMatrix* result = new omxDataMatrix(inMat);
	result->transpose();
	return result;
}

omxDataMatrix* omxMatrixInvert(omxDataMatrix inMat)
{
	int ipiv[inMat.rows];
	int l = 0;
/* TODO: Does this interact with transposition? */
	omxDataMatrix* result = new omxDataMatrix(inMat);
	F77_CALL(dgetrf)(&(result->cols), &(result->rows), result->data, &(result->cols), ipiv, &l);
	return result;
}

omxDataMatrix* omxElementPower(omxDataMatrix inMat, omxDataMatrix power)
{
	omxDataMatrix* result = new omxDataMatrix(inMat);
	
	if(power.rows == 1 && power.cols ==1) {
		for(int j = 0; j < result->cols * result->rows; j++) {
			result->data[j] = pow(result->data[j], power.data[0]);
		}
	} else if(power.rows == inMat.rows && power.cols == inMat.cols) {
		for(int j = 0; j < result->rows; j++) 
			for(int k = 0; k < result->cols; k++)
				result->setElement(j, k, pow(result->element(j,k), power.element(j,k)));
	} else {
		error("Non-conformable matrices in elementPower.");
	}
};

omxDataMatrix* omxMatrixMult(omxDataMatrix preMul, omxDataMatrix postMul) {
	static double zero = 0.0;
	static double one = 1.0;
	
	/* Conformability Check! */
	if(preMul.lagging != postMul.leading) {
		error("Non-conformable matrices in Matrix Multiply.");
	}
	
	omxDataMatrix* newMat;
	newMat = new omxDataMatrix(preMul.rows, postMul.cols, TRUE);
	
	F77_CALL(dgemm)(&(preMul.majority), &(postMul.majority), &(preMul.leading), &(postMul.lagging), &(preMul.lagging), &one, preMul.data, &(postMul.leading), postMul.data, &(postMul.cols), &zero, newMat->data, &(newMat->leading));
	
	return newMat;
};

omxDataMatrix* omxMatrixDot(omxDataMatrix* preDot, omxDataMatrix* postDot) {};
omxDataMatrix* omxKroneckerProd(omxDataMatrix* preMul, omxDataMatrix* postMul) {};
omxDataMatrix* omxQuadraticProd(omxDataMatrix* preMul, omxDataMatrix* postMul) {};
omxDataMatrix* omxElementDivide(omxDataMatrix* inMat, omxDataMatrix* divisor) {};
omxDataMatrix* omxMatrixAdd(omxDataMatrix* inMat, omxDataMatrix* addend) {};
omxDataMatrix* omxMatrixSubtract(omxDataMatrix* inMat, omxDataMatrix* subtrahend) {};
omxDataMatrix* omxUnaryMinus(omxDataMatrix* inMat) {};
omxDataMatrix* omxMatrixHorizCat(omxDataMatrix* matList, double numArgs) {};
omxDataMatrix* omxMatrixVertCat(omxDataMatrix* matList, double numArgs) {};
