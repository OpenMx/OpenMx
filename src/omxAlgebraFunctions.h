#ifndef _OMX_ALGEBRA_FUNCTIONS_
#define _OMX_ALGEBRA_FUNCTIONS_

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDataMatrix.h"

/* Functional Wrappers */
omxDataMatrix* omxMatrixInvert(omxDataMatrix* inMat);
omxDataMatrix* omxMatrixTranspose(omxDataMatrix* inMat);
omxDataMatrix* omxElementPower(omxDataMatrix* inMat, omxDataMatrix* power);
omxDataMatrix* omxMatrixMult(omxDataMatrix* preMul, omxDataMatrix* postMul);
omxDataMatrix* omxMatrixDot(omxDataMatrix* preDot, omxDataMatrix* postDot);
omxDataMatrix* omxKroneckerProd(omxDataMatrix* preMul, omxDataMatrix* postMul);
omxDataMatrix* omxQuadraticProd(omxDataMatrix* preMul, omxDataMatrix* postMul);
omxDataMatrix* omxElementDivide(omxDataMatrix* inMat, omxDataMatrix* divisor);
omxDataMatrix* omxMatrixAdd(omxDataMatrix* inMat, omxDataMatrix* addend);
omxDataMatrix* omxMatrixSubtract(omxDataMatrix* inMat, omxDataMatrix* subtrahend);
omxDataMatrix* omxUnaryMinus(omxDataMatrix* inMat);
omxDataMatrix* omxMatrixHorizCat(omxDataMatrix* matList, double numArgs);
omxDataMatrix* omxMatrixVertCat(omxDataMatrix* matList, double numArgs);

#endif