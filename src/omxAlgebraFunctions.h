#ifndef _OMX_ALGEBRA_FUNCTIONS_
#define _OMX_ALGEBRA_FUNCTIONS_

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxMatrix.h"

/* Functional Wrappers */
omxMatrix* omxMatrixInvert(omxMatrix* inMat);
omxMatrix* omxMatrixTranspose(omxMatrix* inMat);
omxMatrix* omxElementPower(omxMatrix* inMat, omxMatrix* power);
omxMatrix* omxMatrixMult(omxMatrix* preMul, omxMatrix* postMul);
omxMatrix* omxMatrixDot(omxMatrix* preDot, omxMatrix* postDot);
omxMatrix* omxKroneckerProd(omxMatrix* preMul, omxMatrix* postMul);
omxMatrix* omxQuadraticProd(omxMatrix* preMul, omxMatrix* postMul);
omxMatrix* omxElementDivide(omxMatrix* inMat, omxMatrix* divisor);
omxMatrix* omxMatrixAdd(omxMatrix* inMat, omxMatrix* addend);
omxMatrix* omxMatrixSubtract(omxMatrix* inMat, omxMatrix* subtrahend);
omxMatrix* omxUnaryMinus(omxMatrix* inMat);
omxMatrix* omxMatrixHorizCat(omxMatrix* matList, double numArgs);
omxMatrix* omxMatrixVertCat(omxMatrix* matList, double numArgs);

#endif