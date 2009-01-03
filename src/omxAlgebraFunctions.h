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
void omxMatrixInvert(omxMatrix* inMat, omxMatrix* result);
void omxMatrixTranspose(omxMatrix* inMat, omxMatrix* result);
void omxElementPower(omxMatrix* inMat, omxMatrix* power, omxMatrix* result);
void omxMatrixMult(omxMatrix* preMul, omxMatrix* postMul, omxMatrix* result);
void omxMatrixDot(omxMatrix* preDot, omxMatrix* postDot, omxMatrix* result);
void omxKroneckerProd(omxMatrix* preMul, omxMatrix* postMul, omxMatrix* result);
void omxQuadraticProd(omxMatrix* preMul, omxMatrix* postMul, omxMatrix* result);
void omxElementDivide(omxMatrix* inMat, omxMatrix* divisor, omxMatrix* result);
void omxMatrixAdd(omxMatrix* inMat, omxMatrix* addend, omxMatrix* result);
void omxMatrixSubtract(omxMatrix* inMat, omxMatrix* subtrahend, omxMatrix* result);
void omxUnaryMinus(omxMatrix* inMat, omxMatrix* result);
void omxMatrixHorizCat(omxMatrix* matList, double numArgs, omxMatrix* result);
void omxMatrixVertCat(omxMatrix* matList, double numArgs, omxMatrix* result);

#endif