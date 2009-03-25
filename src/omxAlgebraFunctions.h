/*
 *  Copyright 2007-2009 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

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
void omxMatrixHorizCat(omxMatrix** matList, double numArgs, omxMatrix* result);
void omxMatrixVertCat(omxMatrix** matList, double numArgs, omxMatrix* result);

#endif
