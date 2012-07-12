/*
 *  Copyright 2007-2012 The OpenMx Project
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

#ifndef _OMX_ALGEBRA_FUNCTIONS_
#define _OMX_ALGEBRA_FUNCTIONS_

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <stdlib.h>
#include "omxMatrix.h"

/* Helpers */
void omxStandardizeCovMatrix(omxMatrix* cov, double* corList, double* weights); // TODO: Convert to Algebra?
void checkIncreasing(omxMatrix* om, int column);

/* Functional Wrappers */
void omxMatrixInvert(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxUnaryNegation(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxBinaryOr(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxBinaryAnd(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxBinaryGreaterThan(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxBinaryLessThan(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxBinaryApproxEquals(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixTranspose(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementPower(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixMult(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixElementMult(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxKroneckerProd(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxKroneckerPower(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxQuadraticProd(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementDivide(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixAdd(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixSubtract(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxUnaryMinus(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixHorizCat(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixVertCat(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixDeterminant(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixTrace(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixTotalSum(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixTotalProduct(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixArithmeticMean(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixMaximum(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixMinimum(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixAbsolute(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixDiagonal(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixFromDiagonal(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementCosine(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementCosh(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementSine(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementSinh(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementTangent(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementTanh(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementExponent(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementNaturalLog(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxElementSquareRoot(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixExtract(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixVech(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMatrixVechs(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxMultivariateNormalIntegration(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxAllIntegrationNorms(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxSequenceGenerator(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxRowVectorize(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxColVectorize(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxImaginaryEigenvectors(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxRealEigenvectors(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxRealEigenvalues(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxImaginaryEigenvalues(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxSelectRows(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxSelectCols(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxSelectRowsAndCols(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxAddOwnTranspose(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxCovToCor(omxMatrix** matList, int numArgs, omxMatrix* result);
void omxCholesky(omxMatrix** matList, int numArgs, omxMatrix* result);
#endif
