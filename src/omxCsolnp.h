/*
 *  Copyright 2007-2012 The OpenMx Project
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

#ifndef _OMX_CSOLNP_SPECIFIC_H
#define _OMX_CSOLNP_SPECIFIC_H

#include "omxMatrix.h"
#include "matrix.h"

struct Matrix fillMatrix(int cols, int rows, double* array);

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc, int verbose);

//void omxNPSOLConfidenceIntervals(double *f, double *x, double *g, double *R, int ciMaxIterations);

Param_Obj solnp( Matrix solPars, double (*solFun)( Matrix),  Matrix solEqB, Matrix (*solEqBFun)(Matrix), Matrix (*myineqFun)( Matrix) , Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

double csolnpObjectiveFunction(Matrix myPars);
struct Matrix csolnpEqualityFunction(Matrix myEqBFun_arg);
//struct Matrix csolnpEqB(Matrix* EqB_arg[]);
//struct Matrix csolnpEqB();
//struct Matrix csolnpIneqUB(Matrix* IneqUB_arg[]);
//struct Matrix csolnpIneqUB();
//struct Matrix csolnpIneqLB(Matrix* IneqLB_arg[]);
//struct Matrix csolnpIneqLB();
struct Matrix csolnpIneqFun(Matrix myPars);

#endif // #define _OMX_CSOLNP_SPECIFIC_H
