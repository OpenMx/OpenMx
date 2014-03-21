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

//typedef double (*solFun_t)(struct Matrix myPars, int verbose);
//typedef Matrix (*solEqBFun_t)(struct Matrix myPars, int verbose);
//typedef Matrix (*solIneqFun_t)(struct Matrix myPars, int verbose);

struct Matrix fillMatrix(int cols, int rows, double* array);

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc, int *inform_out, int *iter_out, FreeVarGroup *freeVarGroup, int verbose);

void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose);

Param_Obj solnp(Matrix solPars, double (*solFun)(Matrix, int),
                Matrix solEqB, Matrix (*solEqBFun)(int),  Matrix (*myineqFun)(int),
                Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,
                Matrix solctrl, bool debugToggle, int verbose);

double csolnpObjectiveFunction(Matrix myPars, int verbose);
double csolnpLimitObjectiveFunction(Matrix myPars, int verbose);
struct Matrix csolnpEqualityFunction(int verbose);
//struct Matrix csolnpEqualityFunction(Matrix myEqBFun_arg, int verbose);

struct Matrix csolnpIneqFun(int verbose);
//struct Matrix csolnpIneqFun(Matrix myPars, int verbose);


#endif // #define _OMX_CSOLNP_SPECIFIC_H
