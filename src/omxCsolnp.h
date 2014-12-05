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
#include <Eigen/Core>

struct CSOLNP {
    
    int flag, flag_L, flag_U, index_flag_L, index_flag_U, flag_NormgZ, flag_step, minr_rec;
    Matrix LB;
    Matrix UB;
    Matrix ineqLB;
    Matrix ineqUB;
    Matrix ind;
    Matrix resP;
    double resLambda;
    Matrix resHessv;
    Matrix resY;
    Matrix sx_Matrix;
    int mode_val;
    int* mode;
    //    CSOLNP() {}
    //    CSOLNP(double _EMPTY, int _flag, int _flag_L, int _flag_U, int _index_flag_L, int _index_flag_U, int _flag_NormgZ, int _flag_step, int _minr_rec, )
    
    typedef double (*solFun_t)(Matrix, int*, int);
    typedef Matrix (*solEqBFun_t)(int);
    typedef Matrix (*myineqFun_t)(int);
    
    Param_Obj solnp(Matrix solPars, solFun_t solFun, Matrix solEqB, solEqBFun_t solEqBFun, myineqFun_t myineqFun, Matrix solLB, Matrix solUB, Matrix solIneqUB, Matrix solIneqLB, Matrix solctrl, bool debugToggle, int verbose);
    template <typename T1>
    Matrix subnp(Matrix pars, double (*solFun)(Matrix, int*, int), Matrix (*solEqBFun)(int), Matrix (*myineqFun)(int), Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale, Eigen::ArrayBase<T1> &ctrl, int verbose);
};

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc, int *inform_out,
                     FreeVarGroup *freeVarGroup, int verbose, double *hessOut,
                     double tolerance);

void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose,
                                  double tolerance);

double csolnpObjectiveFunction(Matrix myPars, int* mode, int verbose);
double csolnpLimitObjectiveFunction(Matrix myPars, int* mode, int verbose);
struct Matrix csolnpEqualityFunction(int verbose);
//struct Matrix csolnpEqualityFunction(Matrix myEqBFun_arg, int verbose);

struct Matrix csolnpIneqFun(int verbose);
//struct Matrix csolnpIneqFun(Matrix myPars, int verbose);

void CSOLNPOpt_majIter(const char *optionValue);

void CSOLNPOpt_minIter(const char *optionValue);

void CSOLNPOpt_FuncPrecision(const char *optionValue);

#endif // #define _OMX_CSOLNP_SPECIFIC_H
