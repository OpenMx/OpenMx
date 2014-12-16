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

struct CSOLNPFit {
	Eigen::VectorXd solLB;
	Eigen::VectorXd solUB;

	Eigen::VectorXd solIneqLB;
	Eigen::VectorXd solIneqUB;

	Eigen::VectorXd equality;
	Eigen::VectorXd inequality;

	virtual double solFun(double *myPars, int* mode, int verbose) = 0;
	virtual void solEqBFun(int verbose) = 0;
	virtual void myineqFun(int verbose) = 0;
};

struct CSOLNP {

    int flag, flag_L, flag_U, index_flag_L, index_flag_U, flag_NormgZ, flag_step, minr_rec;
    Matrix LB;
    Matrix UB;
    Matrix resP;
    double resLambda;
    Matrix resHessv;
    Matrix resY;
    Matrix sx_Matrix;
    int mode_val;
    int* mode;
	CSOLNPFit &fit;

	CSOLNP(CSOLNPFit &_fit) : fit(_fit) {};
    
	enum Control {
		ControlRho=0,
		ControlMajorLimit,
		ControlMinorLimit,
		ControlFuncPrecision,
		ControlTolerance,
		NumControl
	};

	Param_Obj solnp(double *pars, const Eigen::Array<double, NumControl, 1> &solctrl, int verbose);

    Matrix subnp(Matrix pars, Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,
		 const Eigen::Array<double, 4, 1> &ctrl, int verbose);

	enum indParam {
		indNumParam=0,
		indHasGradient,
		indHasHessian,
		indHasIneq,
		indHasJacobianIneq,
		indHasEq,
		indHasJacobianEq,
		indVectorLength  // must be last
	};

	Eigen::Array<double, int(indVectorLength), 1> ind;
};

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc, int *inform_out,
                     FreeVarGroup *freeVarGroup, int verbose, double *hessOut,
                     double tolerance);

void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose,
                                  double tolerance);

void CSOLNPOpt_majIter(const char *optionValue);

void CSOLNPOpt_minIter(const char *optionValue);

void CSOLNPOpt_FuncPrecision(const char *optionValue);

#endif // #define _OMX_CSOLNP_SPECIFIC_H
