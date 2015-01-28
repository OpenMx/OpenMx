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

#include <ctype.h>
#include <limits>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "glue.h"
#include "omxImportFrontendState.h"
#include "omxCsolnp.h"
#include "omxBuffer.h"

// move to Global? TODO
static int majIter = 400;
static int minIter = 800;
static double funcPrecision = 1.0e-7;

void omxCSOLNP(double *est, GradientOptimizerContext &go)
{
	go.optName = "CSOLNP";
	go.ControlMajorLimit = majIter;
	go.ControlMinorLimit = minIter;
	go.ControlFuncPrecision = funcPrecision;
	solnp(est, go);
}

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc,
                     int *inform_out, FreeVarGroup *freeVarGroup,
                     int verbose, double *hessOut, double tolerance)
{
	GradientOptimizerContext rf(verbose);
	rf.fc = fc;
	rf.fitMatrix = fitMatrix;
	rf.ControlTolerance = std::isfinite(tolerance)? tolerance : 1.0e-9;
	omxCSOLNP(fc->est, rf);
    
    *inform_out = rf.informOut;
    
    if (rf.gradOut.size()) {
	    int n = int(freeVarGroup->vars.size());
	    fc->grad = rf.gradOut.tail(n);
	    Eigen::Map< Eigen::MatrixXd > hess(hessOut, n, n);
	    hess = rf.hessOut.bottomRightCorner(n, n);
    }
}

void CSOLNPOpt_majIter(const char *optionValue)
{
    majIter = atoi(optionValue);
}

void CSOLNPOpt_minIter(const char *optionValue)
{
    minIter = atoi(optionValue);
}

void CSOLNPOpt_FuncPrecision(const char *optionValue)
{
    funcPrecision = atof(optionValue);
}
