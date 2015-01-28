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

// Mostly duplicated code in omxNPSOLConfidenceIntervals
// needs to be refactored so there is only 1 copy of CI
// code that can use whatever optimizer is provided.
void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *opt, int verbose, double tolerance)
{
    const int ciMaxIterations = Global->ciMaxIterations;
    FitContext fc(opt, opt->varGroup);
    fc.createChildren();
    FreeVarGroup *freeVarGroup = fc.varGroup;
    
    const int n = int(freeVarGroup->vars.size());
    
    if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }
    
    const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
    
    for(int i = 0; i < (int) Global->intervalList.size(); i++) {
        omxConfidenceInterval *currentCI = Global->intervalList[i];
        
        const char *matName = "anonymous matrix";
        if (currentCI->matrix->name) {
            matName = currentCI->matrix->name;
        }
        
        currentCI->lbound += opt->fit;          // Convert from offsets to targets
        currentCI->ubound += opt->fit;          // Convert from offsets to targets
        
	for (int lower=0; lower <= 1; ++lower) {
		if (lower  && !std::isfinite(currentCI->lbound)) continue;
		if (!lower && !std::isfinite(currentCI->ubound)) continue;

		memcpy(fc.est, opt->est, n * sizeof(double)); // Reset to previous optimum
        
		int tries = 0;
		int inform = -1;
		double bestFit = std::numeric_limits<double>::max();
            
		while (inform!= 0 && ++tries <= ciMaxIterations) {
			Global->checkpointMessage(opt, opt->est, "%s[%d, %d] %s CI (try %d)",
						  matName, currentCI->row + 1, currentCI->col + 1,
						  lower? "lower" : "upper", tries);

			ConfidenceIntervalFit cif(0);
			cif.optName = "CSOLNP";
			cif.fc = &fc;
			cif.fitMatrix = fitMatrix;
			cif.currentInterval = i;
			cif.calcLower = lower;
			cif.ControlMajorLimit = majIter;
			cif.ControlMinorLimit = minIter;
			cif.ControlFuncPrecision = funcPrecision;
			cif.ControlTolerance = std::isfinite(tolerance)? tolerance : 1.0e-16;
			solnp(fc.est, cif);
                
			fc.copyParamToModel();
			ComputeFit("ComputeGD", fitMatrix, FF_COMPUTE_FIT, &fc);
			const double fitOut = fc.fit;

			if (fitOut < bestFit) {
				double val = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
				if (lower) currentCI->min = val;
				else       currentCI->max = val;
				bestFit = fitOut;
			}

			inform = cif.informOut;
			if (lower) currentCI->lCode = inform;
			else       currentCI->uCode = inform;
			if(verbose>=1) { mxLog("inform(%d,%d) is: %d", i, lower, inform);}
			if(inform == 0) break;

			bool jitter = TRUE;
			for(int j = 0; j < n; j++) {
				if(fabs(fc.est[j] - opt->est[j]) > objDiff) {
					jitter = FALSE;
					break;
				}
			}
			if(jitter) {
				for(int j = 0; j < n; j++) {
					double sign = 2 * (tries % 2) - 1;
					fc.est[j] = opt->est[j] + sign * objDiff * tries;
				}
			}
		}
        }
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
