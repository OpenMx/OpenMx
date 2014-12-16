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

static const char* anonMatrix = "anonymous matrix";
static int majIter = 400;
static int minIter = 800;
static double funcPrecision = 1.0e-7;

struct RegularFit : CSOLNPFit {
	FitContext *fc;
	omxMatrix *fitMatrix;

	RegularFit(FitContext *fc, omxMatrix *fmat) : fc(fc), fitMatrix(fmat) {
		FreeVarGroup *varGroup = fc->varGroup;
		solLB.resize(fc->numParam);
		solUB.resize(fc->numParam);
		for(int index = 0; index < int(fc->numParam); index++) {
			solLB[index] = varGroup->vars[index]->lbound;
			solUB[index] = varGroup->vars[index]->ubound;
		}

		setupIneqConstraintBounds();
	};

	void setupIneqConstraintBounds()
	{
		omxState *globalState = fc->state;
		int eqn = 0;
		for(int j = 0; j < globalState->numConstraints; j++) {
			if (globalState->conList[j].opCode == omxConstraint::EQUALITY) {
				eqn += globalState->conList[j].size;
			}
		}
		int nineqn = globalState->ncnln - eqn;
		equality.resize(eqn);
		inequality.resize(nineqn);
		if (nineqn == 0) return;

		solIneqLB.derived().resize(nineqn);
		solIneqUB.derived().resize(nineqn);
        
		int cur=0;
		for(int j = 0; j < globalState->numConstraints; j++) {
			omxConstraint &con = globalState->conList[j];
			if (con.size == 0 || con.opCode == omxConstraint::EQUALITY) continue;

			double lb, ub;
			if (con.opCode == omxConstraint::LESS_THAN) {
				lb = NEG_INF;
				ub = -0.0;
			} else {
				lb = 0.0;
				ub = INF;
			}

			for (int en=0; en < con.size; ++en) {
				solIneqLB[cur+en] = lb;
				solIneqUB[cur+en] = ub;
			}

			cur += con.size;
		}
	};

	virtual double solFun(double *myPars, int* mode, int verbose)
	{
		fc->iterations += 1;   // ought to be major iterations only
    
		if (fc->est != myPars) memcpy(fc->est, myPars, sizeof(double) * fc->numParam);
		fc->copyParamToModel();
    
		ComputeFit(fitMatrix, FF_COMPUTE_FIT, fc);
    
		if (!std::isfinite(fc->fit) || isErrorRaised()) {
			*mode = -1;
		}
    
		return fc->fit;
	};
	virtual void solEqBFun(int verbose) {
		const int eq_n = (int) equality.size();
		omxState *globalState = fc->state;
    
		if (verbose >= 3) {
			mxLog("Starting csolnpEqualityFunction %d/%d.",
			      eq_n, globalState->numConstraints);
		}
    
		if (!eq_n) return;

		int cur = 0;
		for(int j = 0; j < globalState->numConstraints; j++) {
			omxConstraint &con = globalState->conList[j];
			if (con.opCode != omxConstraint::EQUALITY) continue;

			omxRecompute(con.result, fc);
			for(int k = 0; k < globalState->conList[j].size; k++) {
				equality[cur] = con.result->data[k];
				++cur;
			}
		}
	};
	virtual void myineqFun(int verbose) {
		const int ineq_n = (int) inequality.size();
		omxState *globalState = fc->state;
    
		if (verbose >= 3) {
			mxLog("Starting csolnpInequalityFunction %d/%d.",
			      ineq_n, globalState->numConstraints);
		}
    
		if (!ineq_n) return;

		int cur = 0;
		for (int j = 0; j < globalState->numConstraints; j++) {
			omxConstraint &con = globalState->conList[j];
			if (con.opCode == omxConstraint::EQUALITY) continue;

			omxRecompute(con.result, fc);
			for(int k = 0; k < globalState->conList[j].size; k++){
				inequality[cur] = con.result->data[k];
				++cur;
			}
		}
	};
};

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc,
                     int *inform_out, FreeVarGroup *freeVarGroup,
                     int verbose, double *hessOut, double tolerance)
{
    fc->grad.resize(fc->numParam);
    
    int n = int(freeVarGroup->vars.size());
    
    Eigen::Map< Eigen::VectorXd > myPars(fc->est, n);
    
    Eigen::Array<double, 5, 1> myControl;
    myControl[0] = 1.0;
    myControl[1] = majIter;
    myControl[2] = minIter;
    myControl[3] = funcPrecision;
    myControl[4] = std::isfinite(tolerance)? tolerance : 1.0e-9;
    
    if(OMX_DEBUG) {
        mxLog("--------------------------");
        mxLog("Starting Values (%d) are:", n);
    }
    for(int k = 0; k < n; k++) {
        if((myPars[k] == 0.0)) {
            myPars[k] += 0.1;
        }
        if(OMX_DEBUG) { mxLog("%d: %.8f", k, myPars[k]); }
    }
    
    if(OMX_DEBUG) {
        mxLog("--------------------------");
    }
    
    RegularFit rf(fc, fitMatrix);
    solnp(myPars.data(), rf, myControl, verbose);
    
    fc->fit = rf.fitOut;
    if (verbose >= 1) {
        mxLog("final objective value is: \n");
        mxLog("%2f", fc->fit);
    }
    
    if (verbose>= 1){
        mxLog("final myPars value is: \n");
        for (int i = 0; i < myPars.size(); i++) mxLog("%f", myPars[i]);
    }
    
    *inform_out = rf.informOut;
    
    if (rf.gradOut.size()) {
	    fc->grad = rf.gradOut.head(n);
	    Eigen::Map< Eigen::MatrixXd > hess(hessOut, n, n);
	    hess = rf.hessOut.topLeftCorner(n, n);
    }
    
    fc->copyParamToModel();
}

struct ConfidenceIntervalFit : RegularFit {
	typedef RegularFit super;
	int currentInterval;

	ConfidenceIntervalFit(FitContext *fc, omxMatrix *fmat, int curInt) :
		super(fc, fmat), currentInterval(curInt) {};

	virtual double solFun(double *myPars, int* mode, int verbose)
	{
		//double* f = NULL;
		if (verbose >= 3) {
			mxLog("myPars inside obj is: ");
			for (int i = 0; i < int(fc->numParam); i++)
				mxLog("%f", myPars[i]);
		}
    
		fc->fit = super::solFun(myPars, mode, verbose);
    
		omxConfidenceInterval *oCI = Global->intervalList[currentInterval];
    
		omxRecompute(oCI->matrix, fc);
    
		double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);
    
		if(verbose >= 2) {
			mxLog("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.",
			      oCI->lbound, oCI->ubound, fc->fit, CIElement);
		}
    
		/* Catch boundary-passing condition */
		if(std::isnan(CIElement) || std::isinf(CIElement)) {
			fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			return fc->fit;
		}
    
		if(oCI->calcLower) {
			double diff = oCI->lbound - fc->fit;		// Offset - likelihood
			fc->fit = diff * diff + CIElement;
			// Minimize element for lower bound.
		} else {
			double diff = oCI->ubound - fc->fit;			// Offset - likelihood
			fc->fit = diff * diff - CIElement;
			// Maximize element for upper bound.
		}
    
		if(verbose >= 2) {
			mxLog("Interval fit function in previous iteration was calculated to be %f.", fc->fit);
		}
    
		return fc->fit;
	};
};

// Mostly duplicated code in omxNPSOLConfidenceIntervals
// needs to be refactored so there is only 1 copy of CI
// code that can use whatever optimizer is provided.
void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *opt, int verbose, double tolerance)
{
    int ciMaxIterations = Global->ciMaxIterations;
    // Will fail if we re-enter after an exception
    //if (NPSOL_fitMatrix) Rf_error("NPSOL is not reentrant");
    
    FitContext fc(opt, opt->varGroup);
    FreeVarGroup *freeVarGroup = fc.varGroup;
    
    int n = int(freeVarGroup->vars.size());
    double f = opt->fit;
    omxBuffer< double > gradient(n);
    omxBuffer< double > hessian(n * n);
    
    Eigen::Array<double, 5, 1> myControl;
    myControl[0] = 1.0;
    myControl[1] = majIter;
    myControl[2] = minIter;
    myControl[3] = funcPrecision;
    myControl[4] = std::isfinite(tolerance)? tolerance : 1.0e-16;

    /* Set up actual run */
    
    if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }
    
    const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
    
    for(int i = 0; i < (int) Global->intervalList.size(); i++) {
        omxConfidenceInterval *currentCI = Global->intervalList[i];
        
        const char *matName = anonMatrix;
        if (currentCI->matrix->name) {
            matName = currentCI->matrix->name;
        }
        
        Global->checkpointMessage(opt, opt->est, "%s[%d, %d] begin lower interval",
                                  matName, currentCI->row + 1, currentCI->col + 1);
        
        memcpy(fc.est, opt->est, n * sizeof(double)); // Reset to previous optimum
	Eigen::Map< Eigen::VectorXd > myPars(fc.est, n);
        
        currentCI->lbound += opt->fit;          // Convert from offsets to targets
        currentCI->ubound += opt->fit;          // Convert from offsets to targets
        
        if (std::isfinite(currentCI->lbound))
        {
		int inform = -1;
            // Number of times to keep trying.
            int cycles = ciMaxIterations;
            
            double value = INF;
            
            while(inform!= 0 && cycles > 0) {
                /* Find lower limit */
                currentCI->calcLower = TRUE;
		ConfidenceIntervalFit cif(&fc, fitMatrix, i);
                solnp(fc.est, cif, myControl, verbose);
                
                f = cif.fitOut;
                inform = cif.informOut;
                
                if(verbose>=1) { mxLog("inform_lower is: %d", inform);}
                
                currentCI->lCode = inform;
                
                if(f < value) {
                    currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    
                    value = f;
                }
                
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of lower interval %d failed: Bad inform value of %d",
                          i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opt->est[j]) > objDiff) {
                            jitter = FALSE;
                            if(OMX_DEBUG) {mxLog("are u here?????");}
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opt->est[j] + objDiff;
                        }
                    }
                }
            }
        }
        
        if (std::isfinite(currentCI->ubound)) {
            currentCI->calcLower = FALSE;
            Global->checkpointMessage(opt, opt->est, "%s[%d, %d] begin upper interval",
                                      matName, currentCI->row + 1, currentCI->col + 1);
            
            
            memcpy(fc.est, opt->est, n * sizeof(double)); // Reset to previous optimum
            
            /* Reset for the upper bound */
            double value = INF;
            int inform = -1;
            int cycles = ciMaxIterations;
            
            while(inform != 0 && cycles >= 0) {
                /* Find upper limit */
                currentCI->calcLower = FALSE;
		ConfidenceIntervalFit cif(&fc, fitMatrix, i);
                solnp(myPars.data(), cif, myControl, verbose);
                
                f = cif.fitOut;
                inform = cif.informOut;

                if(verbose >= 1) { mxLog("inform_upper is: %d", inform);}
                currentCI->uCode = inform;
                
                if(f < value) {
                    currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
                }
                
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of upper interval %d failed: Bad inform value of %d",
                          i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opt->est[j]) > objDiff){
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opt->est[j] + objDiff;
                        }
                    }
                }
            }
            
            if(OMX_DEBUG) {mxLog("Found Upper bound %d.", i);}
            
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
