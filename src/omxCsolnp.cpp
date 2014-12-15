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
#include "matrix.h"
#include "omxCsolnp.h"
#include "omxBuffer.h"

static const char* anonMatrix = "anonymous matrix";
static int majIter = 400;
static int minIter = 800;
static double funcPrecision = 1.0e-7;

/* NPSOL-related functions */
//************************* npsol ****************************//
//int solnp(Matrix solPars, double (*solFun)(Matrix),  Matrix solEqB,  Matrix (*solEqBFun)( Matrix),  Matrix (*solEqBStartFun)(Matrix),  Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

static omxMatrix *GLOB_fitMatrix = NULL;
static FitContext *GLOB_fc = NULL;
static int CSOLNP_currentInterval = -1;

template <typename T1>
static void setupIneqConstraintBounds(FitContext *fc, Eigen::MatrixBase<T1> &solIneqLB, Eigen::MatrixBase<T1> &solIneqUB)
{
	omxState *globalState = fc->state;
        int eqn = 0;
        for(int j = 0; j < globalState->numConstraints; j++) {
		if (globalState->conList[j].opCode == omxConstraint::EQUALITY) {
			eqn += globalState->conList[j].size;
		}
        }
        int nineqn = globalState->ncnln - eqn;
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
}

struct RegularFit : CSOLNPFit {
	virtual double solFun(Matrix myPars, int* mode, int verbose) {
		omxMatrix* fitMatrix = GLOB_fitMatrix;
    
		GLOB_fc->iterations += 1;   // ought to be major iterations only
    
		memcpy(GLOB_fc->est, myPars.t, sizeof(double) * myPars.cols);
		GLOB_fc->copyParamToModel();
    
		ComputeFit(fitMatrix, FF_COMPUTE_FIT, GLOB_fc);
    
		if (!std::isfinite(GLOB_fc->fit) || isErrorRaised()) {
			*mode = -1;
		}
    
		return GLOB_fc->fit;
	};
	virtual Matrix solEqBFun(int verbose) {
		int i, j, k, eq_n = 0;
		int l = 0;
		double EMPTY = -999999.0;
		Matrix myEqBFun;
		omxState *globalState = GLOB_fc->state;
    
		if (verbose >= 3) mxLog("Starting csolnpEqualityFunction.");
    
		for(j = 0; j < globalState->numConstraints; j++) {
			if (globalState->conList[j].opCode == 1) {
				eq_n += globalState->conList[j].size;
			}
		}
    
		if (verbose >= 3) {
			mxLog("no.of constraints is: %d.", globalState->numConstraints);
			mxLog("neq is: %d.", eq_n);
		}
    
		if (eq_n) {
			myEqBFun = fill(eq_n, 1, EMPTY);
        
			for(j = 0; j < globalState->numConstraints; j++) {
				if (globalState->conList[j].opCode == 1) {
					if (verbose >= 3) {
						mxLog("result is: %2f", globalState->conList[j].result->data[0]);
					}
					omxRecompute(globalState->conList[j].result, GLOB_fc);
					if (verbose >= 3) {
						mxLog("%.16f", globalState->conList[j].result->data[0]);
						mxLog("size is: %d", globalState->conList[j].size);
					}
				}
				for(k = 0; k < globalState->conList[j].size; k++){
					M(myEqBFun,l,0) = globalState->conList[j].result->data[k];
					l = l + 1;
				}
			}
		}
		if (verbose >= 3) {
			mxLog("myEqBFun is: ");
			for(i = 0; i < myEqBFun.cols; i++)
				{   mxLog("%f", myEqBFun.t[i]);}
		}
		return myEqBFun;
	};
	virtual Matrix myineqFun(int verbose) {
		int j, k, ineq_n = 0;
		int l = 0;
		double EMPTY = -999999.0;
		Matrix myIneqFun;
		omxState *globalState = GLOB_fc->state;
    
		if (verbose >= 3) mxLog("Starting csolnpIneqFun.");
    
		for(j = 0; j < globalState->numConstraints; j++) {
			if ((globalState->conList[j].opCode == 0) || (globalState->conList[j].opCode == 2))
				{
					ineq_n += globalState->conList[j].size;
				}
		}
    
		if (verbose >= 3) {
			mxLog("no.of constraints is: %d.", globalState->numConstraints); putchar('\n');
			mxLog("ineq_n is: %d.", ineq_n); putchar('\n');
		}
    
		if (ineq_n == 0)
			{
				myIneqFun = fill(1, 1, EMPTY);
			}
		else
			{
				myIneqFun = fill(ineq_n, 1, EMPTY);
        
				for(j = 0; j < globalState->numConstraints; j++) {
					if ((globalState->conList[j].opCode == 0) || globalState->conList[j].opCode == 2) {
						omxRecompute(globalState->conList[j].result, GLOB_fc);
					}
					for(k = 0; k < globalState->conList[j].size; k++){
						M(myIneqFun,l,0) = globalState->conList[j].result->data[k];
						l = l + 1;
                
					}
				}
			}
    
		return myIneqFun;
	};
};

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc,
                     int *inform_out, FreeVarGroup *freeVarGroup,
                     int verbose, double *hessOut, double tolerance)
{
    freeMatrices(); // maybe left overs from an aborted optimization attempt
    
    GLOB_fitMatrix = fitMatrix;
    GLOB_fc = fc;
    
    fc->grad.resize(fc->numParam);
    
    int k;
    int inform = 0;
    
    //double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
    omxState *globalState = fc->state;
    const int ncnln = globalState->ncnln;
    int n = int(freeVarGroup->vars.size());
    
    Param_Obj p_obj;
    Matrix param_hess;
    
    Eigen::Map< Eigen::VectorXd > myPars(fc->est, n);
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = majIter;
    M(myControl,2,0) = minIter;
    M(myControl,3,0) = funcPrecision;
    M(myControl,4,0) = std::isfinite(tolerance)? tolerance : 1.0e-9;
    M(myControl,5,0) = 0.0;
    
    bool myDEBUG = false;
    /* Set up actual run */
    
    Eigen::VectorXd EsolIneqLB;
    Eigen::VectorXd EsolIneqUB;
    setupIneqConstraintBounds(fc, EsolIneqLB, EsolIneqUB);
    Matrix solIneqLB(EsolIneqLB);
    Matrix solIneqUB(EsolIneqUB);

    Eigen::VectorXd bl(n);
    Eigen::VectorXd bu(n);
    for(int index = 0; index < n; index++) {
        bl[index] = freeVarGroup->vars[index]->lbound;
        bu[index] = freeVarGroup->vars[index]->ubound;
    }
    Matrix blvar(bl);
    Matrix buvar(bu);
    
    /* Initialize Starting Values */
    if(OMX_DEBUG) {
        mxLog("--------------------------");
        mxLog("Starting Values (%d) are:", n);
    }
    for(k = 0; k < n; k++) {
        if((myPars[k] == 0.0)) {
            myPars[k] += 0.1;
        }
        if(OMX_DEBUG) { mxLog("%d: %.8f", k, myPars[k]); }
    }
    
    if(OMX_DEBUG) {
        mxLog("--------------------------");
    }
    
    
    RegularFit rf;
    CSOLNP solnpContext(rf);
    p_obj = solnpContext.solnp(myPars, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
    
    
    fc->fit = p_obj.objValue;
    if (verbose >= 1) {
        mxLog("final objective value is: \n");
        mxLog("%2f", fc->fit);
    }
    param_hess = p_obj.parameter;
    
    subset(param_hess, 0, 0, n-1, myPars);
    
    if (verbose>= 1){
        mxLog("final myPars value is: \n");
        for (int i = 0; i < myPars.size(); i++) mxLog("%f", myPars[i]);
    }
    
    Matrix inform_m = subset(param_hess, 0, param_hess.cols-1, param_hess.cols-1);
    
    inform = M(inform_m, 0, 0);
    
    // Mahsa, remove this conditional
    if (ncnln == 0) {
	    subset(param_hess, 0, myPars.size() + (myPars.size()*myPars.size()), param_hess.cols-2, fc->grad);
    
	    Eigen::Map< Eigen::VectorXd > hessVec(hessOut, myPars.size() * myPars.size());
	    subset(param_hess, 0, n, param_hess.cols - myPars.size() - 2, hessVec);
    }
    
    GLOB_fc->copyParamToModel();
    
    *inform_out = inform;
    
    GLOB_fitMatrix = NULL;
    GLOB_fc = NULL;
    freeMatrices();
}

struct ConfidenceIntervalFit : RegularFit {
	typedef RegularFit super;
	virtual double solFun(Matrix myPars, int* mode, int verbose) {
		//double* f = NULL;
		if (verbose >= 3) {
			mxLog("myPars inside obj is: ");
			for (int i = 0; i < myPars.cols; i++)
				mxLog("%f", myPars.t[i]);
		}
    
		GLOB_fc->fit = super::solFun(myPars, mode, verbose);
    
		omxConfidenceInterval *oCI = Global->intervalList[CSOLNP_currentInterval];
    
		omxRecompute(oCI->matrix, GLOB_fc);
    
		double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);
    
		if(verbose >= 2) {
			mxLog("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.",
			      oCI->lbound, oCI->ubound, GLOB_fc->fit, CIElement);
		}
    
		/* Catch boundary-passing condition */
		if(std::isnan(CIElement) || std::isinf(CIElement)) {
			GLOB_fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
			return GLOB_fc->fit;
		}
    
		if(oCI->calcLower) {
			double diff = oCI->lbound - GLOB_fc->fit;		// Offset - likelihood
			GLOB_fc->fit = diff * diff + CIElement;
			// Minimize element for lower bound.
		} else {
			double diff = oCI->ubound - GLOB_fc->fit;			// Offset - likelihood
			GLOB_fc->fit = diff * diff - CIElement;
			// Maximize element for upper bound.
		}
    
		if(verbose >= 2) {
			mxLog("Interval fit function in previous iteration was calculated to be %f.", GLOB_fc->fit);
		}
    
		return GLOB_fc->fit;
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
    
    GLOB_fitMatrix = fitMatrix;
    GLOB_fc = &fc;
    
    FreeVarGroup *freeVarGroup = fc.varGroup;
    
    int inform;
    
    int n = int(freeVarGroup->vars.size());
    double f = opt->fit;
    omxBuffer< double > gradient(n);
    omxBuffer< double > hessian(n * n);
    
    Param_Obj p_obj_conf;
    Matrix param_hess;
    Matrix myhess = fill(n*n, 1, (double)0.0);
    Matrix mygrad;
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = majIter;
    M(myControl,2,0) = minIter;
    M(myControl,3,0) = funcPrecision;
    M(myControl,4,0) = std::isfinite(tolerance)? tolerance : 1.0e-16;
    M(myControl,5,0) = 0.0;
    
    bool myDEBUG = false;
    /* Set up actual run */
    
    Eigen::VectorXd EsolIneqLB;
    Eigen::VectorXd EsolIneqUB;
    setupIneqConstraintBounds(&fc, EsolIneqLB, EsolIneqUB);
    Matrix solIneqLB(EsolIneqLB);
    Matrix solIneqUB(EsolIneqUB);
    
    Eigen::VectorXd bl(n);
    Eigen::VectorXd bu(n);
    for(int index = 0; index < n; index++) {
        bl[index] = freeVarGroup->vars[index]->lbound;
        bu[index] = freeVarGroup->vars[index]->ubound;
    }
    Matrix blvar(bl);
    Matrix buvar(bu);
    
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
        CSOLNP_currentInterval = i;
        
        
        currentCI->lbound += opt->fit;          // Convert from offsets to targets
        currentCI->ubound += opt->fit;          // Convert from offsets to targets
        
        if (std::isfinite(currentCI->lbound))
        {
            /* Set up for the lower bound */
            inform = -1;
            // Number of times to keep trying.
            int cycles = ciMaxIterations;
            
            double value = INF;
            
            while(inform!= 0 && cycles > 0) {
                /* Find lower limit */
                currentCI->calcLower = TRUE;
		ConfidenceIntervalFit cif;
                CSOLNP solnpContext1(cif);
                p_obj_conf = solnpContext1.solnp(myPars, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
                
                f = p_obj_conf.objValue;
                
                subset(p_obj_conf.parameter, 0, 0, n-1, myPars);
                myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.size() - 2);
                
                mygrad = subset(p_obj_conf.parameter, 0, myPars.size() + (myPars.size()*myPars.size()), p_obj_conf.parameter.cols-2);
                
                
                Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
                
                inform = M(inform_m, 0, 0);
                
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
            inform = -1;
            int cycles = ciMaxIterations;
            
            while(inform != 0 && cycles >= 0) {
                /* Find upper limit */
                currentCI->calcLower = FALSE;
		ConfidenceIntervalFit cif;
                CSOLNP solnpContext1(cif);
                p_obj_conf = solnpContext1.solnp(myPars, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
                
                f = p_obj_conf.objValue;
                
		subset(p_obj_conf.parameter, 0, 0, n-1, myPars);
                myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.size() - 2);
                
                mygrad = subset(p_obj_conf.parameter, 0, myPars.size() + (myPars.size()*myPars.size()), p_obj_conf.parameter.cols-2);
                
                Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
                
                
                inform = M(inform_m, 0, 0);
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
    
    GLOB_fc = NULL;
    GLOB_fitMatrix = NULL;
    CSOLNP_currentInterval = -1;
    freeMatrices();
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
