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


//****** Objective Function *********//
double csolnpObjectiveFunction(Matrix myPars, int* mode, int verbose)
{
	omxMatrix* fitMatrix = GLOB_fitMatrix;
    
	GLOB_fc->iterations += 1;   // ought to be major iterations only

	memcpy(GLOB_fc->est, myPars.t, sizeof(double) * myPars.cols);
	GLOB_fc->copyParamToModel();

	ComputeFit(fitMatrix, FF_COMPUTE_FIT, GLOB_fc);

	if (!std::isfinite(GLOB_fc->fit) || isErrorRaised()) {
        *mode = -1;
	}
    
	return GLOB_fc->fit;
}


/* Objective function for confidence interval limit finding.
 * Replaces the standard objective function when finding confidence intervals. */
double csolnpLimitObjectiveFunction(Matrix myPars, int* mode, int verbose)
{
    //double* f = NULL;
	if (verbose >= 3) {
		mxLog("myPars inside obj is: ");
        for (int i = 0; i < myPars.cols; i++)
            mxLog("%f", myPars.t[i]);
	}
    
    GLOB_fc->fit = csolnpObjectiveFunction(myPars, mode, verbose);
    
    omxConfidenceInterval *oCI = Global->intervalList[CSOLNP_currentInterval];
    
    omxRecompute(oCI->matrix, FF_COMPUTE_FIT, GLOB_fc);
    
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
}


/* (Non)Linear Constraint Functions */
Matrix csolnpEqualityFunction(int verbose)
{
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
    
    if (eq_n == 0)
    {
        myEqBFun = fill(1, 1, EMPTY);
    }
    else {
	    myEqBFun = fill(eq_n, 1, EMPTY);
	    
	    for(j = 0; j < globalState->numConstraints; j++) {
		    if (globalState->conList[j].opCode == 1) {
			    if (verbose >= 3) {
				    mxLog("result is: %2f", globalState->conList[j].result->data[0]);
			    }
			    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, GLOB_fc);
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
}


Matrix csolnpIneqFun(int verbose)
{
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
		    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, GLOB_fc);
	    }
            for(k = 0; k < globalState->conList[j].size; k++){
                M(myIneqFun,l,0) = globalState->conList[j].result->data[k];
                l = l + 1;
                
            }
        }
    }
    
    return myIneqFun;
}

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc,
                     int *inform_out, FreeVarGroup *freeVarGroup,
                     int verbose, double *hessOut, double tolerance)

{
	freeMatrices();
    freeMatrices_l(); // maybe left overs from an aborted optimization attempt
    
	GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
    double *x = fc->est;
    fc->grad.resize(fc->numParam);
    double *g = fc->grad.data();
    
    
    int k;
    int inform = 0;
    
    //double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
    omxState *globalState = fc->state;
    int ncnln = globalState->ncnln;
    int n = int(freeVarGroup->vars.size());
    
    double EMPTY = -999999.0;
    
    Param_Obj p_obj;
    Matrix param_hess;
    Matrix myhess = fill(n*n, 1, (double)0.0);
    Matrix mygrad;
    Matrix solIneqLB;
    Matrix solIneqUB;
    Matrix solEqB;
    
    Matrix myPars = fillMatrix(n, 1, fc->est);
    
    double (*solFun)(struct Matrix myPars, int* mode, int verbose);
    solFun = &csolnpObjectiveFunction;
    Matrix (*solEqBFun)(int verbose);
    solEqBFun = &csolnpEqualityFunction;
    Matrix (*solIneqFun)(int verbose);
    solIneqFun = &csolnpIneqFun;
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = majIter;
    M(myControl,2,0) = minIter;
    M(myControl,3,0) = funcPrecision;
    M(myControl,4,0) = std::isfinite(tolerance)? tolerance : 1.0e-9;
    M(myControl,5,0) = 0.0;
    
    bool myDEBUG = false;
    /* Set up actual run */
    
    /* needs treatment*/
    if (ncnln == 0)
    {
        solIneqLB = fill(1, 1, EMPTY);
        solIneqUB = fill(1, 1, EMPTY);
        solEqB = fill(1, 1, EMPTY);
    }
    else{
        int j;
        int nineqn;
        int eqn = 0;
        for(j = 0; j < globalState->numConstraints; j++) {
            if (globalState->conList[j].opCode == 1)
            {
                eqn += globalState->conList[j].size;
            }
        }
        if (eqn == ncnln) nineqn = 1;
        else nineqn = ncnln - eqn;
        
        solIneqLB = fill(nineqn, 1, EMPTY);
        solIneqUB = fill(nineqn, 1, EMPTY);
	    if (eqn == 0) {
		    solEqB = fill(1, 1, EMPTY);
	    } else {
		    solEqB = fill(eqn, 1, EMPTY);
	    }
        
	    omxProcessConstraintsCsolnp(fc, &solIneqLB, &solIneqUB, &solEqB);

        if (verbose == 2) {
            mxLog("solIneqLB is: ");
            for (int i = 0; i < solIneqLB.cols; i++){mxLog("%f", solIneqLB.t[i]);}
            mxLog("solIneqUB is: ");
            for (int i = 0; i < solIneqUB.cols; i++){mxLog("%f", solIneqUB.t[i]);}
            mxLog("solEqB is: ");
            for (int i = 0; i < solEqB.cols; i++){mxLog("%f", solEqB.t[i]);}
        }
        
    }

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
        if((M(myPars, k, 0) == 0.0)) {
            M(myPars, k, 0) += 0.1;
        }
        if(OMX_DEBUG) { mxLog("%d: %.8f", k, M(myPars, k, 0)); }
    }

    if(OMX_DEBUG) {
        mxLog("--------------------------");
    }
        
    
    p_obj = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
    
    
    fc->fit = p_obj.objValue;
    if (verbose >= 1) {
	    mxLog("final objective value is: \n");
	    mxLog("%2f", fc->fit);
    }
    param_hess = p_obj.parameter;
    
    int i;
    myPars = subset(param_hess, 0, 0, n-1);
    
    if (verbose>= 1){
	    mxLog("final myPars value is: \n");
        for (i = 0; i < myPars.cols; i++) mxLog("%f", myPars.t[i]);
    }
    myhess = subset(param_hess, 0, n, param_hess.cols - myPars.cols - 2);
    
    Matrix inform_m = subset(param_hess, 0, param_hess.cols-1, param_hess.cols-1);
    
    inform = M(inform_m, 0, 0);
    
    if (verbose >= 2){
	    mxLog("myhess is: \n");
        for (i = 0; i < myhess.cols; i++)
	    mxLog("%f", myhess.t[i]);
    }
    
    mygrad = subset(param_hess, 0, myPars.cols + (myPars.cols*myPars.cols), param_hess.cols-2);
    
    for (i = 0; i < myPars.cols * myPars.cols; i++){
        hessOut[i] = myhess.t[i];
    }
    
    for (i = 0; i < myPars.cols; i++){
        x[i] = myPars.t[i];
        g[i] = mygrad.t[i];
    }
    
    GLOB_fc->copyParamToModel();
    
    *inform_out = inform;
    
    GLOB_fitMatrix = NULL;
    GLOB_fc = NULL;
    freeMatrices_l();
}



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
    omxState *globalState = opt->state;
    int ncnln = globalState->ncnln;
    
    double f = opt->fit;
    omxBuffer< double > gradient(n);
    omxBuffer< double > hessian(n * n);
    
    /* CSOLNP Arguments */
    double EMPTY = -999999.0;
    
    Param_Obj p_obj_conf;
    Matrix param_hess;
    Matrix myhess = fill(n*n, 1, (double)0.0, FALSE);
    Matrix mygrad;
    Matrix solIneqLB;
    Matrix solIneqUB;
    Matrix solEqB;
    
    Matrix myPars = fillMatrix(n, 1, opt->est, FALSE);
    double (*solFun)(struct Matrix myPars, int* mode, int verbose);
    solFun = &csolnpLimitObjectiveFunction;
    Matrix (*solEqBFun)(int verbose);
    solEqBFun = &csolnpEqualityFunction;
    Matrix (*solIneqFun)(int verbose);
    solIneqFun = &csolnpIneqFun;
    
    
    /* Set boundaries and widths. */
    /* Allocate arrays */
    std::vector<double> blBuf(n+ncnln);
    std::vector<double> buBuf(n+ncnln);
    double *bl = blBuf.data();
    double *bu = buBuf.data();
    
    
    struct Matrix myControl = fill(6,1,(double)0.0, FALSE);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = majIter;
    M(myControl,2,0) = minIter;
    M(myControl,3,0) = funcPrecision;
    M(myControl,4,0) = std::isfinite(tolerance)? tolerance : 1.0e-16;
    M(myControl,5,0) = 0.0;
    
    bool myDEBUG = false;
    /* Set up actual run */
    
    /* needs treatment*/
    if (ncnln == 0)
    {
        solIneqLB = fill(1, 1, EMPTY, FALSE);
        solIneqUB = fill(1, 1, EMPTY, FALSE);
        solEqB = fill(1, 1, EMPTY, FALSE);
    }
    else{
        int j;
        int nineqn;
        int eqn = 0;
        for(j = 0; j < globalState->numConstraints; j++) {
            if (globalState->conList[j].opCode == 1)
            {
                eqn += globalState->conList[j].size;
            }
        }
        if (eqn == ncnln) nineqn = 1;
        else nineqn = ncnln - eqn;
        
        solIneqLB = fill(nineqn, 1, EMPTY, FALSE);
        solIneqUB = fill(nineqn, 1, EMPTY, FALSE);
	    if (eqn == 0) {
		    solEqB = fill(1, 1, EMPTY, FALSE);
	    } else {
		    solEqB = fill(eqn, 1, EMPTY, FALSE);
	    }
        
	    omxProcessConstraintsCsolnp(opt, &solIneqLB, &solIneqUB, &solEqB);
        if (verbose == 2) {
            printf("solIneqLB is: ");
            print(solIneqLB); putchar('\n');
            printf("solIneqUB is: ");
            print(solIneqUB); putchar('\n');
            printf("solEqB is: ");
            print(solEqB); putchar('\n');
        }
    }
    
    omxSetupBoundsAndConstraints(opt, bl, bu);
    Matrix blvar = fillMatrix(n, 1, bl, FALSE);
    Matrix buvar = fillMatrix(n, 1, bu, FALSE);
    
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
        myPars = fillMatrix(n, 1, opt->est, FALSE);
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
                    p_obj_conf = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
            
                    f = p_obj_conf.objValue;
                        
                    myPars = subset(p_obj_conf.parameter, 0, 0, n-1);
                    myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.cols - 2);
            
                    mygrad = subset(p_obj_conf.parameter, 0, myPars.cols + (myPars.cols*myPars.cols), p_obj_conf.parameter.cols-2);
            
            
                    Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
            
                    inform = M(inform_m, 0, 0);
            
                    if(verbose>=1) { mxLog("inform_lower is: %d", inform);}
                    
                    currentCI->lCode = inform;
            
                    if(f < value) {
                        currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                
                        value = f;
			memcpy(fc.est, myPars.t, sizeof(double) * myPars.cols);
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
        myPars = fillMatrix(n, 1, opt->est, FALSE);
        
        
        /* Reset for the upper bound */
        double value = INF;
        inform = -1;
        int cycles = ciMaxIterations;

        while(inform != 0 && cycles >= 0) {
            /* Find upper limit */
            currentCI->calcLower = FALSE;
            p_obj_conf = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
            
            f = p_obj_conf.objValue;
            
            myPars = subset(p_obj_conf.parameter, 0, 0, n-1);
            myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.cols - 2);
            
            mygrad = subset(p_obj_conf.parameter, 0, myPars.cols + (myPars.cols*myPars.cols), p_obj_conf.parameter.cols-2);
            
            Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
            
            
            inform = M(inform_m, 0, 0);
            if(verbose >= 1) { mxLog("inform_upper is: %d", inform);}
            currentCI->uCode = inform;
            
            if(f < value) {
                currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);

                value = f;
		memcpy(fc.est, myPars.t, sizeof(double) * myPars.cols);
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
    freeMatrices_l();
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
