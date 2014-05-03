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

static const char* anonMatrix = "anonymous matrix";

/* NPSOL-related functions */
//************************* npsol ****************************//
//int solnp(Matrix solPars, double (*solFun)(Matrix),  Matrix solEqB,  Matrix (*solEqBFun)( Matrix),  Matrix (*solEqBStartFun)(Matrix),  Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

static omxMatrix *GLOB_fitMatrix = NULL;
static FitContext *GLOB_fc = NULL;
static int CSOLNP_currentInterval = -1;

Matrix fillMatrix(int cols, int rows, double* array)
{
    Matrix t = new_matrix(cols, rows);
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++) {
			M(t,j,i)=array[j];
		}
	}
	return t;
}


//****** Objective Function *********//
double csolnpObjectiveFunction(Matrix myPars, int verbose)
{
	if(OMX_DEBUG) {mxLog("Starting Objective Run.");}
    
	omxMatrix* fitMatrix = GLOB_fitMatrix;
    
	R_CheckUserInterrupt();

	GLOB_fc->iterations += 1;   // ought to be major iterations only

	memcpy(GLOB_fc->est, myPars.t, sizeof(double) * myPars.cols);
	GLOB_fc->copyParamToModel(globalState);

	ComputeFit(fitMatrix, FF_COMPUTE_FIT, GLOB_fc);

	if (!std::isfinite(fitMatrix->data[0])) {
		GLOB_fc->fit = 1e24;
	}
    
	if(verbose >= 1) {
		mxLog("Fit function value is: %.32f", fitMatrix->data[0]);
	}

	return GLOB_fc->fit;
}


/* Objective function for confidence interval limit finding.
 * Replaces the standard objective function when finding confidence intervals. */
double csolnpLimitObjectiveFunction(Matrix myPars, int verbose)
{
    //double* f = NULL;
	if (verbose >= 3) {
		mxLog("myPars inside obj is: ");
        for (int i = 0; i < myPars.cols; i++)
            mxLog("%f", myPars.t[i]);
	}
    
    GLOB_fc->fit = csolnpObjectiveFunction(myPars, verbose);
    
    omxConfidenceInterval *oCI = &(Global->intervalList[CSOLNP_currentInterval]);
    
    omxRecompute(oCI->matrix);
    
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
			    omxRecompute(globalState->conList[j].result);
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
            if ((globalState->conList[j].opCode == 0) || globalState->conList[j].opCode == 2)
            {   omxRecompute(globalState->conList[j].result);}
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
	freeMatrices(); // maybe left overs from an aborted optimization attempt
    
	GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
    double *x = fc->est;
    fc->grad.resize(fc->numParam);
    double *g = fc->grad.data();
    
    
    int k;
    int inform = 0;
    
    //double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
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
    
    double (*solFun)(struct Matrix myPars, int verbose);
    solFun = &csolnpObjectiveFunction;
    Matrix (*solEqBFun)(int verbose);
    solEqBFun = &csolnpEqualityFunction;
    Matrix (*solIneqFun)(int verbose);
    solIneqFun = &csolnpIneqFun;
    
    /* Set boundaries and widths. */
    
    std::vector<double> blBuf(n+ncnln);
    std::vector<double> buBuf(n+ncnln);
    double *bl = blBuf.data();
    double *bu = buBuf.data();
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = 400.0;
    M(myControl,2,0) = 800.0;
    M(myControl,3,0) = 1.0e-7;
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
        
        omxProcessConstraintsCsolnp(&solIneqLB, &solIneqUB, &solEqB);

        if (verbose == 2) {
            mxLog("solIneqLB is: ");
            for (int i = 0; i < solIneqLB.cols; i++){mxLog("%f", solIneqLB.t[i]);}
            mxLog("solIneqUB is: ");
            for (int i = 0; i < solIneqUB.cols; i++){mxLog("%f", solIneqUB.t[i]);}
            mxLog("solEqB is: ");
            for (int i = 0; i < solEqB.cols; i++){mxLog("%f", solEqB.t[i]);}
        }
        
    }
    omxSetupBoundsAndConstraints(freeVarGroup, bl, bu);
    
    Matrix blvar = fillMatrix(n, 1, bl);
    Matrix buvar = fillMatrix(n, 1, bu);
        
    /* Initialize Starting Values */
    if(OMX_VERBOSE) {
        mxLog("--------------------------");
        mxLog("Starting Values (%d) are:", n);
    }
    for(k = 0; k < n; k++) {
        if((M(myPars, k, 0) == 0.0)) {
            M(myPars, k, 0) += 0.1;
        }
        if(OMX_VERBOSE) { mxLog("%d: %.8f", k, M(myPars, k, 0)); }
    }

    if(OMX_DEBUG) {
        mxLog("--------------------------");
        mxLog("Setting up optimizer...");
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
    
    GLOB_fc->copyParamToModel(globalState);
    
    *inform_out = inform;
    
    GLOB_fitMatrix = NULL;
    GLOB_fc = NULL;
    freeMatrices();
}



// Mostly duplicated code in omxNPSOLConfidenceIntervals
// needs to be refactored so there is only 1 copy of CI
// code that can use whatever optimizer is provided.
void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose, double tolerance)
{
	int ciMaxIterations = Global->ciMaxIterations;
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) Rf_error("NPSOL is not reentrant");
    
    GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
	FreeVarGroup *freeVarGroup = fc->varGroup;
    
    int inform;
    
    int n = int(freeVarGroup->vars.size());
    int ncnln = globalState->ncnln;
    
    double optimum = fc->fit;
    
    double *optimalValues = fc->est;
    
    double f = optimum;
    std::vector< double > x(n, *optimalValues);
    std::vector< double > gradient(n);
    std::vector< double > hessian(n * n);
    
    /* CSOLNP Arguments */
    double EMPTY = -999999.0;
    
    Param_Obj p_obj_conf;
    Matrix param_hess;
    Matrix myhess = fill(n*n, 1, (double)0.0);
    Matrix mygrad;
    Matrix solIneqLB;
    Matrix solIneqUB;
    Matrix solEqB;
    
    Matrix myPars = fillMatrix(n, 1, fc->est);
    double (*solFun)(struct Matrix myPars, int verbose);
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
    
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = 400.0;
    M(myControl,2,0) = 800.0;
    M(myControl,3,0) = 1.0e-7;
    M(myControl,4,0) = std::isfinite(tolerance)? tolerance : 1.0e-16;
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
        
        omxProcessConstraintsCsolnp(&solIneqLB, &solIneqUB, &solEqB);
        if (verbose == 2) {
            printf("solIneqLB is: ");
            print(solIneqLB); putchar('\n');
            printf("solIneqUB is: ");
            print(solIneqUB); putchar('\n');
            printf("solEqB is: ");
            print(solEqB); putchar('\n');
        }
    }
    
    omxSetupBoundsAndConstraints(freeVarGroup, bl, bu);
    Matrix blvar = fillMatrix(n, 1, bl);
    Matrix buvar = fillMatrix(n, 1, bu);
    
    if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }
    
    
    if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals.");
        mxLog("numIntervals is: %d", Global->numIntervals);
    }
    
    const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?

    for(int i = 0; i < Global->numIntervals; i++) {
        omxConfidenceInterval *currentCI = &(Global->intervalList[i]);
        
	const char *matName = anonMatrix;
	if (currentCI->matrix->name) {
		matName = currentCI->matrix->name;
	}

	Global->checkpointMessage(fc, fc->est, "%s[%d, %d] begin lower interval",
				  matName, currentCI->row + 1, currentCI->col + 1);
        
        memcpy(x.data(), optimalValues, n * sizeof(double)); // Reset to previous optimum
        myPars = fillMatrix(n, 1, x.data());
        CSOLNP_currentInterval = i;
        
        
        currentCI->lbound += optimum;          // Convert from offsets to targets
        currentCI->ubound += optimum;          // Convert from offsets to targets
        
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
                        for (int ii = 0; ii < myPars.cols; ii++){
                            x.data()[ii] = myPars.t[ii];
                        }
                    }
            
                    if(inform != 0 && OMX_DEBUG) {
                        mxLog("Calculation of lower interval %d failed: Bad inform value of %d",
                              i, inform);
                    }
                    cycles--;
                    if(inform != 0) {
                        unsigned int jitter = TRUE;
                        for(int j = 0; j < n; j++) {
                            if(fabs(x[j] - optimalValues[j]) > objDiff) {
                                jitter = FALSE;
                                if(OMX_DEBUG) {mxLog("are u here?????");}
                                break;
                            }
                        }
                        if(jitter) {
                            for(int j = 0; j < n; j++) {
                                x[j] = optimalValues[j] + objDiff;
                            }
                        }
                    }
                }
            }
        
	if (std::isfinite(currentCI->ubound)) {
            currentCI->calcLower = FALSE;
	Global->checkpointMessage(fc, fc->est, "%s[%d, %d] begin upper interval",
				  matName, currentCI->row + 1, currentCI->col + 1);

        
        memcpy(x.data(), optimalValues, n * sizeof(double));
        myPars = fillMatrix(n, 1, x.data());
        
        
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
				for (int ii = 0; ii < myPars.cols; ii++){
        			x.data()[ii] = myPars.t[ii];
                }
            }
            
            if(inform != 0 && OMX_DEBUG) {
                mxLog("Calculation of upper interval %d failed: Bad inform value of %d",
                      i, inform);
            }
            cycles--;
            if(inform != 0) {
                unsigned int jitter = TRUE;
                for(int j = 0; j < n; j++) {
                    if(fabs(x[j] - optimalValues[j]) > objDiff){
                        jitter = FALSE;
                        break;
                    }
                }
                if(jitter) {
                    for(int j = 0; j < n; j++) {
                        x[j] = optimalValues[j] + objDiff;
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
