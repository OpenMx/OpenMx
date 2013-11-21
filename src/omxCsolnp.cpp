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
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
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
	    
	unsigned short int checkpointNow = FALSE;
    
	if(OMX_DEBUG) {printf("Starting Objective Run.");}
    
	omxMatrix* fitMatrix = GLOB_fitMatrix;
    
	omxResetStatus(globalState);						// Clear Error State recursively
    
    /* Interruptible? */
	R_CheckUserInterrupt();
    /* This allows for abitrary repopulation of the free parameters.
     * Typically, the default is for repopulateFun to be NULL,
     * and then handleFreeVarList is invoked */
    
	GLOB_fc->copyParamToModel(globalState, myPars.t);
    
    
    omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, GLOB_fc);

    
    if (!R_FINITE(fitMatrix->data[0])) {
        omxRaiseErrorf(globalState, "Fit function returned %g at iteration %d.%d",
                       fitMatrix->data[0], globalState->majorIteration, globalState->minorIteration);
    }
    if(isErrorRaised(globalState))
    {
        if(OMX_DEBUG)
        {    mxLog("Error status reported.");}
    }
    if (isErrorRaised(globalState) || !R_FINITE(fitMatrix->data[0]))
    {    GLOB_fc->fit = 1e24;}
    else {  GLOB_fc->fit = fitMatrix->data[0];}
    
    
	if(verbose >= 1) {
		mxLog("Fit function value is: %.32f\n", fitMatrix->data[0]);
	}
    
	if(checkpointNow && globalState->numCheckpoints != 0) {	// If it's a new major iteration
		omxSaveCheckpoint(myPars.t, GLOB_fc->fit, FALSE);		// Check about saving a checkpoint
	}
	return GLOB_fc->fit;
}


/* Objective function for confidence interval limit finding.
 * Replaces the standard objective function when finding confidence intervals. */
double csolnpLimitObjectiveFunction(Matrix myPars, int verbose)
{
    //double* f = NULL;
	if (verbose >= 3) {
		printf("myPars inside obj is: ");
		print(myPars); putchar('\n');
	}
    
    if(OMX_VERBOSE) mxLog("Calculating interval %d, %s boundary:", CSOLNP_currentInterval, (Global->intervalList[CSOLNP_currentInterval].calcLower?"lower":"upper"));
    
    GLOB_fc->fit = csolnpObjectiveFunction(myPars, verbose);
    
    omxConfidenceInterval *oCI = &(Global->intervalList[CSOLNP_currentInterval]);
    
    omxRecompute(oCI->matrix);
    
    double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);
    
    if(OMX_DEBUG) {
        mxLog("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.",
              oCI->lbound, oCI->ubound, GLOB_fc->fit, CIElement);
    }
    
    /* Catch boundary-passing condition */
    if(isnan(CIElement) || isinf(CIElement)) {
        omxRaiseError(globalState, -1,
                      "Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
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
    
    if(OMX_DEBUG) {
        mxLog("Interval fit function in previous iteration was calculated to be %f.", GLOB_fc->fit);
    }
    
    return GLOB_fc->fit;
}


/* (Non)Linear Constraint Functions */
Matrix csolnpEqualityFunction(Matrix myPars, int verbose)
{
	int j, k, eq_n = 0;
    int l = 0;
    double EMPTY = -999999.0;
    Matrix myEqBFun;
    
    if (verbose) mxLog("Starting csolnpEqualityFunction.");
    
    for(j = 0; j < globalState->numConstraints; j++) {
	    if (globalState->conList[j].opCode == 1) {
		    eq_n += globalState->conList[j].size;
	    }
    }
    
    if (verbose >= 1) {
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
			    if (verbose >= 1) {
				    mxLog("result is: %2f", globalState->conList[j].result->data[0]);
			    }
			    omxRecompute(globalState->conList[j].result);
			    if (verbose >= 1) {
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
    if (verbose >= 1) {
	    printf("myEqBFun is: ");
	    print(myEqBFun);
	    putchar('\n');
    }
    return myEqBFun;
}


Matrix csolnpIneqFun(Matrix myPars, int verbose)
{
   	int j, k, ineq_n = 0;
    int l = 0;
    double EMPTY = -999999.0;
    Matrix myIneqFun;
    
    if (verbose) mxLog("Starting csolnpIneqFun.");
        
	for(j = 0; j < globalState->numConstraints; j++) {
		if ((globalState->conList[j].opCode == 0) || (globalState->conList[j].opCode == 2))
        {
            ineq_n += globalState->conList[j].size;
        }
    }
    
	if (verbose >= 1) {
		printf("no.of constraints is: %d.", globalState->numConstraints); putchar('\n');
		printf("ineq_n is: %d.", ineq_n); putchar('\n');
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
                     int *inform_out, int *iter_out, bool useGradient, FreeVarGroup *freeVarGroup,
                     int verbose)

{
	freeMatrices(); // maybe left overs from an aborted optimization attempt
    
	GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
    double *x = fc->est;
    double *g = fc->grad;
    
    
    int k, eq_n, ineq_n, iter = -1;
    int inform = 0;
    
    double *bl=NULL, *bu=NULL;
    
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
    int ncnln = globalState->ncnln;
    int n = int(freeVarGroup->vars.size());
    
    double EMPTY = -999999.0;
    int j;
    
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
    Matrix (*solEqBFun)(struct Matrix myPars, int verbose);
    solEqBFun = &csolnpEqualityFunction;
    Matrix (*solIneqFun)(struct Matrix myPars, int verbose);
    solIneqFun = &csolnpIneqFun;
    
    /* Set boundaries and widths. */
    
    /* Allocate arrays */
    bl      = Calloc ( n + ncnln, double );
    bu      = Calloc ( n + ncnln, double );
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = 400.0;
    M(myControl,2,0) = 800.0;
    M(myControl,3,0) = 1.0e-7;
    M(myControl,4,0) = 1.0e-8;
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
        
    /* Initialize Starting Values */
    if(OMX_VERBOSE) {
        printf("--------------------------");
        printf("Starting Values (%d) are:", n);
    }
    for(k = 0; k < n; k++) {
        if((M(myPars, k, 0) == 0.0)) {
            M(myPars, k, 0) += 0.1;
        }
        if(OMX_VERBOSE) { printf("%d: %f", k, M(myPars, k, 0)); }
    }
    if(OMX_DEBUG) {
        printf("--------------------------");
        printf("Setting up optimizer...");
    }
        
    
    p_obj = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
    
    
    fc->fit = p_obj.objValue;
    if (verbose >= 1) {
	    printf("final objective value is: \n");
	    printf("%2f", fc->fit); putchar('\n');
    }
    param_hess = p_obj.parameter;
    
    int i;
    myPars = subset(param_hess, 0, 0, n-1);
    
    if (verbose>= 1){
	    printf("final myPars value is: \n");
	    print(myPars); putchar('\n');
    }
    myhess = subset(param_hess, 0, n, param_hess.cols - myPars.cols - 2);
    
    
    Matrix inform_m = subset(param_hess, 0, param_hess.cols-1, param_hess.cols-1);
    
    inform = M(inform_m, 0, 0);
    
    if (verbose >= 2){
	    printf("myhess is: \n");
	    print(myhess); putchar('\n');
    }
    
    double Hess[myhess.cols];
    
    for (i = 0; i < myhess.cols; i++)
    {
        Hess[i] = myhess.t[i];
      
    }
    
    mygrad = subset(param_hess, 0, myPars.cols + (myPars.cols*myPars.cols), param_hess.cols-2);
    
    
    memcpy(fc->hess, Hess, sizeof(double) * n * n);

    for (i = 0; i < myPars.cols; i++){
        x[i] = myPars.t[i];
    }
    
    omxSaveCheckpoint(x, fc->fit, TRUE);
    
    GLOB_fc->copyParamToModel(globalState);
    
    *inform_out = inform;
    *iter_out   = iter;
    
    GLOB_fitMatrix = NULL;
    GLOB_fc = NULL;
    freeMatrices();
}



void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose)
{
	int ciMaxIterations = Global->ciMaxIterations;
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) error("NPSOL is not reentrant");
    
    GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
	FreeVarGroup *freeVarGroup = fitMatrix->fitFunction->freeVarGroup;
    
    int eq_n, ineq_n;
    double inform;
    double *bl=NULL, *bu=NULL;
    
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
    int j;
    
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
    Matrix (*solEqBFun)(struct Matrix myPars, int verbose);
    solEqBFun = &csolnpEqualityFunction;
    Matrix (*solIneqFun)(struct Matrix myPars, int verbose);
    solIneqFun = &csolnpIneqFun;
    int nclin = 0;
    
    
    /* Set boundaries and widths. */
    /* Allocate arrays */
    bl      = Calloc ( n + ncnln, double );
    bu      = Calloc (n + ncnln, double );
    
    
    struct Matrix myControl = fill(6,1,(double)0.0);
    M(myControl,0,0) = 1.0;
    M(myControl,1,0) = 400.0;
    M(myControl,2,0) = 800.0;
    M(myControl,3,0) = 1.0e-7;
    M(myControl,4,0) = 1.0e-16;
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
    
    int count = 0;
    for(int i = 0; i < Global->numIntervals; i++) {
        if (i == 0)
        {count++;}
        omxConfidenceInterval *currentCI = &(Global->intervalList[i]);
        
        int msgLength = 45;
        
        if (currentCI->matrix->name == NULL) {
            msgLength += strlen(anonMatrix);
        } else {
            msgLength += strlen(currentCI->matrix->name);
        }
        
        char *message = Calloc(msgLength, char);
        
        if (currentCI->matrix->name == NULL) {
            snprintf(message, msgLength, "%s[%d, %d] begin lower interval",
                     anonMatrix, currentCI->row + 1, currentCI->col + 1);
        } else {
            snprintf(message, msgLength, "%s[%d, %d] begin lower interval",
                     currentCI->matrix->name, currentCI->row + 1, currentCI->col + 1);
        }
        
        omxWriteCheckpointMessage(message);
        
        memcpy(x.data(), optimalValues, n * sizeof(double)); // Reset to previous optimum
        myPars = fillMatrix(n, 1, x.data());
        CSOLNP_currentInterval = i;
        
        
        currentCI->lbound += optimum;          // Convert from offsets to targets
        currentCI->ubound += optimum;          // Convert from offsets to targets
        
        /* Set up for the lower bound */
        inform = -1;
        // Number of times to keep trying.
        
        
        int cycles = ciMaxIterations;
        
        double value = INF;
        double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
        
        while(inform!= 0 && cycles > 0) {
            /* Find lower limit */
            currentCI->calcLower = TRUE;
            
            p_obj_conf = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
            
            f = p_obj_conf.objValue;
            
            myPars = subset(p_obj_conf.parameter, 0, 0, n-1);
            myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.cols - 2);
            
            mygrad = subset(p_obj_conf.parameter, 0, myPars.cols + (myPars.cols*myPars.cols), p_obj_conf.parameter.cols-2);
            
            
            Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
            
            inform = M(inform_m, 1, 0);
            
            if(verbose>=1) { mxLog("inform_lower is: %f", inform);}
            currentCI->lCode = inform;
            
            
            if(f < value) {
                currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                
                value = f;
				for (int ii = 0; ii < myPars.cols; ii++){
        			x.data()[ii] = myPars.t[ii];
    			}
                omxSaveCheckpoint(x.data(), f, TRUE);
            }
            
            if(inform != 0 && OMX_DEBUG) {
                mxLog("Calculation of lower interval %d failed: Bad inform value of %f",
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
        
        if(OMX_DEBUG) { mxLog("Found lower bound %d.  Seeking upper.", i); }
        // TODO: Repopulate original optimizer state in between CI calculations
        
        if (currentCI->matrix->name == NULL) {
            snprintf(message, msgLength, "%s[%d, %d] begin upper interval",
                     anonMatrix, currentCI->row + 1, currentCI->col + 1);
        } else {
            snprintf(message, msgLength, "%s[%d, %d] begin upper interval",
                     currentCI->matrix->name, currentCI->row + 1, currentCI->col + 1);
        }
        
        omxWriteCheckpointMessage(message);
        
        Free(message);
        
        memcpy(x.data(), optimalValues, n * sizeof(double));
        myPars = fillMatrix(n, 1, x.data());
        
        
        /* Reset for the upper bound */
        value = INF;
        inform = -1;
        cycles = ciMaxIterations;
        if(verbose >= 1) { mxLog("cycles_upper is: %d", cycles); }
        while(inform != 0 && cycles >= 0) {
            /* Find upper limit */
            currentCI->calcLower = FALSE;
            p_obj_conf = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
            
            f = p_obj_conf.objValue;
            
            myPars = subset(p_obj_conf.parameter, 0, 0, n-1);
            myhess = subset(p_obj_conf.parameter, 0, n, p_obj_conf.parameter.cols - myPars.cols - 2);
            
            mygrad = subset(p_obj_conf.parameter, 0, myPars.cols + (myPars.cols*myPars.cols), p_obj_conf.parameter.cols-2);
            
            Matrix inform_m = subset(p_obj_conf.parameter, 0, p_obj_conf.parameter.cols-1, p_obj_conf.parameter.cols-1);
            
            
            inform = M(inform_m, 1, 0);
            if(verbose >= 1) { mxLog("inform_upper is: %f", inform);}
            currentCI->uCode = inform;
            
            if(f < value) {
                currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                
                value = f;
				for (int ii = 0; ii < myPars.cols; ii++){
        			x.data()[ii] = myPars.t[ii];
                }
                
                omxSaveCheckpoint(x.data(), f, TRUE);
            }
            
            if(inform != 0 && OMX_DEBUG) {
                mxLog("Calculation of upper interval %d failed: Bad inform value of %f",
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
    
	GLOB_fc = NULL;
	GLOB_fitMatrix = NULL;
	CSOLNP_currentInterval = -1;
}
