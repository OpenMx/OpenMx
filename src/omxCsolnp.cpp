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


/* NPSOL-related functions */
//************************* npsol ****************************//
//int solnp(Matrix solPars, double (*solFun)(Matrix),  Matrix solEqB,  Matrix (*solEqBFun)( Matrix),  Matrix (*solEqBStartFun)(Matrix),  Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

//extern void F77_SUB(npoptn)(char* string, int length);

static omxMatrix *GLOB_fitMatrix = NULL;
static FitContext *GLOB_fc = NULL;

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
	if (verbose) {
		printf("myPars inside obj is: ");
		print(myPars); putchar('\n');
	}
    
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
    
    
    if (isErrorRaised(globalState) || !R_FINITE(fitMatrix->data[0])) {
        GLOB_fc->fit = 1e24;
    } else {
        GLOB_fc->fit = fitMatrix->data[0];
    }
    
	if(verbose) {
		mxLog("Fit function value is: %.32f\n", fitMatrix->data[0]);
	}
    
	if(checkpointNow && globalState->numCheckpoints != 0) {	// If it's a new major iteration
		omxSaveCheckpoint(myPars.t, GLOB_fc->fit, FALSE);		// Check about saving a checkpoint
	}
	return GLOB_fc->fit;
}


/* Objective function for confidence interval limit finding.
 * Replaces the standard objective function when finding confidence intervals. */

/* (Non)Linear Constraint Functions */
Matrix csolnpEqualityFunction(Matrix myPars, int verbose)
{
	int j, k, eq_n = 0;
    int l = 0;
    double EMPTY = -999999.0;
    Matrix myEqBFun;
    
    if (verbose) mxLog("Starting csolnpEqualityFunction.");
        
    GLOB_fc->copyParamToModel(globalState, myPars.t);
    
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

    GLOB_fc->copyParamToModel(globalState, myPars.t);
    
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
	GLOB_fitMatrix = fitMatrix;
	GLOB_fc = fc;
    
    double *x = fc->est;
	double *g = fc->grad;
    
    double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *w=NULL; //  *g, *R, *cJac,
    
    int k, ldA, ldJ, ldR, leniw, lenw, eq_n, ineq_n;
    
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
    int ncnln = globalState->ncnln;
    int n = freeVarGroup->vars.size();
    
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
    
    solFun_t solFun = &csolnpObjectiveFunction;
    
    solEqBFun_t solEqBFun = &csolnpEqualityFunction;
    
    solIneqFun_t solIneqFun = &csolnpIneqFun;
    
    /* Set boundaries and widths. */
    
    /* Allocate arrays */
    bl      = Calloc ( n, double );
    bu      = Calloc (n, double );
    
	if (verbose >= 2) {
		for (int i = 0; i < n; i++) {
			printf("bl is: ");
			printf("%2f", bl[i]); putchar('\n');
		}
	}
    
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
        int eqn = 0;
        for(j = 0; j < globalState->numConstraints; j++) {
            if (globalState->conList[j].opCode == 1)
            {
                eqn += globalState->conList[j].size;
            }
        }
        int nineqn = ncnln - eqn;
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
    
    //Matrix blvar = fill(1, 1, EMPTY);
    //Matrix buvar = fill(1, 1, EMPTY);
    
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
    
    /* if (globalState->numConstraints == 0)
     {
     solIneqLB = fill(1,1,-999999.0);
     solIneqUB = fill(1,1,-999999.0);
     }*/
    
    // Matrix EqBFunValue = solEqBFun(myPars);
    // Matrix EqBStartFunValue = solEqBStartFun(myPars);
    if(verbose == 2) { printf("myPars is: ");
        print(myPars); putchar('\n');
        printf("3rd call is: ");
        printf("%2f", solFun(myPars, 0)); putchar('\n');
        printf("solEqB is: ");
        print(solEqB); putchar('\n');
        printf("solEqBFun is: ");
        print(solEqBFun(myPars, 0)); putchar('\n');
        printf("solIneqFun is: ");
        print(solIneqFun(myPars, 0)); putchar('\n');
        printf("blvar is: ");
        print(blvar); putchar('\n');
        printf("buvar is: ");
        print(buvar); putchar('\n');
        printf("solIneqUB is: ");
        print(solIneqUB); putchar('\n');
        printf("solIneqLB is: ");
        print(solIneqLB); putchar('\n');
    }
    
    p_obj = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG, verbose);
    
    
    fc->fit = *p_obj.objValue;
    if (verbose >= 1) {
	    printf("final objective value is: \n");
	    printf("%2f", fc->fit); putchar('\n');
    }
    param_hess = *p_obj.parameter;
    
    int i;
    myPars = subset(param_hess, 0, 0, n-1);
    if (verbose>= 1){    
	    printf("final myPars value is: \n");
	    print(myPars); putchar('\n');
    }
    myhess = subset(param_hess, 0, n, param_hess.cols - myPars.cols - 1);
    
    if (verbose >= 2){
	    printf("myhess is: \n");
	    print(myhess); putchar('\n');
    }

    double Hess[myhess.cols];
    
    for (i = 0; i < myhess.cols; i++)
    {
        Hess[i] = myhess.t[i];
    }

    mygrad = subset(param_hess, 0, myPars.cols + (myPars.cols*myPars.cols), param_hess.cols-1);
    
    memcpy(fc->hess, Hess, sizeof(double) * n * n);
        
    for (i = 0; i < myPars.cols; i++){
        x[i] = myPars.t[i];
    }
    
    omxSaveCheckpoint(x, fc->fit, TRUE);
    
    fc->copyParamToModel(globalState);
    GLOB_fitMatrix = NULL;
    GLOB_fc = NULL;
}
