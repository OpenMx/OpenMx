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
#include <stdbool.h>
#include "omxState.h"
#include "omxGlobalState.h"
#include "omxNPSOLSpecific.h"
#include "omxOptimizer.h"
#include "omxMatrix.h"
#include "npsolWrap.h"
#include "omxImportFrontendState.h"
#include "subnp.h"

//#include "matrix.h"

/* NPSOL-specific globals */
const double NPSOL_BIGBND = 1e20;
const double NEG_INF = -2e20;
const double INF = 2e20;

const char* anonMatrix = "anonymous matrix";



/* NPSOL-related functions */
//************************* npsol ****************************//
//int solnp(Matrix solPars, double (*solFun)(Matrix),  Matrix solEqB,  Matrix (*solEqBFun)( Matrix),  Matrix (*solEqBStartFun)(Matrix),  Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

//extern void F77_SUB(npoptn)(char* string, int length);



Matrix fillMatrix(int cols, int rows, double* array)
{
    Matrix t = new_matrix(cols, rows);
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++) {
            printf("array is: \n");
            printf("%2f", array[j]); putchar('\n');
			M(t,j,i)=array[j];
		}
	}
	return t;
}

//****** Objective Function *********//
double csolnpObjectiveFunction(Matrix myPars)
{
    printf("myPars inside obj is: \n");
    print(myPars); putchar('\n');
    
	unsigned short int checkpointNow = FALSE;
    
	if(OMX_DEBUG) {mxLog("Starting Objective Run.\n");}
    
	omxMatrix* fitMatrix = globalState->fitMatrix;
    printf("fitMatrix is: \n");
    printf("%2f", fitMatrix->data[0]); putchar('\n');

	omxResetStatus(globalState);						// Clear Error State recursively
    printf("fitMatrix is: \n");
    printf("%2f", globalState->fitMatrix->data[0]); putchar('\n');

	/* Interruptible? */
	R_CheckUserInterrupt();
    /* This allows for abitrary repopulation of the free parameters.
     * Typically, the default is for repopulateFun to be NULL,
     * and then handleFreeVarList is invoked */
    
	if (fitMatrix->fitFunction->repopulateFun != NULL) {
		fitMatrix->fitFunction->repopulateFun(fitMatrix->fitFunction, myPars.t, myPars.cols);
        printf("fitMatrix->fitFunction is: \n");
        //printf("%2f", fitMatrix->fitFunction); putchar('\n');

	} else {
		handleFreeVarList(globalState, myPars.t, myPars.cols);
        printf("fitMatrix is: \n");
        printf("%2f", fitMatrix->data[0]); putchar('\n');
    }
    
	if (globalState->analyticGradients == 0 || globalState->currentInterval >= 0)
    {
        omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, NULL);
        printf("fitMatrix inside important if is: \n");
        printf("%2f", fitMatrix->data[0]); putchar('\n');

    }
    
	omxExamineFitOutput(globalState, fitMatrix);
    
    double *ObjectiveValue;
    double doubleValue = 0.0;
    ObjectiveValue = &doubleValue;
	*ObjectiveValue = fitMatrix->data[0];
	if(OMX_VERBOSE) {
		mxLog("Fit function value is: %.32f \n", fitMatrix->data[0]);
	}
    
	if(OMX_DEBUG) { mxLog("-======================================================-\n"); }
    
	if(checkpointNow && globalState->numCheckpoints != 0) {	// If it's a new major iteration
		omxSaveCheckpoint(globalState, myPars.t, ObjectiveValue, FALSE);		// Check about saving a checkpoint
	}
    return *ObjectiveValue;
    
}


/* Objective function for confidence interval limit finding.
 * Replaces the standard objective function when finding confidence intervals. */

//********************* npsolLimitObjectiveFunction is deleted ***********************//

/* (Non)Linear Constraint Functions */
Matrix csolnpEqualityFunction(Matrix myPars)
{
	int j, k, eq_n = 0;
    int l = 0;
    double EMPTY = -999999.0;
    Matrix myEqBFun;
    
    mxLog("Starting csolnpEqualityFunction.\n");
    printf("myPars is: \n");
    print(myPars); putchar('\n');
	handleFreeVarList(globalState, myPars.t, myPars.cols);
    printf("myPars is: \n");
    print(myPars); putchar('\n');
	for(j = 0; j < globalState->numConstraints; j++) {
		if (globalState->conList[j].opCode == 1)
        {
            eq_n += globalState->conList[j].size;
        }
    }
    
    mxLog("no.of constraints is: %d.\n", globalState->numConstraints);
    mxLog("neq is: %d.\n", eq_n);
    
    if (eq_n == 0)
    {
        myEqBFun = fill(1, 1, EMPTY);
    }
    else
    {
        myEqBFun = fill(eq_n, 1, EMPTY);
        
        for(j = 0; j < globalState->numConstraints; j++) {
            if (globalState->conList[j].opCode == 1)
            {   printf("result is: \n");
                printf("%2f", globalState->conList[j].result->data[0]); putchar('\n');
                    omxRecompute(globalState->conList[j].result);                printf("%.16f", globalState->conList[j].result->data[0]);putchar('\n');
                printf("size is: \n");
                printf("%d", globalState->conList[j].size); putchar('\n');}
                for(k = 0; k < globalState->conList[j].size; k++){
                    M(myEqBFun,l,0) = globalState->conList[j].result->data[k];
                    l = l + 1;
                }
        }
    }
    printf("myEqBFun is: \n");
    print(myEqBFun); putchar('\n');
    return myEqBFun;
}


Matrix csolnpIneqFun(Matrix myPars)
{
   	int j, k, ineq_n = 0;
    int l = 0;
    double EMPTY = -999999.0;
    Matrix myIneqFun;
    
    mxLog("Starting csolnpIneqFun.\n");
	handleFreeVarList(globalState, myPars.t, myPars.cols);
    
	for(j = 0; j < globalState->numConstraints; j++) {
		if ((globalState->conList[j].opCode == 0) || (globalState->conList[j].opCode == 2))
        {
            ineq_n += globalState->conList[j].size;
        }
    }
    
    mxLog("no.of constraints is: %d.\n", globalState->numConstraints);
    mxLog("ineq_n is: %d.\n", ineq_n);
    
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

void omxInvokeNPSOL(double *f, double *x, double *g, double *R, int disableOptimizer) {
    
    double *A=NULL, *bl=NULL, *bu=NULL, *c=NULL, *clambda=NULL, *w=NULL; //  *g, *R, *cJac,
    
    int k, ldA, ldJ, ldR, inform, iter, leniw, lenw, eq_n, ineq_n;
    
    double *cJac = NULL;    // Hessian (Approx) and Jacobian
    
    int *iw = NULL;
    
    int *istate = NULL;                 // Current state of constraints (0 = no, 1 = lower, 2 = upper, 3 = both (equality))
    
    int nctotl, nlinwid, nlnwid;    // Helpful side variables.
    
    int nclin = globalState->nclin;
    int ncnln = globalState->ncnln;
    int n = globalState->numFreeParams;
    
    double EMPTY = -999999.0;
    int j;
    
   // for(j = 0; j < globalState->numConstraints; j++) {
   //		if (globalState->conList[j].opCode == 1){ eq_n++;}
  //  }
    
    Matrix solIneqLB;
    Matrix solIneqUB;
    Matrix solEqB;

    Matrix myPars = fillMatrix(n, 1, x);
    
    double (*solFun)(struct Matrix myPars);
    solFun = &csolnpObjectiveFunction;

    Matrix (*solEqBFun)(struct Matrix myPars);
    solEqBFun = &csolnpEqualityFunction;
    
    Matrix (*solIneqFun)(struct Matrix myPars);
    solIneqFun = &csolnpIneqFun;

    
    if(n == 0) {            // Special Case for the evaluation-only condition
        
        if(OMX_DEBUG) { mxLog("No free parameters.  Avoiding Optimizer Entirely.\n"); }
        //int mode = 0, nstate = -1;
        *f = 0;
        x = NULL;
        g = NULL;
        
        if(globalState->fitMatrix != NULL) {
            *f = solFun(myPars);
            printf("1st call \n");
        };
        globalState->numIntervals = 0;  // No intervals if there's no free params
        inform = 0;
        iter = 0;
        
        for(int index = 0; index < globalState->numMats; index++) {
            omxMarkDirty(globalState->matrixList[index]);
        }
        for(int index = 0; index < globalState->numAlgs; index++) {
            omxMarkDirty(globalState->algebraList[index]);
        }
        omxStateNextEvaluation(globalState);    // Advance for a final evaluation.
        
    } else {
        
        /* Set boundaries and widths. */
        
               /* Allocate arrays */
        int i;
        bl      = (double*) R_alloc ( n, sizeof ( double ) );
        bu      = (double*) R_alloc (n, sizeof ( double ) );
        for (i = 0; i < n; i++)
        {
            printf("bl is: \n");
            printf("%2f", bl[i]); putchar('\n');
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
            int eqn;
            for(j = 0; j < globalState->numConstraints; j++) {
                if (globalState->conList[j].opCode == 1)
                {
                    eqn += globalState->conList[j].size;
                }
            }
            int nineqn = ncnln - eqn;
            solIneqLB = fill(nineqn, 1, EMPTY);
            solIneqUB = fill(nineqn, 1, EMPTY);
            solEqB = fill(eqn, 1, EMPTY);
        omxProcessConstraintsCsolnp(&solIneqLB, &solIneqUB, &solEqB);
        printf("solIneqLB is: \n");
        print(solIneqLB); putchar('\n');
        printf("solIneqUB is: \n");
        print(solIneqUB); putchar('\n');
        printf("solEqB is: \n");
        print(solEqB); putchar('\n');
        }
        omxSetupBoundsAndConstraints(bl, bu, n, nclin);
        Matrix blvar = fillMatrix(n, 1, bl);
		Matrix buvar = fillMatrix(n, 1, bu);
        
        //Matrix blvar = fill(1, 1, EMPTY);
        //Matrix buvar = fill(1, 1, EMPTY);
        
        /* Initialize Starting Values */
        if(OMX_VERBOSE) {
            mxLog("--------------------------\n");
            mxLog("Starting Values (%d) are:\n", n);
        }
        for(k = 0; k < n; k++) {
            if((M(myPars, k, 0) == 0.0) && !disableOptimizer) {
                M(myPars, k, 0) += 0.1;
                markFreeVarDependencies(globalState, k);
            }
            if(OMX_VERBOSE) { mxLog("%d: %f\n", k, M(myPars, k, 0)); }
        }
        if(OMX_DEBUG) {
            mxLog("--------------------------\n");
            mxLog("Setting up optimizer...");
        }
        
        //Matrix myPars = fillMatrix(n, 1, x);
        
        /*  F77_CALL(npsol)
         (   int *n,                 -- Number of variables
         int *nclin,             -- Number of linear constraints
         int *ncnln,             -- Number of nonlinear constraints
         int *ldA,               -- Row dimension of A (Linear Constraints)
         int *ldJ,               -- Row dimension of cJac (Jacobian)
         int *ldR,               -- Row dimension of R (Hessian)
         double *A,              -- Linear Constraints Array A (in Column-major order)
         double *bl,             -- Lower Bounds Array (at least n + nclin + ncnln long)
         double *bu,             -- Upper Bounds Array (at least n + nclin + ncnln long)
         function funcon,        -- Nonlinear constraint function
         function funobj,        -- Objective function
         int *inform,            -- Used to report state.  Need not be initialized.
         int *iter,              -- Used to report number of major iterations performed.  Need not be initialized.
         int *istate,            -- Initial State.  Need not be initialized unless using Warm Start.
         double *c,              -- Array of length ncnln.  Need not be initialized.  Reports nonlinear constraints at final iteration.
         double *cJac,           -- Array of Row-length ldJ.  Unused if ncnln = 0. Generally need not be initialized.
         double *clambda,        -- Array of length n+nclin+ncnln.  Need not be initialized unless using Warm Start. Reports final QP multipliers.
         double *f,              -- Used to report final objective value.  Need not be initialized.
         double *g,              -- Array of length n. Used to report final objective gradient.  Need not be initialized.
         double *R,              -- Array of length ldR.  Need not be intialized unless using Warm Start.
         double *x,              -- Array of length n.  Contains initial solution estimate.
         int *iw,                -- Array of length leniw. Need not be initialized.  Provides workspace.
         int *leniw,             -- Length of iw.  Must be at least 3n + nclin + ncnln.
         double *w,              -- Array of length lenw. Need not be initialized.  Provides workspace.
         int *lenw               -- Length of w.  Must be at least 2n^2 + n*nclin + 2*n*ncnln + 20*n + 11*nclin +21*ncnln
         )
         
         bl, bu, istate, and clambda are all length n+nclin+ncnln.
         First n elements refer to the vars, in order.
         Next nclin elements refer to bounds on Ax
         Last ncnln elements refer to bounds on c(x)
         
         All arrays must be in column-major order.
         
         */
        
        if(OMX_DEBUG) {
            mxLog("Set.\n");
        }
        
        if (disableOptimizer) {
            int mode = 0, nstate = -1;
            if(globalState->fitMatrix != NULL) {
                *f = solFun(myPars);
                printf("2nd call \n");
            };
            
            inform = 0;
            iter = 0;
            
            omxStateNextEvaluation(globalState);    // Advance for a final evaluation.
        } else {
            
           /* if (globalState->numConstraints == 0)
            {
                solIneqLB = fill(1,1,-999999.0);
                solIneqUB = fill(1,1,-999999.0);
            }*/

            // Matrix EqBFunValue = solEqBFun(myPars);
            // Matrix EqBStartFunValue = solEqBStartFun(myPars);
            if(OMX_DEBUG) { printf("myPars is: \n");
                print(myPars); putchar('\n');
                printf("3rd call is: \n");
                printf("%2f", solFun(myPars)); putchar('\n');
                printf("solEqB is: \n");
                print(solEqB); putchar('\n');
                printf("solEqBFun is: \n");
                print(solEqBFun(myPars)); putchar('\n');
                printf("solIneqFun is: \n");
                print(solIneqFun(myPars)); putchar('\n');
                printf("blvar is: \n");
                print(blvar); putchar('\n');
                printf("buvar is: \n");
                print(buvar); putchar('\n');
                printf("solIneqUB is: \n");
                print(solIneqUB); putchar('\n');
                printf("solIneqLB is: \n");
                print(solIneqLB); putchar('\n');
            }
        
			myPars = solnp(myPars, solFun, solEqB, solEqBFun, solIneqFun, blvar, buvar, solIneqUB, solIneqLB, myControl, myDEBUG);
        }
        
        if(OMX_DEBUG) { printf("myPars's final value is: \n");
						print(myPars);
						mxLog("Final Objective Value is: %f.\n", solFun(myPars)); 
					}
        
        omxSaveCheckpoint(globalState, myPars.t, f, TRUE);
        
        handleFreeVarList(globalState, myPars.t, myPars.cols);
        
    } // END OF PERFORM OPTIMIZATION CASE
    
    globalState->inform = inform;
    globalState->iter = iter;
    
}


/* static void
 friendlyStringToLogical(const char *key, const char *str, int *out)
 {
 int understood = FALSE;
 int newVal;
 if (matchCaseInsensitive(str, "Yes")) {
 understood = TRUE;
 newVal = 1;
 } else if (matchCaseInsensitive(str, "No")) {
 understood = TRUE;
 newVal = 0;
 } else if (isdigit(str[0]) && (atoi(str) == 1 || atoi(str) == 0)) {
 understood = TRUE;
 newVal = atoi(str);
 }
 if (!understood) {
 warning("Expecting 'Yes' or 'No' for '%s' but got '%s', ignoring", key, str);
 return;
 }
 if(OMX_DEBUG) { mxLog("%s=%d\n", key, newVal); }
 *out = newVal;
 }
 
 void omxSetNPSOLOpts(SEXP options, int *numHessians, int *calculateStdErrors,
 int *ciMaxIterations, int *disableOptimizer, int *numThreads,
 int *analyticGradients, int numFreeParams) {
 
 char optionCharArray[250] = "";			// For setting options
 int numOptions = length(options);
 SEXP optionNames;
 PROTECT(optionNames = GET_NAMES(options));
 for(int i = 0; i < numOptions; i++) {
 const char *nextOptionName = CHAR(STRING_ELT(optionNames, i));
 const char *nextOptionValue = STRING_VALUE(VECTOR_ELT(options, i));
 if(matchCaseInsensitive(nextOptionName, "Calculate Hessian")) {
 if(OMX_DEBUG) { mxLog("Found hessian option... Value: %s. ", nextOptionValue);};
 if(!matchCaseInsensitive(nextOptionValue, "No")) {
 if(OMX_DEBUG) { mxLog("Enabling explicit hessian calculation.\n");}
 if (numFreeParams > 0) {
 *numHessians = 1;
 }
 }
 } else if(matchCaseInsensitive(nextOptionName, "Standard Errors")) {
 friendlyStringToLogical(nextOptionName, nextOptionValue, calculateStdErrors);
 if (*calculateStdErrors == TRUE && numFreeParams > 0) {
 *numHessians = 1;
 }
 } else if(matchCaseInsensitive(nextOptionName, "CI Max Iterations")) {
 int newvalue = atoi(nextOptionValue);
 if (newvalue > 0) *ciMaxIterations = newvalue;
 } else if(matchCaseInsensitive(nextOptionName, "useOptimizer")) {
 if(OMX_DEBUG) { mxLog("Found useOptimizer option...");};
 if(matchCaseInsensitive(nextOptionValue, "No")) {
 if(OMX_DEBUG) { mxLog("Disabling optimization.\n");}
 *disableOptimizer = 1;
 }
 } else if(matchCaseInsensitive(nextOptionName, "Analytic Gradients")) {
 friendlyStringToLogical(nextOptionName, nextOptionValue, analyticGradients);
 } else if(matchCaseInsensitive(nextOptionName, "Number of Threads")) {
 *numThreads = atoi(nextOptionValue);
 if(OMX_DEBUG) { mxLog("Found Number of Threads option (# = %d)...\n", *numThreads);};
 } else {
 sprintf(optionCharArray, "%s %s", nextOptionName, nextOptionValue);
 //F77_CALL(npoptn)(optionCharArray, strlen(optionCharArray));
 if(OMX_DEBUG) { mxLog("Option %s \n", optionCharArray); }
 }
 }
 UNPROTECT(1); // optionNames
 }
 */
