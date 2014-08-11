//
//  nloptc.c
//  
//
//  Created by Mahsa Zahery on 7/22/14.
//
//

#include <stdio.h>
#include <cstdlib>
#include <ctype.h>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxMatrix.h"
#include "glue.h"
#include "omxBuffer.h"
#include "nloptcpp.h"
#include <R_ext/Applic.h>
#include "omxImportFrontendState.h"

#ifdef HAS_NLOPT
#include <nlopt.h>

static const char* anonMatrix = "anonymous matrix";
static omxMatrix * nlopt_sann_fitMatrix = NULL;
static FitContext * nlopt_sann_fc = NULL;
static int nlopt_sann_currentInterval = -1;

typedef struct ex_struct
{
    SEXP R_fcall;    /* function */
    SEXP R_gcall;    /* gradient */
    SEXP R_env;      /* where to evaluate the calls */
    double* ndeps;   /* tolerances for numerical derivatives */
    double fnscale;  /* scaling for objective */
    double* parscale;/* scaling for parameters */
    int usebounds;
    double* lower, *upper;
    SEXP names;	     /* names for par */
} ex_struct, *exStruct;

static double fn(int n, double *par, void *ex)
{
    omxMatrix* fitMatrix = nlopt_sann_fitMatrix;
    
	nlopt_sann_fc->iterations += 1;   // ought to be major iterations only
    
	memcpy(nlopt_sann_fc->est, par, sizeof(double) * n);
	nlopt_sann_fc->copyParamToModel();
    
	ComputeFit(fitMatrix, FF_COMPUTE_FIT, nlopt_sann_fc);
    
	if (!std::isfinite(nlopt_sann_fc->fit) || isErrorRaised()) {
        nlopt_sann_fc->fit = 1e24;
	}
    return nlopt_sann_fc->fit;
}

double nloptObjectiveFunction(unsigned n, const double* x, double* grad, void* f_data)
{
	omxMatrix* fitMatrix = nlopt_sann_fitMatrix;
    
	nlopt_sann_fc->iterations += 1;   // ought to be major iterations only
    
	memcpy(nlopt_sann_fc->est, x, sizeof(double) * n);
	nlopt_sann_fc->copyParamToModel();
    
	ComputeFit(fitMatrix, FF_COMPUTE_FIT, nlopt_sann_fc);
    
	if (!std::isfinite(nlopt_sann_fc->fit) || isErrorRaised()) {
        nlopt_sann_fc->fit = 1e24;
	}
    
    return nlopt_sann_fc->fit;
}

void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
    int j, k, l = 0;
    
    omxState *globalState = nlopt_sann_fc->state;

    for(j = 0; j < globalState->numConstraints; j++) {
        if ((globalState->conList[j].opCode == 0) || globalState->conList[j].opCode == 2) {
		    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc);
	    }
        for(k = 0; k < globalState->conList[j].size; k++){
            result[l] = globalState->conList[j].result->data[k];
            l = l + 1;
        }
    }
}


void nloptEqualityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)

{
    omxState *globalState = nlopt_sann_fc->state;
    int j, k, l = 0;
    
    for(j = 0; j < globalState->numConstraints; j++) {
        if (globalState->conList[j].opCode == 1) {
            omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc);
        }
        for(k = 0; k < globalState->conList[j].size; k++){
            result[l] = globalState->conList[j].result->data[k];
            l = l + 1;
        }
    }
}

void omxInvokeNLOPTorSANN(omxMatrix *fitMatrix, FitContext *fc,
                    int *inform_out, FreeVarGroup *freeVarGroup,
                    int verbose, double *hessOut, double tolerance)
{
    char alg[] = "NLOPT algorithms"; //For now, if another algorithm needs to be used, this variable should be changed.
    /* possible algorithm options are:
     - MMA
     - SLSQP
     - ISRES
     - COBYLA
     - DIRECT
     - BOBYQA
     - csr2lm
     - lbfgs
     - mlsl
     - neldermead
     - newuoa
    */
	nlopt_sann_fitMatrix = fitMatrix;
    nlopt_sann_fc = fc;
	double *x = fc->est;
	fc->grad.resize(fc->numParam);
	double *g = fc->grad.data();
    int n = int(freeVarGroup->vars.size());
    omxState *globalState = fc->state;
    int ncnln = globalState->ncnln;
    int nineq = 0;
    int neq = 0;

    //Define samin parameters
    double Fmin;
    exStruct ex;
    ex_struct zz;
    ex = (exStruct)malloc(sizeof(zz));
    ex->R_gcall = R_NilValue;
    
    if (strcmp(alg, "simulated annealing") == 0)
    {
        if (ncnln == 0)
        {
            samin(n, x, &Fmin, fn, 10000, 10, 10, 1, ex);
            /* arg 5 is maxit: It gives the total number of function evaluations. Defaults to 10000.
               arg 6 is tmax: The number of function evaluations at each temperature. Defaults to 10.
               arg7 is temp: The starting temperature for the cooling schedule. Defaults to 10.
               arg 8 is trace: Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information.
             */
            nlopt_sann_fc->fit = Fmin;
            nlopt_sann_fc->copyParamToModel();
            nlopt_sann_fitMatrix = NULL;
            nlopt_sann_fc = NULL;
            return;
        }
        else
            Rf_error("change the algorithm, SANN cannot be used for constrained problems.\n");
    }

    else{
    /* For constrained problems NLOPT should be used */    
    /* Set up lower and upper bounds */
    std::vector<double> blBuf(n);
    std::vector<double> buBuf(n);
    double *bl = blBuf.data();
    double *bu = buBuf.data();
    omxSetupBoundsAndConstraints(fc, bl, bu);

    /* Initialize Starting Values */
    if(OMX_DEBUG) {
        mxLog("--------------------------");
        mxLog("Starting Values (%d) are:", n);
    }
    for(int k = 0; k < n; k++) {
        if((x[k] == 0.0)) {
            x[k] += 0.1;
        }
        if(OMX_DEBUG) { mxLog("%d: %f", k, x[k]); }
    }
    if(OMX_DEBUG) {
        mxLog("--------------------------");
    }
    
    /* Set up nlopt variables */
    nlopt_opt opt, local_opt;
    double opt_f;
    opt = nlopt_create(NLOPT_AUGLAG, n); // Main algorithm
    local_opt = nlopt_create(NLOPT_LN_COBYLA, n); // Subsidiary algorithm
    
    /* Subsidiary algorithm could be one of the followings:
     
     DIRECT and DIRECT-L
     - NLOPT_GN_DIRECT
     - NLOPT_GN_DIRECT_L
     - NLOPT_GLOBAL_DIRECT_L_RAND
     - NLOPT_GLOBAL_DIRECT_NOSCAL
     - NLOPT_GLOBAL_DIRECT_L_NOSCAL
     - NLOPT_GLOBAL_DIRECT_L_RAND_NOSCAL
     - NLOPT_GN_ORIG_DIRECT
     - NLOPT_GN_ORIG_DIRECT_L
     
     Controlled Random Search (CRS) with local mutation
     - NLOPT_GN_CRS2_LM
     
     MLSL (Multi-Level Single-Linkage)
     - NLOPT_G_MLSL_LDS
     - NLOPT_G_MLSL
     
     StoGO
     - NLOPT_GD_STOGO
     - NLOPT_GD_STOGO_RAND
     
     ISRES (Improved Stochastic Ranking Evolution Strategy)
     - NLOPT_GN_ISRES
     
     ESCH (evolutionary algorithm)
     - NLOPT_GN_ESCH
     
     COBYLA (Constrained Optimization BY Linear Approximations)
     - NLOPT_LN_COBYLA
     
     BOBYQA
     - NLOPT_LN_BOBYQA
     
     NEWUOA + bound constraints
     - NLOPT_LN_NEWUOA
     - NLOPT_LN_NEWUOA_BOUND
     
     PRAXIS (PRincipal AXIS)
     - NLOPT_LN_PRAXIS
     
     Nelder-Mead Simplex
     - NLOPT_LN_NELDERMEAD
     
     Sbplx (based on Subplex)
     - NLOPT_LN_SBPLX
     
     MMA (Method of Moving Asymptotes) and CCSA
     - NLOPT_LD_MMA
     - NLOPT_LD_CCSAQ
     
     SLSQP
     - NLOPT_LD_SLSQP
     
     Preconditioned truncated Newton
     - NLOPT_LD_TNEWTON_PRECOND_RESTART
     - NLOPT_LD_TNEWTON_PRECOND
     - NLOPT_LD_TNEWTON_RESTART
     - NLOPT_LD_TNEWTON
     
     Shifted limited-memory variable-metric
     - NLOPT_LD_VAR2
     - NLOPT_LD_VAR1
    */
        
    nlopt_set_local_optimizer(opt, local_opt);
    nlopt_set_lower_bounds(opt, bl);
    nlopt_set_upper_bounds(opt, bu);
    nlopt_set_xtol_rel(opt, 1e-15);
    nlopt_set_xtol_rel(local_opt, 1e-15);

    if (ncnln == 0)
    {
        nlopt_set_min_objective(opt, nloptObjectiveFunction, NULL);
    }
    
    else{
        /* Distinguish between equality and inequality constraints */
        for(int j = 0; j < globalState->numConstraints; j++) {
            if (globalState->conList[j].opCode == 1)
            {
                neq += globalState->conList[j].size;
            }
        }
        if (neq == ncnln) nineq = 0;
        else nineq = ncnln - neq;

        if (nlopt_get_algorithm(local_opt) == NLOPT_LN_COBYLA)
        {
            if(nineq > 0){
                nlopt_set_min_objective(opt, nloptObjectiveFunction, NULL);
                double tolineq[nineq];
                std::fill_n(tolineq, nineq, 1e-15);
                double *tol_ineq = tolineq;
                nlopt_add_inequality_mconstraint(opt, nineq, nloptInequalityFunction, NULL, tol_ineq);
            }
            if (neq > 0){
            Rf_error("COBYLA cannot handle equality constraints");
            }
        }
        else
        {
            if(nineq > 0){
                double tolineq[nineq];
                std::fill_n(tolineq, nineq, 1e-15);
                double *tol_ineq = tolineq;
                nlopt_add_inequality_mconstraint(opt, nineq, nloptInequalityFunction, NULL, tol_ineq);
            }
            
            if (neq > 0){
                double toleq[neq];
                std::fill_n(toleq, neq, 1e-15);
                double *tol_eq = toleq;
                nlopt_add_equality_mconstraint(opt, neq, nloptEqualityFunction, NULL, tol_eq);
            }
        }
        nlopt_set_min_objective(opt, nloptObjectiveFunction, NULL);
    }
    
    if(nlopt_optimize(opt, x, &opt_f) < 0)
    {   
        printf("nlopt failed! the status code is: %d\n",  nlopt_optimize(opt, x, &opt_f));
    }
    else{
        printf("found minimum at\n");
        for (int i = 0; i < n; i++){
            printf("%g\n", x[i]);
        }
        printf("objective value is: \n"); printf("%g", opt_f);
    }

    
    nlopt_sann_fc->fit = Fmin;
	nlopt_sann_fc->copyParamToModel();
    
    nlopt_sann_fitMatrix = NULL;
    nlopt_sann_fc = NULL;

    nlopt_destroy(local_opt);
    nlopt_destroy(opt);
    }
}

#endif
