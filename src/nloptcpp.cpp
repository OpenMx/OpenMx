#ifdef HAS_NLOPT

#include <stdio.h>
#include <nlopt.h>
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
#include <float.h>
#include <math.h>

static const char* anonMatrix = "anonymous matrix";
static omxMatrix * nlopt_sann_fitMatrix = NULL;
static FitContext * nlopt_sann_fc = NULL;
static int nlopt_sann_currentInterval = -1;
static bool nlopt_sann_useGradient;
static FitContext * nlopt_sann_fc_hhpx = NULL;
static FitContext * nlopt_sann_fc_hhmx = NULL;

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

void nloptSetupBounds(FreeVarGroup *freeVarGroup, double * bl, double * bu)
{
	size_t n = freeVarGroup->vars.size();
    omxState *globalState = nlopt_sann_fc->state;

	/* Set min and max limits */
	for(size_t index = 0; index < n; index++) {
		bl[index] = freeVarGroup->vars[index]->lbound;
		bu[index] = freeVarGroup->vars[index]->ubound;
	}
    int index = n;
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {		// Nonlinear constraints:
		if(OMX_DEBUG) { mxLog("Constraint %d: ", constraintIndex);}
		switch(globalState->conList[constraintIndex].opCode) {
            case 0:									// Less than: Must be strictly less than 0.
                if(OMX_DEBUG) { mxLog("Bounded at (-INF, 0).");}
                for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
                    bl[index] = NEG_INF;
                    bu[index] = -0.0;
                    index++;
                }
                break;
            case 1:									// Equal: Must be roughly equal to 0.
                if(OMX_DEBUG) { mxLog("Bounded at (-0, 0).");}
                for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
                    bl[index] = -0.0;
                    bu[index] = 0.0;
                    index++;
                }
                break;
            case 2:									// Greater than: Must be strictly greater than 0.
                if(OMX_DEBUG) { mxLog("Bounded at (0, INF).");}
                for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
                    if(OMX_DEBUG) { mxLog("\tBounds set for constraint %d.%d.", constraintIndex, offset);}
                    bl[index] = 0.0;
                    bu[index] = INF;
                    index++;
                }
                break;
            default:
                if(OMX_DEBUG) { mxLog("Bounded at (-INF, INF).");}
                for(int offset = 0; offset < globalState->conList[constraintIndex].size; offset++) {
                    bl[index] = NEG_INF;
                    bu[index] = INF;
                    index++;
                }
                break;
		}
	}
    
}

double fn(int n, double *par, void *ex)
{
    omxMatrix* fitMatrix = nlopt_sann_fitMatrix;
    
	nlopt_sann_fc->iterations += 1;   // ought to be major iterations only
    
	memcpy(nlopt_sann_fc->est, par, sizeof(double) * n);
	nlopt_sann_fc->copyParamToModel();
    
	ComputeFit(fitMatrix, FF_COMPUTE_FIT, nlopt_sann_fc);
    
	if (!std::isfinite(nlopt_sann_fc->fit) || isErrorRaised()) {
        nlopt_sann_fc->fit = 1e24;
	}
    
    printf("nlopt_sann_fc is what: %g \n", nlopt_sann_fc->fit);
    return nlopt_sann_fc->fit;
}

double nloptObjectiveFunction(unsigned n, const double* x, double* grad, void* f_data)
{
    // calculate gradient
    if (grad != NULL)
    {
        omxMatrix* fitMatrix_hhpx = nlopt_sann_fitMatrix;
        omxMatrix* fitMatrix_hhmx = nlopt_sann_fitMatrix;
        double fit_hhpx, fit_hhmx;
        double* hh = (double *)malloc(sizeof(double) * n);
        double* hhpx = (double *)malloc(sizeof(double) * n);
        double* hhmx = (double *)malloc(sizeof(double) * n);
        double* gr = (double *)malloc(sizeof(double) * n);
        double heps = pow(DBL_EPSILON, (double)1/3);
        
        for (int i = 0; i < n; i++){
            hh[i] = 0;
            hhpx[i] = 0;
            hhmx[i] = 0;
        }
        
        for (int i = 0; i < n; i++)
        {
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
        }
        for (int i = 0; i < n; i++)
        {
            hh[i] = heps;
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
            memcpy(nlopt_sann_fc_hhpx->est, hhpx, sizeof(double) * n);
            nlopt_sann_fc_hhpx->copyParamToModel();
            ComputeFit(fitMatrix_hhpx, FF_COMPUTE_FIT, nlopt_sann_fc_hhpx);
            if (!std::isfinite(nlopt_sann_fc_hhpx->fit) || isErrorRaised()) {
                nlopt_sann_fc_hhpx->fit = 1e24;
            }
            fit_hhpx = nlopt_sann_fc_hhpx->fit;
            memcpy(nlopt_sann_fc_hhmx->est, hhmx, sizeof(double) * n);
            nlopt_sann_fc_hhmx->copyParamToModel();
            ComputeFit(fitMatrix_hhmx, FF_COMPUTE_FIT, nlopt_sann_fc_hhmx);
            if (!std::isfinite(nlopt_sann_fc_hhmx->fit) || isErrorRaised()) {
                nlopt_sann_fc_hhmx->fit = 1e24;
            }
            
            fit_hhmx = nlopt_sann_fc_hhmx->fit;
            gr[i] = (fit_hhpx - fit_hhmx) / (2*heps);
            grad[i] = gr[i];
            hh[i] = 0;
            hhmx[i] = x[i];
            hhpx[i] = x[i];
        }
    }
    
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

double nloptLimitObjectiveFunction(unsigned n, const double* x, double* grad, void* f_data)
{
    
    nlopt_sann_fc->fit = nloptObjectiveFunction(n, x, grad, NULL);
    
    omxConfidenceInterval *oCI = Global->intervalList[nlopt_sann_currentInterval];
    
    omxRecompute(oCI->matrix, FF_COMPUTE_FIT, nlopt_sann_fc);
    
    double CIElement = omxMatrixElement(oCI->matrix, oCI->row, oCI->col);
    
    if(OMX_DEBUG) {
        mxLog("Finding Confidence Interval Likelihoods: lbound is %f, ubound is %f, estimate likelihood is %f, and element current value is %f.",
              oCI->lbound, oCI->ubound, nlopt_sann_fc->fit, CIElement);
    }
    
    /* Catch boundary-passing condition */
    if(std::isnan(CIElement) || std::isinf(CIElement)) {
	    nlopt_sann_fc->recordIterationError("Confidence interval is in a range that is currently incalculable. Add constraints to keep the value in the region where it can be calculated.");
        return nlopt_sann_fc->fit;
    }
    
    if(oCI->calcLower) {
        double diff = oCI->lbound - nlopt_sann_fc->fit;		// Offset - likelihood
        nlopt_sann_fc->fit = diff * diff + CIElement;
        // Minimize element for lower bound.
    } else {
        double diff = oCI->ubound - nlopt_sann_fc->fit;			// Offset - likelihood
        nlopt_sann_fc->fit = diff * diff - CIElement;
        // Maximize element for upper bound.
    }
    
    if(OMX_DEBUG) {
        mxLog("Interval fit function in previous iteration was calculated to be %f.", nlopt_sann_fc->fit);
    }
    
    return nlopt_sann_fc->fit;
}

void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
    int j, k, l = 0;
    
    omxState *globalState = nlopt_sann_fc->state;
    
    int nc = globalState->ncnln;
    
    // calculate Jacobian
    if (grad != NULL)
    {
        double* sub_result = (double *)malloc(sizeof(double) * m);
        double* result_hhpx = (double *)malloc(sizeof(double) * m);
        double* result_hhmx = (double *)malloc(sizeof(double) * m);
        double* res_jac = (double *)malloc(sizeof(double) * n);
        double* col_res = (double *)malloc(sizeof(double) * n);
        double* hh = (double *)malloc(sizeof(double) * n);
        double* hhpx = (double *)malloc(sizeof(double) * n);
        double* hhmx = (double *)malloc(sizeof(double) * n);
        double heps = pow(DBL_EPSILON, (double)1/3);
        
        for (int i = 0; i < n; i++){
            
            hh[i] = 0;
            hhpx[i] = 0;
            hhmx[i] = 0;
        }
        
        for (int i = 0; i < n; i++)
        {
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
        }
        for (int i = 0; i < n; i++)
        {
            l = 0;
            hh[i] = heps;
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
            
            memcpy(nlopt_sann_fc_hhpx->est, hhpx, sizeof(double) * n);
            nlopt_sann_fc_hhpx->copyParamToModel();
            
            for(j = 0; j < globalState->numConstraints; j++) {
                if ((globalState->conList[j].opCode == 0) || globalState->conList[j].opCode == 2) {
                    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc_hhpx);
                }
                for(k = 0; k < globalState->conList[j].size; k++){
                    result_hhpx[l] = globalState->conList[j].result->data[k];
                    l = l + 1;
                }
            }
            
            l = 0;
            
            memcpy(nlopt_sann_fc_hhmx->est, hhmx, sizeof(double) * n);
            nlopt_sann_fc_hhmx->copyParamToModel();
            
            for(j = 0; j < globalState->numConstraints; j++) {
                if ((globalState->conList[j].opCode == 0) || globalState->conList[j].opCode == 2) {
                    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc_hhmx);
                }
                for(k = 0; k < globalState->conList[j].size; k++){
                    result_hhmx[l] = globalState->conList[j].result->data[k];
                    l = l + 1;
                }
            }
            
            for (int ind = 0; ind < m; ind++){
                sub_result[ind] = (result_hhpx[ind] - result_hhmx[ind]) / (2*heps);
                grad[i + ind * n] = sub_result[ind];
            }
            
            hh[i] = 0;
            hhpx[i] = x[i];
            hhmx[i] = x[i];
        }
    }
    
    l = 0;
    
    memcpy(nlopt_sann_fc->est, x, sizeof(double) * n);
    nlopt_sann_fc->copyParamToModel();
    
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


void nloptEqualityFunction(unsigned m, double* result, unsigned n, const double* x, double* grad, void* h_data)

{
    int j, k, l = 0;
    
    omxState *globalState = nlopt_sann_fc->state;
    
    // calculate Jacobian
    if (grad != NULL)
    {
        double* sub_result = (double *)malloc(sizeof(double) * m);
        double* result_hhpx = (double *)malloc(sizeof(double) * m);
        double* result_hhmx = (double *)malloc(sizeof(double) * m);
        double* res_jac = (double *)malloc(sizeof(double) * n);
        double* col_res = (double *)malloc(sizeof(double) * n);
        double* hh = (double *)malloc(sizeof(double) * n);
        double* hhpx = (double *)malloc(sizeof(double) * n);
        double* hhmx = (double *)malloc(sizeof(double) * n);
        double* gr = (double *)malloc(sizeof(double) * n);
        double heps = pow(DBL_EPSILON, (double)1/3);
        
        for (int i = 0; i < n; i++){
            hh[i] = 0;
            hhpx[i] = 0;
            hhmx[i] = 0;
        }
        
        for (int i = 0; i < n; i++)
        {
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
        }
        for (int i = 0; i < n; i++)
        {
            l = 0;
            hh[i] = heps;
            hhpx[i] = x[i] + hh[i];
            hhmx[i] = x[i] - hh[i];
            memcpy(nlopt_sann_fc_hhpx->est, hhpx, sizeof(double) * n);
            nlopt_sann_fc_hhpx->copyParamToModel();
            
            for(j = 0; j < globalState->numConstraints; j++) {
                if (globalState->conList[j].opCode == 1) {
                    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc_hhpx);
                }
                for(k = 0; k < globalState->conList[j].size; k++){
                    result_hhpx[l] = globalState->conList[j].result->data[k];
                    l = l + 1;
                }
            }
            
            l = 0;
            
            memcpy(nlopt_sann_fc_hhmx->est, hhmx, sizeof(double) * n);
            nlopt_sann_fc_hhmx->copyParamToModel();
            
            for(j = 0; j < globalState->numConstraints; j++) {
                if (globalState->conList[j].opCode == 1) {
                    omxRecompute(globalState->conList[j].result, FF_COMPUTE_FIT, nlopt_sann_fc_hhmx);
                }
                for(k = 0; k < globalState->conList[j].size; k++){
                    result_hhmx[l] = globalState->conList[j].result->data[k];
                    l = l + 1;
                }
            }
            
            for (int ind = 0; ind < m; ind++){
                sub_result[ind] = (result_hhpx[ind] - result_hhmx[ind]) / (2*heps);
                grad[i + ind * n] = sub_result[ind];
            }
            
            hh[i] = 0;
            hhpx[i] = x[i];
            hhmx[i] = x[i];
        }
    }
    
    l = 0;
    memcpy(nlopt_sann_fc->est, x, sizeof(double) * n);
    nlopt_sann_fc->copyParamToModel();
    
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
    omxState *globalState = nlopt_sann_fc->state;
	double *x = fc->est;
	fc->grad.resize(fc->numParam);
	double *g = fc->grad.data();
    int n = int(freeVarGroup->vars.size());
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
             arg 7 is temp: The starting temperature for the cooling schedule. Defaults to 10.
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
        std::vector<double> blBuf(n + ncnln);
        std::vector<double> buBuf(n + ncnln);
        double *bl = blBuf.data();
        double *bu = buBuf.data();
        nloptSetupBounds(freeVarGroup, bl, bu);
        
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
        nlopt_opt opt; //local_opt;
        double opt_f;
        opt = nlopt_create(NLOPT_LD_SLSQP, n); // Main algorithm
        //local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm
        nlopt_sann_useGradient = TRUE;
        if (nlopt_sann_useGradient == TRUE)
        {
            nlopt_sann_fc_hhpx = fc;
            nlopt_sann_fc_hhmx = fc;
        }
        
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
        
        //nlopt_set_local_optimizer(opt, local_opt);
        nlopt_set_lower_bounds(opt, bl);
        nlopt_set_upper_bounds(opt, bu);
        nlopt_set_xtol_rel(opt, 1e-20);
        //nlopt_set_xtol_rel(local_opt, 1e-15);
        
        if (ncnln == 0)
        {
            nlopt_set_min_objective(opt, nloptObjectiveFunction, NULL);
        }
        
        else{
            nlopt_set_min_objective(opt, nloptObjectiveFunction, NULL);
            /* Distinguish between equality and inequality constraints */
            for(int j = 0; j < globalState->numConstraints; j++) {
                if (globalState->conList[j].opCode == 1)
                {
                    neq += globalState->conList[j].size;
                }
            }
            if (neq == ncnln) nineq = 0;
            else nineq = ncnln - neq;
            
            if (nlopt_get_algorithm(opt) == NLOPT_LN_COBYLA)
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
                    //double tol_eq = 1e-15;
                    nlopt_add_equality_mconstraint(opt, neq, nloptEqualityFunction, NULL, tol_eq);
                }
            }
        }
        
        if (nlopt_optimize(opt, x, &opt_f) < 0)
        {
            printf("nlopt failed! the status code is: %d\n",  nlopt_optimize(opt, x, &opt_f));
        }
        printf("found minimum at\n");
        for (int i = 0; i < n; i++){
            printf("%g\n", x[i]);
        }
        printf("objective value is: \n"); printf("%g", opt_f);
        
        
        nlopt_sann_fc->fit = opt_f;
        nlopt_sann_fc->copyParamToModel();
        
        nlopt_sann_fitMatrix = NULL;
        nlopt_sann_fc = NULL;
        
        //nlopt_destroy(local_opt);
        nlopt_destroy(opt);
    }
}

void omxNLOPTorSANNConfidenceIntervals(omxMatrix *fitMatrix, FitContext *opti, double tolerance)
{
    FitContext fc(opti, opti->varGroup);
    
	int ciMaxIterations = Global->ciMaxIterations;
	// Will fail if we re-enter after an exception
	//if (NPSOL_fitMatrix) Rf_error("NPSOL is not reentrant");
	nlopt_sann_fitMatrix = fitMatrix;
	nlopt_sann_fc = &fc;
	FreeVarGroup *freeVarGroup = fc.varGroup;
    
    int n = int(freeVarGroup->vars.size());
    omxState *globalState = nlopt_sann_fc->state;
    int ncnln = globalState->ncnln;
    int nineq = 0;
    int neq = 0;
    
    double f = opti->fit;
    omxBuffer< double > gradient(n);
    omxBuffer< double > hessian(n * n);
    
    /* Set boundaries and widths. */
    /* Allocate arrays */
    std::vector<double> blBuf(n+ncnln);
    std::vector<double> buBuf(n+ncnln);
    double *bl = blBuf.data();
    double *bu = buBuf.data();
    int inform;
    
    nloptSetupBounds(freeVarGroup, bl, bu);
    
    if(OMX_DEBUG) { mxLog("Calculating likelihood-based confidence intervals."); }
    
	const double objDiff = 1.e-4;     // TODO : Use function precision to determine CI jitter?
    
    for(int i = 0; i < (int) Global->intervalList.size(); i++) {
        
		omxConfidenceInterval *currentCI = Global->intervalList[i];
        
		const char *matName = anonMatrix;
		if (currentCI->matrix->name) {
			matName = currentCI->matrix->name;
		}
		Global->checkpointMessage(opti, opti->est, "%s[%d, %d] begin lower interval",
                                  matName, currentCI->row + 1, currentCI->col + 1);
        
        
		memcpy(fc.est, opti->est, n * sizeof(double)); // Reset to previous optimum
        double *x = opti->est;
        nlopt_sann_currentInterval = i;
        
        currentCI->lbound += opti->fit;          // Convert from offsets to targets
        currentCI->ubound += opti->fit;          // Convert from offsets to targets
        
        nlopt_opt opt;
        double opt_f = opti->fit;
        opt = nlopt_create(NLOPT_LD_SLSQP, n);
        nlopt_sann_useGradient = TRUE;
        if (nlopt_sann_useGradient == TRUE)
        {
            nlopt_sann_fc_hhpx = &fc;
            nlopt_sann_fc_hhmx = &fc;
        }
        nlopt_set_lower_bounds(opt, bl);
        nlopt_set_upper_bounds(opt, bu);
        nlopt_set_xtol_rel(opt, 1e-20);
        
        if (ncnln == 0)
        {
            nlopt_set_min_objective(opt, nloptLimitObjectiveFunction, NULL);
        }
        else{
            nlopt_set_min_objective(opt, nloptLimitObjectiveFunction, NULL);
            /* Distinguish between equality and inequality constraints */
            for(int j = 0; j < globalState->numConstraints; j++) {
                if (globalState->conList[j].opCode == 1)
                {
                    neq += globalState->conList[j].size;
                }
            }
            if (neq == ncnln) nineq = 0;
            else nineq = ncnln - neq;
            
            if (nlopt_get_algorithm(opt) == NLOPT_LN_COBYLA)
            {
                if(nineq > 0){
                    nlopt_set_min_objective(opt, nloptLimitObjectiveFunction, NULL);
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
                    //double tol_eq = 1e-15;
                    nlopt_add_equality_mconstraint(opt, neq, nloptEqualityFunction, NULL, tol_eq);
                }
            }
        }
        
        
	    if (std::isfinite(currentCI->lbound)) {
            /* Set up for the lower bound */
            inform = -1;
            // Number of times to keep trying.
            int cycles = ciMaxIterations;
            double value = INF;
            while(inform != 0 && cycles > 0) {
                /* Find lower limit */
                currentCI->calcLower = TRUE;
                if (nlopt_optimize(opt, x, &opt_f) == -1)   inform = 6;
                else inform = 1;
                f = opt_f;
                currentCI->lCode = inform;
                if(f < value) {
                    currentCI->min = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
                    memcpy(fc.est, x, sizeof(double) * n);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opti->est[j]) > objDiff) {
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opti->est[j] + objDiff;
                        }
                    }
                }
            }
	    }
        
	    if (std::isfinite(currentCI->ubound)) {
            currentCI->calcLower = FALSE;
            Global->checkpointMessage(opti, opti->est, "%s[%d, %d] begin upper interval",
                                      matName, currentCI->row + 1, currentCI->col + 1);
            
            memcpy(fc.est, opti->est, n * sizeof(double)); // Reset to previous optimum
            x = opti->est;
            /* Reset for the upper bound */
            double value = INF;
            inform = -1;
            double cycles = ciMaxIterations;
            
            while(inform != 0 && cycles >= 0) {
                /* Find upper limit */
                currentCI->calcLower = FALSE;
                
                if (nlopt_optimize(opt, x, &opt_f) == -1)   inform = 6;
                else inform = 1;
                f = opt_f;
                currentCI->uCode = inform;
                if(f < value) {
                    currentCI->max = omxMatrixElement(currentCI->matrix, currentCI->row, currentCI->col);
                    value = f;
                    memcpy(fc.est, x, sizeof(double) * n);
                }
                
                if(inform != 0 && OMX_DEBUG) {
                    mxLog("Calculation of upper interval %d failed: Bad inform value of %d",
                          i, inform);
                }
                cycles--;
                if(inform != 0) {
                    unsigned int jitter = TRUE;
                    for(int j = 0; j < n; j++) {
                        if(fabs(fc.est[j] - opti->est[j]) > objDiff){
                            jitter = FALSE;
                            break;
                        }
                    }
                    if(jitter) {
                        for(int j = 0; j < n; j++) {
                            fc.est[j] = opti->est[j] + objDiff;
                        }
                    }
                }
            }
            if(OMX_DEBUG) {mxLog("Found Upper bound %d.", i);}
        }
	}
    
	nlopt_sann_fc = NULL;
    nlopt_sann_fitMatrix = NULL;
	nlopt_sann_currentInterval = -1;
}

#endif
