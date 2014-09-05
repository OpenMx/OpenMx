//
//  nloptcpp.h
//  
//
//  Created by Mahsa Zahery on 7/23/14.
//
//
#include <stdio.h>

void nloptSetupBounds(FreeVarGroup *freeVarGroup, double * bl, double * bu);
double nloptObjectiveFunction(int* n, double* x, double* f);
void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);
void nloptEqualityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* h_data);
void omxInvokeNLOPTorSANN(omxMatrix *fitMatrix, FitContext *fc, int *inform_out, FreeVarGroup *freeVarGroup, int verbose, double *hessOut, double tolerance);
double fn(int n, double *par, void *ex);
void omxNLOPTorSANNConfidenceIntervals(omxMatrix *fitMatrix, FitContext *opt, double tolerance);
double nloptLimitObjectiveFunction(unsigned n, const double* x, double* grad, void* f_data);