//
//  File.c
//  
//
//  Created by Mahsa Zahery on 7/23/14.
//
//

double nloptObjectiveFunction(int* n, double* x, double* f);
void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);
void nloptEqualityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);
void omxInvokeNLOPTorSANN(omxMatrix *fitMatrix, FitContext *fc, int *inform_out, FreeVarGroup *freeVarGroup, int verbose, double *hessOut, double tolerance);
