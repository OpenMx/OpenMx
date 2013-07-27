#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H

extern Matrix result_csolnpEqBStartFun, result_csolnpEqB, result_ineqLB, result_ineqUB,  result_ineqVal;
int matchCaseInsensitive(const char *source, const char *target);
int omxProcessMxDataEntities(SEXP data);
int omxProcessMxMatrixEntities(SEXP matList);
void omxProcessMxAlgebraEntities(SEXP algList);
int omxProcessMxExpectationEntities(SEXP expList);
int omxCompleteMxExpectationEntities();
void omxInitialMatrixAlgebraCompute();
void omxProcessFitFunction(SEXP fitFunction);
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, int n);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);
void omxProcessConstraintsCsolnp(struct Matrix *lb_ineq, struct Matrix *ub_ineq, struct Matrix *eqb);
//Matrix** omxProcessConstraints(SEXP constraints);
void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin);
void setupIneqLess(Matrix bl_ineqless, Matrix bu_ineqless, int size);
void setupIneqGreater(struct Matrix* bl_ineqmore, struct Matrix* bu_ineqmore, int size);
void setupEqB(Matrix bl_eqb, Matrix bu_eqb, int size);
#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
