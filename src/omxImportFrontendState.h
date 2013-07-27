#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H

extern Matrix result_csolnpEqBStartFun, result_csolnpEqB, result_ineqLB, result_ineqUB,  result_ineqVal;
int matchCaseInsensitive(const char *source, const char *target);
void omxProcessMxDataEntities(SEXP data);
void omxProcessMxMatrixEntities(SEXP matList);
void omxProcessMxAlgebraEntities(SEXP algList);
void omxProcessMxExpectationEntities(SEXP expList);
void omxCompleteMxExpectationEntities();
void omxProcessMxFitFunction(SEXP algList);
void omxCompleteMxFitFunction(SEXP algList);
void omxProcessMxComputeEntities(SEXP computeList);
void omxInitialMatrixAlgebraCompute();
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);

void omxProcessConstraintsCsolnp(struct Matrix *lb_ineq, struct Matrix *ub_ineq, struct Matrix *eqb);
void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin);
void setupIneqLess(Matrix bl_ineqless, Matrix bu_ineqless, int size);
void setupIneqGreater(struct Matrix* bl_ineqmore, struct Matrix* bu_ineqmore, int size);
void setupEqB(Matrix bl_eqb, Matrix bu_eqb, int size);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
