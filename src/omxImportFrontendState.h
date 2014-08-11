#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H

extern Matrix result_csolnpEqBStartFun, result_csolnpEqB, result_ineqLB, result_ineqUB,  result_ineqVal;
int matchCaseInsensitive(const char *source, const char *target);
void omxProcessMxDataEntities(SEXP data);
void omxInitialMatrixAlgebraCompute(omxState *state, FitContext *fc);
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues);

void omxProcessConstraintsCsolnp(FitContext *fc, struct Matrix *lb_ineq, struct Matrix *ub_ineq, struct Matrix *eqb);
void omxSetupBoundsAndConstraints(FitContext *fc, double * bl, double * bu);
void setupIneqLess(Matrix bl_ineqless, Matrix bu_ineqless, int size);
void setupIneqGreater(struct Matrix* bl_ineqmore, struct Matrix* bu_ineqmore, int size);
void setupEqB(Matrix bl_eqb, Matrix bu_eqb, int size);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
