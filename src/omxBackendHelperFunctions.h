#ifndef _OMX_BACKEND_HELPER_FUNCTIONS_H
#define _OMX_BACKEND_HELPER_FUNCTIONS_H


int matchCaseInsensitive(const char *source, int lenSource, const char *target);
int omxProcessMxDataEntities(SEXP data);
int omxProcessMxMatrixEntities(SEXP matList);
int omxProcessMxAlgebraEntities(SEXP algList);
int omxInitialMatrixAlgebraCompute();
int omxProcessObjectiveFunction(SEXP objective, int *n);
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, int n);

#endif // #define _OMX_BACKEND_HELPER_FUNCTIONS_H
