#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H


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
void omxProcessFreeVarList(SEXP fgNames, SEXP varList);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
