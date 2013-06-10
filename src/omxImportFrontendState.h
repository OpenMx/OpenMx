#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H


int matchCaseInsensitive(const char *source, const char *target);
void omxProcessMxDataEntities(SEXP data);
void omxProcessMxMatrixEntities(SEXP matList);
void omxProcessMxAlgebraEntities(SEXP algList);
void omxProcessMxExpectationEntities(SEXP expList);
void omxCompleteMxExpectationEntities();
void omxInitialMatrixAlgebraCompute();
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);
void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
