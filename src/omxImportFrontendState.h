#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H


int matchCaseInsensitive(const char *source, const char *target);
int omxProcessMxDataEntities(SEXP data);
int omxProcessMxMatrixEntities(SEXP matList);
void omxProcessMxAlgebraEntities(SEXP algList);
int omxProcessMxExpectationEntities(SEXP expList);
int omxCompleteMxExpectationEntities();
void omxInitialMatrixAlgebraCompute();
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, int n);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);
void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
