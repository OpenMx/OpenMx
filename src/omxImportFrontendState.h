#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H


int matchCaseInsensitive(const char *source, int lenSource, const char *target);
int omxProcessMxDataEntities(SEXP data);
int omxProcessMxMatrixEntities(SEXP matList);
int omxProcessMxAlgebraEntities(SEXP algList);
int omxProcessMxExpectationEntities(SEXP expList);
int omxCompleteMxExpectationEntities();
int omxInitialMatrixAlgebraCompute();
int omxProcessFitFunction(SEXP fitFunction);
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, int n);
void omxProcessConfidenceIntervals(SEXP intervalList);
void omxProcessConstraints(SEXP constraints);
void omxSetupBoundsAndConstraints(double * bl, double * bu, int n, int nclin);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
