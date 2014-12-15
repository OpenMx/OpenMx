#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H

int matchCaseInsensitive(const char *source, const char *target);
void omxProcessMxDataEntities(SEXP data);
void omxInitialMatrixAlgebraCompute(omxState *state, FitContext *fc);
void omxProcessCheckpointOptions(SEXP checkpointList);
void omxProcessFreeVarList(SEXP varList, std::vector<double> *startingValues);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
