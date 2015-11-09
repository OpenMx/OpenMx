#ifndef _OMX_IMPORT_FRONTEND_STATE_H
#define _OMX_IMPORT_FRONTEND_STATE_H

int matchCaseInsensitive(const char *source, const char *target);
void omxInitialMatrixAlgebraCompute(omxState *state, FitContext *fc);
void omxProcessCheckpointOptions(SEXP checkpointList);

#endif // #define _OMX_IMPORT_FRONTEND_STATE_H
