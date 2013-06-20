#ifndef _OMX_EXPORT_BACKEND_STATE_H
#define _OMX_EXPORT_BACKEND_STATE_H

#include "types.h"
#include "omxState.h"

void omxExportResults(omxState *currentState, MxRList *out);
void omxPopulateFitFunction(omxMatrix *om, MxRList *result);
void omxPopulateConfidenceIntervals(SEXP intervals, SEXP intervalCodes);

#endif // #define _OMX_EXPORT_BACKEND_STATE_H
