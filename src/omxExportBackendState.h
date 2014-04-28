#ifndef _OMX_EXPORT_BACKEND_STATE_H
#define _OMX_EXPORT_BACKEND_STATE_H

#include "omxDefines.h"
#include "omxState.h"

void omxExportResults(omxState *currentState, MxRList *out);
void omxPopulateFitFunction(omxMatrix *om, MxRList *result);

#endif // #define _OMX_EXPORT_BACKEND_STATE_H
