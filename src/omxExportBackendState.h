#ifndef _OMX_EXPORT_BACKEND_STATE_H
#define _OMX_EXPORT_BACKEND_STATE_H

#include "types.h"
#include "omxState.h"

void omxFinalAlgebraCalculation(omxState* currentState, SEXP matrices, SEXP algebras, SEXP expectations);
void omxPopulateFitFunction(omxMatrix *om, MxRList *result);
void omxPopulateConfidenceIntervals(omxState* currentState, SEXP intervals, SEXP intervalCodes);

#endif // #define _OMX_EXPORT_BACKEND_STATE_H
