#ifndef _OMX_EXPORT_BACKEND_STATE_H
#define _OMX_EXPORT_BACKEND_STATE_H

#include "types.h"
#include "omxState.h"

void omxFinalAlgebraCalculation(omxState* currentState, SEXP matrices, SEXP algebras, SEXP expectations);
void omxPopulateFitFunction(omxState *currentState, MxRList *result);
void omxPopulateHessians(int numHessians, omxMatrix* currentFit, 
	SEXP calculatedHessian, SEXP stdErrors, int calculateStdErrors, int n);
void omxPopulateConfidenceIntervals(omxState* currentState, SEXP intervals, SEXP intervalCodes);

#endif // #define _OMX_EXPORT_BACKEND_STATE_H
