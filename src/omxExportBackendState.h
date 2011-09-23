#ifndef _OMX_EXPORT_BACKEND_STATE_H
#define _OMX_EXPORT_BACKEND_STATE_H

void omxFinalAlgebraCalculation(SEXP matrices, SEXP algebras);
void omxPopulateObjectiveFunction(int numReturns, SEXP *ans, SEXP *names);
void omxPopulateHessians(int numHessians, omxMatrix* currentObjective, 
	SEXP calculatedHessian, SEXP stdErrors, int calculateStdErrors, int n);
void omxPopulateConfidenceIntervals(SEXP intervals, SEXP intervalCodes);

#endif // #define _OMX_EXPORT_BACKEND_STATE_H
