#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

SEXP omx_get_rpf_names();
void ba81Gradient(omxExpectation* oo, double *out);
int ba81ExpectationHasGradients(omxExpectation* oo);
double ba81ComputeFit(omxExpectation* oo);

omxRListElement *
ba81EAP(omxExpectation *oo, int *numReturns);


#endif
