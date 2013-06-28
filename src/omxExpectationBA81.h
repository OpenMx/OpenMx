#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

SEXP omx_get_rpf_names();
void ba81Gradient(omxExpectation* oo, double *out);
double ba81ComputeFit(omxExpectation* oo, int want, double *gradient);

omxRListElement *
ba81EAP(omxExpectation *oo, int *numReturns);


#endif
