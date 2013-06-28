#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

SEXP omx_get_rpf_names();
void ba81Gradient(omxExpectation* oo, double *out);
int ba81ExpectationHasGradients(omxExpectation* oo);
double ba81ComputeFit(omxExpectation* oo);
void ba81EAP(omxExpectation *oo, omxRListElement *out);

#endif
