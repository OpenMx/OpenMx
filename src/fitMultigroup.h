#ifndef _expectationMultigroup_H_
#define _expectationMultigroup_H_

omxFitFunction *initFitMultigroup();
void omxMultigroupAdd(omxFitFunction *oo, omxFitFunction *grp);
void mgSetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg);

#endif
