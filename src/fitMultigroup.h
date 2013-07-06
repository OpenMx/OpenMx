#ifndef _expectationMultigroup_H_
#define _expectationMultigroup_H_

void initFitMultigroup(omxFitFunction *off);
void omxMultigroupAdd(omxFitFunction *oo, omxFitFunction *grp);
void mgSetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg);

#endif
