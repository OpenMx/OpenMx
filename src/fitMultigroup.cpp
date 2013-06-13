#include "omxExpectation.h"
#include "omxOptimizer.h"
#include "fitMultigroup.h"
#include <vector>

// http://openmx.psyc.virginia.edu/issue/2013/01/multigroup-fit-function

struct FitMultigroup {
	std::vector< int > fits;  // store pointers or index numbers? TODO
	bool checkedRepopulate;
	FitMultigroup() : checkedRepopulate(0) {}
};

static void mgDestroy(omxFitFunction* oo)
{
	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	delete mg;
}

static void checkRepopulate(omxFitFunction* oo)
{
	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	omxState *os = oo->matrix->currentState;
	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = os->algebraList[mg->fits[ex]];
		omxFitFunction *ff = f1->fitFunction;
		if (!ff || ff->repopulateFun == handleFreeVarList) continue;
		error("Cannot add %s to multigroup fit", f1->name);
	}
}

static void mgCompute(omxFitFunction* oo, int ffcompute, double* grad)
{
	omxMatrix *fitMatrix  = oo->matrix;
	omxState *os = fitMatrix->currentState;
	fitMatrix->data[0] = 0;

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;

	if (!mg->checkedRepopulate) {
		checkRepopulate(oo);
		mg->checkedRepopulate = TRUE;
	}

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = os->algebraList[mg->fits[ex]];
		if (f1->fitFunction) {
			// possibly invalidate gradients TODO
			omxFitFunctionCompute(f1->fitFunction, ffcompute, grad);
		} else {
			// invalidate gradients TODO
			omxRecompute(f1);

			// This should really be checked elsewhere. TODO
			if(f1->rows != 1 || f1->cols != 1) {
				error("%s algebra %d does not evaluate to a 1x1 matrix", oo->fitType, ex);
			}
		}
		fitMatrix->data[0] += f1->data[0];
	}
	if(OMX_DEBUG) { Rprintf("Fit Function sum of %d groups is %f.\n", mg->fits.size(), fitMatrix->data[0]); }
}

void initFitMultigroup(omxFitFunction *oo)
{
	oo->expectation = NULL;  // don't care about this
	oo->computeFun = mgCompute;
	oo->destructFun = mgDestroy;
	oo->repopulateFun = handleFreeVarList;

	FitMultigroup *mg = new FitMultigroup;
	oo->argStruct = mg;

	SEXP rObj = oo->rObj;
	if (!rObj) return;

	int myIndex = oo->matrix->matrixNumber;

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("groups")));
	int *fits = INTEGER(slotValue);
	for(int gx = 0; gx < length(slotValue); gx++) {
		if (fits[gx] == myIndex) error("Cannot add multigroup to itself");
		mg->fits.push_back(fits[gx]);
	}
}

void omxMultigroupAdd(omxFitFunction *oo, omxFitFunction *grp)
{
	if (oo->initFun != initFitMultigroup) error("%s is not the multigroup fit", oo->fitType);
	if (!oo->initialized) error("Fit %p not initialized", oo);

	int myIndex = ~oo->matrix->matrixNumber;
	int grpIndex = ~grp->matrix->matrixNumber;
	if (grpIndex == myIndex) error("Cannot add multigroup to itself");

	omxState *os = oo->matrix->currentState;
	if (os->algebraList[grpIndex] != grp->matrix) {
		error("Attempt to add group %d missing from algebraList", grpIndex);
	}

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	mg->fits.push_back(grpIndex);
	//addFreeVarDependency(oo->matrix->currentState, oo->matrix, grp->matrix);
}
