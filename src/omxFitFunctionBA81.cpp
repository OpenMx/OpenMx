/*
 * This is throw-away proof-of-concept code that will likely be
 * replaced by something else.
 *
 * JNP 2012Dec10
 */

#include "omxFitFunction.h"
#include "omxExpectationBA81.h"

static const char *NAME = "FitFunctionBA81";

typedef struct {

	omxData *data;

} omxBA81State;


static void ba81Destroy(omxFitFunction *oo) {
	if(OMX_DEBUG) {
		Rprintf("Freeing %s function.\n", NAME);
	}
	//omxBA81State *mml = (omxBA81State *) oo->argStruct;
	// nothing to do yet
}

// TODO: Don't trample the Expectation/FitFunction separation.

static omxRListElement *ba81SetFinalReturns(omxFitFunction *off, int *numReturns)
{
	return ba81EAP(off->expectation, numReturns);
}

static void ba81Compute(omxFitFunction *oo, int want, double *gradient, double *hessian) {
	if(OMX_DEBUG_MML) {Rprintf("Beginning %s Computation.\n", NAME);}

	omxExpectation* expectation = oo->expectation;
  
	oo->matrix->data[0] = ba81ComputeFit(expectation, want, gradient, hessian);
}

void omxInitFitFunctionBA81(omxFitFunction* oo, SEXP rObj) {
	//omxExpectation* expectation = oo->expectation;

	//omxState* currentState = oo->matrix->currentState;
	
	if(OMX_DEBUG) {
	  Rprintf("Initializing %s.\n", NAME);
	}
	
	//omxBA81State *newObj = (omxBA81State*) R_alloc(1, sizeof(omxBA81State));
	
	//newObj->data = oo->expectation->data;

	omxExpectationCompute(oo->expectation, COMPUTE_EXPECT_INITIALIZE);

	oo->computeFun = ba81Compute;
	oo->setFinalReturns = ba81SetFinalReturns;
	oo->destructFun = ba81Destroy;
}
