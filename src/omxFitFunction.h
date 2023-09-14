/*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/***********************************************************
*
*  omxFitFunction.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxFitFunction class
*   omxFitFunction objects hold necessary information to simplify
* 	fit function calculation.
*
**********************************************************/

#ifndef u_OMXFITFUNCTION_H_
#define u_OMXFITFUNCTION_H_

#include <functional>
#include "omxDefines.h"
#include <R_ext/Rdynload.h>

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxData.h"
#include "omxState.h"
#include "omxExpectation.h"

class omxFitFunction {
 private:
  bool applyPenalty;
  std::vector< Penalty* > penalties;
protected:
 public:
  void subCompute(int want, FitContext *fc); // used by multigroup, hidden markov
  void recompute(int want, FitContext *fc) { subCompute(want, fc); } // used by omxRecompute

	SEXP rObj;
	omxExpectation* expectation;

	// This is always a pointer to a static string.
	// We do not need to allocate or free it.
	const char* fitType;

	omxMatrix* matrix;
  double scale;
	bool initialized;
	bool hessianAvailable;
	FitStatisticUnits units;
	bool canDuplicate;
	bool openmpUser; // can decide this in omxAlgebraPreeval
  int verbose;

	int derivCount;
	std::vector<int> gradMap;
	std::vector<int> missingGrad;

	omxFitFunction() : rObj(0), expectation(0), scale(1), initialized(false),
		hessianAvailable(false), units(FIT_UNITS_UNINITIALIZED), canDuplicate(false),
    openmpUser(false), verbose(0), derivCount(0) {};
	virtual ~omxFitFunction() {};
	virtual omxFitFunction *initMorph();
	virtual void init()=0;
	void compute(int ffcompute, FitContext *fc);
	virtual void compute2(int ffcompute, FitContext *fc)=0;
	virtual void invalidateCache() {};
	virtual void traverse(std::function<void(omxMatrix*)> fn);

	// addOutput should only be used for returning global results
	virtual void addOutput(MxRList *out) {};

	// populateAttr should be used for returning results specific to fit functions or expectations
	virtual void populateAttr(SEXP algebra) {};

  void buildGradMap(FitContext *fc, std::vector<const char *> &names, bool strict);
  void invalidateGradient(FitContext *fc);
	void setUnitsFromName(const char *name);
	const char *name() const { return matrix->name(); }
  void connectPenalties();
	virtual void sufficientDerivs2Grad() {};
};

/* Initialize and Destroy */
void omxFillMatrixFromMxFitFunction(omxMatrix* om, int matrixNumber, SEXP rObj);

omxFitFunction *omxChangeFitType(omxFitFunction *oo, const char *fitType);

void omxCompleteFitFunction(omxMatrix *om);

	void omxGetFitFunctionStandardErrors(omxFitFunction *oo);					// Get Standard Errors

/* FitFunction-specific implementations of matrix functions */
	void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* targetState);

omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* state, const char* slotName);
omxMatrix *omxNewMatrixFromSlotOrAnon(SEXP rObj, omxState* currentState, const char* slotName,
																			int rows, int cols);

omxFitFunction *omxInitFIMLFitFunction();
omxFitFunction *omxInitAlgebraFitFunction();
omxFitFunction *omxInitWLSFitFunction();
omxFitFunction *omxInitRowFitFunction();
omxFitFunction *omxInitMLFitFunction();
omxFitFunction *omxInitRFitFunction();
omxFitFunction *omxInitFitFunctionBA81();
omxFitFunction *InitMarkovFF();
omxFitFunction *omxInitGREMLFitFunction();
omxFitFunction *InitFellnerFitFunction();
omxFitFunction *ssMLFitInit();
omxFitFunction *GRMFIMLFitInit();

void ba81SetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg);

void ComputeFit(const char *callerName, omxMatrix *fitMat, int want, FitContext *fc);
void loglikelihoodCIFun(omxFitFunction* oo, int ffcompute, FitContext *fc);

const char *fitUnitsToName(FitStatisticUnits units);
bool fitUnitsIsChiSq(FitStatisticUnits units);
SEXP makeFitUnitsFactor(SEXP obj);

#endif /* u_OMXFITFUNCTION_H_ */
