/*
 *  Copyright 2007-2020 by the individuals mentioned in the source code history
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/***********************************************************
*
*  omxFitFunction.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	FitFunction objects are a subclass of data matrix that evaluates
*   itself anew at each iteration, so that any changes to
*   free parameters can be incorporated into the update.
*   // Question: Should FitFunction be a ``subtype'' of
*   // omxAlgebra or a separate beast entirely?
*
**********************************************************/

#include "omxFitFunction.h"
#include "fitMultigroup.h"
#include "Compute.h"
#include "finiteDifferences.h"
#include "EnableWarnings.h"

typedef struct omxFitFunctionTableEntry omxFitFunctionTableEntry;

struct omxFitFunctionTableEntry {

	char name[32];
	omxFitFunction *(*allocate)();

};

omxFitFunction *omxFitFunction::initMorph()
{
	init();
	return this;
}

void omxFitFunction::buildGradMap(FitContext *fc, std::vector<const char *> &names,
                                  bool strict)
{
  if (fc->getNumFree() == -1) mxThrow("Forgot to call fc->calcNumFree");

  std::vector<bool> haveGrad(fc->getNumFree(), false);
  derivCount = 0;
  int numNames = names.size();
	gradMap.resize(numNames);
  auto &fim = fc->freeToIndexMap;

	for (int nx=0; nx < numNames; ++nx) {
		auto it = fim.find(names[nx]);
    if (it == fim.end()) {
      gradMap[nx] = -1;
      if (strict) {
        mxThrow("Fit function '%s' has a derivative entry for unrecognized "
                "parameter '%s'. If this is not an mistake and you "
                "have merely fixed this parameter then you can "
                "use the strict=FALSE argument to mxFitFunction "
                "to turn off this precautionary check", name(), names[nx]);
      }
    } else {
      gradMap[nx] = it->second;
      haveGrad[it->second] = true;
      ++derivCount;
      if (verbose) {
        mxLog("%s: name '%s' mapped to free parameter %d",
              matrix->name(), names[nx], it->second);
      }
    }
	}

  missingGrad.clear();
  missingGrad.reserve(fc->getNumFree() - derivCount);
  for (int fx=0; fx < fc->getNumFree(); ++fx) {
    if (haveGrad[fx]) continue;
    missingGrad.push_back(fx);
  }
}

void omxFitFunction::invalidateGradient(FitContext *fc)
{
  if (derivCount == 0) {
    fc->gradZ.setConstant(NA_REAL);
  } else {
    for (int ix : missingGrad) fc->gradZ[ix] = NA_REAL;
  }
}

static const omxFitFunctionTableEntry omxFitFunctionSymbolTable[] = {
	{"MxFitFunctionAlgebra", 			&omxInitAlgebraFitFunction},
	{"MxFitFunctionWLS",				&omxInitWLSFitFunction},
	{"MxFitFunctionRow", 				&omxInitRowFitFunction},
	{"MxFitFunctionML", 				&omxInitMLFitFunction},
	{"imxFitFunctionFIML", &omxInitFIMLFitFunction},
	{"MxFitFunctionR",					&omxInitRFitFunction},
	{"MxFitFunctionMultigroup", &initFitMultigroup},
	{"MxFitFunctionGREML", &omxInitGREMLFitFunction},
	{"imxFitFunctionFellner", &InitFellnerFitFunction},
	{"imxFitFunctionBA81", &omxInitFitFunctionBA81},
	{"imxFitFunciontStateSpace", &ssMLFitInit},
	{"imxFitFunciontHiddenMarkov", &InitMarkovFF},
	{"imxFitFunciontGRMFIML", &GRMFIMLFitInit},
};

void omxFitFunction::setUnitsFromName(const char *name)
{
	if (strEQ(name, "-2lnL")) {
		units = FIT_UNITS_MINUS2LL;
  } else if (strEQ(name, "r'Wr")) {
    units = FIT_UNITS_SQUARED_RESIDUAL;
	} else {
		Rf_warning("Unknown units '%s' passed to fit function '%s'",
			   name, matrix->name());
		units = FIT_UNITS_UNKNOWN;
	}
}

bool fitUnitsIsChiSq(FitStatisticUnits units)
{
	return units == FIT_UNITS_MINUS2LL || units == FIT_UNITS_SQUARED_RESIDUAL_CHISQ;
}

static const char* FitUnitNames[] = { "?", "Pr", "-2lnL", "r'Wr", "r'Wr" };

SEXP makeFitUnitsFactor(SEXP obj)
{
	return makeFactor(obj, OMX_STATIC_ARRAY_SIZE(FitUnitNames), FitUnitNames);
}

const char *fitUnitsToName(FitStatisticUnits units)
{
	switch (units) {
	case FIT_UNITS_UNINITIALIZED: return "";
	case FIT_UNITS_UNKNOWN:
	case FIT_UNITS_PROBABILITY:
	case FIT_UNITS_MINUS2LL:
	case FIT_UNITS_SQUARED_RESIDUAL:
	case FIT_UNITS_SQUARED_RESIDUAL_CHISQ:
		return FitUnitNames[units-1];
	default: mxThrow("Don't know how to stringify units %d", units);
	}
}

void omxDuplicateFitMatrix(omxMatrix *tgt, const omxMatrix *src, omxState* newState) {

	if(tgt == NULL || src == NULL) return;

	omxFitFunction *ff = src->fitFunction;
	if(ff == NULL) return;

	omxFillMatrixFromMxFitFunction(tgt, src->matrixNumber, ff->rObj);
}

static void omxFitFunctionComputeAuto(omxFitFunction *ff, int want, FitContext *fc)
{
	if (!ff->initialized) return;

	if (!fc->ciobj) {
		ff->compute(want, fc);
	} else {
    if (fitUnitsIsChiSq(ff->units)) {
      fc->ciobj->evalFit(ff, want, fc);
    } else {
      mxThrow("Confidence intervals are not supported for units %s",
              fitUnitsToName(ff->units));
    }
	}

	if (fc) fc->wanted |= want;
}

// Can replace with
//   fc->withoutCIobjective([&](){ ComputeFit("CI", fitMat, FF_COMPUTE_FIT, fc); });
// ?
void omxFitFunctionCompute(omxFitFunction *off, int want, FitContext *fc)
{
       if (!off->initialized) return;

       off->compute(want, fc);
       if (fc) fc->wanted |= want;
}

double totalLogLikelihood(omxMatrix *fitMat)
{
	if (fitMat->rows != 1) {
		omxFitFunction *ff = fitMat->fitFunction;
		if (ff->units == FIT_UNITS_PROBABILITY) {
			// NOTE: Floating-point addition is not
			// associative. If we compute this in parallel
			// then we introduce non-determinancy.
			double sum = 0;
			for(int i = 0; i < fitMat->rows; i++) {
				sum += log(omxVectorElement(fitMat, i));
			}
			if (!Global->rowLikelihoodsWarning) {
				Rf_warning("%s does not evaluate to a 1x1 matrix. Fixing model by adding "
					   "mxAlgebra(-2*sum(log(%s)), 'm2ll'), mxFitFunctionAlgebra('m2ll')",
					   fitMat->name(), fitMat->name());
				Global->rowLikelihoodsWarning = true;
			}
			return sum * Global->llScale;
		} else {
			omxRaiseErrorf("%s of type %s returned %d values instead of 1, not sure how to proceed",
				       fitMat->name(), ff->fitType, fitMat->rows);
			return nan("unknown");
		}
	} else {
		return fitMat->data[0];
	}
}

static void numericalGradientApprox(omxFitFunction *ff, FitContext *fc, bool haveFreshFit)
{
  if (isErrorRaised()) return;

  double fitSave = fc->fit;
  const int numFree = fc->getNumFree();

  if (!fc->numericalGradTool) {
    // allow option customization TODO
    fc->numericalGradTool =
      std::unique_ptr< AutoTune<JacobianGadget> >(new AutoTune<JacobianGadget>("numericalGradTool"));
    fc->numericalGradTool->setWork(std::unique_ptr<JacobianGadget>(new JacobianGadget(numFree)));
    fc->numericalGradTool->setMaxThreads(fc->numOptimizerThreads());
  }
  auto &ngt = fc->numericalGradTool->work();

  if (ngt.needRefFit() && !haveFreshFit) {
    ComputeFit("gradient", ff->matrix, FF_COMPUTE_FIT, fc);
  }

  Eigen::ArrayXd ref(1);
  ref[0] = fc->fit;
  Eigen::Map< Eigen::RowVectorXd > gradOut(fc->gradZ.data(), fc->gradZ.size());

	(*fc->numericalGradTool)([&](double *myPars, int thrId, Eigen::Ref<Eigen::ArrayXd> result)->void{
			FitContext *fc2 = thrId >= 0? fc->childList[thrId] : fc;
			Eigen::Map< Eigen::VectorXd > Est(myPars, fc2->numParam);
			// Only 1 parameter is different so we could
			// update only that parameter instead of all
			// of them.
			fc2->setEstFromOptimizer(myPars);
			ComputeFit("gradient", fc2->lookupDuplicate(ff->matrix), FF_COMPUTE_FIT, fc2);
			double fit = fc2->fit;
			if (fc2->outsideFeasibleSet()) fit = nan("infeasible");
      result[0] = fit;
      }, ref, [&fc](){ return fc->getCurrentFree(); }, true, gradOut);

  robustifyInplace(fc->gradZ);

  fc->fit = fitSave;
}

void ComputeFit(const char *callerName, omxMatrix *fitMat, int want, FitContext *fc)
{
	fc->incrComputeCount();
	fc->skippedRows = 0;
	omxFitFunction *ff = fitMat->fitFunction;
	if (ff) {
		if (want & FF_COMPUTE_GRADIENT) fc->initGrad();
		omxFitFunctionComputeAuto(ff, want, fc);
	} else {
		if (want != FF_COMPUTE_FIT) mxThrow("Only fit is available");
		if (fc->ciobj) mxThrow("CIs cannot be computed for unitless algebra");
		omxRecompute(fitMat, fc);
	}
	if (ff) {
		if (want & FF_COMPUTE_FIT) {
			fc->fit = totalLogLikelihood(fitMat);
			if (std::isfinite(fc->fit)) {
				fc->resetIterationError();
			}
			Global->checkpointPostfit(callerName, fc, false);
			if (OMX_DEBUG) {
				mxLog("%s: completed evaluation, fit=%.12g skippedRows=%d",
				      fitMat->name(), fc->fit, fc->skippedRows);
			}
		}
		if (want & FF_COMPUTE_GRADIENT) {
      if (!Global->analyticGradients) fc->gradZ.setConstant(NA_REAL);
      if (Global->NPSOL_HACK==0 && !fc->gradZ.allFinite()) {
        numericalGradientApprox(ff, fc, want & FF_COMPUTE_FIT);
      }
    }
	}
}

static omxFitFunction *omxNewInternalFitFunction(omxState* os, const char *fitType,
						 omxExpectation *expect, omxMatrix *matrix, bool rowLik)
{
	omxFitFunction *obj = 0;

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if(strcmp(fitType, entry->name) == 0) {
			obj = entry->allocate();
			obj->fitType = entry->name;
			break;
		}
	}
	if (!obj) mxThrow("omxNewInternalFitFunction: cannot find '%s'", fitType);

	if (!matrix) {
		obj->matrix = omxInitMatrix(1, 1, TRUE, os);
		obj->matrix->hasMatrixNumber = TRUE;
		obj->matrix->matrixNumber = ~os->algebraList.size();
		os->algebraList.push_back(obj->matrix);
	} else {
		obj->matrix = matrix;
	}

	obj->matrix->fitFunction = obj;

	obj->expectation = expect;

	if (rowLik && expect && expect->data) {
		omxData *dat = expect->data;
		omxResizeMatrix(matrix, dat->nrows(), 1);
	} else {
		omxResizeMatrix(matrix, 1, 1);
		matrix->data[0] = NA_REAL;
	}

	return obj;
}

void omxFillMatrixFromMxFitFunction(omxMatrix* om, int matrixNumber, SEXP rObj)
{
	om->hasMatrixNumber = TRUE;
	om->matrixNumber = matrixNumber;

	ProtectedSEXP fitFunctionClass(STRING_ELT(Rf_getAttrib(rObj, R_ClassSymbol), 0));
	const char *fitType = CHAR(fitFunctionClass);

	omxExpectation *expect = NULL;
	ProtectedSEXP slotValue(R_do_slot(rObj, Rf_install("expectation")));
	if (Rf_length(slotValue) == 1) {
		int expNumber = Rf_asInteger(slotValue);
		if(expNumber != NA_INTEGER) {
			expect = omxExpectationFromIndex(expNumber, om->currentState);
		}
	}

	bool rowLik = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));

	omxFitFunction *ff =
		omxNewInternalFitFunction(om->currentState, fitType, expect, om, rowLik);
	ff->rObj = rObj;
}

omxFitFunction *omxChangeFitType(omxFitFunction *oo, const char *fitType)
{
	if (oo->initialized) {
		mxThrow("%s: cannot omxChangeFitType from %s to %s; already initialized",
			 oo->matrix->name(), oo->fitType, fitType);
	}

	for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxFitFunctionSymbolTable); fx++) {
		const omxFitFunctionTableEntry *entry = omxFitFunctionSymbolTable + fx;
		if (strEQ(fitType, entry->name)) {
			auto *newObj = entry->allocate();
			newObj->rObj = oo->rObj;
			newObj->expectation = oo->expectation;
			newObj->fitType = entry->name;
			newObj->matrix = oo->matrix;
			newObj->units = oo->units;
			oo->matrix = 0;
			newObj->matrix->fitFunction = newObj;
			delete oo;
			// Need to call initMorph again? Probably never have 2 levels
			// of specialization?
			newObj->init();
			return newObj;
		}
	}

	mxThrow("Cannot find fit type '%s'", fitType);
}

void omxCompleteFitFunction(omxMatrix *om)
{
	omxFitFunction *obj = om->fitFunction;
	if (obj->initialized) return;

	int depth = Global->mpi->getDepth();

	if (obj->expectation) {
		omxCompleteExpectation(obj->expectation);
	}

	obj = obj->initMorph();

	if (Global->mpi->getDepth() != depth) mxThrow("%s improperly used the R protect stack", om->name());

	obj->initialized = TRUE;
}

/* Helper functions */
omxMatrix* omxNewMatrixFromSlot(SEXP rObj, omxState* currentState, const char* slotName) {
	SEXP slotValue;
	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install(slotName)));
	omxMatrix* newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
	return newMatrix; // NULL when length(slot)==0
}

omxMatrix *omxNewMatrixFromSlotOrAnon(SEXP rObj, omxState* currentState, const char* slotName,
																			int rows, int cols)
{
	ProtectedSEXP slotValue(R_do_slot(rObj, Rf_install(slotName)));
	omxMatrix *newMatrix;
	if (Rf_length(slotValue) == 0) {
		newMatrix = omxInitMatrix(rows, cols, currentState);
	} else {
		newMatrix = omxMatrixLookupFromState1(slotValue, currentState);
		if (newMatrix->rows != rows || newMatrix->cols != cols) {
			mxThrow("Matrix '%s' must be dimension %dx%d instead of %dx%d",
							slotName, rows, cols, newMatrix->rows, newMatrix->cols);
		}
	}
	return newMatrix;
}

void omxFitFunction::traverse(std::function<void(omxMatrix*)> &fn)
{
	fn(matrix);
}
