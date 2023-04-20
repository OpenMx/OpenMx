/*
 *  Copyright (C) 2016 Joshua N. Pritikin <jpritikin@pobox.com>
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

#include "omxDefines.h"
#include "omxExpectation.h"
#include "fitMultigroup.h"
#include "omxExportBackendState.h"
#include "Compute.h"

// http://openmx.ssri.psu.edu/issue/2013/01/multigroup-fit-function

struct FitMultigroup : omxFitFunction {
	std::vector< FreeVarGroup* > varGroups;
	std::vector< omxMatrix* > fits;
	int verbose;

	virtual void init() override;
	virtual void compute2(int ffcompute, FitContext *fc) override;
	virtual void addOutput(MxRList *out) override;
	virtual void traverse(std::function<void(omxMatrix*)> fn) override;
};

void FitMultigroup::compute2(int want, FitContext *fc)
{
	omxMatrix *fitMatrix = matrix;
	double fit = 0;
	double mac = 0;
  double totalWeight = 0;

	FitMultigroup *mg = (FitMultigroup*) this;

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
		if (f1->fitFunction) {
      totalWeight += f1->fitFunction->scale;
    } else {
      totalWeight += 1;
    }
  }

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
    double subScale = 1;
		if (f1->fitFunction) {
      subScale = f1->fitFunction->scale;
			f1->fitFunction->subCompute(want, fc);
			if (want & FF_COMPUTE_MAXABSCHANGE) {
				mac = std::max(fc->mac, mac);
			}
			if (want & FF_COMPUTE_PREOPTIMIZE) {
				if (units == FIT_UNITS_UNINITIALIZED) {
					units = f1->fitFunction->units;
        } else if (units == FIT_UNITS_ANY) {
          // don't check anything
				} else if (units != f1->fitFunction->units) {
					mxThrow("%s: cannot combine units %s and %s (from %s)",
						matrix->name(), fitUnitsToName(units),
						fitUnitsToName(f1->fitFunction->units), f1->name());
				}
			}
		} else {
			omxRecompute(f1, fc);
      if (want & FF_COMPUTE_GRADIENT) invalidateGradient(fc);
		}
		if (want & FF_COMPUTE_FIT) {
			if(f1->rows != 1 || f1->cols != 1) {
				omxRaiseErrorf("%s[%d]: %s of type %s does not evaluate to a 1x1 matrix",
					       fitMatrix->name(), (int)ex, f1->name(), f1->fitFunction->fitType);
			}
      if (units != FIT_UNITS_SQUARED_RESIDUAL &&
          units != FIT_UNITS_SQUARED_RESIDUAL_CHISQ) {
        // Need to restrict to WLS only because we expose the fit value
        // in the mxAlgebra API. Although combining groups using proportions
        // improves numerical precision, it's a hazard for backward compatibility.
        // For example, custom mxConstraints will not be scaled correctly if we enable
        // this for all models.
        // TODO: Make an option?
        totalWeight = 1;
      }
      double prop = subScale / totalWeight;
			fit += prop * f1->data[0];
			if (mg->verbose >= 1) {
        mxLog("%s: %s fit=%.6g prop=%f", fitMatrix->name(), f1->name(), f1->data[0], prop);
      }
		}
	}

	if (fc) fc->mac = mac;

	if (want & FF_COMPUTE_FIT) {
    scale = totalWeight;
		fitMatrix->data[0] = fit;
		if (mg->verbose >= 1) { mxLog("%s: fit=%.6g scale=%f", fitMatrix->name(), fit, totalWeight); }
	}
}

void FitMultigroup::addOutput(MxRList *out)
{
	FitMultigroup *mg = this;

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
		if (!f1->fitFunction) continue;
		omxPopulateFitFunction(f1, out);
	}
}

omxFitFunction *initFitMultigroup()
{ return new FitMultigroup; }

void FitMultigroup::init()
{
	auto *oo = this;
	FitMultigroup *mg =this;

	SEXP rObj = oo->rObj;
	if (!rObj) return;

	if (mg->fits.size()) return; // hack to prevent double initialization, remove TOOD

	oo->units = FIT_UNITS_UNINITIALIZED;
	oo->hessianAvailable = TRUE;
	oo->canDuplicate = true;

	omxState *os = oo->matrix->currentState;

	ProtectedSEXP Rverb(R_do_slot(rObj, Rf_install("verbose")));
	mg->verbose = Rf_asInteger(Rverb);

	ProtectedSEXP Rgroups(R_do_slot(rObj, Rf_install("groups")));
	int *fits = INTEGER(Rgroups);
	for(int gx = 0; gx < Rf_length(Rgroups); gx++) {
		if (isErrorRaised()) break;
		omxMatrix *mat;
		if (fits[gx] >= 0) {
			mat = os->algebraList[fits[gx]];
		} else {
			mxThrow("Can only add algebra and fitfunction");
		}
		if (mat == oo->matrix) mxThrow("Cannot add multigroup to itself");
		mg->fits.push_back(mat);
		if (mat->fitFunction) {
			omxCompleteFitFunction(mat);
			oo->hessianAvailable = (oo->hessianAvailable && mat->fitFunction->hessianAvailable);
		} else {
			oo->hessianAvailable = FALSE;
		}
	}
}

void FitMultigroup::traverse(std::function<void(omxMatrix*)> fn)
{
	fn(matrix);
	for (auto &f1 : fits) {
		fn(f1);
	}
}

/* TODO
void omxMultigroupAdd(omxFitFunction *oo, omxFitFunction *fit)
{
	if (oo->initFun != initFitMultigroup) mxThrow("%s is not the multigroup fit", oo->fitType);
	if (!oo->initialized) mxThrow("Fit not initialized", oo);

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	mg->fits.push_back(fit->matrix);
	//addFreeVarDependency(oo->matrix->currentState, oo->matrix, fit->matrix);
}
*/
