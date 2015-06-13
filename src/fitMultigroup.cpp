/*
  Copyright 2012-2014 Joshua Nathaniel Pritikin and contributors

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "omxDefines.h"
#include "omxExpectation.h"
#include "fitMultigroup.h"
#include "omxExportBackendState.h"

// http://openmx.psyc.virginia.edu/issue/2013/01/multigroup-fit-function

struct FitMultigroup {
	std::vector< FreeVarGroup* > varGroups;
	std::vector< omxMatrix* > fits;
	int verbose;
};

static void mgDestroy(omxFitFunction* oo)
{
	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	delete mg;
}

static void mgCompute(omxFitFunction* oo, int want, FitContext *fc)
{
	omxMatrix *fitMatrix  = oo->matrix;
	double fit = 0;
	double mac = 0;

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
		if (f1->fitFunction) {
			omxFitFunctionCompute(f1->fitFunction, want, fc);
			if (want & FF_COMPUTE_MAXABSCHANGE) {
				mac = std::max(fc->mac, mac);
			}
		} else {
			omxRecompute(f1, fc);
		}
		if (want & FF_COMPUTE_FIT) {
			if(f1->rows != 1 || f1->cols != 1) {
				omxRaiseErrorf("%s[%d]: %s of type %s does not evaluate to a 1x1 matrix",
					       fitMatrix->name(), (int)ex, f1->name(), f1->fitFunction->fitType);
			}
			fit += f1->data[0];
			if (mg->verbose >= 1) { mxLog("%s: %s fit=%f", fitMatrix->name(), f1->name(), f1->data[0]); }
		}
	}

	if (fc) fc->mac = mac;

	if (want & FF_COMPUTE_FIT) {
		fitMatrix->data[0] = fit;
		if (mg->verbose >= 1) { mxLog("%s: fit=%f", fitMatrix->name(), fit); }
	}
}

void mgSetFreeVarGroup(omxFitFunction *oo, FreeVarGroup *fvg)
{
	if (!oo->argStruct) initFitMultigroup(oo); // ugh TODO

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;

	if (!mg->fits.size()) {
		mg->varGroups.push_back(fvg);
	} else {
		for (size_t ex=0; ex < mg->fits.size(); ex++) {
			omxMatrix *f1 = mg->fits[ex];
			if (!f1->fitFunction) {  // simple algebra
				oo->freeVarGroup = fvg;
				continue;
			}
			setFreeVarGroup(f1->fitFunction, fvg);
			oo->freeVarGroup = f1->fitFunction->freeVarGroup;
		}
	}
}

void mgAddOutput(omxFitFunction* oo, MxRList *out)
{
	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
		if (!f1->fitFunction) continue;
		omxPopulateFitFunction(f1, out);
	}
}

void initFitMultigroup(omxFitFunction *oo)
{
	oo->expectation = NULL;  // don't care about this
	oo->computeFun = mgCompute;
	oo->destructFun = mgDestroy;
	oo->setVarGroup = mgSetFreeVarGroup;
	oo->addOutput = mgAddOutput;

	if (!oo->argStruct) oo->argStruct = new FitMultigroup;
	FitMultigroup *mg = (FitMultigroup *) oo->argStruct;

	SEXP rObj = oo->rObj;
	if (!rObj) return;

	if (mg->fits.size()) return; // hack to prevent double initialization, remove TOOD

	oo->units = FIT_UNITS_UNINITIALIZED;
	oo->gradientAvailable = TRUE;
	oo->hessianAvailable = TRUE;

	omxState *os = oo->matrix->currentState;

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	mg->verbose = Rf_asInteger(slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("groups")));
	int *fits = INTEGER(slotValue);
	for(int gx = 0; gx < Rf_length(slotValue); gx++) {
		if (isErrorRaised()) break;
		omxMatrix *mat;
		if (fits[gx] >= 0) {
			mat = os->algebraList[fits[gx]];
		} else {
			Rf_error("Can only add algebra and fitfunction");
		}
		if (mat == oo->matrix) Rf_error("Cannot add multigroup to itself");
		mg->fits.push_back(mat);
		if (mat->fitFunction) {
			for (size_t vg=0; vg < mg->varGroups.size(); ++vg) {
				setFreeVarGroup(mat->fitFunction, mg->varGroups[vg]);
				oo->freeVarGroup = mat->fitFunction->freeVarGroup;
			}
			omxCompleteFitFunction(mat);
			oo->gradientAvailable = (oo->gradientAvailable && mat->fitFunction->gradientAvailable);
			oo->hessianAvailable = (oo->hessianAvailable && mat->fitFunction->hessianAvailable);
			if (oo->units == FIT_UNITS_UNINITIALIZED) {
				oo->units = mat->fitFunction->units;
				oo->ciFun = mat->fitFunction->ciFun;
			} else if (oo->units != mat->fitFunction->units) {
				Rf_error("%s: cannot combine units %d and %d (from %s)",
					 oo->matrix->name(),
					 oo->units, mat->fitFunction->units, mat->name());
			}
		} else {
			oo->gradientAvailable = FALSE;
			oo->hessianAvailable = FALSE;
		}
	}
}

/* TODO
void omxMultigroupAdd(omxFitFunction *oo, omxFitFunction *fit)
{
	if (oo->initFun != initFitMultigroup) Rf_error("%s is not the multigroup fit", oo->fitType);
	if (!oo->initialized) Rf_error("Fit not initialized", oo);

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	mg->fits.push_back(fit->matrix);
	//addFreeVarDependency(oo->matrix->currentState, oo->matrix, fit->matrix);
}
*/
