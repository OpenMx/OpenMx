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
};

static void mgDestroy(omxFitFunction* oo)
{
	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;
	delete mg;
}

static void mgCompute(omxFitFunction* oo, int ffcompute, FitContext *fc)
{
	omxMatrix *fitMatrix  = oo->matrix;
	fitMatrix->data[0] = 0;
	double mac = 0;

	FitMultigroup *mg = (FitMultigroup*) oo->argStruct;

	for (size_t ex=0; ex < mg->fits.size(); ex++) {
		omxMatrix* f1 = mg->fits[ex];
		if (f1->fitFunction) {
			omxFitFunctionCompute(f1->fitFunction, ffcompute, fc);
			if (ffcompute & FF_COMPUTE_MAXABSCHANGE) {
				mac = std::max(fc->mac, mac);
			}
			if (OMX_DEBUG) { mxLog("mg fit %s %d", f1->name, ffcompute); }
		} else {
			omxRecompute(f1);

			// This should really be checked elsewhere. TODO
			if(f1->rows != 1 || f1->cols != 1) {
				Rf_error("%s algebra %d does not evaluate to a 1x1 matrix", oo->fitType, ex);
			}
		}
		fitMatrix->data[0] += f1->data[0];
	}
	if (fc) fc->mac = mac;
	if(OMX_DEBUG) { mxLog("Fit Function sum of %lu groups is %f.", mg->fits.size(), fitMatrix->data[0]); }
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

	oo->gradientAvailable = TRUE;
	oo->hessianAvailable = TRUE;
	oo->parametersHaveFlavor = TRUE;

	omxState *os = oo->matrix->currentState;

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("groups")));
	int *fits = INTEGER(slotValue);
	for(int gx = 0; gx < Rf_length(slotValue); gx++) {
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
			oo->parametersHaveFlavor = (oo->parametersHaveFlavor && mat->fitFunction->parametersHaveFlavor);
		} else {
			// TODO derivs for algebra
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
