 /*
 *  Copyright 2007-2017 The OpenMx Project
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

#include "omxFitFunction.h"
#include "EnableWarnings.h"

namespace MarkovFF {

	struct state : omxFitFunction {
		//std::vector< FreeVarGroup* > varGroups;
		std::vector< omxMatrix* > components;
		int verbose;
		omxMatrix *initial;
		omxMatrix *transition;
		FitStatisticUnits componentUnits;

		virtual void init();
		virtual void compute(int ffcompute, FitContext *fc);
	};

	void state::compute(int want, FitContext *fc)
	{
		state *st = (state*) this;
		auto *oo = this;

		omxExpectationCompute(fc, expectation, NULL);

		for (auto c1 : components) {
			if (c1->fitFunction) {
				omxFitFunctionCompute(c1->fitFunction, want, fc);
			} else {
				omxRecompute(c1, fc);
			}
		}
		if (!(want & FF_COMPUTE_FIT)) return;

		int nrow = components[0]->rows;
		for (auto c1 : components) {
			if (c1->rows != nrow) {
				Rf_error("%s: component '%s' has %d rows but component '%s' has %d rows",
					 oo->name(), components[0]->name(), nrow, c1->name(), c1->rows);
			}
		}

		EigenVectorAdaptor Einitial(st->initial);
		Eigen::VectorXd expect(Einitial);
		Eigen::VectorXd tp(components.size());
		double lp=0;
		for (int rx=0; rx < nrow; ++rx) {
			for (int cx=0; cx < int(components.size()); ++cx) {
				EigenVectorAdaptor Ecomp(components[cx]);
				tp[cx] = Ecomp[rx];
			}
			if (st->verbose >= 4) mxPrintMat("tp", tp);
			if (st->transition) {
				if (expectation->loadDefVars(rx)) {
					omxExpectationCompute(fc, expectation, NULL);
				}
				EigenMatrixAdaptor Etransition(st->transition);
				expect = (Etransition * expect).eval();
			}
			expect = tp.array() * expect.array();
			double rowp = expect.sum();
			expect /= rowp;
			lp += log(rowp);
		}
		oo->matrix->data[0] = Global->llScale * lp;
		if (st->verbose >= 2) mxLog("%s: fit=%f", oo->name(), lp);
	}

	// Does vector=TRUE mean something sensible? Mixture of mixtures?
	void state::init()
	{
		auto *oo = this;
		auto *ms = this;
		if (!oo->expectation) { Rf_error("%s requires an expectation", oo->fitType); }

		oo->units = FIT_UNITS_MINUS2LL;
		oo->canDuplicate = true;

		omxState *currentState = oo->matrix->currentState;
		const char *myex = "MxExpectationHiddenMarkov";
		if (!expectation || !strEQ(expectation->expType, myex))
			Rf_error("%s must be paired with an %s", oo->name(), myex);

		ProtectedSEXP Rverbose(R_do_slot(oo->rObj, Rf_install("verbose")));
		ms->verbose = Rf_asInteger(Rverbose);

		ProtectedSEXP Rcomponents(R_do_slot(oo->rObj, Rf_install("components")));
		int nc = Rf_length(Rcomponents);
		int *cvec = INTEGER(Rcomponents);
		componentUnits = FIT_UNITS_UNINITIALIZED;
		for (int cx=0; cx < nc; ++cx) {
			omxMatrix *fmat = currentState->algebraList[ cvec[cx] ];
			if (fmat->fitFunction) {
				omxCompleteFitFunction(fmat);
				auto ff = fmat->fitFunction;
				if (ff->units != FIT_UNITS_PROBABILITY) {
					omxRaiseErrorf("%s: component %s must be in probability units",
						       oo->name(), ff->name());
					return;
				}
				if (componentUnits == FIT_UNITS_UNINITIALIZED) {
					componentUnits = ff->units;
				} else if (ff->units != componentUnits) {
					omxRaiseErrorf("%s: components with heterogenous units %s and %s in same mixture",
						       oo->name(), fitUnitsToName(ff->units), fitUnitsToName(componentUnits));
				}
			}
			ms->components.push_back(fmat);
		}
		if (componentUnits == FIT_UNITS_UNINITIALIZED) componentUnits = FIT_UNITS_PROBABILITY;

		ms->initial = expectation->getComponent("initial");
		ms->transition = expectation->getComponent("transition");
	}
};

omxFitFunction *InitMarkovFF()
{ return new MarkovFF::state; }

