 /*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
				mxThrow("%s: component '%s' has %d rows but component '%s' has %d rows",
					 oo->name(), components[0]->name(), nrow, c1->name(), c1->rows);
			}
		}

		Eigen::VectorXd expect;
		Eigen::VectorXd rowResult;
		int numC = components.size();
		Eigen::VectorXd tp(numC);
		double lp=0;
		for (int rx=0; rx < nrow; ++rx) {
			if (expectation->loadDefVars(rx) || rx == 0) {
				omxExpectationCompute(fc, expectation, NULL);
				if (!st->transition || rx == 0) {
					EigenVectorAdaptor Einitial(st->initial);
					expect = Einitial;
					if (expect.rows() != numC || expect.cols() != 1) {
						omxRaiseErrorf("%s: initial prob matrix must be %dx%d not %dx%d",
							       name(), numC, 1, expect.rows(), expect.cols());
						return;
					}
				}
				if (st->transition && (st->transition->rows != numC || st->transition->cols != numC)) {
					omxRaiseErrorf("%s: transition prob matrix must be %dx%d not %dx%d",
						       name(), numC, numC, st->transition->rows, st->transition->cols);
					return;
				}
			}
			for (int cx=0; cx < int(components.size()); ++cx) {
				EigenVectorAdaptor Ecomp(components[cx]);
				tp[cx] = Ecomp[rx];
			}
			if (st->verbose >= 4) {
				mxPrintMat("tp", tp);
			}
			if (st->transition) {
				EigenMatrixAdaptor Etransition(st->transition);
				expect = (Etransition * expect).eval();
			}
			rowResult = tp.array() * expect.array();
			double rowp = rowResult.sum();
			rowResult /= rowp;
			lp += log(rowp);
			if (st->transition) expect = rowResult;
		}
		oo->matrix->data[0] = Global->llScale * lp;
		if (st->verbose >= 2) mxLog("%s: fit=%f", oo->name(), lp);
	}

	// Does vector=TRUE mean something sensible? Mixture of mixtures?
	void state::init()
	{
		auto *oo = this;
		auto *ms = this;
		if (!oo->expectation) { mxThrow("%s requires an expectation", oo->fitType); }

		oo->units = FIT_UNITS_MINUS2LL;
		oo->canDuplicate = true;

		omxState *currentState = oo->matrix->currentState;
		const char *myex1 = "MxExpectationHiddenMarkov";
		const char *myex2 = "MxExpectationMixture";
		if (!expectation || !(strEQ(expectation->name, myex1) ||
				      strEQ(expectation->name, myex2)))
			mxThrow("%s must be paired with %s or %s", oo->name(), myex1, myex2);

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
}

omxFitFunction *InitMarkovFF()
{ return new MarkovFF::state; }

