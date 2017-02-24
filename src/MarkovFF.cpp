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

	struct state {
		//std::vector< FreeVarGroup* > varGroups;
		std::vector< omxMatrix* > components;
		int verbose;
		omxMatrix *initial;
		omxMatrix *transition;
	};

	static void dtor(omxFitFunction* oo)
	{
		state *mg = (state*) oo->argStruct;
		delete mg;
	}

	static void compute(omxFitFunction* oo, int want, FitContext *fc)
	{
		state *st = (state*) oo->argStruct;

		omxExpectation* expectation = oo->expectation;
		omxExpectationCompute(fc, expectation, NULL);

		auto components = st->components;
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

		EigenMatrixAdaptor Einitial(st->initial);
		Eigen::VectorXd expect = Einitial;
		expect /= expect.sum();
		if (st->verbose >= 3) mxPrintMat("expect", expect);

		if (st->transition) {
			EigenArrayAdaptor Etransition(st->transition);
			Eigen::ArrayXd v = Etransition.colwise().sum();
			Etransition.colwise() /= v;
			if (st->verbose >= 3) mxPrintMat("transition", Etransition);
		}

		Eigen::VectorXd tp(components.size());
		double lp=0;
		double rowp=1;
		for (int rx=0; rx < nrow; ++rx) {
			for (int cx=0; cx < int(components.size()); ++cx) {
				EigenVectorAdaptor Ecomp(components[cx]);
				tp[cx] = Ecomp[rx];
			}
			if (st->verbose >= 4) mxPrintMat("tp", tp);
			if (st->transition) {
				EigenMatrixAdaptor Etransition(st->transition);
				expect = (Etransition * expect).eval();
			}
			expect = tp.array() * expect.array();
			rowp = expect.sum();
			expect /= rowp;
			lp += log(rowp);
		}
		oo->matrix->data[0] = Global->llScale * lp;
		if (st->verbose >= 2) mxLog("%s: fit=%f", oo->name(), lp);
	}

};

void InitMarkovFF(omxFitFunction* oo)
{
	if (!oo->expectation) { Rf_error("%s requires an expectation", oo->fitType); }

	oo->units = FIT_UNITS_MINUS2LL;
	oo->destructFun = MarkovFF::dtor;
	oo->computeFun = MarkovFF::compute;
	oo->canDuplicate = true;

	if (oo->argStruct) Rf_error("double initialization");

	omxState *currentState = oo->matrix->currentState;
	omxExpectation *expectation = oo->expectation;
	const char *myex = "MxExpectationHiddenMarkov";
	if (!expectation || !strEQ(expectation->expType, myex))
		Rf_error("%s must be paired with an %s", oo->name(), myex);

	auto ms = new MarkovFF::state;
	oo->argStruct = ms;

	ProtectedSEXP Rverbose(R_do_slot(oo->rObj, Rf_install("verbose")));
	ms->verbose = Rf_asInteger(Rverbose);

	ProtectedSEXP Rcomponents(R_do_slot(oo->rObj, Rf_install("components")));
	int nc = Rf_length(Rcomponents);
	int *cvec = INTEGER(Rcomponents);
	for (int cx=0; cx < nc; ++cx) {
		omxMatrix *fmat = currentState->algebraList[ cvec[cx] ];
		auto ff = fmat->fitFunction;
		if (ff) {
			omxCompleteFitFunction(fmat);
			if (ff->units != FIT_UNITS_MINUS2LL) {
				omxRaiseErrorf("%s: component %s must be in likelihood units",
					       oo->name(), ff->name());
				return;
			}
		}
		ms->components.push_back(fmat);
	}

	ms->initial = expectation->getComponent("initial");
	ms->transition = expectation->getComponent("transition");
}
