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

#include "omxExpectation.h"
#include "EnableWarnings.h"

class MarkovExpectation : public omxExpectation {
public:
	std::vector< omxExpectation* > components;
	omxMatrix *initial;  // mixture at time=0
	omxMatrix *transition;
	int verbose;
	
	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual omxMatrix *getComponent(const char*);
};

omxExpectation *InitHiddenMarkovExpectation()
{ return new MarkovExpectation; }

void MarkovExpectation::init()
{
	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose);

	ProtectedSEXP Rcomponents(R_do_slot(rObj, Rf_install("components")));
	int *cvec = INTEGER(Rcomponents);
	int nc = Rf_length(Rcomponents);
	for (int cx=0; cx < nc; ++cx) {
		components.push_back(omxExpectationFromIndex(cvec[cx], currentState));
	}

	initial = omxNewMatrixFromSlot(rObj, currentState, "T0");
	transition = omxNewMatrixFromSlot(rObj, currentState, "G");

	if (data) Rf_warning("'%s' does not require data (ignored)", name);
}

void MarkovExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	for (auto c1 : components) {
		c1->compute(fc, what, how);
	}

	omxRecompute(initial, fc);
	if (transition)
		omxRecompute(transition, fc);
}

omxMatrix *MarkovExpectation::getComponent(const char* component)
{
	omxMatrix *retval = 0;

	if (strEQ("initial", component)) {
		retval = initial;
	} else if (strEQ("transition", component)) {
		retval = transition;
	}
	return retval;
}
