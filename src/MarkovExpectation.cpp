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

#include "omxExpectation.h"
#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>
#include <RcppEigenWrap.h>
#include "EnableWarnings.h"

class MarkovExpectation : public omxExpectation {
	typedef omxExpectation super;
public:
	enum ScaleType { SCALE_SOFTMAX, SCALE_SUM, SCALE_NONE };

	std::vector< omxExpectation* > components;
	omxMatrix *initial;
	omxMatrix *transition;
	unsigned initialV;
	unsigned transitionV;
	int verbose;
	ScaleType scale;
	omxMatrix *scaledInitial;
	omxMatrix *scaledTransition;
	const bool isMixtureInterface;

	MarkovExpectation(omxState *st, int num, bool _isMixtureInterface)
		: super(st, num), initialV(0), transitionV(0),
			isMixtureInterface(_isMixtureInterface) {};
	virtual ~MarkovExpectation();
	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual omxMatrix *getComponent(const char*);
	virtual void populateAttr(SEXP expectation);
};

omxExpectation *InitHiddenMarkovExpectation(omxState *st, int num)
{ return new MarkovExpectation(st, num, false); }

omxExpectation *InitMixtureExpectation(omxState *st, int num)
{ return new MarkovExpectation(st, num, true); }

MarkovExpectation::~MarkovExpectation()
{
	omxFreeMatrix(scaledInitial);
	omxFreeMatrix(scaledTransition);
}

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

	if (isMixtureInterface) {
		initial = omxNewMatrixFromSlot(rObj, currentState, "weights");
		transition = 0;
	} else {
		initial = omxNewMatrixFromSlot(rObj, currentState, "initial");
		transition = omxNewMatrixFromSlot(rObj, currentState, "transition");
	}

	ProtectedSEXP Rscale(R_do_slot(rObj, Rf_install("scale")));
	auto scaleName = CHAR(STRING_ELT(Rscale, 0));
	if (strEQ(scaleName, "softmax")) {
		scale = SCALE_SOFTMAX;
	} else if (strEQ(scaleName, "sum")) {
		scale = SCALE_SUM;
	} else if (strEQ(scaleName, "none")) {
		scale = SCALE_NONE;
	} else {
		mxThrow("%s: unknown scale '%s'", name, scaleName);
	}

	scaledInitial = omxInitMatrix(1, 1, TRUE, currentState);
	scaledTransition = 0;
	if (transition) {
		scaledTransition = omxInitMatrix(1, 1, TRUE, currentState);
	}
}

void MarkovExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	if (fc) {
		for (auto c1 : components) {
			c1->compute(fc, what, how);
		}
	}

	omxRecompute(initial, fc);
	if (initialV != omxGetMatrixVersion(initial)) {
		omxCopyMatrix(scaledInitial, initial);
		EigenVectorAdaptor Ei(scaledInitial);
		if (scale == SCALE_SOFTMAX) Ei.derived() = Ei.array().exp();
		if (scale != SCALE_NONE) {
			Ei /= Ei.sum();
		}
		if (verbose >= 2) mxPrintMat("initial", Ei);
		initialV = omxGetMatrixVersion(initial);
	}

	if (transition) {
		omxRecompute(transition, fc);
		if (transitionV != omxGetMatrixVersion(transition)) {
			omxCopyMatrix(scaledTransition, transition);
			EigenArrayAdaptor Et(scaledTransition);
			if (scale == SCALE_SOFTMAX) Et.derived() = Et.array().exp();
			if (scale != SCALE_NONE) {
				Eigen::ArrayXd v = Et.colwise().sum();
				Et.rowwise() /= v.transpose();
			}
			if (verbose >= 2) mxPrintMat("transition", Et);
			transitionV = omxGetMatrixVersion(transition);
		}
	}
}

void MarkovExpectation::populateAttr(SEXP robj)
{
	compute(0, 0, 0); // needed? TODO

	MxRList out;

	EigenVectorAdaptor Ei(scaledInitial);
	const char *initialName = isMixtureInterface? "weights" : "initial";
	out.add(initialName, Rcpp::wrap(Ei));

	if (scaledTransition) {
		EigenMatrixAdaptor Et(scaledTransition);
		out.add("transition", Rcpp::wrap(Et));
	}

	Rf_setAttrib(robj, Rf_install("output"), out.asR());
}

omxMatrix *MarkovExpectation::getComponent(const char* component)
{
	omxMatrix *retval = 0;

	if (strEQ("initial", component)) {
		retval = scaledInitial;
	} else if (strEQ("transition", component)) {
		retval = scaledTransition;
	}
	return retval;
}
