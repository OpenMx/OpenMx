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

/***********************************************************
* 
*  omxExpectation.h
*
*  Created: Timothy R. Brick 	Date: 2009-02-17
*
*	Contains header information for the omxExpectation class
*   omxExpectation objects hold necessary information to simplify
* 	Expectation function calculation.
*
**********************************************************/

#ifndef _OMXEXPECTATION_H_
#define _OMXEXPECTATION_H_

#include "omxDefines.h"
#include <R_ext/Rdynload.h> 

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxData.h"
#include "omxState.h"

/* Expectation structure itself */
class omxExpectation {					// An Expectation
	typedef omxExpectation base;
	int *dataColumnsPtr;
	std::vector<const char *> dataColumnNames;
	std::vector< omxThresholdColumn > thresholds;  // size() == numDataColumns

 public:
	int numDataColumns;
	SEXP rObj;
        const char* expType;   // pointer to a static string, no need to allocate or free
	omxData* data;
	omxMatrix *thresholdsMat;
	int numOrdinal;  // number of thresholds with matrix != 0
	/* Replication of some of the structures from Matrix */
	unsigned short isComplete;													// Whether or not this expectation has been initialize
	omxState* currentState;
	int expNum;  // index in omxState's vector

	// omxExpectation should not need to know about free variables.
	FreeVarGroup *freeVarGroup; // TODO remove
	const char *name;

	bool canDuplicate;
	bool dynamicDataSource;

	omxExpectation(omxState *state) :
		dataColumnsPtr(0), numDataColumns(0), rObj(0), expType(0),
		data(0), thresholdsMat(0), numOrdinal(0), isComplete(false), currentState(state),
		expNum(0), freeVarGroup(0), name(0), canDuplicate(false), dynamicDataSource(false) {};
	virtual ~omxExpectation() {};
	virtual void init() {};
	virtual void compute(FitContext *fc, const char *what, const char *how) = 0;
	virtual void print();
	virtual void populateAttr(SEXP expectation) {};

	// getComponent & mutate probably take encapsulation a little too seriously.
	// The Fit function should probably just include the structure definition
	// for the expectation and access fields directly or through object methods.
	virtual omxMatrix *getComponent(const char*) { return 0; }
	virtual void mutate(const char*, omxMatrix*) {};
	virtual void invalidateCache() {};
	virtual void generateData(FitContext *fc, MxRList &out);
	virtual int numSummaryStats();
	virtual void asVector1(FitContext *fc, int row, Eigen::Ref<Eigen::VectorXd> out);
	template <typename T> void asVector(FitContext *fc, int row, Eigen::MatrixBase<T> &out) {
		asVector1(fc, row, out.derived());
	}

	virtual bool usesDataColumnNames() const { return true; }
	void loadFromR();
	bool loadDefVars(int row);

	void saveDataColumnsInfo(SEXP vec) {
		numDataColumns = Rf_length(vec);
		dataColumnsPtr = INTEGER(vec);
	}

	typedef Eigen::Matrix<int, Eigen::Dynamic, 1> DataColumnIndexVector;
	virtual const Eigen::Map<DataColumnIndexVector> getDataColumns();
	virtual const std::vector<const char *> &getDataColumnNames() const;
	virtual void getExogenousPredictors(std::vector<int> &out) {};
	virtual std::vector< omxThresholdColumn > &getThresholdInfo();

	void loadThresholds();
};

	void omxCompleteExpectation(omxExpectation *ox);

	void omxFreeExpectationArgs(omxExpectation* Expectation);					// Frees all args
omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os);
	omxExpectation* omxNewIncompleteExpectation(SEXP mxobj, int expNum, omxState* os);
	

/* Expectation-specific implementations of matrix functions */
void omxExpectationRecompute(FitContext *fc, omxExpectation *ox);
void omxExpectationCompute(FitContext *fc, omxExpectation *ox, const char *what, const char *how);

static inline void omxExpectationCompute(FitContext *fc, omxExpectation *ox, const char *what)
{ omxExpectationCompute(fc, ox, what, NULL); }

	omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState);
	
	void omxExpectationPrint(omxExpectation *source, char* d);					// Pretty-print a (small-form) expectation
	
omxMatrix* omxGetExpectationComponent(omxExpectation *ox, const char* component);
	
void omxSetExpectationComponent(omxExpectation *ox, const char* component, omxMatrix *om);

omxExpectation *omxInitNormalExpectation(omxState *);
omxExpectation *omxInitLISRELExpectation(omxState *);
omxExpectation *omxInitStateSpaceExpectation(omxState *);
omxExpectation *omxInitRAMExpectation(omxState *);
omxExpectation *omxInitExpectationBA81(omxState *);
omxExpectation *omxInitGREMLExpectation(omxState *);
omxExpectation *InitHiddenMarkovExpectation(omxState *);
omxExpectation *InitMixtureExpectation(omxState *);
omxExpectation *povRAMExpectationInit(omxState *);

void complainAboutMissingMeans(omxExpectation *off);

void normalToStdVector(omxMatrix *cov, omxMatrix *mean, omxMatrix *slope, omxMatrix *thr,
		       int numOrdinal, std::vector< omxThresholdColumn > &ti,
		       Eigen::Ref<Eigen::VectorXd> out);

#endif /* _OMXEXPECTATION_H_ */
