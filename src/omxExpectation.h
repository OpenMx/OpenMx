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
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxMatrix.h"
#include "omxAlgebra.h"
#include "omxData.h"
#include "omxState.h"

/* Expectation structure itself */
class omxExpectation {					// An Expectation
	typedef omxExpectation base;
	int *dataColumnsPtr;
	std::vector< omxThresholdColumn > thresholds;

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
	int expNum;

	// omxExpectation should not need to know about free variables.
	FreeVarGroup *freeVarGroup; // TODO remove
	const char *name;

	bool canDuplicate;
	bool dynamicDataSource;

	omxExpectation() : dataColumnsPtr(0), numDataColumns(0), rObj(0), expType(0),
		data(0), thresholdsMat(0), numOrdinal(0), isComplete(false), currentState(0),
		expNum(0), freeVarGroup(0), name(0), canDuplicate(false), dynamicDataSource(false) {};
	virtual ~omxExpectation() {};
	virtual void init() {};
	virtual void compute(FitContext *fc, const char *what, const char *how) = 0;
	virtual void print();
	virtual void populateAttr(SEXP expectation) {};

	virtual bool hasRowWeights() { return false; };
	virtual int getNumRows() { return 0; }
	virtual double *getRowWeights() { return 0; }
	virtual void setRowWeights(double *) {};
	
	// getComponent & mutate probably take encapsulation a little too seriously.
	// The Fit function should probably just include the structure definition
	// for the expectation and access fields directly or through object methods.
	virtual omxMatrix *getComponent(const char*) { return 0; }
	virtual void mutate(const char*, omxMatrix*) {};

	void loadFromR();
	bool loadDefVars(int row);

	void saveDataColumnsInfo(SEXP vec) {
		numDataColumns = Rf_length(vec);
		dataColumnsPtr = INTEGER(vec);
	}

	typedef Eigen::Matrix<int, Eigen::Dynamic, 1> DataColumnType;
	virtual const Eigen::Map<DataColumnType> getDataColumns();
	virtual std::vector< omxThresholdColumn > &getThresholdInfo();

	void loadThresholds(int numCols, int *thresholdColumn, int *thresholdNumber);
};

omxExpectation *
omxNewInternalExpectation(const char *expType, omxState* os);

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

omxExpectation *omxInitNormalExpectation();
omxExpectation *omxInitLISRELExpectation();
omxExpectation *omxInitStateSpaceExpectation();
omxExpectation *omxInitRAMExpectation();
omxExpectation *omxInitExpectationBA81();
omxExpectation *omxInitGREMLExpectation();
omxExpectation *InitHiddenMarkovExpectation();

void complainAboutMissingMeans(omxExpectation *off);

#endif /* _OMXEXPECTATION_H_ */
