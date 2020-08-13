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
	omxMatrix *thresholdsMat;
	double *discreteSpecPtr;
  bool discreteCheckCount;
  bool _connectedToData;
  const Eigen::Map<Eigen::MatrixXd> getDiscreteSpec()
  { const Eigen::Map<Eigen::MatrixXd> ds(discreteSpecPtr, 2, discreteMat->cols); return ds; }
	omxMatrix *discreteMat;
	std::vector< Eigen::VectorXd > discreteCache;
	std::vector< omxThresholdColumn > thresholds;  // size() == numDataColumns

	void loadThresholds();
protected:
  void setConnectedToData(bool _to);
  bool getConnectedToData() const { return _connectedToData; }

 public:
	int numDataColumns;
	SEXP rObj;
	const char *name;   // pointer to a static string, no need to allocate or free
	omxData* data;
	int numOrdinal;  // number of thresholds with matrix != 0
	/* Replication of some of the structures from Matrix */
	unsigned short isComplete;													// Whether or not this expectation has been initialize
	omxState* currentState;
	bool isClone() const;
	int expNum;  // index in omxState's vector

	// omxExpectation should not need to know about free variables.
	FreeVarGroup *freeVarGroup; // TODO remove

	bool canDuplicate;
	bool dynamicDataSource;

	omxExpectation(omxState *state, int num) :
		dataColumnsPtr(0), thresholdsMat(0),
		discreteSpecPtr(0), _connectedToData(false), discreteMat(0),
    numDataColumns(0), rObj(0), name(0), data(0), numOrdinal(0),
    isComplete(false), currentState(state),
		expNum(num), freeVarGroup(0), canDuplicate(false), dynamicDataSource(false) {};
	virtual ~omxExpectation() {};
	virtual void init() {};
  virtual void connectToData();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual void print();
	virtual void populateAttr(SEXP expectation) {};

  void populateNormalAttr(SEXP robj, MxRList &out);

	// getComponent & mutate probably take encapsulation a little too seriously.
	// The Fit function should probably just include the structure definition
	// for the expectation and access fields directly or through object methods.
	virtual omxMatrix *getComponent(const char*) { return 0; }
	virtual void mutate(const char*, omxMatrix*) {};
	virtual void invalidateCache();
	virtual void generateData(FitContext *fc, MxRList &out);
	virtual int numSummaryStats();
	virtual void asVector1(FitContext *fc, int row, Eigen::Ref<Eigen::VectorXd> out);
	template <typename T> void asVector(FitContext *fc, int row, Eigen::MatrixBase<T> &out) {
		asVector1(fc, row, out.derived());
	}
  Eigen::MatrixXd buildThresholdMatrix();

	virtual bool usesDataColumnNames() const { return true; }
	void loadDataColFromR();
	void loadThresholdFromR();
	bool loadDefVars(int row);
	void loadFakeDefVars();

	void saveDataColumnsInfo(SEXP vec) {
		numDataColumns = Rf_length(vec);
		dataColumnsPtr = INTEGER(vec);
	}

	typedef Eigen::Matrix<int, Eigen::Dynamic, 1> DataColumnIndexVector;
	virtual const Eigen::Map<DataColumnIndexVector> getDataColumns();
	virtual const std::vector<const char *> &getDataColumnNames() const;
	virtual void getExogenousPredictors(std::vector<int> &out) {};
	virtual std::vector< omxThresholdColumn > &getThresholdInfo();
	double getThreshold(int r, int c);
};

	void omxCompleteExpectation(omxExpectation *ox);

	void omxFreeExpectationArgs(omxExpectation* Expectation);					// Frees all args
omxExpectation* omxExpectationFromIndex(int expIndex, omxState* os);
	omxExpectation* omxNewIncompleteExpectation(SEXP mxobj, int expNum, omxState* os);


/* Expectation-specific implementations of matrix functions */
static inline void omxExpectationCompute(FitContext *fc, omxExpectation *ox,
																				 const char *what, const char *how)
{
	if (!ox) return;
	ox->compute(fc, what, how);
}

static inline void omxExpectationCompute(FitContext *fc, omxExpectation *ox, const char *what)
{ omxExpectationCompute(fc, ox, what, NULL); }

static inline void omxExpectationCompute(FitContext *fc, omxExpectation *ox)
{ omxExpectationCompute(fc, ox, NULL); }

	omxExpectation* omxDuplicateExpectation(const omxExpectation *src, omxState* newState);

	void omxExpectationPrint(omxExpectation *source, char* d);					// Pretty-print a (small-form) expectation

omxMatrix* omxGetExpectationComponent(omxExpectation *ox, const char* component);

void omxSetExpectationComponent(omxExpectation *ox, const char* component, omxMatrix *om);

omxExpectation *omxInitNormalExpectation(omxState *, int num);
omxExpectation *omxInitLISRELExpectation(omxState *, int num);
omxExpectation *omxInitStateSpaceExpectation(omxState *, int num);
omxExpectation *omxInitRAMExpectation(omxState *, int num);
omxExpectation *omxInitExpectationBA81(omxState *, int num);
omxExpectation *omxInitGREMLExpectation(omxState *, int num);
omxExpectation *InitHiddenMarkovExpectation(omxState *, int num);
omxExpectation *InitMixtureExpectation(omxState *, int num);

void complainAboutMissingMeans(omxExpectation *off);

// NOTE: omxThresholdColumn.dataColumn is ignored because
// we assume that summary data is permuted into the
// expectation order as part of the summarization process.
template <typename T>
void normalToStdVector(omxMatrix *cov, omxMatrix *mean, omxMatrix *slope, T Eth,
		       std::vector< omxThresholdColumn > &ti, Eigen::Ref<Eigen::VectorXd> out)
{
	// order of elements: (c.f. lav_model_wls, lavaan 0.6-2)
	// 1. thresholds + means (interleaved)
	// 2. slopes (if any, columnwise per exo)
	// 3. variances (continuous indicators only)
	// 4. covariances; not correlations (lower triangle)

	EigenMatrixAdaptor Ecov(cov);
	if (ti.size() == 0) {
		int dx = 0;
		if (mean) {
			EigenVectorAdaptor Emean(mean);
			for (int rx=0; rx < cov->cols; ++rx) {
				out[dx++] = Emean(rx);
			}
		}
		if (slope) {
			EigenMatrixAdaptor Eslope(slope);
			for (int cx=0; cx < Eslope.cols(); ++cx) {
				for (int rx=0; rx < Eslope.rows(); ++rx) {
					out[dx++] = Eslope(rx,cx);
				}
			}
		}
		for (int cx=0; cx < cov->cols; ++cx) {
			out[dx++] = Ecov(cx,cx);
		}
		for (int cx=0; cx < cov->cols-1; ++cx) {
			for (int rx=cx+1; rx < cov->rows; ++rx) {
				out[dx++] = Ecov(rx,cx);
			}
		}
		return;
	}
	if (!mean) mxThrow("ordinal indicators and no mean vector");

	EigenVectorAdaptor Emean(mean);
	Eigen::VectorXd sdTmp(1.0/Ecov.diagonal().array().sqrt());
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> sd(Emean.size());
	sd.setIdentity();

	int dx = 0;
  for (int cx=0; cx < int(ti.size()); ++cx) {
    auto &th = ti[cx];
		for (int t1=0; t1 < th.numThresholds; ++t1) {
			double sd1 = sdTmp[cx];
			out[dx++] = (Eth(t1, cx) - Emean[cx]) * sd1;
			sd.diagonal()[cx] = sd1;
		}
		if (!th.numThresholds) {
			out[dx++] = Emean[cx];
		}
	}

	if (slope) {
		EigenMatrixAdaptor Eslope(slope);
		for (int cx=0; cx < Eslope.cols(); ++cx) {
			for (int rx=0; rx < Eslope.rows(); ++rx) {
				out[dx++] = Eslope(rx,cx);
			}
		}
	}

	Eigen::MatrixXd stdCov(sd * Ecov * sd);

	for (int cx=0; cx < cov->cols; ++cx) {
		if (ti[cx].numThresholds) continue;
		out[dx++] = stdCov(cx,cx);
	}

	for (int cx=0; cx < cov->cols-1; ++cx) {
		for (int rx=cx+1; rx < cov->rows; ++rx) {
			out[dx++] = stdCov(rx,cx);
		}
	}
}

#endif /* _OMXEXPECTATION_H_ */
