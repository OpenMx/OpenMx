/*
 *  Copyright 2007-2016 The OpenMx Project
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
 *
 */
 
#ifndef _OMXFIMLFITFUNCTION_H_
#define _OMXFIMLFITFUNCTION_H_

#include "omxFitFunction.h"
#include "omxSadmvnWrapper.h"

typedef struct omxFIMLRowOutput {  // Output object for each row of estimation.  Mirrors the Mx1 output vector
	double Minus2LL;		// Minus 2 Log Likelihood
	double Mahalanobis;		// Square root of Mahalanobis distance Q = (row - means) %*% solve(sigma) %*% (row - means)
	double Z;				// Z-score of Mahalnobis.  This is ((Q/n )^(1/3) - 1 + 2/(9n)) (9n/2)^(.5)
	double Nrows;			// Number of rows in the data set--Unclear why this is returned
	double Ncols;			// Number of columns in the (censored) data set
	int finalMissed;		// Whether or not likelihood was calculable in the final case
	int modelNumber;		// Not used
} omxFIMLRowOutput;

enum JointStrategy {
	JOINT_AUTO,
	JOINT_CONDCONT,
	JOINT_CONDORD,
	JOINT_OLD
};

struct omxFIMLFitFunction {
	omxFIMLFitFunction *parent;
	int rowwiseParallel;
	omxMatrix* cov;				// Covariance Matrix
	omxMatrix* means;			// Vector of means
	omxData* data;				// The data

	omxMatrix* rowLikelihoods;     // The row-by-row likelihoods
	int returnRowLikelihoods;   // Whether or not to return row-by-row likelihoods
	int populateRowDiagnostics; // Whether or not to populated the row-by-row likelihoods back to R
	omxContiguousData contiguous;		// TODO

	std::vector<bool> isOrdinal;
	int numOrdinal;
	int numContinuous;
	OrdinalLikelihood ol;

	bool inUse;

	enum JointStrategy jointStrat;

	// performance counters
	int expectationComputeCount;
	int conditionCount;
	int invertCount;
	int ordSetupCount;
	int ordDensityCount;
	int contDensityCount;

	// --- old stuff below

	/* Structures for JointFIMLFitFunction */
	omxMatrix* contRow;		    // Memory reserved for continuous data row
	omxMatrix* ordCov;	    	// Memory reserved for ordinal covariance matrix
	omxMatrix* ordMeans;		// Memory reserved for ordinal column means    
    omxMatrix* ordContCov;      // Memory reserved for ordinal/continuous covariance
	omxMatrix* halfCov;         // Memory reserved for computations
    omxMatrix* reduceCov;       // Memory reserved for computations
    
	/* Reserved memory for faster calculation */
	omxMatrix* smallRow;		// Memory reserved for operations on each data row
	omxMatrix* smallCov;		// Memory reserved for operations on covariance matrix
	omxMatrix* smallMeans;		// Memory reserved for operations on the means matrix

	omxMatrix* RCX;				// Memory reserved for computationxs
		
	std::vector<int> indexVector;
	std::vector<int> identicalDefs;
	std::vector<int> identicalMissingness;
	std::vector<int> identicalRows;
	std::vector<bool> rowCompareInfo;
};

omxRListElement* omxSetFinalReturnsFIMLFitFunction(omxFitFunction *oo, int *numReturns);

void omxPopulateFIMLFitFunction(omxFitFunction *oo, SEXP algebra);
void omxInitFIMLFitFunction(omxFitFunction* oo, SEXP rObj);

bool omxFIMLSingleIterationJoint(FitContext *fc, omxFitFunction *localobj,
				 omxMatrix* output,
				 int rowbegin, int rowcount);

class mvnByRow {
	omxFitFunction *localobj;
 public:
	omxFIMLFitFunction *ofo;
	omxFIMLFitFunction *shared_ofo;
	omxExpectation* expectation;
	omxData *data;
	OrdinalLikelihood &ol;
	std::vector<int> &indexVector;
	std::vector<int> &identicalRows;
	int row;
	int lastrow;
	bool firstRow;
	omxMatrix *thresholdsMat;
	std::vector< omxThresholdColumn > &thresholdCols;
	omxMatrix *cov;
	omxMatrix *means;
	FitContext *fc;
	const Eigen::Map<Eigen::VectorXi> dataColumns;
	omxMatrix *rowLikelihoods;
	bool returnRowLikelihoods;
	std::vector<bool> &isOrdinal;
	int numOrdinal;
	int numContinuous;
	EigenVectorAdaptor jointMeans;
	EigenMatrixAdaptor jointCov;
	omxFIMLFitFunction *ofiml;

	int rowOrdinal;
	int rowContinuous;
	Eigen::VectorXd cDataBuf;
	Eigen::VectorXi iDataBuf;
	Eigen::VectorXi ordColBuf;
	std::vector<bool> isMissing;

	mvnByRow(FitContext *_fc, omxFitFunction *_localobj,
		 omxFIMLFitFunction *_ofiml, int rowbegin, int rowcount)
	:
	ofo((omxFIMLFitFunction*) _localobj->argStruct),
		shared_ofo(ofo->parent? ofo->parent : ofo),
		expectation(_localobj->expectation),
		ol(ofo->ol),
		indexVector(shared_ofo->indexVector),
		identicalRows(shared_ofo->identicalRows),
		thresholdCols(expectation->thresholds),
		dataColumns(expectation->dataColumnsPtr, expectation->numDataColumns),
		isOrdinal(_ofiml->isOrdinal),
		jointMeans(ofo->means),
		jointCov(ofo->cov)
	{
		data = ofo->data;
		ol.attach(dataColumns, data, expectation->thresholdsMat, expectation->thresholds);
		row = rowbegin;
		lastrow = rowbegin + rowcount;
		firstRow = true;
		thresholdsMat = expectation->thresholdsMat;
		cov = ofo->cov;
		means = ofo->means;
		fc = _fc;
		ofiml = _ofiml;
		rowLikelihoods = ofiml->rowLikelihoods;
		returnRowLikelihoods = ofiml->returnRowLikelihoods;
		localobj = _localobj;
		omxSetMatrixElement(localobj->matrix, 0, 0, 0.0);
		numOrdinal = ofiml->numOrdinal;
		numContinuous = ofiml->numContinuous;
		cDataBuf.resize(numContinuous);
		iDataBuf.resize(numOrdinal);
		ordColBuf.resize(numOrdinal);
		isMissing.resize(dataColumns.size());

		if (row > 0) {
			int prevIdentical = identicalRows[row - 1];
			row += (prevIdentical - 1);
		}
	};

	void loadRow(int r1)
	{
		rowOrdinal = 0;
		rowContinuous = 0;
		for(int j = 0; j < dataColumns.size(); j++) {
			int var = dataColumns[j];
			if (isOrdinal[j]) {
				int value = omxIntDataElement(data, r1, var);
				isMissing[j] = value == NA_INTEGER;
				if (!isMissing[j]) {
					ordColBuf[rowOrdinal] = j;
					iDataBuf[rowOrdinal++] = value;
				}
			} else {
				double value = omxDoubleDataElement(data, r1, var);
				isMissing[j] = std::isnan(value);
				if (!isMissing[j]) cDataBuf[rowContinuous++] = value;
			}
		}
	}

	void record(double lik)
	{
		if (returnRowLikelihoods) Rf_error("oops");

		EigenVectorAdaptor rl(localobj->matrix);
		rl[0] += log(lik);
	}

	void recordRow(double rowLik)
	{
		int numIdentical = identicalRows[row];
		if (returnRowLikelihoods) {
			EigenVectorAdaptor rl(rowLikelihoods);
			for(int nid = 0; nid < numIdentical; nid++) {
				int to = indexVector[row+nid];
				rl[to] = rowLik;
			}
		} else {
			EigenVectorAdaptor rl(localobj->matrix);
			rl[0] += numIdentical * log(rowLik);
		}
		row += numIdentical;
		firstRow = false;
	}
};

struct condContByRow : mvnByRow {
	typedef mvnByRow super;
	condContByRow(FitContext *_fc, omxFitFunction *_localobj,
		     omxFIMLFitFunction *_ofiml, int rowbegin, int rowcount)
		: super(_fc, _localobj, _ofiml, rowbegin, rowcount) {};
	bool eval();
};

struct oldByRow : mvnByRow {
	typedef mvnByRow super;
	oldByRow(FitContext *_fc, omxFitFunction *_localobj,
		     omxFIMLFitFunction *_ofiml, int rowbegin, int rowcount)
		: super(_fc, _localobj, _ofiml, rowbegin, rowcount) {};
	bool eval();
};

struct condOrdByRow : mvnByRow {
	typedef mvnByRow super;
	condOrdByRow(FitContext *_fc, omxFitFunction *_localobj,
		   omxFIMLFitFunction *_ofiml, int rowbegin, int rowcount)
		: super(_fc, _localobj, _ofiml, rowbegin, rowcount) {};
	bool eval();
};

#endif /* _OMXFIMLFITFUNCTION_H_ */
