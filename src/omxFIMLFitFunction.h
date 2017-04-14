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
 *
 */
 
#ifndef _OMXFIMLFITFUNCTION_H_
#define _OMXFIMLFITFUNCTION_H_

#include "omxFitFunction.h"
#include "omxSadmvnWrapper.h"
#include "Compute.h"

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

#if !_OPENMP && OMX_DEBUG
#define OMX_DEBUG_FIML_STATS 1
#else
#define OMX_DEBUG_FIML_STATS 0
#endif

#if OMX_DEBUG_FIML_STATS
#define INCR_COUNTER(x) ofiml->x##Count += 1
#else
#define INCR_COUNTER(x)
#endif

struct sufficientSet {
	int                  start;
	int                  length;
	int                  rows;   // accounting for row weights
	Eigen::MatrixXd      dataCov;
	Eigen::VectorXd      dataMean;
};

class omxFIMLFitFunction : public omxFitFunction {
	bool builtCache;
 public:
	const int ELAPSED_HISTORY_SIZE = 5;

	omxFIMLFitFunction *parent;
	int rowwiseParallel;
	omxMatrix* cov;				// Covariance Matrix
	omxMatrix* means;			// Vector of means
	omxData* data;				// The data

	omxMatrix* rowLikelihoods;     // The row-by-row likelihoods
	int returnRowLikelihoods;   // Whether or not to return row-by-row likelihoods
	int populateRowDiagnostics; // Whether or not to populated the row-by-row likelihoods back to R

	int skippedRows;
	int origStateId;
	int curParallelism;
	int rowBegin;
	int rowCount;
	int curElapsed;
	std::vector<nanotime_t> elapsed;
	nanotime_t getMedianElapsedTime();

	std::vector<bool> isOrdinal;
	int numOrdinal;
	int numContinuous;
	OrdinalLikelihood ol;

	int verbose;
	bool inUse;
	std::vector<int> indexVector;
	std::vector<bool> sameAsPrevious;
	std::vector<bool> continuousMissingSame;
	std::vector<bool> continuousSame;
	std::vector<bool> missingSameContinuousSame;
	std::vector<bool> missingSameOrdinalSame;
	std::vector<bool> ordinalMissingSame;
	std::vector<bool> missingSame;
	std::vector<bool> ordinalSame;
	bool useSufficientSets;
	std::vector<sufficientSet> sufficientSets;

	enum JointStrategy jointStrat;

	// performance counters
	int expectationComputeCount;
	int conditionMeanCount;
	int conditionCovCount;
	int invertCount;
	int ordSetupCount;
	int ordDensityCount;
	int contDensityCount;

	virtual ~omxFIMLFitFunction();
	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
	virtual void populateAttr(SEXP algebra);
	virtual void invalidateCache();

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
		
	std::vector<int> identicalDefs;
	std::vector<int> identicalMissingness;
	std::vector<int> identicalRows;
};

void omxPopulateFIMLFitFunction(omxFitFunction *oo, SEXP algebra);
void omxInitFIMLFitFunction(omxFitFunction* oo, SEXP rObj);

class mvnByRow {
	omxFitFunction *localobj;
 public:
	omxFIMLFitFunction *ofo;
	omxFIMLFitFunction *shared_ofo;
	omxExpectation* expectation;
	omxData *data;
	OrdinalLikelihood &ol;
	std::vector<int> &indexVector;
	std::vector<bool> &sameAsPrevious;
	int row;
	int lastrow;
	bool firstRow;
	omxMatrix *thresholdsMat;
	const std::vector< omxThresholdColumn > &thresholdCols;
	omxMatrix *cov;
	omxMatrix *means;
	FitContext *fc;
	const Eigen::Map<Eigen::VectorXi> dataColumns;
	omxMatrix *rowLikelihoods;
	bool returnRowLikelihoods;
	std::vector<bool> &isOrdinal;
	int numOrdinal;
	int numContinuous;
	omxFIMLFitFunction *ofiml;
	omxFIMLFitFunction *parent;
	int sortedRow;  // it's really the unsorted row (row in the original data); rename TODO
	bool useSufficientSets;
	double *rowWeight;

	int rowOrdinal;
	int rowContinuous;
	Eigen::VectorXd cDataBuf;
	Eigen::VectorXi iDataBuf;
	Eigen::VectorXi ordColBuf;
	std::vector<bool> isMissing;
	int verbose;
	nanotime_t startTime;

	struct subsetOp {
		std::vector<bool> &isOrdinal;
		std::vector<bool> &isMissing;
		bool wantOrdinal;
	subsetOp(std::vector<bool> &_isOrdinal,
		 std::vector<bool> &_isMissing) : isOrdinal(_isOrdinal), isMissing(_isMissing) {};
		// true to include
		bool operator()(int gx) { return !((wantOrdinal ^ isOrdinal[gx]) || isMissing[gx]); };
	} op;

	mvnByRow(FitContext *_fc, omxFitFunction *_localobj,
		 omxFIMLFitFunction *_parent, omxFIMLFitFunction *_ofiml)
	:
	ofo((omxFIMLFitFunction*) _localobj),
		shared_ofo(ofo->parent? ofo->parent : ofo),
		expectation(_localobj->expectation),
		ol(ofo->ol),
		indexVector(shared_ofo->indexVector),
		sameAsPrevious(shared_ofo->sameAsPrevious),
		thresholdCols(expectation->getThresholdInfo()),
		dataColumns(expectation->getDataColumns()),
		isOrdinal(_ofiml->isOrdinal),
		op(isOrdinal, isMissing)
	{
		data = ofo->data;
		ol.attach(dataColumns, data, expectation->thresholdsMat, expectation->getThresholdInfo());
		row = ofo->rowBegin;
		lastrow = ofo->rowBegin + ofo->rowCount;
		firstRow = true;
		thresholdsMat = expectation->thresholdsMat;
		cov = ofo->cov;
		means = ofo->means;
		fc = _fc;
		ofiml = _ofiml;
		parent = _parent;
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
		useSufficientSets = ofiml->useSufficientSets;
		rowWeight = data->getWeightColumn();
		verbose = ofiml->verbose;

		if (fc->isClone()) {  // rowwise parallel
			startTime = get_nanotime();
		}

		if (row > 0) {
			while (row < lastrow && sameAsPrevious[row]) row += 1;
		}
	};

	~mvnByRow() {
		if (fc->isClone()) {  // rowwise parallel
			double el1 = get_nanotime() - startTime;
			ofo->elapsed[shared_ofo->curElapsed] = el1;
			if (verbose >= 3) mxLog("%d--%d %.2fms", ofo->rowBegin, ofo->rowCount, el1/1000000.0);
		}
	};

	bool loadRow()
	{
		mxLogSetCurrentRow(row);
		sortedRow = indexVector[row];
		rowOrdinal = 0;
		rowContinuous = 0;
		for(int j = 0; j < dataColumns.size(); j++) {
			int var = dataColumns[j];
			if (isOrdinal[j]) {
				int value = omxIntDataElement(data, sortedRow, var);
				isMissing[j] = value == NA_INTEGER;
				if (!isMissing[j]) {
					ordColBuf[rowOrdinal] = j;
					iDataBuf[rowOrdinal++] = value;
				}
			} else {
				double value = omxDoubleDataElement(data, sortedRow, var);
				isMissing[j] = std::isnan(value);
				if (!isMissing[j]) cDataBuf[rowContinuous++] = value;
			}
			//mxLog("col %d datacol %d ordinal=%d missing=%d", j, var, int(isOrdinal[j]), int(isMissing[j]));
		}
		//mxLog("rowOrdinal %d rowContinuous %d", rowOrdinal, rowContinuous);

		bool numVarsFilled = expectation->loadDefVars(sortedRow);
		if (numVarsFilled || firstRow) {
			omxExpectationCompute(fc, expectation, NULL);
			INCR_COUNTER(expectationCompute);

			if (rowOrdinal) {
				omxRecompute(thresholdsMat, fc);
				for (int jx=0; jx < rowOrdinal; jx++) {
					int j = ordColBuf[jx];
					if (!thresholdsIncreasing(thresholdsMat, thresholdCols[j].column,
								  thresholdCols[j].numThresholds, fc)) return false;
				}
			}
		}
		return true;
	}

	void record(double logLik, int rows)
	{
		if (returnRowLikelihoods) Rf_error("oops");
		if (!std::isfinite(logLik)) {
			ofiml->skippedRows += rows;
		} else {
			EigenVectorAdaptor rl(localobj->matrix);
			//mxLog("%g += record(%g)", rl[0], lik);
			rl[0] += logLik;
		}
		firstRow = false;
		row += rows;
	}

	void skipRow()
	{
		int oldRow = row;
		if (returnRowLikelihoods) {
			EigenVectorAdaptor rl(rowLikelihoods);
			double rowLik = 0.0;
			rl[sortedRow] = rowLik;
			row += 1;
			while (row < data->rows && sameAsPrevious[row]) {
				int index = indexVector[row];
				rl[index] = rowLik;
				row += 1;
			}
		} else {
			row += 1;
			while (row < data->rows && sameAsPrevious[row]) {
				row += 1;
			}
		}
		ofiml->skippedRows += row - oldRow;
		firstRow = false;
	}

	void recordRow(double contLogLik, double ordLik)
	{
		if (ordLik == 0.0 || !std::isfinite(contLogLik)) {
			skipRow();
			return;
		}
		if (OMX_DEBUG_ROWS(sortedRow)) {
			mxLog("%d/%d ordLik %g contLogLik %g", row, sortedRow, ordLik, contLogLik);
		}
		if (returnRowLikelihoods) {
			EigenVectorAdaptor rl(rowLikelihoods);
			double rowLik = exp(contLogLik) * ordLik;
			if (rowWeight && rowWeight[sortedRow] != 1.0) rowLik = pow(rowLik, rowWeight[sortedRow]);
			rl[sortedRow] = rowLik;
			row += 1;
			while (row < data->rows && sameAsPrevious[row]) {
				rl[ indexVector[row] ] = rowLik;
				row += 1;
			}
		} else {
			EigenVectorAdaptor rl(localobj->matrix);
			double rowLogLik = contLogLik + log(ordLik);
			if (rowWeight) rowLogLik *= rowWeight[sortedRow];
			rl[0] += rowLogLik;
			row += 1;
			while (row < data->rows && sameAsPrevious[row]) {
				rl[0] += rowLogLik;
				row += 1;
			}
		}
		firstRow = false;
	}

	void recordRowOld(double rowLik) // TODO remove
	{
		auto &identicalRows = shared_ofo->identicalRows;
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

	void reportBadOrdLik(int loc)
	{
		if (fc) fc->recordIterationError("Ordinal covariance is not positive definite "
						 "in data '%s' row %d (loc%d)", data->name, 1+sortedRow, loc);
		if (verbose >= 1) ol.log();
	}

	template <typename T1, typename T2, typename T3>
	void reportBadContRow(const Eigen::MatrixBase<T1> &cdata, const Eigen::MatrixBase<T2> &resid,
			      const Eigen::MatrixBase<T3> &icov)
	{
		if (cdata.size() > 50) {
			if (fc) fc->recordIterationError("In data '%s' row %d continuous variables are too"
							 " far from the model implied distribution",
							 data->name, 1+sortedRow);
		} else {
			std::string empty = std::string("");
			std::string buf;
			buf += mxStringifyMatrix("data", cdata, empty);
			buf += mxStringifyMatrix("resid", resid, empty);
			buf += mxStringifyMatrix("covariance", icov, empty);
			if (fc) fc->recordIterationError("In data '%s' row %d continuous variables are too"
							 " far from the model implied distribution. Details:\n%s",
							 data->name, 1+sortedRow, buf.c_str());
		}
	}

	template <typename T1>
	void reportBadContLik(int loc, const Eigen::MatrixBase<T1> &badCov)
	{
		if (badCov.rows() > 50) {
			if (fc) fc->recordIterationError("The continuous part of the model implied covariance (loc%d) "
							 "is not positive definite in data '%s' row %d",
							 loc, data->name, 1+sortedRow);
		} else {
			std::string empty = std::string("");
			std::string buf = mxStringifyMatrix("covariance", badCov, empty);
			if (fc) fc->recordIterationError("The continuous part of the model implied covariance (loc%d) "
							 "is not positive definite in data '%s' row %d. Detail:\n%s",
							 loc, data->name, 1+sortedRow, buf.c_str());
		}
	}
};

struct condContByRow : mvnByRow {
	typedef mvnByRow super;
	condContByRow(FitContext *_fc, omxFitFunction *_localobj,
		      omxFIMLFitFunction *_parent, omxFIMLFitFunction *_ofiml)
		: super(_fc, _localobj, _parent, _ofiml) {};
	bool eval();
};

struct oldByRow : mvnByRow {
	typedef mvnByRow super;
	oldByRow(FitContext *_fc, omxFitFunction *_localobj,
		 omxFIMLFitFunction *_parent, omxFIMLFitFunction *_ofiml)
		: super(_fc, _localobj, _parent, _ofiml) {};
	bool eval();
};

struct condOrdByRow : mvnByRow {
	typedef mvnByRow super;
	condOrdByRow(FitContext *_fc, omxFitFunction *_localobj,
		     omxFIMLFitFunction *_parent, omxFIMLFitFunction *_ofiml)
		: super(_fc, _localobj, _parent, _ofiml) {};
	bool eval();
};

#endif /* _OMXFIMLFITFUNCTION_H_ */
