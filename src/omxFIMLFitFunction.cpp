/*
 *  Copyright 2007-2016 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxDefines.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"

#ifdef SHADOW_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#endif

template <typename T1, typename T2, typename T3, typename T4>
void upperRightCovariance(const Eigen::MatrixBase<T1> &gcov,
			  T2 filterTest, T3 includeTest,
			  Eigen::MatrixBase<T4> &cov)
{
	for (int gcx=0, cx=0; gcx < gcov.cols(); gcx++) {
		if (filterTest(gcx) || !includeTest(gcx)) continue;
		for (int grx=0, rx=0; grx < gcov.rows(); grx++) {
			if (filterTest(grx) || includeTest(grx)) continue;
			cov(rx,cx) = gcov(grx, gcx);
			rx += 1;
		}
		cx += 1;
	}
}

bool condOrdByRow::eval()
{
	using Eigen::Map;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Eigen::VectorXi;
	using Eigen::ArrayXi;

	SimpCholesky< MatrixXd >  covDecomp;
	double ordLik = 1.0;
	double contLik = 1.0;

	VectorXd contMean;
	MatrixXd contCov;
	VectorXd ordMean;
	MatrixXd ordCov;
	VectorXd xi;
	MatrixXd U11;
	MatrixXd invOrdCov;
		
	int ssx = parent->sufficientSets.size() + 1;
	if (useSufficientSets && parent->sufficientSets.size()) {
		auto &sufficientSets = parent->sufficientSets;
		sufficientSet ssRef;
		ssRef.start = row;
		ssx = std::lower_bound(sufficientSets.begin(), sufficientSets.end(), ssRef,
				       [](const sufficientSet &s1, const sufficientSet &s2)
				       {return s1.start < s2.start; }) - sufficientSets.begin();
		if (ssx < int(parent->sufficientSets.size())) {
			auto &cur = sufficientSets[ssx];
			if (row > cur.start && row <= cur.start + cur.length - 1) row = cur.start + cur.length;
			//mxLog("row %d ssx %d start %d len %d ", row, ssx, cur.start, cur.length);
		}
		if (ssx > 0) {
			auto &prev = sufficientSets[ssx-1];
			if (row <= prev.start + prev.length - 1) row = prev.start + prev.length;
		}
	}

	if (row >= lastrow) {
#pragma atomic
		ofiml->unnecessaryThr += 1;
	}

	while(row < lastrow) {
		if (!loadRow()) return true;
		Map< VectorXd > cData(cDataBuf.data(), rowContinuous);
		Map< VectorXi > iData(iDataBuf.data(), rowOrdinal);
		EigenVectorAdaptor jointMeans(ofo->means);
		EigenMatrixAdaptor jointCov(ofo->cov);

		if (rowOrdinal == 0) {
			ordLik = 1.0;
		} else {
			if (!parent->ordinalMissingSame[row] || firstRow) {
				op.wantOrdinal = true;
				subsetNormalDist(jointMeans, jointCov, op, rowOrdinal, ordMean, ordCov);
				INCR_COUNTER(ordSetup);
				ol.setCovariance(ordCov, fc);
				Map< ArrayXi > ordColumns(ordColBuf.data(), rowOrdinal);
				ol.setColumns(ordColumns);
				ol.setMean(ordMean);
			}
			if (!parent->ordinalSame[row] || firstRow) {
				ordLik = ol.likelihood(sortedRow);
				INCR_COUNTER(ordDensity);
			}

			if (rowContinuous) {
				if (!parent->ordinalSame[row] || firstRow) {
					std::vector< omxThresholdColumn > &colInfo = expectation->thresholds;
					EigenMatrixAdaptor tMat(thresholdsMat);
					VectorXd uThresh(rowOrdinal);
					VectorXd lThresh(rowOrdinal);
					for(int jj=0; jj < rowOrdinal; jj++) {
						int var = dataColumns[ ordColBuf[jj] ];
						if (OMX_DEBUG && !omxDataColumnIsFactor(data, var)) {
							Rf_error("Must be a factor");
						}
						int pick = omxIntDataElement(data, sortedRow, var) - 1;
						if (OMX_DEBUG && (pick < 0 || pick > colInfo[var].numThresholds)) {
							Rf_error("Out of range");
						}
						int tcol = colInfo[var].column;
						if (pick == 0) {
							lThresh[jj] = -std::numeric_limits<double>::infinity();
							uThresh[jj] = (tMat(pick, tcol) - ordMean[jj]);
						} else if (pick == colInfo[var].numThresholds) {
							lThresh[jj] = (tMat(pick-1, tcol) - ordMean[jj]);
							uThresh[jj] = std::numeric_limits<double>::infinity();
						} else {
							lThresh[jj] = (tMat(pick-1, tcol) - ordMean[jj]);
							uThresh[jj] = (tMat(pick, tcol) - ordMean[jj]);
						}
					}

					if (!_mtmvnorm(ordLik, ordCov, lThresh, uThresh, xi, U11)) return true;
					U11 = U11.selfadjointView<Eigen::Upper>();
				}
				if (!parent->ordinalMissingSame[row] || firstRow) {
					invOrdCov = ordCov;
					if (InvertSymmetricPosDef(invOrdCov, 'L')) Rf_error("Non-positive definite");
					invOrdCov = invOrdCov.selfadjointView<Eigen::Lower>();
				}
				if (!parent->missingSameOrdinalSame[row] || firstRow) {
					// Aitken (1934) "Note on Selection from a Multivariate Normal Population"
					// Or Johnson/Kotz (1972), p.70
					MatrixXd V22;  //cont
					op.wantOrdinal = false;
					subsetNormalDist(jointMeans, jointCov, op, rowContinuous, contMean, V22);
					MatrixXd V12(rowOrdinal, rowContinuous);
					upperRightCovariance(jointCov, [&](int xx){ return isMissing[xx]; },
							     [&](int xx){ return !isOrdinal[xx]; }, V12);
					INCR_COUNTER(conditionCov);
					contCov = (V22 - V12.transpose() * (invOrdCov -
									    invOrdCov.selfadjointView<Eigen::Lower>() * U11 *
									    invOrdCov.selfadjointView<Eigen::Lower>()) * V12);
					INCR_COUNTER(invert);
					covDecomp.compute(contCov);
					if (covDecomp.info() != Eigen::Success ||
					    !(covDecomp.vectorD().array() > 0.0).all()) return true;
					covDecomp.refreshInverse();
					INCR_COUNTER(conditionMean);
					contMean += xi.transpose() * invOrdCov.selfadjointView<Eigen::Lower>() * V12;
				}
			}
		}

		if (rowContinuous) {
			if (!rowOrdinal && (!parent->continuousMissingSame[row] || firstRow)) {
				op.wantOrdinal = false;
				subsetNormalDist(jointMeans, jointCov, op, rowContinuous, contMean, contCov);
				INCR_COUNTER(invert);
				covDecomp.compute(contCov);
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) return true;
				covDecomp.refreshInverse();
			}

			const MatrixXd &iV = covDecomp.getInverse();

			if (row < data->rows-1 && parent->missingSameOrdinalSame[row+1] &&
			    useSufficientSets && ssx < (int)parent->sufficientSets.size()) {
				INCR_COUNTER(contDensity);
				sufficientSet &ss = parent->sufficientSets[ssx++];
				if (ss.start != row) Rf_error("oops");
				Eigen::VectorXd resid = ss.dataMean - contMean;
				//mxPrintMat("dataCov", ss.dataCov);
				//mxPrintMat("contMean", contMean);
				//mxPrintMat("dataMean", ss.dataMean);
				//mxPrintMat("resid", resid);
				double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
				double tr1 = trace_prod(iV, ss.dataCov);
				double logDet = covDecomp.log_determinant();
				double cterm = M_LN_2PI * ss.dataMean.size();
				//mxLog("[%d] iqf %f tr1 %f logDet %f cterm %f", ssx, iqf, tr1, logDet, cterm);
				double ll = ss.length * (iqf + logDet + cterm) + (ss.length-1) * tr1;
				record(-0.5 * ll + ss.length * log(ordLik));
				contLik = 1.0;
				row += ss.length;
				continue;
			}

			INCR_COUNTER(contDensity);
			VectorXd resid = cData - contMean;
			double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
			double cterm = M_LN_2PI * resid.size();
			double logDet = covDecomp.log_determinant();
			//mxLog("[%d] cont %f %f %f", sortedRow, iqf, cterm, logDet);
			contLik = exp(-0.5 * (iqf + cterm + logDet));
		} else { contLik = 1.0; }

		recordRow(ordLik * contLik);
	}

	return false;
}

bool condContByRow::eval()
{
	using Eigen::ArrayXi;
	using Eigen::VectorXi;
	using Eigen::VectorXd;
	using Eigen::MatrixXd;
	using Eigen::Map;

	MatrixXd ordCov;
	MatrixXd contCov;
	VectorXd contMean;
	VectorXd ordMean;
	MatrixXd ordAdj;
	SimpCholesky< Eigen::MatrixXd >  covDecomp;
	double contLik = 1.0;
	double ordLik = 1.0;

	while(row < lastrow) {
		if (!loadRow()) return true;
		Map< VectorXd > cData(cDataBuf.data(), rowContinuous);
		Map< VectorXi > iData(iDataBuf.data(), rowOrdinal);
		EigenVectorAdaptor jointMeans(ofo->means);
		EigenMatrixAdaptor jointCov(ofo->cov);

		bool newOrdCov = false;
		if (rowOrdinal && rowContinuous) {
			if (!parent->continuousMissingSame[row] || firstRow) {
				INCR_COUNTER(invert);
				op.wantOrdinal = false;
				subsetCovariance(jointCov, op, rowContinuous, contCov);
				covDecomp.compute(contCov);
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) return true;
				covDecomp.refreshInverse();
			}
			if (!parent->missingSame[row] || firstRow) {
				MatrixXd V12(rowOrdinal, rowContinuous); // avoid allocation TODO
				upperRightCovariance(jointCov, [&](int xx){ return isMissing[xx]; },
						     [&](int xx){ return !isOrdinal[xx]; }, V12);
				const Eigen::MatrixXd &icontCov = covDecomp.getInverse();
				ordAdj = V12 * icontCov.selfadjointView<Eigen::Lower>();

				op.wantOrdinal = true;
				subsetCovariance(jointCov, op, rowOrdinal, ordCov);
				ordCov -= ordAdj * V12.transpose();
				newOrdCov = true;
				INCR_COUNTER(conditionCov);
			}
			if (!parent->missingSameContinuousSame[row] || firstRow) {
				op.wantOrdinal = false;
				subsetVector(jointMeans, op, rowContinuous, contMean);
				op.wantOrdinal = true;
				subsetVector(jointMeans, op, rowOrdinal, ordMean);
				ordMean += ordAdj * (cData - contMean);
				INCR_COUNTER(conditionMean);
			}
		} else if (rowOrdinal) {
			if (!parent->ordinalMissingSame[row] || firstRow) {
				op.wantOrdinal = true;
				subsetNormalDist(jointMeans, jointCov, op, rowOrdinal, ordMean, ordCov);
				newOrdCov = true;
			}
		} else if (rowContinuous) {
			if (!parent->continuousMissingSame[row] || firstRow) {
				op.wantOrdinal = false;
				subsetNormalDist(jointMeans, jointCov, op, rowContinuous, contMean, contCov);
				INCR_COUNTER(invert);
				covDecomp.compute(contCov);
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) return true;
				covDecomp.refreshInverse();
			}
		}

		if (rowContinuous) {
			if (!parent->continuousSame[row] || firstRow) {
				INCR_COUNTER(contDensity);
				Eigen::VectorXd resid = cData - contMean;
				const Eigen::MatrixXd &iV = covDecomp.getInverse();
				double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
				double cterm = M_LN_2PI * resid.size();
				double logDet = covDecomp.log_determinant();
				contLik = exp(-0.5 * (iqf + cterm + logDet));
			}
		} else {
			contLik = 1.0;
		}

		if (rowOrdinal) {
			if (newOrdCov) {
				INCR_COUNTER(ordSetup);
				ol.setCovariance(ordCov, fc);
			}
			Eigen::Map< Eigen::ArrayXi > ordColumns(ordColBuf.data(), rowOrdinal);
			ol.setColumns(ordColumns);
			ol.setMean(ordMean);

			INCR_COUNTER(ordDensity);
			ordLik = ol.likelihood(sortedRow);
			//mxLog("[%d] %.5g", sortedRow, log(ordLik));

			if (ordLik == 0.0) {
				if (fc) fc->recordIterationError("Improper value detected by integration routine "
								 "in data row %d: Most likely the maximum number of "
								 "ordinal variables (20) has been exceeded.  \n"
								 " Also check that expected covariance matrix is not "
								 "positive-definite", sortedRow);
				if(OMX_DEBUG) {mxLog("Improper input to sadmvn in row likelihood.  Skipping Row.");}
				return true;
			}
		} else {
			ordLik = 1.0;
		}

		recordRow(contLik * ordLik);
	}
	return false;
}

static void omxDestroyFIMLFitFunction(omxFitFunction *off)
{
	if(OMX_DEBUG) { mxLog("Destroying FIML fit function object."); }
	omxFIMLFitFunction *argStruct = (omxFIMLFitFunction*) (off->argStruct);

	omxFreeMatrix(argStruct->smallMeans);
	omxFreeMatrix(argStruct->ordMeans);
	omxFreeMatrix(argStruct->contRow);
	omxFreeMatrix(argStruct->ordCov);
	omxFreeMatrix(argStruct->ordContCov);
	omxFreeMatrix(argStruct->halfCov);
	omxFreeMatrix(argStruct->reduceCov);

	omxFreeMatrix(argStruct->smallRow);
	omxFreeMatrix(argStruct->smallCov);
	omxFreeMatrix(argStruct->RCX);
	omxFreeMatrix(argStruct->rowLikelihoods);
	delete argStruct;
}

static void omxPopulateFIMLAttributes(omxFitFunction *off, SEXP algebra)
{
	if(OMX_DEBUG) { mxLog("Populating FIML Attributes."); }
	omxFIMLFitFunction *argStruct = ((omxFIMLFitFunction*)off->argStruct);
	SEXP expCovExt, expMeanExt, rowLikelihoodsExt;
	omxMatrix *expCovInt, *expMeanInt;

	omxExpectationCompute(NULL, off->expectation, NULL);
	expCovInt = argStruct->cov;
	expMeanInt = argStruct->means;
	
	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL && expMeanInt->rows > 0  && expMeanInt->cols > 0) {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, 0, 0));		
	}

	Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
	
	if(argStruct->populateRowDiagnostics){
		omxMatrix *rowLikelihoodsInt = argStruct->rowLikelihoods;
		Rf_protect(rowLikelihoodsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		for(int row = 0; row < rowLikelihoodsInt->rows; row++)
			REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);
		Rf_setAttrib(algebra, Rf_install("likelihoods"), rowLikelihoodsExt);
	}

	if (OMX_DEBUG_FIML_STATS) {
		MxRList count;
		count.add("expectation", Rf_ScalarInteger(argStruct->expectationComputeCount));
		count.add("conditionMean", Rf_ScalarInteger(argStruct->conditionMeanCount));
		count.add("conditionCov", Rf_ScalarInteger(argStruct->conditionCovCount));
		count.add("invert", Rf_ScalarInteger(argStruct->invertCount));
		count.add("ordSetup", Rf_ScalarInteger(argStruct->ordSetupCount));
		count.add("ordDensity", Rf_ScalarInteger(argStruct->ordDensityCount));
		count.add("contDensity", Rf_ScalarInteger(argStruct->contDensityCount));
		Rf_setAttrib(algebra, Rf_install("stats"), count.asR());
	}
}

struct FIMLCompare {
	omxData *data;
	omxExpectation *ex;
	std::vector<bool> ordinal;
	bool ordinalFirst;

	FIMLCompare(omxExpectation *_ex, bool _ordinalFirst) {
		ex = _ex;
		ordinalFirst = _ordinalFirst;
		data = ex->data;

		auto dc = ex->getDataColumns();
		ordinal.resize(dc.size());
		for (int cx=0; cx < dc.size(); ++cx) {
			ordinal[cx] = omxDataColumnIsFactor(data, cx);
		}
	}

	bool compareDataPart(bool part, int la, int ra, bool &mismatch) const
	{
		mismatch = true;
		auto dc = ex->getDataColumns();
		for (int cx=0; cx < dc.size(); ++cx) {
			if (part ^ ordinalFirst ^ ordinal[cx]) continue;
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, la, col);
			// if rm is not matching then result is undefined
			if (lm) continue;
			double lv = omxDoubleDataElement(data, la, col);
			double rv = omxDoubleDataElement(data, ra, col);
			if (lv == rv) continue;
			return lv < rv;
		}

		mismatch = false;
		return false;
	}

	bool isAllMissingnessPart(bool part, int la) const
	{
		auto dc = ex->getDataColumns();
		for (int cx=0; cx < dc.size(); ++cx) {
			if (part ^ ordinalFirst ^ ordinal[cx]) continue;
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, la, col);
			if (!lm) return false;
		}
		return true;
	}

	bool compareMissingnessPart(bool part, int la, int ra, bool &mismatch) const
	{
		mismatch = true;
		auto dc = ex->getDataColumns();
		for (int cx=0; cx < dc.size(); ++cx) {
			if (part ^ ordinalFirst ^ ordinal[cx]) continue;
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, la, col);
			bool rm = omxDataElementMissing(data, ra, col);
			if (lm == rm) continue;
			return lm < rm;
		}

		mismatch = false;
		return false;
	}

	bool compareMissingnessAndDataPart(bool part, int la, int ra, bool &mismatch) const
	{
		int got = compareMissingnessPart(part, la, ra, mismatch);
		if (mismatch) return got;
		got = compareDataPart(part, la, ra, mismatch);
		if (mismatch) return got;
		return false;
	}

	bool compareData(int la, int ra, bool &mismatch) const
	{
		int got = compareDataPart(false, la, ra, mismatch);
		if (mismatch) return got;
		return compareDataPart(true, la, ra, mismatch);
	}

	bool compareMissingness(int la, int ra, bool &mismatch) const
	{
		int got = compareMissingnessPart(false, la, ra, mismatch);
		if (mismatch) return got;
		return compareMissingnessPart(true, la, ra, mismatch);
	}

	bool compareAllDefVars(int la, int ra, bool &mismatch) const
	{
		mismatch = true;

		for (size_t k=0; k < data->defVars.size(); ++k) {
			int col = data->defVars[k].column;
			double lv = omxDoubleDataElement(data, la, col);
			double rv = omxDoubleDataElement(data, ra, col);
			if (doubleEQ(lv, rv)) continue;
			return lv < rv;
		}

		mismatch = false;
		return false;
	}

	bool operator()(int la, int ra) const
	{
		bool mismatch;
		bool got = compareAllDefVars(la, ra, mismatch);
		if (mismatch) return got;
		if (!ordinalFirst) {
			got = compareMissingnessPart(false, la, ra, mismatch);
			if (mismatch) return got;
			got = compareMissingnessPart(true, la, ra, mismatch);
			if (mismatch) return got;
			got = compareDataPart(false,la, ra, mismatch);
			if (mismatch) return got;
			got = compareDataPart(true,la, ra, mismatch);
			if (mismatch) return got;
		} else {
			got = compareMissingnessPart(false, la, ra, mismatch);
			if (mismatch) return got;
			got = compareDataPart(false,la, ra, mismatch);
			if (mismatch) return got;
			got = compareMissingnessPart(true, la, ra, mismatch);
			if (mismatch) return got;
			got = compareDataPart(true,la, ra, mismatch);
			if (mismatch) return got;
		}
		return false;
	}
};

static void recordGap(int rx, int &prev, std::vector<int> &identical)
{
	int gap = rx - prev;
	for (int gx=0; gx < gap; ++gx) identical[prev + gx] = gap - gx;
	prev = rx;
}

static void loadSufficientSet(omxFitFunction *off, int from, sufficientSet &ss)
{
	omxExpectation *ex = off->expectation;
	omxFIMLFitFunction *ofiml = ((omxFIMLFitFunction*)off->argStruct);
	auto& indexVector = ofiml->indexVector;
	omxData *data = ofiml->data;
	std::vector<bool> &isOrdinal = ofiml->isOrdinal;

	auto dc = ex->getDataColumns();
	int perRow = 0;
	for (int cx=0; cx < dc.size(); ++cx) {
		if (isOrdinal[cx]) continue;
		int col = dc[cx];
		perRow += !omxDataElementMissing(data, indexVector[from], col);
	}
	if (perRow == 0) return;

	Eigen::VectorXd dvec(perRow * ss.length);
	for (int row=0; row < ss.length; ++row) {
		int sortedRow = indexVector[from + row];
		for (int cx=0,dx=0; cx < dc.size(); ++cx) {
			if (isOrdinal[cx]) continue;
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, sortedRow, col);
			if (lm) continue;
			if (dx >= perRow) Rf_error("oops");
			dvec[row * perRow + dx] = omxDoubleDataElement(data, sortedRow, col);
			dx += 1;
		}
	}

	computeMeanCov(dvec, perRow, ss.dataMean, ss.dataCov);
}

static void addSufficientSet(omxFitFunction *off, int from, int to)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);
	//mxLog("ss from %d to %d length %d", from, to, 1 + to - from);
	sufficientSet ss1;
	ss1.start = from;
	ss1.length = 1 + to - from;
	loadSufficientSet(off, from, ss1);
	ofiml->sufficientSets.push_back(ss1);
}

static void sortData(omxFitFunction *off)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);
	auto& indexVector = ofiml->indexVector;
	indexVector.clear();
	ofiml->sufficientSets.clear();
	omxData *data = ofiml->data;
	indexVector.reserve(data->rows);
	for (int rx=0; rx < data->rows; ++rx) indexVector.push_back(rx);
	ofiml->sameAsPrevious.assign(data->rows, false);

	FIMLCompare cmp(off->expectation, ofiml->jointStrat == JOINT_CONDORD);

	if (data->needSort) {
		if (OMX_DEBUG) mxLog("sort %s for %s", data->name, off->name());
		//if (ofiml->jointStrat == JOINT_OLD) cmp.old = true;
		//cmp.old=true;
		std::sort(indexVector.begin(), indexVector.end(), cmp);
		//data->omxPrintData("sorted", 1000, indexVector.data());
	}

	switch (ofiml->jointStrat) {
	case JOINT_OLD:{
		auto& identicalDefs = ofiml->identicalDefs;
		auto& identicalMissingness = ofiml->identicalMissingness;
		auto& identicalRows = ofiml->identicalRows;
		identicalDefs.resize(data->rows);
		identicalMissingness.resize(data->rows);
		identicalRows.resize(data->rows);

		int prevDefs=0;
		int prevMissingness=0;
		int prevRows=0;
		for (int rx=1; rx < data->rows; ++rx) {
			bool mismatch;
			cmp.compareAllDefVars(indexVector[prevDefs], indexVector[rx], mismatch);
			if (mismatch) recordGap(rx, prevDefs, identicalDefs);
			cmp.compareMissingness(indexVector[prevMissingness], indexVector[rx], mismatch);
			if (mismatch) {
				recordGap(rx, prevMissingness, identicalMissingness);
			}
			cmp.compareData(indexVector[prevRows], indexVector[rx], mismatch);
			if (mismatch) recordGap(rx, prevRows, identicalRows);
		}
		recordGap(data->rows - 1, prevDefs, identicalDefs);
		recordGap(data->rows - 1, prevMissingness, identicalMissingness);
		recordGap(data->rows - 1, prevRows, identicalRows);
		identicalRows[data->rows - 1] = 1;
		if (0) {
			Eigen::Map< Eigen::VectorXi > m1(identicalDefs.data(), data->rows);
			Eigen::Map< Eigen::VectorXi > m2(identicalMissingness.data(), data->rows);
			Eigen::Map< Eigen::VectorXi > m3(identicalRows.data(), data->rows);
			mxPrintMat("m1", m1);
			mxPrintMat("m2", m2);
			mxPrintMat("m3", m3);
		}
		break;
	};
	case JOINT_AUTO:
	case JOINT_CONDORD:
	case JOINT_CONDCONT:{
		cmp.ordinalFirst = true;
		ofiml->ordinalMissingSame.assign(data->rows, false);
		ofiml->continuousMissingSame.assign(data->rows, false);
		ofiml->missingSameOrdinalSame.assign(data->rows, false);
		ofiml->missingSameContinuousSame.assign(data->rows, false);
		ofiml->continuousSame.assign(data->rows, false);
		ofiml->missingSame.assign(data->rows, false);
		ofiml->ordinalSame.assign(data->rows, false);
		int prevSS = -1;
		for (int rx=1; rx < data->rows; ++rx) {
			bool m1;
			cmp.compareAllDefVars(indexVector[rx-1], indexVector[rx], m1);
			bool m2;
			cmp.compareMissingnessPart(false, indexVector[rx-1], indexVector[rx], m2);
			if (!m1 && !m2) ofiml->ordinalMissingSame[rx] = true;
			bool m3;
			cmp.compareMissingnessPart(true, indexVector[rx-1], indexVector[rx], m3);
			if (!m1 && !m3) ofiml->continuousMissingSame[rx] = true;
			bool m4;
			cmp.compareMissingness(indexVector[rx-1], indexVector[rx], m4);
			if (!m1 && !m4) ofiml->missingSame[rx] = true;
			bool m7;
			cmp.compareDataPart(false, indexVector[rx-1], indexVector[rx], m7);
			if (!m1 && !m2 && !m7) ofiml->ordinalSame[rx] = true;
			if (!m1 && !m4 && !m7) {
				ofiml->missingSameOrdinalSame[rx] = true;
				if (prevSS == -1 && !cmp.isAllMissingnessPart(true, indexVector[rx-1])) {
					prevSS = rx-1;
				}
			} else {
				if (prevSS != -1) {
					addSufficientSet(off, prevSS, rx-1);
					prevSS = -1;
				}
			}
			bool m5;
			cmp.compareDataPart(true, indexVector[rx-1], indexVector[rx], m5);
			if (!m1 && !m3 && !m5) ofiml->continuousSame[rx] = true;
			if (!m1 && !m4 && !m5) ofiml->missingSameContinuousSame[rx] = true;

			if (m1 || m4) continue;
			bool m6;
			cmp.compareData(indexVector[rx-1], indexVector[rx], m6);
			if (m6) continue;
			ofiml->sameAsPrevious[rx] = true;
		}
		if (prevSS != -1) addSufficientSet(off, prevSS, data->rows-1);
		break;}
	}
}

static bool dispatchByRow(FitContext *_fc, omxFitFunction *_localobj,
			  omxFIMLFitFunction *parent, omxFIMLFitFunction *ofiml, int rowbegin, int rowcount)
{
	switch (ofiml->jointStrat) {
	case JOINT_OLD:{
		oldByRow batch(_fc, _localobj, parent, ofiml, rowbegin, rowcount);
		return batch.eval();
	}
	case JOINT_CONDORD:{
		condOrdByRow batch(_fc, _localobj, parent, ofiml, rowbegin, rowcount);
		return batch.eval();
	}
	case JOINT_AUTO:
	case JOINT_CONDCONT:{
		condContByRow batch(_fc, _localobj, parent, ofiml, rowbegin, rowcount);
		return batch.eval();
	}
	default: Rf_error("oops");
	}
}

static void CallFIMLFitFunction(omxFitFunction *off, int want, FitContext *fc)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);

	if (want & FF_COMPUTE_INITIAL_FIT) return;
	if (want & FF_COMPUTE_PREOPTIMIZE) {
		ofiml->inUse = true;
		if (fc->isClone()) {
			omxMatrix *pfitMat = fc->getParentState()->getMatrixFromIndex(off->matrix);
			ofiml->parent = (omxFIMLFitFunction*) pfitMat->fitFunction->argStruct;
		} else {
			off->openmpUser = ofiml->rowwiseParallel != 0;
			if (!ofiml->indexVector.size()) sortData(off);
		}
		return;
	}
	if (want & FF_COMPUTE_FINAL_FIT && !ofiml->inUse) return;

	if(OMX_DEBUG) mxLog("%s: joint FIML; openmpUser=%d", off->name(), off->openmpUser);

	omxMatrix* fitMatrix  = off->matrix;
	int numChildren = fc? fc->childList.size() : 0;

	omxMatrix *means	= ofiml->means;
	omxData* data           = ofiml->data;                            //  read-only

	omxFIMLFitFunction *parent = ofiml->parent? ofiml->parent : ofiml;
	int returnRowLikelihoods = ofiml->returnRowLikelihoods;   //  read-only
	omxExpectation* expectation = off->expectation;
	if (!means) {
		if (want & FF_COMPUTE_FINAL_FIT) return;
		complainAboutMissingMeans(expectation);
		return;
	}

	bool failed = false;

	int parallelism = (numChildren == 0 || !off->openmpUser) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}
	parallelism -= ofiml->unnecessaryThr;

	//mxLog("par=%d", parallelism);
	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		if (OMX_DEBUG) {
			FitContext *kid = fc->childList[0];
			omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) childFit->argStruct);
			if (!ofo->parent) Rf_error("oops");
		}

#pragma omp parallel for num_threads(parallelism) reduction(||:failed)
		for(int i = 0; i < parallelism; i++) {
			int rowbegin = stride*i;
			int rowcount = stride;
			if (i == parallelism-1) rowcount = data->rows - rowbegin;
			FitContext *kid = fc->childList[i];
			omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			failed |= dispatchByRow(kid, childFit, parent, ofiml, rowbegin, rowcount);
		}
	} else {
		failed |= dispatchByRow(fc, off, parent, ofiml, 0, data->rows);
	}

	if (failed) {
		if(!returnRowLikelihoods) {
			omxSetMatrixElement(off->matrix, 0, 0, NA_REAL);
		} else {
			EigenArrayAdaptor got(off->matrix);
			got.setZero(); // not sure if NA_REAL is safe
		}
		return;
	}

	if(!returnRowLikelihoods) {
		if (parallelism > 1) {
			double sum = 0.0;
			for(int i = 0; i < parallelism; i++) {
				FitContext *kid = fc->childList[i];
				omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
				EigenVectorAdaptor got(childMatrix);
				sum += got[0];
			}
			sum *= -2.0;
			omxSetMatrixElement(off->matrix, 0, 0, sum);
		} else {
			EigenVectorAdaptor got(off->matrix);
			got[0] *= -2.0;
		}
		if(OMX_DEBUG) {mxLog("%s: total likelihood is %3.3f", off->name(), off->matrix->data[0]);}
	} else {
		omxCopyMatrix(off->matrix, ofiml->rowLikelihoods);
		if (OMX_DEBUG) {
			omxPrintMatrix(ofiml->rowLikelihoods, "row likelihoods");
		}
	}
}

void omxInitFIMLFitFunction(omxFitFunction* off)
{
	if(OMX_DEBUG) {
		mxLog("Initializing FIML fit function function.");
	}
	off->canDuplicate = TRUE;

	omxExpectation* expectation = off->expectation;
	if(expectation == NULL) {
		omxRaiseError("FIML cannot fit without model expectations.");
		return;
	}

	omxFIMLFitFunction *newObj = new omxFIMLFitFunction;
	newObj->unnecessaryThr = 0;
	newObj->inUse = false;
	newObj->parent = 0;
	newObj->expectationComputeCount = 0;
	newObj->conditionMeanCount = 0;
	newObj->conditionCovCount = 0;
	newObj->invertCount = 0;
	newObj->ordSetupCount = 0;
	newObj->ordDensityCount = 0;
	newObj->contDensityCount = 0;

	omxMatrix *cov = omxGetExpectationComponent(expectation, "cov");
	if(cov == NULL) { 
		omxRaiseError("No covariance expectation in FIML evaluation.");
		return;
	}

	omxMatrix *means = omxGetExpectationComponent(expectation, "means");
	
    newObj->cov = cov;
    newObj->means = means;
    newObj->smallMeans = NULL;
    newObj->ordMeans   = NULL;
    newObj->contRow    = NULL;
    newObj->ordCov     = NULL;
    newObj->ordContCov = NULL;
    newObj->halfCov    = NULL;
    newObj->reduceCov  = NULL;
    
    off->computeFun = CallFIMLFitFunction;
	
	off->destructFun = omxDestroyFIMLFitFunction;
	off->populateAttrFun = omxPopulateFIMLAttributes;

	if(OMX_DEBUG) {
		mxLog("Accessing data source.");
	}
	newObj->data = off->expectation->data;

	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood option.");
	}
	SEXP rObj = off->rObj;
	newObj->rowwiseParallel = Rf_asLogical(R_do_slot(rObj, Rf_install("rowwiseParallel")));

	const char *jointStrat = CHAR(Rf_asChar(R_do_slot(rObj, Rf_install("jointConditionOn"))));
	if (strEQ(jointStrat, "auto")) {
		newObj->jointStrat = JOINT_AUTO;
	} else if (strEQ(jointStrat, "ordinal")) {
		newObj->jointStrat = JOINT_CONDORD;
	} else if (strEQ(jointStrat, "continuous")) {
		newObj->jointStrat = JOINT_CONDCONT;
	} else if (strEQ(jointStrat, "old")) {
		newObj->jointStrat = JOINT_OLD;
	} else { Rf_error("jointConditionOn '%s'?", jointStrat); }

	newObj->returnRowLikelihoods = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));
	newObj->useSufficientSets = !newObj->returnRowLikelihoods;
	newObj->rowLikelihoods = omxInitMatrix(newObj->data->rows, 1, TRUE, off->matrix->currentState);
	
	
	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood population option.");
	}
	newObj->populateRowDiagnostics = Rf_asInteger(R_do_slot(rObj, Rf_install("rowDiagnostics")));


	if(OMX_DEBUG) {
		mxLog("Accessing variable mapping structure.");
	}
	auto dc = off->expectation->getDataColumns();

	if(OMX_DEBUG) {
		mxLog("Accessing Threshold matrix.");
	}

	std::vector<bool> &isOrdinal = newObj->isOrdinal;
	isOrdinal.resize(dc.size());

	int numOrdinal=0;
	int numContinuous=0;
	for(int j = 0; j < dc.size(); j++) {
		int var = dc[j];
		isOrdinal[j] = omxDataColumnIsFactor(newObj->data, var);
		if (isOrdinal[j]) numOrdinal += 1;
		else numContinuous += 1;
	}
	newObj->numOrdinal = numOrdinal;
	newObj->numContinuous = numContinuous;

	if (0 == numOrdinal && newObj->jointStrat == JOINT_AUTO) {
		newObj->jointStrat = JOINT_CONDORD;
	}

    /* Temporary storage for calculation */
    int covCols = newObj->cov->cols;
	if(OMX_DEBUG){mxLog("Number of columns found is %d", covCols);}
    // int ordCols = omxDataNumFactor(newObj->data);        // Unneeded, since we don't use it.
    // int contCols = omxDataNumNumeric(newObj->data);
    newObj->smallRow = omxInitMatrix(1, covCols, TRUE, off->matrix->currentState);
    newObj->smallCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
    newObj->RCX = omxInitMatrix(1, covCols, TRUE, off->matrix->currentState);
//  newObj->zeros = omxInitMatrix(1, newObj->cov->cols, TRUE, off->matrix->currentState);

    omxCopyMatrix(newObj->smallCov, newObj->cov);          // Will keep its aliased state from here on.
    if (means) {
	    newObj->smallMeans = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
	    omxCopyMatrix(newObj->smallMeans, newObj->means);
	    newObj->ordMeans = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
	    omxCopyMatrix(newObj->ordMeans, newObj->means);
    }
    newObj->contRow = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
    omxCopyMatrix(newObj->contRow, newObj->smallRow );
    newObj->ordCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
    omxCopyMatrix(newObj->ordCov, newObj->cov);

    off->argStruct = (void*)newObj;

    if(numOrdinal > 0) {
        if(OMX_DEBUG) mxLog("Ordinal Data detected");

        newObj->ordContCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->halfCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        omxCopyMatrix(newObj->ordContCov, newObj->cov);
    }
}
