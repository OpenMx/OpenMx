/*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
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

#include <limits>
#include "omxDefines.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include <Eigen/Cholesky>
#include "omxFIMLFitFunction.h"
#include "EnableWarnings.h"

nanotime_t omxFIMLFitFunction::getMedianElapsedTime()
{
	std::sort(elapsed.begin(), elapsed.end());
	return elapsed[elapsed.size() / 2];
}

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

bool condOrdByRow::eval() //<--This is what gets called when all manifest variables are continuous.
{
	using Eigen::Map;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;
	using Eigen::VectorXi;
	using Eigen::ArrayXi;

	SimpCholesky< MatrixXd >  covDecomp;
	double ordLik = 1.0;
	double contLogLik = 0.0;
	const double Scale = Global->llScale;

	VectorXd contMean;
	MatrixXd contCov;
	VectorXd ordMean;
	MatrixXd ordCov;
	VectorXd xi;
	MatrixXd U11;
	MatrixXd invOrdCov;
	int prevRowContinuous = 0;

	int ssx = parent->sufficientSets.size() + 1;
	if (parent->sufficientSets.size()) {
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

	while(row < lastrow) {
		loadRow();
		double iqf = NA_REAL;
		double residSize = NA_REAL;
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
				ordLik = ol.likelihood(fc, sortedRow);
				INCR_COUNTER(ordDensity);
			}

			if (rowContinuous) {
				bool firstContinuous = prevRowContinuous == 0;
				if (!parent->ordinalSame[row] || firstContinuous) {
					std::vector< omxThresholdColumn > &colInfo = expectation->getThresholdInfo();
					VectorXd uThresh(rowOrdinal);
					VectorXd lThresh(rowOrdinal);
					for(int jj=0; jj < rowOrdinal; jj++) {
						int col = ordColBuf[jj];
						int var = dataColumns[col];
						if (OMX_DEBUG && !omxDataColumnIsFactor(data, var)) {
							mxThrow("Must be a factor");
						}
						int pick = omxIntDataElement(data, sortedRow, var);
						if (OMX_DEBUG && (pick < 0 || pick > colInfo[col].numThresholds)) {
							mxThrow("Out of range");
						}
            if (!std::isfinite(ordMean[jj])) { reportBadOrdLik(3); return true; }
						if (pick == 0) {
							lThresh[jj] = -std::numeric_limits<double>::infinity();
							uThresh[jj] = (expectation->getThreshold(pick, col) - ordMean[jj]);
              if (!std::isfinite(uThresh[jj])) { reportBadOrdLik(4); return true; }
						} else if (pick == colInfo[col].numThresholds) {
							lThresh[jj] = (expectation->getThreshold(pick-1, col) - ordMean[jj]);
              if (!std::isfinite(lThresh[jj])) { reportBadOrdLik(5); return true; }
							uThresh[jj] = std::numeric_limits<double>::infinity();
						} else {
							lThresh[jj] = (expectation->getThreshold(pick-1, col) - ordMean[jj]);
              if (!std::isfinite(lThresh[jj])) { reportBadOrdLik(5); return true; }
							uThresh[jj] = (expectation->getThreshold(pick, col) - ordMean[jj]);
              if (!std::isfinite(uThresh[jj])) { reportBadOrdLik(4); return true; }
						}
					}
          if ((lThresh.array() >= uThresh.array()).any()) {
            reportBadOrdLik(6); return true;
          }
					if (ordLik == 0.0 ||
					    !u_mtmvnorm(fc, ordLik, ordCov, lThresh, uThresh, xi, U11)) {
						reportBadOrdLik(1);
						return true;
					}
					U11 = U11.selfadjointView<Eigen::Upper>();
				}
				if (!parent->ordinalMissingSame[row] || firstContinuous) {
					invOrdCov = ordCov;
					if (InvertSymmetricPosDef(invOrdCov, 'L')) {
						reportBadOrdLik(2);
						return true;
					}
					invOrdCov = invOrdCov.selfadjointView<Eigen::Lower>();
				}
				if (!parent->missingSameOrdinalSame[row] || firstContinuous) {
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
					    !(covDecomp.vectorD().array() > 0.0).all()) {
						reportBadContLik(1, contCov);
						return true;
					}
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
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) {
					reportBadContLik(2, contCov);
					return true;
				}
				covDecomp.refreshInverse();
			}

			const MatrixXd &iV = covDecomp.getInverse();

			if (row < rows-1 && parent->missingSameOrdinalSame[row+1] &&
			    ssx < (int)parent->sufficientSets.size()) {
				INCR_COUNTER(contDensity);
				sufficientSet &ss = parent->sufficientSets[ssx++];
				if (ss.start != row) OOPS;
				if (ordLik == 0.0) {
					record(-std::numeric_limits<double>::infinity(), ss.length);
					continue;
				}
				Eigen::VectorXd resid = ss.dataMean - contMean;
				residSize = ss.dataMean.size();
				//mxPrintMat("dataCov", ss.dataCov);
				//mxPrintMat("contMean", contMean);
				//mxPrintMat("dataMean", ss.dataMean);
				//mxPrintMat("resid", resid);
				iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
				double tr1 = trace_prod(iV, ss.dataCov);
				double logDet = 2.0 * covDecomp.log_determinant();
				double cterm = M_LN_2PI * residSize;
				//mxLog("[%d] iqf %f tr1 %f logDet %f cterm %f", ssx, iqf, tr1, logDet, cterm);
				double ll = ss.rows * (iqf + logDet + cterm) + (ss.rows-1) * tr1;
				record(-0.5 * ll + ss.rows * log(ordLik), ss.length);
				contLogLik = 0.0;
				if (want & FF_COMPUTE_GRADIENT){
					if(Global->analyticGradients && ofiml->expectation->canProvideSufficientDerivs){
						ofiml->expectation->provideSufficientDerivs(fc, ofiml->dSigma_dtheta, ofiml->dNu_dtheta);
						for(size_t px=0; px < ofiml->dSigma_dtheta.size(); px++){
							double ssDerivCurr=0; //<--Fit derivative for current parameter for current sufficient set
							Eigen::MatrixXd dSigma_dtheta_curr(ofiml->dSigma_dtheta[0].rows(),ofiml->dSigma_dtheta[0].cols());
							Eigen::VectorXd dNu_dtheta_curr(ofiml->dNu_dtheta[0].size());
							//Use `subsetNormalDist()` to filter dSigma_dtheta[px] & dNu_dtheta[px] for missingness...
							subsetNormalDist(ofiml->dNu_dtheta[px], ofiml->dSigma_dtheta[px], op, rowContinuous, dNu_dtheta_curr, dSigma_dtheta_curr);
							//Do analytic derivs here.
							fc->gradZ[px] += Scale * ssDerivCurr;
						}
					}}
				continue;
			}

			INCR_COUNTER(contDensity);
			VectorXd resid = cData - contMean;
			residSize = resid.size();
			iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
			double cterm = M_LN_2PI * residSize;
			double logDet = 2.0 * covDecomp.log_determinant();
			//mxLog("[%d] cont %f %f %f", sortedRow, iqf, cterm, logDet);
			contLogLik = -0.5 * (iqf + cterm + logDet);
			if (!std::isfinite(contLogLik)) reportBadContRow(cData, resid, contCov);
		} else { contLogLik = 0.0; }

		recordRow(contLogLik, ordLik, iqf, residSize);
		prevRowContinuous = rowContinuous;
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
	double contLogLik = 0.0;
	double ordLik = 1.0;
	bool ordConditioned = false;

	while(row < lastrow) {
		loadRow();
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
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) {
					reportBadContLik(3, contCov);
					return true;
				}
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
			ordConditioned = true;
		} else if (rowOrdinal) {
			if (!parent->ordinalMissingSame[row] || firstRow || ordConditioned) {
				op.wantOrdinal = true;
				subsetNormalDist(jointMeans, jointCov, op, rowOrdinal, ordMean, ordCov);
				newOrdCov = true;
				ordConditioned = false;
			}
		} else if (rowContinuous) {
			if (!parent->continuousMissingSame[row] || firstRow) {
				op.wantOrdinal = false;
				subsetNormalDist(jointMeans, jointCov, op, rowContinuous, contMean, contCov);
				INCR_COUNTER(invert);
				covDecomp.compute(contCov);
				if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) {
					reportBadContLik(4, contCov);
					return true;
				}
				covDecomp.refreshInverse();
			}
		}

		double iqf = NA_REAL;
		double residSize = NA_REAL;

		if (rowContinuous) {
			if (!parent->continuousSame[row] || firstRow) {
				INCR_COUNTER(contDensity);
				Eigen::VectorXd resid = cData - contMean;
				residSize = resid.size();
				const Eigen::MatrixXd &iV = covDecomp.getInverse();
				iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
				double cterm = M_LN_2PI * residSize;
				double logDet = 2.0 * covDecomp.log_determinant();
				contLogLik = -0.5 * (iqf + cterm + logDet);
				if (!std::isfinite(contLogLik)) {
					reportBadContRow(cData, resid, contCov);
				}
			}
		} else {
			contLogLik = 0.0;
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
			ordLik = ol.likelihood(fc, sortedRow);
			//mxLog("[%d] %.5g", sortedRow, log(ordLik));
		} else {
			ordLik = 1.0;
		}

		recordRow(contLogLik, ordLik, iqf, residSize);
	}

	return false;
}

omxFIMLFitFunction::~omxFIMLFitFunction()
{
	if(OMX_DEBUG) { mxLog("Destroying FIML fit function object."); }
	omxFIMLFitFunction *argStruct = this;

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
	omxFreeMatrix(argStruct->otherRowwiseValues);
}

void omxFIMLFitFunction::populateAttr(SEXP algebra)
{
	if(OMX_DEBUG) { mxLog("Populating FIML Attributes."); }
	auto *off = this;
	omxFIMLFitFunction *argStruct = this;
	SEXP expCovExt, expMeanExt, rowLikelihoodsExt, rowObsExt, rowDistExt;
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
		omxMatrix *otherRowwiseValuesInt = argStruct->otherRowwiseValues;
		Rf_protect(rowLikelihoodsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		Rf_protect(rowObsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		Rf_protect(rowDistExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		for(int row = 0; row < rowLikelihoodsInt->rows; row++) {
			REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);
			REAL(rowDistExt)[row] = omxMatrixElement(otherRowwiseValuesInt, row, 0);
			REAL(rowObsExt)[row] = omxMatrixElement(otherRowwiseValuesInt, row, 1);
		}
		Rf_setAttrib(algebra, Rf_install("likelihoods"), rowLikelihoodsExt);
		Rf_setAttrib(algebra, Rf_install("rowDist"), rowDistExt);
		Rf_setAttrib(algebra, Rf_install("rowObs"), rowObsExt);
	}

	const char *jointLabels[] = {
		"auto", "continuous", "ordinal", "old"
	};
	Rf_setAttrib(algebra, Rf_install("jointConditionOn"),
		     makeFactor(Rf_ScalarInteger(1+argStruct->jointStrat),
				OMX_STATIC_ARRAY_SIZE(jointLabels), jointLabels));

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

	FIMLCompare(omxExpectation *u_ex) {
		ex = u_ex;
		ordinalFirst = true;
		data = ex->data;

		auto dc = ex->getDataColumns();
		ordinal.resize(dc.size());
		for (int cx=0; cx < dc.size(); ++cx) {
			ordinal[cx] = omxDataColumnIsFactor(data, dc[cx]);
			//mxLog("%d is ordinal=%d", cx, int(ordinal[cx]));
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

static void loadSufficientSet(omxFitFunction *off, int from, sufficientSet &ss)
{
	omxExpectation *ex = off->expectation;
	omxFIMLFitFunction *ofiml = ((omxFIMLFitFunction*)off);
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
	ss.rows = 0;

	for (int row=0; row < ss.length; ++row) {
		int sortedRow = indexVector[from + row];
		ss.rows += 1;
		for (int cx=0,dx=0; cx < dc.size(); ++cx) {
			if (isOrdinal[cx]) continue;
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, sortedRow, col);
			if (lm) continue;
			if (dx >= perRow) OOPS;
			dvec[row * perRow + dx] = omxDoubleDataElement(data, sortedRow, col);
			dx += 1;
		}
	}

	computeMeanCov(dvec, perRow, ss.dataMean, ss.dataCov);
}

static void addSufficientSet(omxFitFunction *off, int from, int to)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off);
	if (!ofiml->useSufficientSets) return;
	if (ofiml->verbose >= 2) {
    mxLog("%s: addSufficientSet from %d to %d length %d", off->name(), from, to, 1 + to - from);
  }
	omxData *data = ofiml->data;
	double *rowWeight = data->getWeightColumn();
	if (rowWeight) return; // too complex, for now
	sufficientSet ss1;
	ss1.start = from;
	ss1.length = 1 + to - from;
	loadSufficientSet(off, from, ss1);
	ofiml->sufficientSets.push_back(ss1);
}

static void sortData(omxFitFunction *off)
{
	//mxLog("%s: sortData", off->matrix->name());
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off);
	auto &rowMult = ofiml->rowMult;
	auto& indexVector = ofiml->indexVector;
	ofiml->sufficientSets.clear();
	omxData *data = ofiml->data;
	data->recalcRowWeights(rowMult, indexVector);
	int rows = int(indexVector.size());
	ofiml->sameAsPrevious.assign(rows, false);

	FIMLCompare cmp(off->expectation);

	if (ofiml->jointStrat == JOINT_AUTO) {
		cmp.ordinalFirst = true;
		std::sort(indexVector.begin(), indexVector.end(), cmp);

		int numUnique = rows;
		for (int rx=1; rx < rows; ++rx) {
			bool m1;
			cmp.compareAllDefVars(indexVector[rx-1], indexVector[rx], m1);
			bool m2;
			cmp.compareMissingnessPart(false, indexVector[rx-1], indexVector[rx], m2);
			bool m7;
			cmp.compareDataPart(false, indexVector[rx-1], indexVector[rx], m7);
			if (!m1 && !m2 && !m7) --numUnique;
		}
		double rowsPerOrdinalPattern = rows / (double)numUnique;
		// magic numbers from logistic regression
		double prediction = -3.58 + 0.36 * rowsPerOrdinalPattern - 0.06 * ofiml->numContinuous;
		if (prediction > 0.0) {
			ofiml->jointStrat = JOINT_CONDORD;
		} else {
			ofiml->jointStrat = JOINT_CONDCONT;
		}
	}

	cmp.ordinalFirst = ofiml->jointStrat == JOINT_CONDORD;

	if (data->needSort) {
		if (ofiml->verbose >= 1) mxLog("%s: sort %s strategy %d for %s",
                                   ofiml->name(), data->name, ofiml->jointStrat, off->name());
		// Maybe already sorted by JOINT_AUTO, but not a big waste to resort
		std::sort(indexVector.begin(), indexVector.end(), cmp);
		//data->omxPrintData("sorted", 1000, indexVector.data());
	}

	cmp.ordinalFirst = true;
	ofiml->ordinalMissingSame.assign(rows, false);
	ofiml->continuousMissingSame.assign(rows, false);
	ofiml->missingSameOrdinalSame.assign(rows, false);
	ofiml->missingSameContinuousSame.assign(rows, false);
	ofiml->continuousSame.assign(rows, false);
	ofiml->missingSame.assign(rows, false);
	ofiml->ordinalSame.assign(rows, false);
	int prevSS = -1;
	for (int rx=1; rx < rows; ++rx) {
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
	if (prevSS != -1) addSufficientSet(off, prevSS, rows-1);
	if (ofiml->verbose >= 3) {
		mxLog("key: row ordinalMissingSame continuousMissingSame missingSameOrdinalSame missingSameContinuousSame continuousSame missingSame ordinalSame");
		for (int rx=0; rx < rows; ++rx) {
			mxLog("row=%d sortedrow=%d %d %d %d %d %d %d %d", rx, indexVector[rx],
			      bool(ofiml->ordinalMissingSame[rx]),
			      bool(ofiml->continuousMissingSame[rx]),
			      bool(ofiml->missingSameOrdinalSame[rx]),
			      bool(ofiml->missingSameContinuousSame[rx]),
			      bool(ofiml->continuousSame[rx]),
			      bool(ofiml->missingSame[rx]),
			      bool(ofiml->ordinalSame[rx]));
		}
	}
}

static bool dispatchByRow(FitContext *u_fc, omxFitFunction *u_localobj,
			  omxFIMLFitFunction *parent, omxFIMLFitFunction *ofiml, int &want)
{
  if (parent->verbose >= 4) mxLog("%s: jointStrat %d", ofiml->name(), ofiml->jointStrat);
	switch (ofiml->jointStrat) {
	case JOINT_CONDORD:{
		condOrdByRow batch(u_fc, u_localobj, parent, ofiml, want);
		return batch.eval();
	}
	case JOINT_AUTO:
	case JOINT_CONDCONT:{
		condContByRow batch(u_fc, u_localobj, parent, ofiml, want);
		return batch.eval();
	}
	default: OOPS;
	}
}

static omxFitFunction *getChildFitObj(FitContext *fc, omxMatrix *fitMatrix, int i)
{
	FitContext *kid = fc->childList[i];
	omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
	omxFitFunction *childFit = childMatrix->fitFunction;
	return childFit;
}

static omxFIMLFitFunction *getChildFIMLObj(FitContext *fc, omxMatrix *fitMatrix, int i)
{
	omxFitFunction *childFit = getChildFitObj(fc, fitMatrix, i);
	omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) childFit);
	return ofo;
}

static void recalcRowBegin(FitContext *fc, omxMatrix *fitMatrix, int parallelism)
{
	if (parallelism == 1) {
		omxFIMLFitFunction *ff = ((omxFIMLFitFunction*)fitMatrix->fitFunction);
		ff->rowBegin = 0;
	} else {
		int begin = 0;
		for(int i = 0; i < parallelism; i++) {
			omxFIMLFitFunction *ff = getChildFIMLObj(fc, fitMatrix, i);
			//mxLog("%d--%d", begin, ff->rowCount);
			ff->rowBegin = begin;
			begin += ff->rowCount;
		}
	}
}

static void setParallelism(FitContext *fc, omxFIMLFitFunction *parent,
			   omxMatrix *fitMatrix, int parallelism)
{
	int rows = parent->indexVector.size();
	if (parallelism == 1) {
		omxFIMLFitFunction *ff = ((omxFIMLFitFunction*)fitMatrix->fitFunction);
		ff->rowCount = rows;
	} else {
		int stride = (rows / parallelism);

		for(int i = 0; i < parallelism; i++) {
			omxFIMLFitFunction *ff = getChildFIMLObj(fc, fitMatrix, i);
			int rowcount = stride;
			if (i == parallelism-1) rowcount = rows - stride*i;
			ff->rowCount = rowcount;
		}
	}
	recalcRowBegin(fc, fitMatrix, parallelism);
	parent->curParallelism = parallelism;
}

void omxFIMLFitFunction::invalidateCache()
{
	builtCache = false;
	indexVector.clear();
	openmpUser = false;

	rowCount = std::numeric_limits<decltype(rowCount)>::max();
	omxResizeMatrix(rowLikelihoods, data->nrows(), 1);
}

void omxFIMLFitFunction::compute2(int want, FitContext *fc)
{
	auto *off = this;
	omxFIMLFitFunction* ofiml = this;

	if (want & FF_COMPUTE_INITIAL_FIT) return;
	
	int numFree = fc->getNumFree();

	if (!builtCache) {
    if (!fc->isClone()) {
			if (!indexVector.size()) sortData(off);
			openmpUser = rowwiseParallel && fc->permitParallel;
      diagParallel(OMX_DEBUG, "%s: openmpUser = %d", name(), openmpUser);
    } else {
			omxMatrix *pfitMat = fc->getParentState()->getMatrixFromIndex(off->matrix);
      auto *pff = (omxFIMLFitFunction*) pfitMat->fitFunction;
      if (!pff->openmpUser) OOPS;
      parent = pff;
      elapsed.resize(ELAPSED_HISTORY_SIZE);
      curElapsed = NA_INTEGER;
		}
		builtCache = true;
	}

	if (want & FF_COMPUTE_PREOPTIMIZE) {
		inUse = true;
		return;
	}
	if (want & FF_COMPUTE_FINAL_FIT && !ofiml->inUse) return;

	if(OMX_DEBUG) mxLog("%s: joint FIML; openmpUser=%d", off->name(), off->openmpUser);

	omxMatrix* fitMatrix  = off->matrix;

	omxFIMLFitFunction *myParent = ofiml->parent? ofiml->parent : ofiml;
	int rows = int(myParent->indexVector.size());
	if (!means) {
		if (want & FF_COMPUTE_FINAL_FIT) return;
		complainAboutMissingMeans(expectation);
		return;
	}

  if (want & FF_COMPUTE_GRADIENT){
  	if(Global->analyticGradients && off->expectation->canProvideSufficientDerivs){
  		if(dSigma_dtheta.size() != size_t(numFree)){ 
  			dSigma_dtheta.resize(numFree);
  		}
  		if(dNu_dtheta.size() != size_t(numFree)){
  			dNu_dtheta.resize(numFree);
  		}
  		off->expectation->provideSufficientDerivs(fc, dSigma_dtheta, dNu_dtheta);
  	}
  	else{	
  		invalidateGradient(fc);
  	}
  }

	bool failed = false;

	if (fc->childList.size() == 0) {
		setParallelism(fc, myParent, fitMatrix, 1);
	}

	int childStateId = 0;
	if (fc->childList.size()) {
		omxFitFunction *ff = getChildFitObj(fc, fitMatrix, 0);
		childStateId = ff->matrix->currentState->getId();
	}
	if (myParent->curParallelism == 0 || myParent->origStateId != childStateId) {
		int numChildren = fc? fc->childList.size() : 0;
		int parallelism = (numChildren == 0 || !off->openmpUser) ? 1 : numChildren;
		if (OMX_DEBUG_FIML_STATS) parallelism = 1;
		if (parallelism > rows) parallelism = rows;
		setParallelism(fc, myParent, fitMatrix, parallelism);
		myParent->origStateId = childStateId;
		myParent->curElapsed = 0;
	}

	nanotime_t startTime = 0;
	if (ofiml->verbose >= 3) {
		startTime = get_nanotime();
		mxLog("%s: start eval with %d threads", off->name(), myParent->curParallelism);
	}

	ofiml->skippedRows = 0;

	bool reduceParallelism = false;
	if (myParent->curParallelism > 1) {
		if (OMX_DEBUG) {
			omxFIMLFitFunction *ofo = getChildFIMLObj(fc, fitMatrix, 0);
			if (!ofo->parent) OOPS;
		}

		for (int tx=0; tx < myParent->curParallelism; ++tx) {
			omxFIMLFitFunction *ofo = getChildFIMLObj(fc, fitMatrix, tx);
			ofo->skippedRows = 0;
		}
#pragma omp parallel for num_threads(myParent->curParallelism) reduction(||:failed)
		for(int i = 0; i < myParent->curParallelism; i++) {
			FitContext *kid = fc->childList[i];
			omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			try {
				failed |= dispatchByRow(kid, childFit, myParent, ofiml, want);
			} catch (const std::exception& e) {
				omxRaiseErrorf("%s", e.what());
			} catch (...) {
				omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
			}
		}
		for (int tx = 0; tx < myParent->curParallelism; tx++) {
			omxFIMLFitFunction *ofo = getChildFIMLObj(fc, fitMatrix, tx);
			ofiml->skippedRows += ofo->skippedRows;
		}

		myParent->curElapsed = (myParent->curElapsed+1) % ELAPSED_HISTORY_SIZE;
		if (myParent->curElapsed == 0) {
			std::vector<nanotime_t> thrElapsed;
			thrElapsed.reserve(myParent->curParallelism);
			for (int tx=0; tx < myParent->curParallelism; ++tx) {
				omxFIMLFitFunction *ofo = getChildFIMLObj(fc, fitMatrix, tx);
				thrElapsed.push_back(ofo->getMedianElapsedTime());
			}
			auto mm = std::minmax_element(thrElapsed.begin(), thrElapsed.end());
			int minT = mm.first - thrElapsed.begin();
			int maxT = mm.second - thrElapsed.begin();
			if (*mm.first < 1000000) { // 1ms
				reduceParallelism = true;
			} else {
				double imbalance = (*mm.second - *mm.first) / (5.0 * (*mm.second + *mm.first));
				imbalance = std::min(imbalance, 0.2);
				omxFIMLFitFunction *minC = getChildFIMLObj(fc, fitMatrix, minT);
				omxFIMLFitFunction *maxC = getChildFIMLObj(fc, fitMatrix, maxT);
				int toMove = maxC->rowCount * imbalance;
				if (toMove > 0) {
					if (ofiml->verbose >= 3) {
						mxLog("transfer work %d (%.2f%%) from thread %d to %d",
						      toMove, imbalance*100, maxT, minT);
					}
					minC->rowCount += toMove;
					maxC->rowCount -= toMove;
					recalcRowBegin(fc, fitMatrix, myParent->curParallelism);
				}
			}
		}
	} else {
		failed |= dispatchByRow(fc, off, myParent, ofiml, want);
	}

	if (ofiml->verbose >= 3) {
    mxLog("%s: done in %.2fms", off->name(), (get_nanotime() - startTime)/1000000.0);
  }

	if (!returnVector && ofiml->skippedRows == rows) {
		// all rows skipped
		failed = true;
	}
	if (failed) {
		if(!returnVector) {
			omxSetMatrixElement(off->matrix, 0, 0, NA_REAL);
		} else {
			EigenArrayAdaptor got(off->matrix);
			got.setZero(); // not sure if NA_REAL is safe
		}
		return;
	}

	if (!returnVector) {
		if (wantRowLikelihoods) {
			double sum = 0.0;
			EigenVectorAdaptor rl(ofiml->rowLikelihoods);
			sum = rl.array().log().sum();
			omxSetMatrixElement(off->matrix, 0, 0, sum);
		} else if (myParent->curParallelism > 1) {
			double sum = 0.0;
			for(int i = 0; i < myParent->curParallelism; i++) {
				FitContext *kid = fc->childList[i];
				omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
				EigenVectorAdaptor got(childMatrix);
				sum += got[0];
			}
			omxSetMatrixElement(off->matrix, 0, 0, sum);
		}
		fc->skippedRows += ofiml->skippedRows;
		EigenVectorAdaptor got(off->matrix);
		got[0] = addSkippedRowPenalty(got[0], ofiml->skippedRows);
		got[0] *= Global->llScale;
		if(ofiml->verbose + OMX_DEBUG >= 2) {
			mxLog("%s: total likelihood is %3.3f skipped %d",
			      off->name(), off->matrix->data[0], ofiml->skippedRows);
		}
	} else {
		omxCopyMatrix(off->matrix, ofiml->rowLikelihoods);
		if (OMX_DEBUG) {
			omxPrintMatrix(ofiml->rowLikelihoods, "row likelihoods");
		}
	}
	if (reduceParallelism) {
		setParallelism(fc, myParent, fitMatrix, myParent->curParallelism - 1);
		if (ofiml->verbose >= 2) {
			mxLog("reducing number of threads to %d", myParent->curParallelism);
		}
	}
	if (openmpUser && !fc->isClone() && want & FF_COMPUTE_BESTFIT) {
    if (curParallelism == 1) {
      diagParallel(OMX_DEBUG, "%s: rowwiseParallel used %d threads; "
                   "recommend rowwiseParallel=FALSE",
                   name(), curParallelism);
    } else {
      diagParallel(OMX_DEBUG, "%s: rowwiseParallel used %d threads", name(), curParallelism);
    }
  }
}

omxFitFunction *omxInitFIMLFitFunction()
{ return new omxFIMLFitFunction; }

void omxFIMLFitFunction::init()
{
	auto *off = this;
	auto *newObj = this;

	if(OMX_DEBUG) {
		mxLog("Initializing FIML fit function function.");
	}
	off->canDuplicate = TRUE;

	if(expectation == NULL) {
		omxRaiseErrorf("FIML cannot fit without model expectations.");
		return;
	}

	newObj->curParallelism = 0;
	newObj->origStateId = 0;
	newObj->inUse = false;
	newObj->parent = 0;
	newObj->expectationComputeCount = 0;
	newObj->conditionMeanCount = 0;
	newObj->conditionCovCount = 0;
	newObj->invertCount = 0;
	newObj->ordSetupCount = 0;
	newObj->ordDensityCount = 0;
	newObj->contDensityCount = 0;
	newObj->wantRowLikelihoods = false;

	cov = omxGetExpectationComponent(expectation, "cov");
	if(cov == NULL) {
		omxRaiseErrorf("%s: covariance not found in expectation '%s'",
			       name(), expectation->name);
		return;
	}

	means = omxGetExpectationComponent(expectation, "means");

    newObj->smallMeans = NULL;
    newObj->ordMeans   = NULL;
    newObj->contRow    = NULL;
    newObj->ordCov     = NULL;
    newObj->ordContCov = NULL;
    newObj->halfCov    = NULL;
    newObj->reduceCov  = NULL;

	if(OMX_DEBUG) {
		mxLog("Accessing data source.");
	}
	newObj->data = off->expectation->data;
	newObj->rowBegin = 0;

	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood option.");
	}

	{
		ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
		newObj->verbose = Rf_asInteger(Rverbose) + OMX_DEBUG;
	}

	const char *jointStratName = CHAR(Rf_asChar(R_do_slot(rObj, Rf_install("jointConditionOn"))));
	if (strEQ(jointStratName, "auto")) {
		newObj->jointStrat = JOINT_AUTO;
	} else if (strEQ(jointStratName, "ordinal")) {
		newObj->jointStrat = JOINT_CONDORD;
	} else if (strEQ(jointStratName, "continuous")) {
		newObj->jointStrat = JOINT_CONDCONT;
	} else { mxThrow("jointConditionOn '%s'?", jointStratName); }

	returnVector = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));

	units = returnVector? FIT_UNITS_PROBABILITY : FIT_UNITS_MINUS2LL;
	if (returnVector) wantRowLikelihoods = true;

	newObj->rowLikelihoods = omxInitMatrix(newObj->data->nrows(), 1, off->matrix->currentState);
	newObj->otherRowwiseValues = omxInitMatrix(newObj->data->nrows(), 2, off->matrix->currentState);

  invalidateCache();

	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood population option.");
	}
	populateRowDiagnostics = Rf_asInteger(R_do_slot(rObj, Rf_install("rowDiagnostics")));
	if (populateRowDiagnostics) wantRowLikelihoods = true;

	newObj->useSufficientSets = !newObj->wantRowLikelihoods;

	if(OMX_DEBUG) {
		mxLog("Accessing variable mapping structure.");
	}
	auto dc = off->expectation->getDataColumns();

	if(OMX_DEBUG) {
		mxLog("Accessing Threshold matrix.");
	}

	isOrdinal.resize(dc.size());

	numOrdinal=0;
	numContinuous=0;
	for(int j = 0; j < dc.size(); j++) {
		int var = dc[j];
		isOrdinal[j] = omxDataColumnIsFactor(newObj->data, var);
		if (isOrdinal[j]) numOrdinal += 1;
		else numContinuous += 1;
	}

	rowwiseParallel = Rf_asLogical(R_do_slot(rObj, Rf_install("rowwiseParallel")));
  if (rowwiseParallel == NA_INTEGER) {
    rowwiseParallel = numOrdinal >= 10;
		if (verbose >= 1) {
			mxLog("%s: set rowwiseParallel=%d", name(), rowwiseParallel);
		}
  }
	if (newObj->jointStrat == JOINT_AUTO && 0 == numOrdinal) {
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

    if(numOrdinal > 0) {
        if(OMX_DEBUG) mxLog("Ordinal Data detected");

        newObj->ordContCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->halfCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        omxCopyMatrix(newObj->ordContCov, newObj->cov);
    }
}
