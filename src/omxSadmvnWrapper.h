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
 
#ifndef _OMXSADMVNWRAPPER_H
#define _OMXSADMVNWRAPPER_H

#include "omxData.h"
#include "Connectedness.h"
#include <limits>
#include <Rcpp.h>
#include <Eigen/CholmodSupport>
#include <Eigen/Cholesky>
#include <RcppEigenWrap.h>
#include "matrix.h"

#ifdef  __cplusplus
extern "C" {
#endif

void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*,
		     double*, double*, double*, double*, int*, int*);

#ifdef  __cplusplus
}
#endif

#ifdef SHADOW_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#endif

void omxSadmvnWrapper(int numVars, 
	double *corList, double *lThresh, double *uThresh, int *Infin, double *likelihood, int *inform);

using namespace UndirectedGraph;

class OrdinalLikelihood { // rename to mvn cdf ? TODO
 private:
	struct block {
		OrdinalLikelihood *ol;
		Eigen::VectorXd uThresh;
		Eigen::VectorXd lThresh;
		Eigen::VectorXi Infin;
		Eigen::VectorXd mean;
		Eigen::ArrayXd corList;
		std::vector<bool> varMask;
		std::vector<int> varMap;

		template <typename T2>
		void setCorrelation(Eigen::MatrixBase<T2> &corIn)
		{
			varMap.clear();
			for (int vx=0; vx < (int) varMask.size(); ++vx) {
				if (!varMask[vx]) continue;
				varMap.push_back(vx);
			}
			uThresh.resize(varMap.size());
			lThresh.resize(varMap.size());
			Infin.resize(varMap.size());

			corList.resize(triangleLoc1(varMap.size() - 1));
			for(int rx = 1, dr=0; rx < corIn.rows(); rx++) {
				if (!varMask[rx]) continue;
				bool any=false;
				for(int cx = 0, dc=0; cx < rx; cx++) {
					if (!varMask[cx]) continue;
					corList[triangleLoc1(dr) + dc] = corIn(rx, cx);
					dc += 1;
					any = true;
				}
				if (any) dr += 1;
			}
			//mxPrintMat("curList", corList);
		}

		void setZeroMean() {
			mean.resize(varMap.size());
			mean.setZero();
		}

		template <typename T1>
		void setMean(Eigen::MatrixBase<T1> &meanIn)
		{
			mean.resize(varMap.size());
			for (int vx=0, dx=0; vx < meanIn.rows(); ++vx) {
				if (!varMask[vx]) continue;
				mean[dx++] = meanIn[vx];
			}
		}

		inline void loadRow(int row);
		inline double likelihood(int row);
		template <typename T1, typename T2>
		double likelihood(const Eigen::MatrixBase<T1> &lbound, const Eigen::MatrixBase<T2> &ubound);
	};

	Eigen::ArrayXd stddev;
	Eigen::MatrixXd cor;
	std::vector<block> blocks;
	Eigen::VectorXi dataColumns;
	omxData *data;
	omxMatrix *thresholdMat;
	std::vector< omxThresholdColumn > *colInfoPtr;
	Eigen::ArrayXi ordColumns;
 public:
	std::vector<int> itemToThresholdCol;
	omxMatrix *thresholdsMat;
	std::vector<int> numThresholds;

	OrdinalLikelihood()
	{
		data = 0;
		thresholdMat = 0;
		colInfoPtr = 0;
	}

	template <typename T>
	void attach(const Eigen::MatrixBase<T> &dc, omxData *_data, omxMatrix *tMat,
		    std::vector< omxThresholdColumn > &colInfo)
	{
		dataColumns = dc;
		this->data = _data;
		this->thresholdMat = tMat;
		this->colInfoPtr = &colInfo;
	};

	struct CorCmp {
		Eigen::MatrixXd &cor;
		CorCmp(Eigen::MatrixXd *_cor) : cor(*_cor) {};
		bool operator()(int i, int j) { return fabs(cor.data()[i]) > fabs(cor.data()[j]); };
	};

	template <typename T1>
	void setCovariance(Eigen::MatrixBase<T1> &cov, FitContext *fc)
	{
		setCovarianceUnsafe(cov);

		for(int i = 1; i < cov.rows(); i++) {
			for(int j = 0; j < i; j++) {
				double val = cor(i,j);
				if (fabs(val) <= 1.0) continue;
				if (fc) {
					fc->recordIterationError("Found correlation with absolute value"
								 " greater than 1 (r=%.2f)", val);
				} else {
					cov(0,0) = NA_REAL; // Signal disaster
				}
			}
		}
	}

	void setStandardNormal(int dims)
	{
		stddev.resize(dims);
		stddev.setConstant(1.0);
		cor.setIdentity(dims, dims);
		setupCorrelation();
	}

	void setupCorrelation()
	{
		std::vector<int> cells;
		for(int i = 1; i < cor.rows(); i++) {
			for(int j = 0; j < i; j++) {
				if (cor(i,j) == 0.0) continue;
				cells.push_back(&cor.coeffRef(i,j) - &cor.coeffRef(0,0));
			}
		}

		std::sort(cells.begin(), cells.end(), CorCmp(&cor));
		//Eigen::VectorXd ec(cells.size());
		//for (int ex=0; ex < ec.size(); ex++) ec[ex] = cor.data()[cells[ex]];
		//mxPrintMat("ec", ec);

		std::vector<int> region;
		Connectedness::SubgraphType subgraph;
		Connectedness cc(region, subgraph, stddev.size(), false);
		for (int cx=0; cx < (int)cells.size(); ++cx) {
			int offset = cells[cx];
			int col = offset / cor.rows();
			int row = offset % cor.rows();
			if (cc.getSizeIfConnected(row, col) <= Global->maxOrdinalPerBlock) {
				cc.connect(row, col);
			}
		}

		//mxLog("split %d vars into %d blocks", stddev.size(), cc.numSubgraphs());
		blocks.clear();
		blocks.resize(cc.numSubgraphs());
		int bx = 0;
		for (int rx=0; rx < (int)region.size(); ++rx) {
			std::set<int> members;
			if (region[rx] == -1) members.insert(rx);
			else std::swap(members, subgraph[ region[rx] ]);
			if (0 == members.size()) continue;

			blocks[bx].ol = this;
			std::vector<bool> &mask = blocks[bx].varMask;
			mask.assign(cor.cols(), false);
			for (std::set<int>::iterator it=members.begin(); it!=members.end(); ++it) {
				mask[*it] = true;
			}
			blocks[bx].setCorrelation(cor);
			bx += 1;
		}
	}

	// only considers lower triangle
	template <typename T1>
	void setCorrelation(const Eigen::MatrixBase<T1> &corIn)
	{
		stddev.resize(corIn.rows());
		stddev.setConstant(1.0);
		cor = corIn;
		setupCorrelation();
	}

	// only considers lower triangle
	template <typename T1>
	void setCovarianceUnsafe(const Eigen::MatrixBase<T1> &cov)
	{
		stddev = cov.diagonal().array().sqrt();

		cor.resize(cov.rows(), cov.cols());
		for(int i = 1; i < cov.rows(); i++) {
			for(int j = 0; j < i; j++) {
				cor(i,j) = cov(i, j) / (stddev[i] * stddev[j]);
			}
		}

		setupCorrelation();
	}

	void setZeroMean() {
		for (int bx=0; bx < (int)blocks.size(); ++bx) {
			blocks[bx].setZeroMean();
		}
	}

	template <typename T1>
	void setMean(Eigen::MatrixBase<T1> &meanIn)
	{
		for (int bx=0; bx < (int)blocks.size(); ++bx) {
			blocks[bx].setMean(meanIn);
		}
	};

	template <typename T1>
	void setColumns(Eigen::ArrayBase<T1> &colIn)
	{
		ordColumns = colIn;
	};

	inline double likelihood(int row)
	{
		double lk = 1.0;
		for (int bx=0; bx < int(blocks.size()); ++bx) {
			lk *= blocks[bx].likelihood(row);  // log space instead?
		}
		return lk;
	};

	template <typename T1, typename T2>
	double likelihood(const Eigen::MatrixBase<T1> &lbound, const Eigen::MatrixBase<T2> &ubound)
	{
		double lk = 1.0;
		for (int bx=0; bx < int(blocks.size()); ++bx) {
			double l1 = blocks[bx].likelihood(lbound, ubound);
			//mxLog("%g %g", lk, l1);
			lk *= l1;
		}
		return lk;
	};
};

template <typename T1, typename T2>
double OrdinalLikelihood::block::likelihood(const Eigen::MatrixBase<T1> &lbound,
					    const Eigen::MatrixBase<T2> &ubound)
{
	for (int ox=0, vx=0; ox < (int)varMask.size(); ++ox) {
		if (!varMask[ox]) continue;
		double sd = ol->stddev[ox];
		uThresh[vx] = (ubound[ox] - mean[vx]) / sd;
		lThresh[vx] = (lbound[ox] - mean[vx]) / sd;
		Infin[vx] = 2;
		if (!R_finite(lThresh[vx])) {
			Infin[vx] -= 2;
		}
		if (!R_finite(uThresh[vx])) {
			Infin[vx] -= 1;
		}
		vx += 1;
	}
	int inform;
	double ordLik;
	omxSadmvnWrapper(mean.size(), corList.data(),
			 lThresh.data(), uThresh.data(),
			 Infin.data(), &ordLik, &inform);
	// if inform == 1, retry? TODO
	if (inform == 2) {
		return 0.0;
	}
	return ordLik;
}

void OrdinalLikelihood::block::loadRow(int row)
{
	std::vector< omxThresholdColumn > &colInfo = *ol->colInfoPtr;
	Eigen::ArrayXi &ordColumns = ol->ordColumns;
	EigenMatrixAdaptor tMat(ol->thresholdMat);
	for (int ox=0, vx=0; ox < ordColumns.size(); ++ox) {
		if (!varMask[ox]) continue;
		int j = ordColumns[ox];
		int var = ol->dataColumns[j];
		int pick = omxIntDataElement(ol->data, row, var) - 1;
		double sd = ol->stddev[ox];
		int tcol = colInfo[j].column;
		if (pick == 0) {
			lThresh[vx] = -std::numeric_limits<double>::infinity();
			uThresh[vx] = (tMat(pick, tcol) - mean[vx]) / sd;
			Infin[vx] = 0;
		} else if (pick == colInfo[j].numThresholds) {
			lThresh[vx] = (tMat(pick-1, tcol) - mean[vx]) / sd;
			uThresh[vx] = std::numeric_limits<double>::infinity();
			Infin[vx] = 1;
		} else {
			lThresh[vx] = (tMat(pick-1, tcol) - mean[vx]) / sd;
			uThresh[vx] = (tMat(pick, tcol) - mean[vx]) / sd;
			Infin[vx] = 2;
		}
		vx += 1;
	}
}

double OrdinalLikelihood::block::likelihood(int row)
{
	loadRow(row);
	int inform;
	double ordLik;
	omxSadmvnWrapper(varMap.size(), corList.data(),
			 lThresh.data(), uThresh.data(),
			 Infin.data(), &ordLik, &inform);
	if (inform == 1) mxLog("Sadmvn error larger than epsilon");
	if (inform == 2) {
		return 0.0;
	}
	return ordLik;
}

/*
# Dichtefunktion und Verteilung einer multivariate truncated normal
#
# Problem ist die Bestimmung der Randverteilung einer Variablen.
#
# 1. Im bivariaten Fall kann explizit eine Formel angegeben werden (vgl. Arnold (1993))
# 2. Im multivariaten Fall kann ein Integral angegeben werden (vgl. Horrace (2005))
# 3. Bestimmung der Dichtefunktion über das Integral möglich?
# 4. Kann die Verteilungsfunktion pmvnorm() helfen? Kann man dann nach einer Variablen differenzieren?

# Literatur:
#
# Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1, 141-150
# Genz, A. (1993). Comparison of methods for the computation of multivariate normal probabilities. Computing Science and Statistics, 25, 400-405
# Horrace (2005).
# Jack Cartinhour (1990): One-dimensional marginal density functions of a truncated multivariate normal density function
# Communications in Statistics - Theory and Methods, Volume 19, Issue 1 1990 , pages 197 - 203

# Dichtefunktion für Randdichte f(xn) einer Truncated Multivariate Normal Distribution,
# vgl. Jack Cartinhour (1990) "One-dimensional marginal density functions of a truncated multivariate normal density function"
*/

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void _dtmvnorm_marginal(double prob, const Eigen::MatrixBase<T1> &xn, int nn,
			const Eigen::MatrixBase<T2> &sigma,
			const Eigen::MatrixBase<T3> &lower, const Eigen::MatrixBase<T4> &upper,
			Eigen::MatrixBase<T5> &density)
{
	using Eigen::VectorXd;
	using Eigen::MatrixXd;

	for (int dx=0; dx < xn.size(); dx++) {
		if (!(lower[nn] <= xn[dx] && xn[dx] <= upper[nn])) Rf_error("xn out of range");
	}

	if (sigma.rows() == 1) {
		double sd = sqrt(sigma(0,0));
		if (!std::isfinite(prob)) {
			prob = Rf_pnorm5(upper[0], 0, sd, true, false) - Rf_pnorm5(lower[0], 0, sd, true, false);
		}
		for (int dx=0; dx < xn.size(); dx++) {
			density[dx] = Rf_dnorm4(xn[dx], 0, sd, false) / prob;
		}
		return;
	}

	MatrixXd AA = sigma;
	if (InvertSymmetricPosDef(AA, 'L')) Rf_error("Non-positive definite");
	
	MatrixXd A_1;
	struct subset1 {
		int nn;
		subset1(int _nn) : nn(_nn) {};
		bool operator()(int rr) { return rr != nn; };
	} op1(nn);

	subsetCovariance(AA, op1, sigma.cols()-1, A_1);
	if (InvertSymmetricPosDef(A_1, 'L')) Rf_error("Non-positive definite");

	double c_nn = 1.0/sigma(nn, nn);
	VectorXd cc(sigma.rows()-1);
	for (int rx=0, dx=0; rx < sigma.rows(); ++rx) {
		if (rx != nn) cc[dx++] = sigma(rx, nn);
	}

	OrdinalLikelihood ol;
	if (!std::isfinite(prob)) {
		ol.setCovarianceUnsafe(sigma);
		ol.setZeroMean();
		prob = ol.likelihood(lower, upper);
	}

	ol.setCovarianceUnsafe(A_1);
	for (int dx=0; dx < 2; dx++) {
		if (fabs(xn[dx]) == std::numeric_limits<double>::infinity()) {
			density[dx] = 0;
			continue;
		}
		VectorXd lower_n(sigma.rows() - 1);
		VectorXd upper_n(sigma.rows() - 1);
		for (int rx=0, tx=0; rx < sigma.rows(); ++rx) {
			if (rx == nn) continue;
			lower_n[tx] = lower[rx];
			upper_n[tx] = upper[rx];
			++tx;
		}
		VectorXd mu = xn[dx] * cc.array() * c_nn;
		ol.setMean(mu);
		density[dx] = exp(-0.5*xn[dx]*xn[dx]*c_nn) * ol.likelihood(lower_n, upper_n);
	}

	density.array() /= prob * sqrt(M_2PI * sigma(nn, nn));
}

/*
# Computation of the bivariate marginal density F_{q,r}(x_q, x_r) (q != r)
# of truncated multivariate normal distribution 
# following the works of Tallis (1961), Leppard and Tallis (1989)
#
# References:
# Tallis (1961): 
#   "The Moment Generating Function of the Truncated Multi-normal Distribution"
# Leppard and Tallis (1989): 
#   "Evaluation of the Mean and Covariance of the Truncated Multinormal"
# Manjunath B G and Stefan Wilhelm (2009): 
#   "Moments Calculation for the Doubly Truncated Multivariate Normal Distribution"
#
# (n-2) Integral, d.h. zweidimensionale Randdichte in Dimension q und r, 
# da (n-2) Dimensionen rausintegriert werden.
# vgl. Tallis (1961), S.224 und Code Leppard (1989), S.550
#
# f(xq=b[q], xr=b[r])
#
# Attention: Function is not vectorized at the moment!
# Idee: Vektorisieren xq, xr --> die Integration Bounds sind immer verschieden,
#       pmvnorm() kann nicht vektorisiert werden. Sonst spart man schon ein bisschen Overhead.
# Der eigentliche bottleneck ist aber pmvnorm().
# Gibt es Unterschiede bzgl. der verschiedenen Algorithmen GenzBretz() vs. Miwa()?
# pmvnorm(corr=) kann ich verwenden
 */
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void _dtmvnorm_marginal2(double alpha, const Eigen::MatrixBase<T1> &xq, const Eigen::MatrixBase<T2> &xr,
			 int qq, int rr,
			 const Eigen::MatrixBase<T3> &sigma,
			 const Eigen::MatrixBase<T4> &lower, const Eigen::MatrixBase<T5> &upper,
			 Eigen::MatrixBase<T6> &density)
{
	using Eigen::VectorXd;
	using Eigen::Vector4d;
	using Eigen::MatrixXd;
	using Eigen::DiagonalMatrix;

	if (sigma.rows() < 2) Rf_error("Dimension must be >= 2");

	for (int dx=0; dx < xq.size(); dx++) {
		if (!(lower[qq] <= xq[dx] && xq[dx] <= upper[qq])) Rf_error("xq[%d] out of range", dx);
	}
	for (int dx=0; dx < xr.size(); dx++) {
		if (!(lower[rr] <= xr[dx] && xr[dx] <= upper[rr])) Rf_error("xr[%d] out of range", dx);
	}

	OrdinalLikelihood ol;
	if (!std::isfinite(alpha)) {
		ol.setCovarianceUnsafe(sigma);
		ol.setZeroMean();
		alpha = ol.likelihood(lower, upper);
	}
	
	struct subset1 {
		int qq, rr;
		bool flip;
		subset1(int _qq, int _rr) : qq(_qq), rr(_rr), flip(false) {};
		bool operator()(int nn) { return (rr == nn || qq == nn) ^ flip; };
	} op1(qq, rr);

	MatrixXd qrSigma;
	subsetCovariance(sigma, op1, 2, qrSigma);
	SimpCholesky< Eigen::MatrixXd > covDecomp;
	covDecomp.compute(qrSigma);
	covDecomp.refreshInverse();
	for (int dx=0; dx < density.size(); dx++) {
		VectorXd loc(2);
		if (fabs(xq[dx]) == std::numeric_limits<double>::infinity() ||
		    fabs(xr[dx]) == std::numeric_limits<double>::infinity()) {
			density[dx] = 0.0;
			continue;
		}
		loc[0] = xq[dx];
		loc[1] = xr[dx];
		double dist = loc.transpose() * covDecomp.getInverse() * loc;
		density[dx] = exp(-0.5 * (M_LN_2PI * 2 + covDecomp.log_determinant() + dist)) / alpha;
	}

	if (sigma.rows() == 2) return;

	VectorXd SD = sigma.diagonal().array().sqrt();
	VectorXd lowerStd = lower.array() / SD.array();
	VectorXd upperStd = upper.array() / SD.array();
	Vector4d xqStd = xq.array() / SD[qq];
	Vector4d xrStd = xr.array() / SD[rr];
	VectorXd iSDcoef = 1.0/SD.array();
	DiagonalMatrix<double, Eigen::Dynamic> iSD(iSDcoef);
	MatrixXd RR = iSD * sigma * iSD;
	MatrixXd RINV = RR;
	if (InvertSymmetricPosDef(RINV, 'L')) Rf_error("Non-positive definite");
	MatrixXd WW;
	op1.flip = true;
	subsetCovariance(RINV, op1, RINV.cols()-2, WW);
	if (InvertSymmetricPosDef(WW, 'L')) Rf_error("Non-positive definite");
	MatrixXd RQR(WW.rows(), WW.cols());
	for (int cx=0; cx < WW.cols(); ++cx) {
		for (int rx=cx; rx < WW.rows(); ++rx) {
			RQR(rx,cx) = WW(rx,cx) / sqrt(WW(rx,rx) * WW(cx,cx));
		}
	}
	Eigen::Matrix<bool, Eigen::Dynamic, 1> isInf =
		(xqStd.array().abs() == std::numeric_limits<double>::infinity() ||
		 xrStd.array().abs() == std::numeric_limits<double>::infinity());
	MatrixXd AQR(xq.size(), WW.rows());
	MatrixXd BQR(xq.size(), WW.rows());
	for (int ii=0, m2=0; ii < RR.rows(); ++ii) {
		if (ii == qq || ii == rr) continue;
		double BSQR = (RR(qq, ii) - RR(qq, rr) * RR(rr, ii)) / (1 - RR(qq, rr)*RR(qq, rr));
		double BSRQ = (RR(rr, ii) - RR(qq, rr) * RR(qq, ii)) / (1 - RR(qq, rr)*RR(qq, rr));
		double RSRQ = (1 - RR(ii, qq)*RR(ii, qq)) * (1 - RR(qq, rr)*RR(qq, rr));
		// partial correlation coefficient R[r,i] given q
		RSRQ = (RR(ii, rr) - RR(ii, qq) * RR(qq, rr)) / sqrt(RSRQ);

		double denom = sqrt((1 - RR(ii, qq)*RR(ii, qq)) * (1 - RSRQ*RSRQ));
		
		// lower integration bound
		AQR.col(m2) = ((lowerStd[ii] - BSQR * xqStd.array() - BSRQ * xrStd.array()) / denom);
		
		// upper integration bound
		BQR.col(m2) = ((upperStd[ii] - BSQR * xqStd.array() - BSRQ * xrStd.array()) / denom);

		for (int rx=0; rx < isInf.size(); ++rx) {
			if (!isInf[rx]) continue;
			AQR(rx,m2) = -std::numeric_limits<double>::infinity();
			BQR(rx,m2) = std::numeric_limits<double>::infinity();
		}

		m2 += 1;
	}
	
	for (int ii=0; ii < xq.size(); ++ii) {
		if (RQR.rows() == 1) {
			ol.setStandardNormal(1);
		} else {
			ol.setCorrelation(RQR);
		}
		ol.setZeroMean();
		double prob = ol.likelihood(AQR.row(ii), BQR.row(ii));
		density[ii] *= prob;
	}
}

// Translated from tmvtnorm 1.4-10

// Could search for independent blocks but maybe not worth it
// because conditioning on ordinal is not efficient when there
// are a large number of ordinal patterns.

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void _mtmvnorm(double prob, const Eigen::MatrixBase<T1> &sigma,
	       const Eigen::MatrixBase<T2> &lower, const Eigen::MatrixBase<T3> &upper,
	       Eigen::MatrixBase<T4> &xi, Eigen::MatrixBase<T5> &U11)
{
	using Eigen::Vector2d;
	using Eigen::Vector4d;
	using Eigen::VectorXd;
	using Eigen::MatrixXd;
	using Eigen::Matrix;

	int kk = sigma.rows();

	Vector2d marginalBounds;
	Vector2d marginalOut;
	VectorXd F_a(kk);
	VectorXd F_b(kk);
	for (int qq=0; qq < kk; ++qq) {
		marginalBounds[0] = lower[qq];
		marginalBounds[1] = upper[qq];
		_dtmvnorm_marginal(prob, marginalBounds, qq, sigma, lower, upper, marginalOut);
		F_a[qq] = marginalOut[0];
		F_b[qq] = marginalOut[1];
	}
	xi = sigma * (F_a - F_b);

	MatrixXd F2(kk,kk);
	F2.diagonal().setZero();
	Vector4d xq;
	Vector4d xr;
	Vector4d marginal2Out;
	for (int qq=0; qq < kk-1; ++qq) {          //col
		for (int ss=qq+1; ss < kk; ++ss) { //row
			xq << lower[qq], upper[qq], lower[qq], upper[qq];
			xr << lower[ss], lower[ss], upper[ss], upper[ss];
			_dtmvnorm_marginal2(prob, xq, xr, qq, ss, sigma, lower, upper, marginal2Out);
			F2(qq,ss) = (marginal2Out[0] - marginal2Out[1]) - (marginal2Out[2] - marginal2Out[3]);
		}
	}
	F2 = F2.selfadjointView<Eigen::Upper>();
	
	F_a.array() *= lower.array();
	F_b.array() *= upper.array();
	for (int kx=0; kx < kk; ++kx) {
		if (!std::isfinite(F_a[kx])) F_a[kx] = 0.0;
		if (!std::isfinite(F_b[kx])) F_b[kx] = 0.0;
	};
	VectorXd F_ab = (F_a - F_b).array() / sigma.diagonal().array();

	U11.derived().resize(kk,kk);
	for (int ii =0; ii < kk; ++ii) {          //col
		for (int jj =ii; jj < kk; ++jj) { //row
			double sum = 0.0;
			for (int sr=0; sr < kk; ++sr) {  //aka qq
				sum += sigma(ii,sr) * sigma(jj,sr) * F_ab[sr];
				if (jj != sr) {
					double sum2 = 0.0;
					for (int sc=0; sc < kk; ++sc) {
						double tt = sigma(jj,sc) - sigma(sr,sc) * sigma(jj,sr) / sigma(sr,sr);
						sum2 += tt * F2(sr,sc);
					}
					sum2 *= sigma(ii,sr);
					sum += sum2;
				}
			}
			U11(ii,jj) = sigma(ii,jj) + sum;
		}
	}
	U11 -= xi * xi.transpose();
}

#endif 
