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

#ifdef  __cplusplus
extern "C" {
#endif

void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*,
		     double*, double*, double*, double*, int*, int*);

#ifdef  __cplusplus
}
#endif

void omxSadmvnWrapper(int numVars, 
	double *corList, double *lThresh, double *uThresh, int *Infin, double *likelihood, int *inform);

class OrdinalLikelihood {
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

		template <typename T1>
		void setCorrelation(Eigen::MatrixBase<T1> &corIn)
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
				for(int cx = 0, dc=0; cx < rx; cx++) {
					if (!varMask[cx]) continue;
					corList[triangleLoc1(dr) + dc] = corIn(rx, cx);
					dc += 1;
				}
				dr += 1;
			}
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

		template <typename T1>
		double likelihood(int row, Eigen::ArrayBase<T1> &ordIndices);
	};

	Eigen::ArrayXd stddev;
	Eigen::MatrixXd cor;
	std::vector<block> blocks;
	omxMatrix *dataColumns;
	omxData *data;
	omxMatrix *thresholdMat;
	std::vector< omxThresholdColumn > *colInfoPtr;
 public:
	std::vector<int> itemToThresholdCol;
	omxMatrix *thresholdsMat;
	std::vector<int> numThresholds;

	OrdinalLikelihood()
	{
		blocks.resize(1);
		blocks[0].ol = this;
		dataColumns = 0;
		data = 0;
		thresholdMat = 0;
		colInfoPtr = 0;
	}

	void attach(omxMatrix *dc, omxData *data, omxMatrix *tMat,
		    std::vector< omxThresholdColumn > &colInfo)
	{
		dataColumns = dc;
		this->data = data;
		this->thresholdMat = tMat;
		this->colInfoPtr = &colInfo;
	};

	template <typename T1>
	void setCovariance(Eigen::MatrixBase<T1> &cov, FitContext *fc)
	{
		stddev = cov.diagonal().array().sqrt();
		cor.resize(cov.rows(), cov.cols());
		for(int i = 1; i < cov.rows(); i++) {
			for(int j = 0; j < i; j++) {
				double val = cov(i, j) / (stddev[i] * stddev[j]);
				if (fabs(val) > 1.0) {
					if (fc) {
						fc->recordIterationError("Found correlation with absolute value"
									 " greater than 1 (r=%.2f)", val);
					} else {
						val = NA_REAL;  // Signal disaster
					}
				}
				cor(i,j) = val;
			}
		}

		// analyze covariances and split into blocks TODO
		std::vector<bool> &mask = blocks[0].varMask;
		mask.assign(cov.cols(), true);
		blocks[0].setCorrelation(cor);
	}

	template <typename T1>
	void setMean(Eigen::MatrixBase<T1> &meanIn)
	{
		for (int bx=0; bx < (int)blocks.size(); ++bx) {
			blocks[bx].setMean(meanIn);
		}
	};

	template <typename T1>
	double likelihood(int row, Eigen::ArrayBase<T1> &indices)
	{
		double lk = 1.0;
		for (int bx=0; bx < int(blocks.size()); ++bx) {
			lk *= blocks[bx].likelihood(row, indices);  // log space instead?
		}
		return lk;
	};
};

template <typename T1>
double OrdinalLikelihood::block::likelihood(int row, Eigen::ArrayBase<T1> &ordIndices)
{
	// varMask TODO
	std::vector< omxThresholdColumn > &colInfo = *ol->colInfoPtr;
	EigenMatrixAdaptor tMat(ol->thresholdMat);
	for (int ox=0; ox < ordIndices.size(); ++ox) {
		int j = ordIndices[ox];
		int var = omxVectorElement(ol->dataColumns, j);
		int pick = omxIntDataElement(ol->data, row, var) - 1;
		double sd = ol->stddev[ox];
		int tcol = colInfo[j].column;
		if (pick == 0) {
			lThresh[ox] = NA_REAL;
			uThresh[ox] = (tMat(pick, tcol) - mean[ox]) / sd;
			Infin[ox] = 0;
		} else if (pick == colInfo[j].numThresholds) {
			lThresh[ox] = (tMat(pick-1, tcol) - mean[ox]) / sd;
			uThresh[ox] = NA_REAL;
			Infin[ox] = 1;
		} else {
			lThresh[ox] = (tMat(pick-1, tcol) - mean[ox]) / sd;
			uThresh[ox] = (tMat(pick, tcol) - mean[ox]) / sd;
			Infin[ox] = 2;
		}
	}
	int inform;
	double ordLik;
	omxSadmvnWrapper(varMap.size(), corList.data(),
			 lThresh.data(), uThresh.data(),
			 Infin.data(), &ordLik, &inform);
	// if inform == 1, retry? TODO
	if (inform == 2) {
		return 0.0;
	}
	return ordLik;
}

#endif 
