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

using namespace UndirectedGraph;

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

		template <typename T1>
		void setMean(Eigen::MatrixBase<T1> &meanIn)
		{
			mean.resize(varMap.size());
			for (int vx=0, dx=0; vx < meanIn.rows(); ++vx) {
				if (!varMask[vx]) continue;
				mean[dx++] = meanIn[vx];
			}
		}

		inline double likelihood(int row);
		template <typename T1, typename T2>
		double likelihood(Eigen::MatrixBase<T1> &lbound, Eigen::MatrixBase<T2> &ubound);
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
	void attach(const Eigen::MatrixBase<T> &dc, omxData *data, omxMatrix *tMat,
		    std::vector< omxThresholdColumn > &colInfo)
	{
		dataColumns = dc;
		this->data = data;
		this->thresholdMat = tMat;
		this->colInfoPtr = &colInfo;
	};

	struct CorCmp {
		Eigen::MatrixXd &cor;
		CorCmp(Eigen::MatrixXd *cor) : cor(*cor) {};
		bool operator()(int i, int j) { return fabs(cor.data()[i]) > fabs(cor.data()[j]); };
	};

	template <typename T1>
	void setCovariance(Eigen::MatrixBase<T1> &cov, FitContext *fc)
	{
		stddev = cov.diagonal().array().sqrt();

		std::vector<int> cells;
		cor.resize(cov.rows(), cov.cols());
		for(int i = 1; i < cov.rows(); i++) {
			for(int j = 0; j < i; j++) {
				double val = cov(i, j) / (stddev[i] * stddev[j]);
				if (fabs(val) > 1.0) {
					if (fc) {
						fc->recordIterationError("Found correlation with absolute value"
									 " greater than 1 (r=%.2f)", val);
					} else {
						// Signal disaster
						cov(0,0) = NA_REAL;
						val = NA_REAL;
					}
				}
				cor(i,j) = val;
				if (val != 0.0) {
					cells.push_back(&cor.coeffRef(i,j) - &cor.coeffRef(0,0));
				}
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
			mask.assign(cov.cols(), false);
			for (std::set<int>::iterator it=members.begin(); it!=members.end(); ++it) {
				mask[*it] = true;
			}
			blocks[bx].setCorrelation(cor);
			bx += 1;
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
	double likelihood(Eigen::MatrixBase<T1> &lbound, Eigen::MatrixBase<T2> &ubound)
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
double OrdinalLikelihood::block::likelihood(Eigen::MatrixBase<T1> &lbound, Eigen::MatrixBase<T2> &ubound)
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

double OrdinalLikelihood::block::likelihood(int row)
{
	std::vector< omxThresholdColumn > &colInfo = *ol->colInfoPtr;
	Eigen::ArrayXi &ordColumns = ol->ordColumns;
	EigenMatrixAdaptor tMat(ol->thresholdMat);
	for (int ox=0, vx=0; ox < ordColumns.size(); ++ox) {
		if (!varMask[ox]) continue;
		int j = ordColumns[ox];
		int var = omxVectorElement(ol->dataColumns, j);
		int pick = omxIntDataElement(ol->data, row, var) - 1;
		double sd = ol->stddev[ox];
		int tcol = colInfo[j].column;
		if (pick == 0) {
			lThresh[vx] = NA_REAL;
			uThresh[vx] = (tMat(pick, tcol) - mean[vx]) / sd;
			Infin[vx] = 0;
		} else if (pick == colInfo[j].numThresholds) {
			lThresh[vx] = (tMat(pick-1, tcol) - mean[vx]) / sd;
			uThresh[vx] = NA_REAL;
			Infin[vx] = 1;
		} else {
			lThresh[vx] = (tMat(pick-1, tcol) - mean[vx]) / sd;
			uThresh[vx] = (tMat(pick, tcol) - mean[vx]) / sd;
			Infin[vx] = 2;
		}
		vx += 1;
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
