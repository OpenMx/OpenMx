/*
 * Copyright 2012-2017 Joshua Nathaniel Pritikin and contributors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

#include "omxExpectation.h"
#include "ba81quad.h"

enum expectation_type {
	EXPECTATION_AUGMENTED, // E-M
	EXPECTATION_OBSERVED,  // regular
};

template <typename T>
struct BA81Estep {
	void begin(ifaGroup *state);
	void addRow(class ifaGroup *state, int px, int thrId);
	void recordTable(class ifaGroup *state);
	bool hasEnd() { return true; }
};

template <typename T>
struct BA81LatentFixed {
	bool wantSummary() { return false; };
	void normalizeWeights(class ifaGroup *state, T extraData, int px, double weight, int thrid);
	void end(class ifaGroup *state, T extraData) {};
	bool hasEnd() { return false; }
};

template <typename T>
struct BA81LatentSummary {
	bool wantSummary() { return true; };
	void normalizeWeights(class ifaGroup *state, T extraData, int px, double weight, int thrId);
	void end(class ifaGroup *state, T extraData);
	bool hasEnd() { return true; }
};

class BA81Expect : public omxExpectation {
 public:
	virtual ~BA81Expect();
	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual void populateAttr(SEXP expectation);
	virtual omxMatrix *getComponent(const char*);

	class ifaGroup grp;
	int totalOutcomes() { return grp.totalOutcomes; }
	const double *itemSpec(int ix) { return grp.spec[ix]; }
	int numItems() { return grp.numItems(); }
	int getNumUnique() { return (int) grp.rowMap.size(); }
	int itemOutcomes(int ix) { return grp.itemOutcomes[ix]; }

	double LogLargestDouble;       // should be const but need constexpr
	double LargestDouble;          // should be const but need constexpr

	// data characteristics
	double weightSum;                // sum of rowWeight

	// quadrature related
	struct ba81NormalQuad &getQuad() { return grp.quad; }

	// estimation related
	omxMatrix *itemParam;
	double *EitemParam;
	double SmallestPatternLik;
	bool expectedUsed;
	int ElatentVersion;

	omxMatrix *_latentMeanOut;
	omxMatrix *_latentCovOut;
	template <typename Tmean, typename Tcov>
	void getLatentDistribution(FitContext *fc, Eigen::MatrixBase<Tmean> &mean, Eigen::MatrixBase<Tcov> &cov);

	omxMatrix *estLatentMean;
	omxMatrix *estLatentCov;

	unsigned itemParamVersion;
	unsigned latentParamVersion;
	enum expectation_type type;
	int verbose;
	bool debugInternal;
	struct omxFitFunction *fit;  // weak pointer

	BA81Expect() : grp(Global->numThreads, true) {};
	const char *getLatentIncompatible(BA81Expect *other);

	void refreshPatternLikelihood(bool hasFreeLatent);

	virtual bool hasRowWeights() { return true; };
	virtual int getNumRows() { return data->rows; }
	virtual double *getRowWeights() { return grp.rowWeight; }
	virtual void setRowWeights(double *rw) { grp.rowWeight = rw; };
};

template <typename Tmean, typename Tcov>
void BA81Expect::getLatentDistribution(FitContext *fc, Eigen::MatrixBase<Tmean> &mean, Eigen::MatrixBase<Tcov> &cov)
{
	int dim = grp.quad.abilities();
	mean.derived().resize(dim);
	if (!_latentMeanOut) {
		mean.setZero();
	} else {
		omxRecompute(_latentMeanOut, fc);
		memcpy(mean.derived().data(), _latentMeanOut->data, sizeof(double) * dim);
	}
	
	cov.derived().resize(dim, dim);
	if (!_latentCovOut) {
		cov.setIdentity();
	} else {
		omxRecompute(_latentCovOut, fc);
		memcpy(cov.derived().data(), _latentCovOut->data, sizeof(double) * dim * dim);
	}
}

extern const struct rpf *Grpf_model;
extern int Grpf_numModels;

void ba81RefreshQuadrature(omxExpectation* oo);

void ba81AggregateDistributions(std::vector<struct omxExpectation *> &expectation,
				int *version, omxMatrix *meanMat, omxMatrix *covMat);

#endif
