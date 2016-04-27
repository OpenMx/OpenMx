/*
 * Copyright 2012-2016 Joshua Nathaniel Pritikin and contributors
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

template <typename CovType>
struct BA81EstepBase {
	void addRow1(class ifaGroup *state, int px, double *Qweight, double *out);
};

template <typename T, typename CovType>
struct BA81Estep : BA81EstepBase<CovType> {
	std::vector<double> thrExpected;

	void begin(class ifaGroup *state, T extraData);
	void addRow(class ifaGroup *state, T extraData, int px, double *Qweight, int thrId);
	void recordTable(class ifaGroup *state, T extraData);
	bool hasEnd() { return true; }
};

template <typename T>
struct BA81LatentFixed {
	void begin(class ifaGroup *state, T extraData) {}
	void normalizeWeights(class ifaGroup *state, T extraData, int px, double *Qweight, double weight, int thrid);
	void end(class ifaGroup *state, T extraData) {};
	bool hasEnd() { return false; }
};

template <typename T>
struct BA81LatentSummary {
	void begin(class ifaGroup *state, T extraData);
	void normalizeWeights(class ifaGroup *state, T extraData, int px, double *Qweight, double weight, int thrId);
	void end(class ifaGroup *state, T extraData);
	bool hasEnd() { return true; }
};

class BA81Expect {
 public:
	class ifaGroup grp;
	int totalOutcomes() { return grp.totalOutcomes; }
	const double *itemSpec(int ix) { return grp.spec[ix]; }
	int numItems() { return grp.numItems(); }
	int getNumUnique() { return (int) grp.rowMap.size(); }
	int itemOutcomes(int ix) { return grp.itemOutcomes[ix]; }

	const char *name;              // from omxExpectation
	double LogLargestDouble;       // should be const but need constexpr
	double LargestDouble;          // should be const but need constexpr

	// data characteristics
	omxData *data;
	double weightSum;                // sum of rowWeight
	// aggregate distribution of data in quadrature
	std::vector<double> thrDweight;  // quad.weightTableSize * numThreads

	// quadrature related
	struct ba81NormalQuad &getQuad() { return grp.quad; }

	// estimation related
	omxMatrix *itemParam;
	double *EitemParam;
	double SmallestPatternLik;
	double *expected;                     // totalOutcomes * totalQuadPoints (E-step table)
	bool expectedUsed;
	int ElatentVersion;

	omxMatrix *_latentMeanOut;
	omxMatrix *_latentCovOut;
	template <typename Tmean, typename Tcov>
	void getLatentDistribution(FitContext *fc, Eigen::MatrixBase<Tmean> &mean, Eigen::MatrixBase<Tcov> &cov);

	omxMatrix *estLatentMean;
	omxMatrix *estLatentCov;

	int itemParamVersion;
	int latentParamVersion;
	enum expectation_type type;
	int verbose;
	bool debugInternal;
	struct omxFitFunction *fit;  // weak pointer

	BA81Expect() : grp(Global->numThreads, true) {};
	const char *getLatentIncompatible(BA81Expect *other);
};

template <typename Tmean, typename Tcov>
void BA81Expect::getLatentDistribution(FitContext *fc, Eigen::MatrixBase<Tmean> &mean, Eigen::MatrixBase<Tcov> &cov)
{
	mean.derived().resize(grp.maxAbilities);
	if (!_latentMeanOut) {
		mean.setZero();
	} else {
		omxRecompute(_latentMeanOut, fc);
		memcpy(mean.derived().data(), _latentMeanOut->data, sizeof(double) * grp.maxAbilities);
	}
	
	cov.derived().resize(grp.maxAbilities, grp.maxAbilities);
	if (!_latentCovOut) {
		cov.setIdentity();
	} else {
		omxRecompute(_latentCovOut, fc);
		memcpy(cov.derived().data(), _latentCovOut->data, sizeof(double) * grp.maxAbilities * grp.maxAbilities);
	}
}

extern const struct rpf *Grpf_model;
extern int Grpf_numModels;

OMXINLINE static void
gramProduct(double *vec, size_t len, double *out)
{
	int cell = 0;
	for (size_t v1=0; v1 < len; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			out[cell] = vec[v1] * vec[v2];
			++cell;
		}
	}
}

void ba81SetupQuadrature(omxExpectation* oo);

void ba81AggregateDistributions(std::vector<struct omxExpectation *> &expectation,
				int *version, omxMatrix *meanMat, omxMatrix *covMat);

static const double BA81_MIN_VARIANCE = .01;

#endif
