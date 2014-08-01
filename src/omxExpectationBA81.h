/*
  Copyright 2012-2014 Joshua Nathaniel Pritikin and contributors

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _OMX_EXPECTATIONBA81_H_
#define _OMX_EXPECTATIONBA81_H_

#include "omxExpectation.h"
#include "omxOpenmpWrap.h"
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
		omxRecompute(_latentMeanOut, FF_COMPUTE_FIT, fc);
		memcpy(mean.derived().data(), _latentMeanOut->data, sizeof(double) * grp.maxAbilities);
	}
	
	cov.derived().resize(grp.maxAbilities, grp.maxAbilities);
	if (!_latentCovOut) {
		cov.setIdentity();
	} else {
		omxRecompute(_latentCovOut, FF_COMPUTE_FIT, fc);
		memcpy(cov.derived().data(), _latentCovOut->data, sizeof(double) * grp.maxAbilities * grp.maxAbilities);
	}
}

extern const struct rpf *rpf_model;
extern int rpf_numModels;

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
