/*
  Copyright 2012-2017 Joshua Nathaniel Pritikin and contributors

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

#include "ba81quad.h"
#include "EnableWarnings.h"

using namespace ba81quad;

const double ba81NormalQuad::MIN_VARIANCE = 0.01;

ba81NormalQuad::ba81NormalQuad() : numThreads(-1)
{
	setOne(1);
	layers.resize(1, layer(this));
}

ba81NormalQuad::ba81NormalQuad(ba81NormalQuad &quad) : numThreads(-1)
{
	setOne(quad.One);
	layers.resize(quad.layers.size(), layer(this));
	width = quad.width;
	gridSize = quad.gridSize;
	Qpoint = quad.Qpoint;
	hasBifactorStructure = quad.hasBifactorStructure;
	for (size_t lx=0; lx < quad.layers.size(); ++lx) {
		layers[lx].copyStructure(quad.layers[lx]);
	}
}

int ba81NormalQuad::abilities()
{
	int sum=0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		sum += l1.numAbil();
	}
	return sum;
}

void ba81NormalQuad::layer::copyStructure(ba81NormalQuad::layer &orig)
{
	abilitiesMask = orig.abilitiesMask;
	abilitiesMap = orig.abilitiesMap;
	maxDims = orig.maxDims;
	totalQuadPoints = orig.totalQuadPoints;
	weightTableSize = orig.weightTableSize;
	numSpecific = orig.numSpecific;
	primaryDims = orig.primaryDims;
	totalPrimaryPoints = orig.totalPrimaryPoints;
}

// Depends on item parameters, but not latent distribution
void ba81NormalQuad::cacheOutcomeProb(double *param, bool wantLog)
{
	if (layers.size() != 1) mxThrow("layers.size() != 1");

	layer &l1 = layers[0];
	l1.outcomeProbX.resize(l1.totalOutcomes * l1.totalQuadPoints);

#pragma omp parallel for num_threads(numThreads)
	for (int ix=0; ix < l1.numItems(); ix++) {
		const double *ispec = l1.spec[ix];
		int id = ispec[RPF_ISpecID];
		double *iparam = param + l1.paramRows * ix;
		rpf_prob_t prob_fn = wantLog? Glibrpf_model[id].logprob : Glibrpf_model[id].prob;

		Eigen::VectorXi abx(abscissaDim());
		Eigen::VectorXd abscissa(abscissaDim());

		l1.cacheOutcomeProb(ispec, iparam, prob_fn, ix, abx, abscissa);
	}
	// for (size_t lx=0; lx < layers.size(); ++lx) {
	// 	layer &l1 = layers[lx];
	// 	mxPrintMat("op", l1.outcomeProbX);
	// }
}

void ba81NormalQuad::layer::addSummary(ba81NormalQuad::layer &l1)
{
	Dweight.col(0) += l1.Dweight.col(0);
}

void ba81NormalQuad::addSummary(ba81NormalQuad &quad)
{
	allocSummary();
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].prepSummary();
		layers[lx].addSummary(quad.layers[lx]);
	}
}

ifaGroup::ifaGroup(bool _twotier) :
	itemDims(-1), qwidth(6.0), qpoints(49),
	twotier(_twotier), mean(0), cov(0),
	weightColumnName(0), rowWeight(0), freqColumnName(0), rowFreq(0),
	minItemsPerScore(NA_INTEGER), excludedPatterns(-1)
{}

// The idea here is to avoid denormalized values if they are
// enabled (5e-324 vs 2e-308).  It would be bad if results
// changed depending on the denormalization setting.
// Moreover, we don't lose too much even if denormalized
// values are disabled. This mainly affects models with
// more than a thousand items.
const double ifaGroup::SmallestPatternLik = 1e16 * std::numeric_limits<double>::min();  //constexpr

ifaGroup::~ifaGroup()
{
}

void ifaGroup::importSpec(const List &slotValue)
{
	for (int sx=0; sx < slotValue.size(); ++sx) {
		S4 model = slotValue[sx];
		NumericVector s1 = model.slot("spec");
		spec.push_back(s1.begin());
	}

	dataColumns.reserve(spec.size());
	itemOutcomes.reserve(spec.size());

	impliedParamRows = 0;
	totalOutcomes = 0;
	maxOutcomes = 0;
	for (int cx = 0; cx < numItems(); cx++) {
		const double *ispec = spec[cx];
		int id = ispec[RPF_ISpecID];
		int dims = ispec[RPF_ISpecDims];
		if (itemDims == -1) {
			itemDims = dims;
		} else if (dims != itemDims) {
			mxThrow("All items must have the same number of factors (%d != %d)",
				 itemDims, dims);
		}
		int no = ispec[RPF_ISpecOutcomes];
		itemOutcomes.push_back(no);
		maxOutcomes = std::max(maxOutcomes, no);
		totalOutcomes += no;

		int numParam = (*Glibrpf_model[id].numParam)(ispec);
		if (impliedParamRows < numParam)
			impliedParamRows = numParam;
	}
}

void ifaGroup::verifyFactorNames(const List &dimnames, const char *matName)
{
	using ba81quad::strEQ;
	static const char *dimname[] = { "row", "col" };

	if (dimnames.size() != 2) return;

	for (int dx=0; dx < 2; ++dx) {
		RObject d1 = dimnames[dx];
		if (d1.isNULL()) continue;
		CharacterVector names = as<CharacterVector>(d1);
		if (int(factorNames.size()) != names.size()) {
			mxThrow("%s %snames must be length %d",
					 matName, dimname[dx], (int) factorNames.size());
		}
		int nlen = names.size();
		for (int nx=0; nx < nlen; ++nx) {
			const char *name = names[nx];
			if (strEQ(factorNames[nx].c_str(), name)) continue;
			mxThrow("%s %snames[%d] is '%s', does not match factor name '%s'",
					 matName, dimname[dx], 1+nx, name, factorNames[nx].c_str());
		}
	}
}

void ifaGroup::learnMaxAbilities()
{
	int maxAbilities = 0;
	Eigen::ArrayXi loadings(itemDims);
	loadings.setZero();
	for (int cx = 0; cx < numItems(); cx++) {
		for (int dx=0; dx < itemDims; ++dx) {
			if (getItemParam(cx)[dx] != 0) loadings[dx] += 1;
		}
	}
	maxAbilities = (loadings != 0).count();
	if (itemDims != maxAbilities) {
		for (int lx=0; lx < itemDims; ++lx) {
			if (loadings[lx] == 0) mxThrow("Factor %d does not load on any items", 1+lx);
		}
	}
}

void ifaGroup::import(const List &Rlist)
{
	using ba81quad::strEQ;
	CharacterVector argNames = Rlist.attr("names");
	if (Rlist.size() != argNames.size()) {
		mxThrow("All list elements must be named");
	}

	std::vector<const char *> dataColNames;

	paramRows = -1;
	int pmatCols=-1;
	int mips = 1;
	int dataRows = 0;
	NumericVector Rmean;
	NumericMatrix Rcov;

	for (int ax=0; ax < Rlist.size(); ++ax) {
		const char *key = argNames[ax];
		RObject slotValue = Rlist[ax];
		if (strEQ(key, "spec")) {
			importSpec(as<List>(slotValue));
		} else if (strEQ(key, "param")) {
			if (!is<NumericVector>(slotValue)) mxThrow("'param' must be a numeric matrix of item parameters");
			NumericMatrix Rparam = as<NumericMatrix>(slotValue);
			param = Rparam.begin();
			paramRows = Rparam.nrow();
			pmatCols = Rparam.ncol();

			List dimnames = Rparam.attr("dimnames");
			if (dimnames.size() == 2) {
				if (dimnames[0] != R_NilValue) {
					CharacterVector names = dimnames[0];
					factorNames.resize(names.size());
					for (int nx=0; nx < names.size(); ++nx) {
						factorNames[nx] = names[nx];
					}
				}
				if (dimnames[1] != R_NilValue) {
					CharacterVector names = dimnames[1];
					itemNames.resize(names.size());
					for (int nx=0; nx < names.size(); ++nx) {
						itemNames[nx] = names[nx];
					}
				}
			}
		} else if (strEQ(key, "mean")) {
			Rmean = as<NumericVector>(slotValue);
			mean = Rmean.begin();
		} else if (strEQ(key, "cov")) {
			Rcov = as<NumericMatrix>(slotValue);
			cov = Rcov.begin();
		} else if (strEQ(key, "data")) {
			Rdata = as<DataFrame>(slotValue);
			IntegerVector col0 = Rdata[0];
			dataRows = col0.size();

			CharacterVector names = Rdata.attr("names");
			dataColNames.reserve(names.size());
			for (int nx=0; nx < names.size(); ++nx) {
				dataColNames.push_back(names[nx]);
			}
			dataRowNames = Rdata.attr("row.names");
		} else if (strEQ(key, "weightColumn")) {
			weightColumnName = as<const char *>(slotValue);
		} else if (strEQ(key, "freqColumn")) {
			freqColumnName = as<const char *>(slotValue);
		} else if (strEQ(key, "qwidth")) {
			qwidth = as<double>(slotValue);
		} else if (strEQ(key, "qpoints")) {
			qpoints = as<int>(slotValue);
		} else if (strEQ(key, "minItemsPerScore")) {
			mips = as<int>(slotValue);
		} else {
			// ignore
		}
	}

	learnMaxAbilities();

	if (itemDims < (int) factorNames.size())
		factorNames.resize(itemDims);

	if (int(factorNames.size()) < itemDims) {
		factorNames.reserve(itemDims);
		const int SMALLBUF = 24;
		char buf[SMALLBUF];
		while (int(factorNames.size()) < itemDims) {
			snprintf(buf, SMALLBUF, "s%d", int(factorNames.size()) + 1);
			factorNames.push_back(CHAR(Rf_mkChar(buf)));
		}
	}

	if (Rmean.size()) {
		if (Rmean.size() != itemDims) {
			mxThrow("mean must be a vector of length %d (not %d)", itemDims, Rmean.size());
		}

		verifyFactorNames(Rmean.attr("dimnames"), "mean");
	}

	if (Rcov.size()) {
		int nrow = Rcov.nrow();
		int ncol = Rcov.ncol();
		if (nrow != itemDims || ncol != itemDims) {
			mxThrow("cov must be %dx%d matrix", itemDims, itemDims);
		}

		verifyFactorNames(Rcov.attr("dimnames"), "cov");
	}

	setLatentDistribution(mean, cov);

	setMinItemsPerScore(mips);

	if (numItems() != pmatCols) {
		mxThrow("item matrix implies %d items but spec is length %d",
			 pmatCols, numItems());
	}

	if (Rdata.size()) {
		if (itemNames.size() == 0) mxThrow("Item matrix must have colnames");
		for (int ix=0; ix < numItems(); ++ix) {
			bool found=false;
			for (int dc=0; dc < int(dataColNames.size()); ++dc) {
				if (strEQ(itemNames[ix], dataColNames[dc])) {
					IntegerVector col = Rdata[dc];
					if (!Rf_isFactor(col)) {
						mxThrow("Column '%s' is an integer but "
								 "not an ordered factor",
								 dataColNames[dc]);
					}
					dataColumns.push_back(col.begin());
					found=true;
					break;
				}
			}
			if (!found) {
				mxThrow("Cannot find item '%s' in data", itemNames[ix]);
			}
		}
		if (weightColumnName) {
			for (int dc=0; dc < int(dataColNames.size()); ++dc) {
				if (strEQ(weightColumnName, dataColNames[dc])) {
					NumericVector col = Rdata[dc];
					rowWeight = col.begin();
					break;
				}
			}
			if (!rowWeight) {
				mxThrow("Cannot find weight column '%s'", weightColumnName);
			}
		}
		if (freqColumnName) {
			for (int dc=0; dc < int(dataColNames.size()); ++dc) {
				if (strEQ(freqColumnName, dataColNames[dc])) {
					IntegerVector col = Rdata[dc];
					rowFreq = col.begin();
					break;
				}
			}
			if (!rowFreq) {
				mxThrow("Cannot find frequency column '%s'", freqColumnName);
			}
		}
		rowMap.reserve(dataRows);
		for (int rx=0; rx < dataRows; ++rx) rowMap.push_back(rx);
	}

	Eigen::Map< Eigen::ArrayXXd > Eparam(param, paramRows, numItems());
	Eigen::Map< Eigen::VectorXd > meanVec(mean, itemDims);
	Eigen::Map< Eigen::MatrixXd > covMat(cov, itemDims, itemDims);

	quad.setStructure(qwidth, qpoints, Eparam, meanVec, covMat, twotier);
	quad.setupOutcomes(*this);

	if (paramRows < impliedParamRows) {
		mxThrow("At least %d rows are required in the item parameter matrix, only %d found",
			 impliedParamRows, paramRows);
	}

	quad.refresh(meanVec, covMat);
}

void ifaGroup::setLatentDistribution(double *_mean, double *_cov)
{
	if (!_mean) {
		mean = (double *) R_alloc(itemDims, sizeof(double));
		if (itemDims) memset(mean, 0, itemDims * sizeof(double));
	} else {
		mean = _mean;
	}

	if (!_cov) {
		cov = (double *) R_alloc(itemDims * itemDims, sizeof(double));
		Eigen::Map< Eigen::MatrixXd > covMat(cov, itemDims, itemDims);
		covMat.setIdentity();
	} else {
		cov = _cov;
	}
}

void ba81NormalQuad::layer::setupOutcomes(ifaGroup &ig)
{
	dataColumns.clear();
	dataColumns.reserve(numItems());
	totalOutcomes = 0;
	for (int ix=0; ix < numItems(); ++ix) {
		int outcomes = ig.itemOutcomes[ itemsMap[ix] ];
		itemOutcomes.push_back(outcomes);
		cumItemOutcomes.push_back(totalOutcomes);
		totalOutcomes += outcomes;
		dataColumns.push_back(ig.dataColumns[ itemsMap[ix] ]);
	}

	spec = ig.spec;
	paramRows = ig.paramRows;
}

void ba81NormalQuad::setupOutcomes(class ifaGroup &ig)
{ layers[0].setupOutcomes(ig); }

void ifaGroup::setMinItemsPerScore(int mips)
{
	if (numItems() && mips > numItems()) {
		mxThrow("minItemsPerScore (=%d) cannot be larger than the number of items (=%d)",
			 mips, numItems());
	}
	minItemsPerScore = mips;
}

void ifaGroup::buildRowMult()
{
	weightSum = 0.0;
	rowMult.resize(rowMap.size());
	for (int rx=0; rx < int(rowMap.size()); ++rx) {
		double mm = 1.0;
		if (rowWeight) mm *= rowWeight[rx];
		if (rowFreq) mm *= rowFreq[rx];
		weightSum += mm;
		rowMult[rx] = mm;
	}
}

void ifaGroup::buildRowSkip()
{
	rowSkip.assign(rowMap.size(), false);

	if (itemDims == 0) return;

	// Rows with no information about an ability will obtain the
	// prior distribution as an ability estimate. This will
	// throw off multigroup latent distribution estimates.
	for (size_t rx=0; rx < rowMap.size(); rx++) {
		bool hasNA = false;
		std::vector<int> contribution(itemDims);
		for (int ix=0; ix < numItems(); ix++) {
			int pick = dataColumn(ix)[ rowMap[rx] ];
			if (pick == NA_INTEGER) {
				hasNA = true;
				continue;
			}
			const double *ispec = spec[ix];
			int dims = ispec[RPF_ISpecDims];
			double *iparam = getItemParam(ix);
			for (int dx=0; dx < dims; dx++) {
				// assume factor loadings are the first item parameters
				if (iparam[dx] == 0) continue;
				contribution[dx] += 1;
			}
		}
		if (!hasNA) continue;
		if (minItemsPerScore == NA_INTEGER) {
			mxThrow("You have missing data. You must set minItemsPerScore");
		}
		for (int ax=0; ax < itemDims; ++ax) {
			if (contribution[ax] < minItemsPerScore) {
				// We could compute the other scores, but estimation of the
				// latent distribution is in the hot code path. We can reconsider
				// this choice when we try generating scores instead of the
				// score distribution.
				rowSkip[rx] = true;
			}
		}
	}
}

void ba81NormalQuad::layer::allocSummary(int numThreads)
{
	Dweight.resize(weightTableSize, numThreads);
	Dweight.setZero();
}

void ba81NormalQuad::allocSummary()
{
	if (numThreads < 1) mxThrow("numThreads < 1");
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].allocSummary(numThreads);
	}
}

void ba81NormalQuad::layer::prepSummary()
{
	for (int tx=1; tx < Dweight.cols(); ++tx) Dweight.col(0) += Dweight.col(tx);
}

void ba81NormalQuad::prepSummary()
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].prepSummary();
	}
}

void ba81NormalQuad::layer::allocBuffers(int numThreads)
{
	Qweight.resize(weightTableSize, numThreads);

	if (!numSpecific) return;

	thrEi.resize(totalPrimaryPoints, numThreads);
	thrEis.resize(totalPrimaryPoints * numSpecific, numThreads);
}

void ba81NormalQuad::allocBuffers()
{
	if (numThreads < 1) mxThrow("numThreads < 1");
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].allocBuffers(numThreads);
	}
}

void ba81NormalQuad::layer::releaseBuffers()
{
	Qweight.resize(0,0);
	thrEi.resize(0,0);
	thrEis.resize(0,0);
}

void ba81NormalQuad::releaseBuffers()
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].releaseBuffers();
	}
}

void ba81NormalQuad::releaseDerivCoefCache()
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].derivCoef.resize(0,0);
	}
}

void ifaGroup::setGridFineness(double width, int points)
{
	if (std::isfinite(width)) qwidth = width;
	if (points != NA_INTEGER) qpoints = points;
}

void ba81NormalQuad::prepExpectedTable()
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		for (int tx=1; tx < l1.expected.cols(); ++tx) {
			l1.expected.col(0) += l1.expected.col(tx);
		}
	}
}

void ba81NormalQuad::allocEstep()
{
	if (numThreads < 1) mxThrow("numThreads < 1");
	if (layers.size() != 1) mxThrow("layers.size() != 1");
	layer &l1 = layers[0];
	l1.expected.resize(l1.totalOutcomes * l1.totalQuadPoints, numThreads);
	l1.expected.setZero();
}

double ba81NormalQuad::mstepFit()
{
	double ll = 0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		ll += layers[lx].outcomeProbX.transpose().matrix() * layers[lx].expected.col(0).matrix();
	}
	return ll;
}

void ba81NormalQuad::releaseEstep()
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		l1.expected.resize(0,0);
	}
}

void ifaGroup::setFactorNames(std::vector<const char *> &names)
{
	if (int(names.size()) < itemDims) mxThrow("Not enough names");
	factorNames.resize(itemDims);
	for (int fx=0; fx < itemDims; ++fx) factorNames[fx] = names[fx];
}
