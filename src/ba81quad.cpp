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

#include "ba81quad.h"
#include "dmvnorm.h"

static inline void gramProduct(double *vec, size_t len, double *out)
{
	int cell = 0;
	for (size_t v1=0; v1 < len; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			out[cell] = vec[v1] * vec[v2];
			++cell;
		}
	}
}

struct sortAreaHelper {  // could be generalized with a template
	std::vector<double> &target;
	bool operator() (int i,int j) { return target[i] > target[j]; }
	sortAreaHelper(std::vector<double> &tgt) : target(tgt) {}
};

ba81NormalQuad::ba81NormalQuad() :
	quadGridSize(0), maxDims(-1), primaryDims(-1), numSpecific(-1),
	maxAbilities(-1)
{
	setOne(1);
}

void ba81NormalQuad::pointToWhere(const int *quad, double *where, int upto)
{
	for (int dx=0; dx < upto; dx++) {
		where[dx] = Qpoint[quad[dx]];
	}
}

void ba81NormalQuad::decodeLocation(int qx, const int dims, int *quad)
{
	for (int dx=dims-1; dx >= 0; --dx) {
		quad[dx] = qx % quadGridSize;
		qx = qx / quadGridSize;
	}
}

void ba81NormalQuad::setup0()
{
	quadGridSize = 1;
	numSpecific = 0;
	primaryDims = 0;
	maxDims = 1;
	maxAbilities = 0;
	totalQuadPoints = 1;
	totalPrimaryPoints = 1;
	weightTableSize = 1;
	Qpoint.clear();
	Qpoint.reserve(1);
	Qpoint.push_back(0);
	priQarea.clear();
	priQarea.push_back(One);
	wherePrep.clear();
	wherePrep.push_back(0);
}

void ba81NormalQuad::setup(double Qwidth, int Qpoints, double *means,
			   Eigen::MatrixXd &priCov, Eigen::VectorXd &sVar)
{
	quadGridSize = Qpoints;
	numSpecific = sVar.rows() * sVar.cols();
	primaryDims = priCov.rows();
	maxDims = primaryDims + (numSpecific? 1 : 0);
	maxAbilities = primaryDims + numSpecific;

	totalQuadPoints = 1;
	for (int dx=0; dx < maxDims; dx++) {
		totalQuadPoints *= quadGridSize;
	}

	if (int(Qpoint.size()) != quadGridSize) {
		Qpoint.clear();
		Qpoint.reserve(quadGridSize);
		double qgs = quadGridSize-1;
		for (int px=0; px < quadGridSize; ++px) {
			Qpoint.push_back(px * 2 * Qwidth / qgs - Qwidth);
		}
	}

	std::vector<double> tmpWherePrep(totalQuadPoints * maxDims);

	Eigen::VectorXi quad(maxDims);
	for (int qx=0; qx < totalQuadPoints; qx++) {
		double *wh = tmpWherePrep.data() + qx * maxDims;
		decodeLocation(qx, maxDims, quad.data());
		pointToWhere(quad.data(), wh, maxDims);
	}

	totalPrimaryPoints = totalQuadPoints;
	weightTableSize = totalQuadPoints;

	if (numSpecific) {
		totalPrimaryPoints /= quadGridSize;
		speQarea.resize(quadGridSize * numSpecific);
		weightTableSize *= numSpecific;
	}

	std::vector<double> tmpPriQarea;
	tmpPriQarea.reserve(totalPrimaryPoints);
	{
		Eigen::VectorXd where(primaryDims);
		for (int qx=0; qx < totalPrimaryPoints; qx++) {
			decodeLocation(qx, primaryDims, quad.data());
			pointToWhere(quad.data(), where.data(), primaryDims);
			double den = exp(dmvnorm(primaryDims, where.data(), means, priCov.data()));
			tmpPriQarea.push_back(den);
		}
	}

	std::vector<int> priOrder;
	priOrder.reserve(totalPrimaryPoints);
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		priOrder.push_back(qx);
	}
	if (0) {
		sortAreaHelper priCmp(tmpPriQarea);
		std::sort(priOrder.begin(), priOrder.end(), priCmp);
	}

	priQarea.clear();
	priQarea.reserve(totalPrimaryPoints);

	double totalArea = 0;
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		double den = tmpPriQarea[priOrder[qx]];
		priQarea.push_back(den);
		//double prevTotalArea = totalArea;
		totalArea += den;
		// if (totalArea == prevTotalArea) {
		// 	mxLog("%.4g / %.4g = %.4g", den, totalArea, den / totalArea);
		// }
	}

	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		priQarea[qx] *= One;
		priQarea[qx] /= totalArea;
		//mxLog("%.5g,", priQarea[qx]);
	}

	for (int sgroup=0; sgroup < numSpecific; sgroup++) {
		totalArea = 0;
		double mean = means[primaryDims + sgroup];
		double var = sVar(sgroup);
		for (int qx=0; qx < quadGridSize; qx++) {
			double den = exp(dmvnorm(1, &Qpoint[qx], &mean, &var));
			speQarea[sIndex(sgroup, qx)] = den;
			totalArea += den;
		}
		for (int qx=0; qx < quadGridSize; qx++) {
			speQarea[sIndex(sgroup, qx)] /= totalArea;
		}
	}
	//pda(speQarea.data(), numSpecific, quadGridSize);

	for (int sx=0; sx < int(speQarea.size()); ++sx) {
		speQarea[sx] *= One;
	}
	//pda(speQarea.data(), numSpecific, quadGridSize);

	wherePrep.clear();
	wherePrep.reserve(totalQuadPoints * maxDims);

	if (numSpecific == 0) {
		for (int qx=0; qx < totalPrimaryPoints; qx++) {
			int sortq = priOrder[qx] * maxDims;
			for (int dx=0; dx < maxDims; ++dx) {
				wherePrep.push_back(tmpWherePrep[sortq + dx]);
			}
		}
	} else {
		for (int qx=0; qx < totalPrimaryPoints; ++qx) {
			int sortq = priOrder[qx] * quadGridSize;
			for (int sx=0; sx < quadGridSize; ++sx) {
				int base = (sortq + sx) * maxDims;
				for (int dx=0; dx < maxDims; ++dx) {
					wherePrep.push_back(tmpWherePrep[base + dx]);
				}
			}
		}
	}

	// recompute whereGram because the order might have changed
	whereGram.resize(triangleLoc1(maxDims), totalQuadPoints);

	for (int qx=0; qx < totalQuadPoints; qx++) {
		double *wh = wherePrep.data() + qx * maxDims;
		gramProduct(wh, maxDims, &whereGram.coeffRef(0, qx));
	}
}

ifaGroup::ifaGroup(int cores, bool _twotier) : Rdata(NULL),
		numThreads(cores), qwidth(6.0), qpoints(49),
		twotier(_twotier),
		maxAbilities(0),
		numSpecific(0),
		mean(0),
					       cov(0), dataRowNames(0),
	    weightColumnName(0), rowWeight(0),
					       minItemsPerScore(NA_INTEGER),
					       outcomeProb(0), excludedPatterns(-1)
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
	Free(outcomeProb);
}

void ifaGroup::importSpec(SEXP slotValue)
{
	for (int sx=0; sx < Rf_length(slotValue); ++sx) {
		SEXP model = VECTOR_ELT(slotValue, sx);
		if (!OBJECT(model)) {
			Rf_error("Item models must inherit rpf.base");
		}
		SEXP Rspec;
		ScopedProtect p1(Rspec, R_do_slot(model, Rf_install("spec")));
		spec.push_back(REAL(Rspec));
	}

	dataColumns.reserve(spec.size());
	itemOutcomes.reserve(spec.size());
	cumItemOutcomes.reserve(spec.size());

	paramRows = 0;
	totalOutcomes = 0;
	maxItemDims = 0;
	for (int cx = 0; cx < numItems(); cx++) {
		const double *ispec = spec[cx];
		int id = ispec[RPF_ISpecID];
		int dims = ispec[RPF_ISpecDims];
		if (maxItemDims < dims)
			maxItemDims = dims;
		int no = ispec[RPF_ISpecOutcomes];
		itemOutcomes.push_back(no);
		cumItemOutcomes.push_back(totalOutcomes);
		totalOutcomes += no;

		int numParam = (*librpf_model[id].numParam)(ispec);
		if (paramRows < numParam)
			paramRows = numParam;
	}
}

void ifaGroup::verifyFactorNames(SEXP mat, const char *matName)
{
	static const char *dimname[] = { "row", "col" };

	SEXP dimnames;
	Rf_protect(dimnames = Rf_getAttrib(mat, R_DimNamesSymbol));
	if (!Rf_isNull(dimnames) && Rf_length(dimnames) == 2) {
		for (int dx=0; dx < 2; ++dx) {
			SEXP names;
			Rf_protect(names = VECTOR_ELT(dimnames, dx));
			if (!Rf_length(names)) continue;
			if (int(factorNames.size()) != Rf_length(names)) {
				Rf_error("%s %snames must be length %d",
					 matName, dimname[dx], (int) factorNames.size());
			}
			int nlen = Rf_length(names);
			for (int nx=0; nx < nlen; ++nx) {
				const char *name = CHAR(STRING_ELT(names, nx));
				if (strEQ(factorNames[nx], name)) continue;
				Rf_error("%s %snames[%d] is '%s', does not match factor name '%s'",
					 matName, dimname[dx], 1+nx, name, factorNames[nx]);
			}
		}
	}
}

void ifaGroup::learnMaxAbilities()
{
	maxAbilities = 0;
	Eigen::ArrayXi loadings(maxItemDims);
	loadings.setZero();
	for (int cx = 0; cx < numItems(); cx++) {
		for (int dx=0; dx < maxItemDims; ++dx) {
			if (getItemParam(cx)[dx] != 0) loadings[dx] += 1;
		}
	}
	maxAbilities = (loadings != 0).count();
	if (maxItemDims != maxAbilities) {
		for (int lx=0; lx < maxItemDims; ++lx) {
			if (loadings[lx] == 0) Rf_error("Factor %d does not load on any items", 1+lx);
		}
	}
}

void ifaGroup::import(SEXP Rlist)
{
	SEXP argNames;
	Rf_protect(argNames = Rf_getAttrib(Rlist, R_NamesSymbol));
	if (Rf_length(Rlist) != Rf_length(argNames)) {
		Rf_error("All list elements must be named");
	}

	std::vector<const char *> dataColNames;

	int pmatRows=-1, pmatCols=-1;
	int mips = 1;
	int dataRows = 0;
	SEXP Rmean=0, Rcov=0;

	for (int ax=0; ax < Rf_length(Rlist); ++ax) {
		const char *key = R_CHAR(STRING_ELT(argNames, ax));
		SEXP slotValue = VECTOR_ELT(Rlist, ax);
		if (strEQ(key, "spec")) {
			importSpec(slotValue);
		} else if (strEQ(key, "param")) {
			if (!Rf_isReal(slotValue)) Rf_error("'param' must be a numeric matrix of item parameters");
			param = REAL(slotValue);
			getMatrixDims(slotValue, &pmatRows, &pmatCols);

			SEXP dimnames;
			Rf_protect(dimnames = Rf_getAttrib(slotValue, R_DimNamesSymbol));
			if (!Rf_isNull(dimnames) && Rf_length(dimnames) == 2) {
				SEXP names;
				Rf_protect(names = VECTOR_ELT(dimnames, 0));
				int nlen = Rf_length(names);
				factorNames.resize(nlen);
				for (int nx=0; nx < nlen; ++nx) {
					factorNames[nx] = CHAR(STRING_ELT(names, nx));
				}
				Rf_protect(names = VECTOR_ELT(dimnames, 1));
				nlen = Rf_length(names);
				itemNames.resize(nlen);
				for (int nx=0; nx < nlen; ++nx) {
					itemNames[nx] = CHAR(STRING_ELT(names, nx));
				}
			}
		} else if (strEQ(key, "mean")) {
			Rmean = slotValue;
			if (!Rf_isReal(slotValue)) Rf_error("'mean' must be a numeric vector or matrix");
			mean = REAL(slotValue);
		} else if (strEQ(key, "cov")) {
			Rcov = slotValue;
			if (!Rf_isReal(slotValue)) Rf_error("'cov' must be a numeric matrix");
			cov = REAL(slotValue);
		} else if (strEQ(key, "data")) {
			Rdata = slotValue;
			dataRows = Rf_length(VECTOR_ELT(Rdata, 0));

			SEXP names;
			Rf_protect(names = Rf_getAttrib(Rdata, R_NamesSymbol));
			int nlen = Rf_length(names);
			dataColNames.reserve(nlen);
			for (int nx=0; nx < nlen; ++nx) {
				dataColNames.push_back(CHAR(STRING_ELT(names, nx)));
			}
			Rf_protect(dataRowNames = Rf_getAttrib(Rdata, R_RowNamesSymbol));
		} else if (strEQ(key, "weightColumn")) {
			if (Rf_length(slotValue) != 1) {
				Rf_error("You can only have one weightColumn");
			}
			weightColumnName = CHAR(STRING_ELT(slotValue, 0));
		} else if (strEQ(key, "qwidth")) {
			qwidth = Rf_asReal(slotValue);
		} else if (strEQ(key, "qpoints")) {
			qpoints = Rf_asInteger(slotValue);
		} else if (strEQ(key, "minItemsPerScore")) {
			mips = Rf_asInteger(slotValue);
		} else {
			// ignore
		}
	}

	learnMaxAbilities();

	if (maxAbilities < (int) factorNames.size())
		factorNames.resize(maxAbilities);

	if (!factorNames.size()) {
		factorNames.reserve(maxAbilities);
		const int SMALLBUF = 10;
		char buf[SMALLBUF];
		for (int sx=0; sx < maxAbilities; ++sx) {
			snprintf(buf, SMALLBUF, "s%d", sx+1);
			factorNames.push_back(CHAR(Rf_mkChar(buf)));
		}
	}

	if (Rmean) {
		if (Rf_isMatrix(Rmean)) {
			int nrow, ncol;
			getMatrixDims(Rmean, &nrow, &ncol);
			if (!(nrow * ncol == maxAbilities && (nrow==1 || ncol==1))) {
				Rf_error("mean must be a column or row matrix of length %d", maxAbilities);
			}
		} else {
			if (Rf_length(Rmean) != maxAbilities) {
				Rf_error("mean must be a vector of length %d", maxAbilities);
			}
		}

		verifyFactorNames(Rmean, "mean");
	}

	if (Rcov) {
		if (Rf_isMatrix(Rcov)) {
			int nrow, ncol;
			getMatrixDims(Rcov, &nrow, &ncol);
			if (nrow != maxAbilities || ncol != maxAbilities) {
				Rf_error("cov must be %dx%d matrix", maxAbilities, maxAbilities);
			}
		} else {
			if (Rf_length(Rcov) != 1) {
				Rf_error("cov must be %dx%d matrix", maxAbilities, maxAbilities);
			}
		}

		verifyFactorNames(Rcov, "cov");
	}

	setLatentDistribution(maxAbilities, mean, cov);

	setMinItemsPerScore(mips);

	if (numItems() != pmatCols) {
		Rf_error("item matrix implies %d items but spec is length %d",
			 pmatCols, numItems());
	}

	if (Rdata) {
		if (itemNames.size() == 0) Rf_error("Item matrix must have colnames");
		for (int ix=0; ix < numItems(); ++ix) {
			bool found=false;
			for (int dc=0; dc < int(dataColNames.size()); ++dc) {
				if (strEQ(itemNames[ix], dataColNames[dc])) {
					dataColumns.push_back(INTEGER(VECTOR_ELT(Rdata, dc)));
					found=true;
					break;
				}
			}
			if (!found) {
				Rf_error("Cannot find item '%s' in data", itemNames[ix]);
			}
		}
		if (weightColumnName) {
			for (int dc=0; dc < int(dataColNames.size()); ++dc) {
				if (strEQ(weightColumnName, dataColNames[dc])) {
					rowWeight = REAL(VECTOR_ELT(Rdata, dc));
					break;
				}
			}
			if (!rowWeight) {
				Rf_error("Cannot find weight column '%s'", weightColumnName);
			}
		}
		rowMap.reserve(dataRows);
		for (int rx=0; rx < dataRows; ++rx) rowMap.push_back(rx);
	}

	detectTwoTier();
	sanityCheck();

	if (pmatRows < paramRows) {
		Rf_error("At least %d rows are required in the item parameter matrix, only %d found",
			 paramRows, pmatRows);
	}

	Eigen::Map<Eigen::MatrixXd> fullCov(cov, maxAbilities, maxAbilities);
	int dense = maxAbilities - numSpecific;
	Eigen::MatrixXd priCov = fullCov.block(0, 0, dense, dense);
	Eigen::VectorXd sVar = fullCov.diagonal().tail(numSpecific);

	quad.setup(qwidth, qpoints, mean, priCov, sVar);
}

void ifaGroup::setLatentDistribution(int dims, double *_mean, double *_cov)
{
	maxAbilities = dims;
	if (maxAbilities < 0) Rf_error("maxAbilities must be non-negative");

	if (!mean) {
		mean = (double *) R_alloc(maxAbilities, sizeof(double));
		memset(mean, 0, maxAbilities * sizeof(double));
	} else {
		mean = _mean;
	}

	if (!cov) {
		cov = (double *) R_alloc(maxAbilities * maxAbilities, sizeof(double));
		Eigen::Map< Eigen::MatrixXd > covMat(cov, maxAbilities, maxAbilities);
		covMat.setIdentity();
	} else {
		cov = _cov;
	}
}

void ifaGroup::detectTwoTier()
{
	int mlen = maxAbilities;

	if (!twotier || mlen < 3) return;

	std::vector<int> orthogonal;
	if (mlen >= 3) {
		Eigen::Map<Eigen::MatrixXd> Ecov(cov, mlen, mlen);
		Eigen::Matrix<Eigen::DenseIndex, Eigen::Dynamic, 1> numCov((Ecov.array() != 0.0).matrix().colwise().count());
		std::vector<int> candidate;
		for (int fx=0; fx < numCov.rows(); ++fx) {
			if (numCov(fx) == 1) candidate.push_back(fx);
		}
		if (candidate.size() > 1) {
			std::vector<bool> mask(numItems());
			for (int cx=candidate.size() - 1; cx >= 0; --cx) {
				std::vector<bool> loading(numItems());
				for (int ix=0; ix < numItems(); ++ix) {
					loading[ix] = param[ix * paramRows + candidate[cx]] != 0;
				}
				std::vector<bool> overlap(loading.size());
				std::transform(loading.begin(), loading.end(),
					       mask.begin(), overlap.begin(),
					       std::logical_and<bool>());
				if (std::find(overlap.begin(), overlap.end(), true) == overlap.end()) {
					std::transform(loading.begin(), loading.end(),
						       mask.begin(), mask.begin(),
						       std::logical_or<bool>());
					orthogonal.push_back(candidate[cx]);
				}
			}
		}
		std::reverse(orthogonal.begin(), orthogonal.end());
	}
	if (orthogonal.size() == 1) orthogonal.clear();
	if (orthogonal.size() && orthogonal[0] != mlen - int(orthogonal.size())) {
		Rf_error("Independent factors must be given after dense factors");
	}

	numSpecific = orthogonal.size();

	if (numSpecific) {
		Sgroup.assign(numItems(), 0);
		for (int ix=0; ix < numItems(); ix++) {
			for (int dx=orthogonal[0]; dx < maxAbilities; ++dx) {
				if (param[ix * paramRows + dx] != 0) {
					Sgroup[ix] = dx - orthogonal[0];
					continue;
				}
			}
		}
	}
}

void ifaGroup::setMinItemsPerScore(int mips)
{
	if (numItems() && mips > numItems()) {
		Rf_error("minItemsPerScore (=%d) cannot be larger than the number of items (=%d)",
			 mips, numItems());
	}
	minItemsPerScore = mips;
}

void ifaGroup::buildRowSkip()
{
	rowSkip.assign(rowMap.size(), false);

	if (maxAbilities == 0) return;

	// Rows with no information about an ability will obtain the
	// prior distribution as an ability estimate. This will
	// throw off multigroup latent distribution estimates.
	for (size_t rx=0; rx < rowMap.size(); rx++) {
		bool hasNA = false;
		std::vector<int> contribution(maxAbilities);
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
			Rf_error("You have missing data. You must set minItemsPerScore");
		}
		for (int ax=0; ax < maxAbilities; ++ax) {
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

void ifaGroup::sanityCheck() // remove TODO
{
}

void ifaGroup::setGridFineness(double width, int points)
{
	if (std::isfinite(width)) qwidth = width;
	if (points != NA_INTEGER) qpoints = points;
}
