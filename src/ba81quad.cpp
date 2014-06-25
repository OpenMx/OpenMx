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

static inline void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    Rf_protect(matrixDims = Rf_getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    Rf_unprotect(1);
}

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
	One(1.0), quadGridSize(0), maxDims(-1), primaryDims(-1), numSpecific(-1),
	maxAbilities(-1)
{
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
			Qpoint.push_back(Qwidth - px * 2 * Qwidth / qgs);
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
	if (numSpecific) {
		totalPrimaryPoints /= quadGridSize;
		speQarea.resize(quadGridSize * numSpecific);
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
	sortAreaHelper priCmp(tmpPriQarea);
	std::sort(priOrder.begin(), priOrder.end(), priCmp);

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

void ba81NormalQuad::mapDenseSpace(double piece, const double *where,
				   const double *whereGram, double *latentDist)
{
	const int pmax = primaryDims;
	int gx = 0;
	int cx = maxAbilities;
	for (int d1=0; d1 < pmax; d1++) {
		double piece_w1 = piece * where[d1];
		latentDist[d1] += piece_w1;
		for (int d2=0; d2 <= d1; d2++) {
			double piece_cov = piece * whereGram[gx];
			latentDist[cx] += piece_cov;
			++cx; ++gx;
		}
	}
}

void ba81NormalQuad::mapSpecificSpace(int sgroup, double piece, const double *where,
				      const double *whereGram, double *latentDist)
{
	const int pmax = primaryDims;

	int sdim = pmax + sgroup;
	double piece_w1 = piece * where[pmax];
	latentDist[sdim] += piece_w1;

	double piece_var = piece * whereGram[triangleLoc0(pmax)];
	int to = maxAbilities + triangleLoc0(sdim);
	latentDist[to] += piece_var;
}

void ba81NormalQuad::EAP(double *thrDweight, double scalingFactor, double *scorePad)
{
	if (numSpecific == 0) { // use template to handle this branch at compile time? TODO
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			mapDenseSpace(thrDweight[qx], &wherePrep[qx * maxDims],
				      &whereGram.coeffRef(0, qx), scorePad);
		}
	} else {
		int qloc=0;
		for (int qx=0; qx < totalQuadPoints; qx++) {
			const double *whPrep = &wherePrep[qx * maxDims];
			const double *whGram = &whereGram.coeffRef(0, qx);
			mapDenseSpace(thrDweight[qloc], whPrep, whGram, scorePad);
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapSpecificSpace(Sgroup, thrDweight[qloc], whPrep, whGram, scorePad);
				++qloc;
			}
		}
	}

	const int padSize = maxAbilities + triangleLoc1(maxAbilities);
	for (int d1=0; d1 < padSize; d1++) {
		scorePad[d1] *= scalingFactor;
	}

	int cx = maxAbilities;
	for (int a1=0; a1 < primaryDims; ++a1) {
		for (int a2=0; a2 <= a1; ++a2) {
			double ma1 = scorePad[a1];
			double ma2 = scorePad[a2];
			scorePad[cx] -= ma1 * ma2;
			++cx;
		}
	}
	for (int sx=0; sx < numSpecific; sx++) {
		int sdim = primaryDims + sx;
		double ma1 = scorePad[sdim];
		scorePad[maxAbilities + triangleLoc0(sdim)] -= ma1 * ma1;
	}
}

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
		Rf_protect(Rspec = R_do_slot(model, Rf_install("spec")));
		spec.push_back(REAL(Rspec));
	}

	dataColumns.reserve(spec.size());
	itemOutcomes.reserve(spec.size());
	cumItemOutcomes.reserve(spec.size());

	paramRows = 0;
	totalOutcomes = 0;
	for (int cx = 0; cx < numItems(); cx++) {
		const double *ispec = spec[cx];
		int id = ispec[RPF_ISpecID];
		int no = ispec[RPF_ISpecOutcomes];
		itemOutcomes.push_back(no);
		cumItemOutcomes.push_back(totalOutcomes);
		totalOutcomes += no;

		int numParam = (*librpf_model[id].numParam)(ispec);
		if (paramRows < numParam)
			paramRows = numParam;
	}
}

void ifaGroup::import(SEXP Rlist)
{
	SEXP argNames;
	Rf_protect(argNames = Rf_getAttrib(Rlist, R_NamesSymbol));

	std::vector<const char *> dataColNames;

	int mlen = 0;
	int nrow=0, ncol=0; // cov size

	int pmatRows=-1, pmatCols=-1;

	for (int ax=0; ax < Rf_length(Rlist); ++ax) {
		const char *key = R_CHAR(STRING_ELT(argNames, ax));
		SEXP slotValue = VECTOR_ELT(Rlist, ax);
		if (strEQ(key, "spec")) {
			importSpec(slotValue);
		} else if (strEQ(key, "param")) {
			param = REAL(slotValue);
			getMatrixDims(slotValue, &pmatRows, &pmatCols);

			SEXP dimnames;
			Rf_protect(dimnames = Rf_getAttrib(slotValue, R_DimNamesSymbol));
			if (!Rf_isNull(dimnames) && Rf_length(dimnames) == 2) {
				SEXP names;
				Rf_protect(names = VECTOR_ELT(dimnames, 1));
				int nlen = Rf_length(names);
				itemNames.resize(nlen);
				for (int nx=0; nx < nlen; ++nx) {
					itemNames[nx] = CHAR(STRING_ELT(names, nx));
				}
			}
		} else if (strEQ(key, "mean")) {
			mlen = Rf_length(slotValue);
			mean = REAL(slotValue);
		} else if (strEQ(key, "cov")) {
			getMatrixDims(slotValue, &nrow, &ncol);
			if (nrow != ncol) Rf_error("cov must be a square matrix (not %dx%d)", nrow, ncol);
			cov = REAL(slotValue);

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
			}
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
		} else {
			// ignore
		}
	}
	if (mlen != nrow) Rf_error("Mean length %d does not match cov size %d", mlen, nrow);

	if (numItems() != pmatCols) {
		Rf_error("param implies %d items but spec is length %d",
			 pmatCols, numItems());
	}

	if (Rdata) {
		if (itemNames.size() == 0) Rf_error("Item parameter matrix must have colnames");
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
	}

	maxAbilities = mlen;
	detectTwoTier();

	if (pmatRows < paramRows) {
		Rf_error("At least %d rows are required in the item parameter matrix, only %d found",
			 paramRows, pmatRows);
	}
}

void ifaGroup::setLatentDistribution(int dims, double *_mean, double *_cov)
{
	maxAbilities = dims;
	mean = _mean;
	cov = _cov;
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

void ifaGroup::sanityCheck()
{
	for (int ix=0; ix < numItems(); ++ix) {
		const int dims = spec[ix][RPF_ISpecDims];

		int loadings = 0;
		for (int dx=0; dx < dims; ++dx) {
			if (getItemParam(ix)[dx] != 0) loadings += 1;
		}
		if (loadings > maxAbilities) {
			omxRaiseErrorf("Item %d has more factor loadings (%d) than there are factors (%d)",
				       1+ix, loadings, maxAbilities);
		}
	}
}

double ifaGroup::area(int qx, int ix)
{
	if (numSpecific == 0) {
		return quad.priQarea[qx];
	} else {
		int px = qx / quad.quadGridSize;
		int sx = qx % quad.quadGridSize;
		return quad.priQarea[px] * quad.speQarea[sx * quad.numSpecific + Sgroup[ix]];
	}
}

// Depends on item parameters, but not latent distribution
void ifaGroup::ba81OutcomeProb(double *param, bool wantLog)
{
	const int maxDims = quad.maxDims;
	outcomeProb = Realloc(outcomeProb, totalOutcomes * quad.totalQuadPoints, double);

#pragma omp parallel for num_threads(Global->numThreads)
	for (int ix=0; ix < numItems(); ix++) {
		double *qProb = outcomeProb + cumItemOutcomes[ix] * quad.totalQuadPoints;
		const double *ispec = spec[ix];
		int id = ispec[RPF_ISpecID];
		int dims = ispec[RPF_ISpecDims];
		Eigen::VectorXd ptheta(dims);
		double *iparam = param + paramRows * ix;
		rpf_prob_t prob_fn = wantLog? librpf_model[id].logprob : librpf_model[id].prob;

		for (int qx=0; qx < quad.totalQuadPoints; qx++) {
			double *where = quad.wherePrep.data() + qx * maxDims;
			for (int dx=0; dx < dims; dx++) {
				ptheta[dx] = where[std::min(dx, maxDims-1)];
			}

			(*prob_fn)(ispec, iparam, ptheta.data(), qProb);
			qProb += itemOutcomes[ix];
		}
	}
}

void ifaGroup::ba81LikelihoodSlow2(const int px, double *out)
{
	const int totalQuadPoints = quad.totalQuadPoints;
	double *oProb = outcomeProb;
	std::vector<double> &priQarea = quad.priQarea;

	for (int qx=0; qx < totalQuadPoints; ++qx) {
		out[qx] = priQarea[qx];
	}

	const int row = rowMap[px];
	for (int ix=0; ix < numItems(); ix++) {
		int pick = dataColumns[ix][row];
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		for (int qx=0; qx < totalQuadPoints; ++qx) {
			out[qx] *= oProb[pick];
			oProb += itemOutcomes[ix];
		}
	}
}

void ifaGroup::cai2010EiEis(const int px, double *lxk, double *Eis, double *Ei)
{
	double *oProb = outcomeProb;
	const int totalQuadPoints = quad.totalQuadPoints;
	const int totalPrimaryPoints = quad.totalPrimaryPoints;
	const int specificPoints = quad.quadGridSize;
	std::vector<double> &speQarea = quad.speQarea;
	std::vector<double> &priQarea = quad.priQarea;

	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sx=0; sx < specificPoints * numSpecific; sx++) {
			lxk[qloc] = speQarea[sx];
			++qloc;
		}
	}

	const int row = rowMap[px];
	for (int ix=0; ix < numItems(); ix++) {
		int pick = dataColumns[ix][row];
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;
		int Sgroup1 = Sgroup[ix];
		double *out1 = lxk;
		for (int qx=0; qx < quad.totalQuadPoints; qx++) {
			out1[Sgroup1] *= oProb[pick];
			oProb += itemOutcomes[ix];
			out1 += numSpecific;
		}
	}

	for (int qx=0; qx < totalPrimaryPoints * numSpecific; ++qx) Eis[qx] = 0;
	for (int qx=0; qx < totalPrimaryPoints; ++qx) Ei[qx] = priQarea[qx];

	int eisloc = 0;
	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sx=0; sx < specificPoints; sx++) {
			for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
				double piece = lxk[qloc];
				Eis[eisloc + sgroup] += piece;
				++qloc;
			}
		}
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Ei[qx] *= Eis[eisloc + sgroup] * quad.getReciprocalOfOne();
		}
		eisloc += numSpecific;
	}
}
