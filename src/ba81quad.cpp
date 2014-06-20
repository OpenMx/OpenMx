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
	One(1.0), quadGridSize(0)
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
