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

void ba81NormalQuad::setup(double Qwidth, int Qpoints, double *means,
			   Eigen::MatrixXd &priCov, Eigen::VectorXd &sVar)
{
	quadGridSize = Qpoints;
	numSpecific = sVar.rows() * sVar.cols();
	int priDims = priCov.rows();
	int maxDims = priDims + (numSpecific? 1 : 0);

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

	for (int qx=0; qx < totalQuadPoints; qx++) {
		double *wh = tmpWherePrep.data() + qx * maxDims;
		int quad[maxDims];
		decodeLocation(qx, maxDims, quad);
		pointToWhere(quad, wh, maxDims);
	}

	totalPrimaryPoints = totalQuadPoints;
	if (numSpecific) {
		totalPrimaryPoints /= quadGridSize;
		speQarea.resize(quadGridSize * numSpecific);
	}

	std::vector<double> tmpPriQarea;
	tmpPriQarea.reserve(totalPrimaryPoints);
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, quad);
		double where[priDims];
		pointToWhere(quad, where, priDims);
		double den = exp(dmvnorm(priDims, where, means, priCov.data()));
		tmpPriQarea.push_back(den);
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
		double mean = means[priDims + sgroup];
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
}

