/*
  Copyright 2012-2014, 2016 Joshua Nathaniel Pritikin and contributors

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

#ifndef _BA81QUAD_H_
#define _BA81QUAD_H_

#include "glue.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues> 
#include "libifa-rpf.h"
#include "dmvnorm.h"

extern const struct rpf *Glibrpf_model;
extern int Glibrpf_numModels;

inline void gramProduct(double *vec, size_t len, double *out)
{
	int cell = 0;
	for (size_t v1=0; v1 < len; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			out[cell] = vec[v1] * vec[v2];
			++cell;
		}
	}
}

class ba81NormalQuad {
 private:
	template <typename T1>
	static inline void decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out, int dims);

	double width;
	int gridSize;
	std::vector<double> Qpoint;           // gridSize

	// Currently, there is always only one layer.
	class layer {
	private:
		template <typename T1, typename T2, typename T3, typename T4>
		void calcDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
				   Eigen::MatrixBase<T3> &icov, Eigen::MatrixBase<T4> &where, int qx);
		template <typename T1, typename T2, typename T4>
		void calcDerivCoef1(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
				    Eigen::MatrixBase<T4> &where, int qx, int sgroup);
		template <typename T2>
		void mapLatentDeriv(double piece, int qx, Eigen::ArrayBase<T2> &derivOut);
		template <typename T2>
		void mapLatentDerivS(int sgroup, double piece, int qx, int curGroup,
				     Eigen::ArrayBase<T2> &derivOut);
	public:

		class ba81NormalQuad *quad;

		// subset
		std::vector<bool> abilitiesMask;      // true=include, false=exclude
		std::vector<int> abilitiesMap;        // local to global index
		std::vector<bool> itemsMask;          // true=include, false=exclude
		std::vector<int> itemsMap;            // local to global index
		std::vector<int> glItemsMap;          // global to local index
		int numItems() const { return (int) itemsMap.size(); };
		std::vector<int> itemOutcomes;
		std::vector<int> cumItemOutcomes;
		int totalOutcomes;
		std::vector<const int*> dataColumns;

		// quadrature
		int maxDims;                          // integration dimension after dimension reduction
		int totalQuadPoints;                  // gridSize ^ maxDims
		int weightTableSize;                  // dense: totalQuadPoints; 2tier: totalQuadPoints * numSpecific
		Eigen::ArrayXd outcomeProbX;          // outcomes (within item) * totalQuadPoints * item
		Eigen::ArrayXXd expected;             // outcomes (within item) * totalQuadPoints * item * numThreads
		std::vector<double> priQarea;         // totalPrimaryPoints
		Eigen::ArrayXXd Qweight;              // weightTableSize * numThreads (for 1 response pattern)
		Eigen::ArrayXXd Dweight;              // weightTableSize * numThreads (for entire dataset)

		// Cai (2010) dimension reduction
		int numSpecific;
		int primaryDims;
		int totalPrimaryPoints;               // totalQuadPoints except for specific dim
		std::vector<int> Sgroup;              // item's specific group 0..numSpecific-1
		std::vector<double> speQarea;         // gridSize * numSpecific
		Eigen::ArrayXXd thrEi;
		Eigen::ArrayXXd thrEis;

		// Mislevy (1984) deriv stuff
		Eigen::ArrayXXd derivCoef;            // deriv * totalQuadPoints

		layer(class ba81NormalQuad *quad)
			: quad(quad), maxDims(-1),
			totalQuadPoints(-1), weightTableSize(-1),
			numSpecific(-1), primaryDims(-1), totalPrimaryPoints(-1) {};
		int numAbil() const { return (int) abilitiesMap.size(); }
		inline int sIndex(int sx, int qx) {
			//if (sx < 0 || sx >= state->numSpecific) Rf_error("Out of domain");
			//if (qx < 0 || qx >= state->gridSize) Rf_error("Out of domain");
			return qx * numSpecific + sx;
		};
		template <typename T1>
		inline void mapDenseSpace(double piece, const double *where,
					  const double *whereGram, Eigen::ArrayBase<T1> &latentDist);
		template <typename T1>
		inline void mapSpecificSpace(int sgroup, double piece, const double *where,
					     const double *whereGram, Eigen::ArrayBase<T1> &latentDist);
		template <typename T1>
		inline void finalizeLatentDist(const double sampleSize, Eigen::ArrayBase<T1> &scorePad);

		void allocBuffers(int numThreads);
		void releaseBuffers();
		inline double computePatternLik(int thrId, int row);
		inline void prepLatentDist(int thrId);
		template <typename T1, typename T2, typename T3>
		void detectTwoTier(Eigen::ArrayBase<T1> &param,
				   Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
		template <typename T1, typename T2, typename T3, typename T4>
		void globalToLocalDist(Eigen::MatrixBase<T1> &gmean, Eigen::MatrixBase<T2> &gcov,
				       Eigen::MatrixBase<T3> &mean, Eigen::MatrixBase<T4> &cov);
		template <typename T1, typename T2, typename T3>
		void setStructure(Eigen::ArrayBase<T1> &param,
				  Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
		template <typename T1, typename T2>
		inline void refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov);
		inline void addToExpected(int thrId, int px);
		template <typename T1, typename T2, typename T3>
		void mstepIter(int ix, Eigen::MatrixBase<T1> &abx,
			       Eigen::MatrixBase<T2> &abscissa, T3 &op);
		inline void weightBy(int thrId, double weight);
		inline void weightByAndSummarize(int thrId, double weight);
		template <typename T1, typename T2>
		int cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov);
		template <typename T1, typename T2>
		void pointToGlobalAbscissa(int qx, Eigen::MatrixBase<T1> &abx,
					   Eigen::MatrixBase<T2> &abscissa);
		template <typename T1, typename T2>
		void pointToLocalAbscissa(int qx, Eigen::MatrixBase<T1> &abx,
					  Eigen::MatrixBase<T2> &abscissa);
		template <typename T, typename T1, typename T2, typename T3>
		void computeRowDeriv(int thrId, Eigen::MatrixBase<T2> &abx,
				     Eigen::MatrixBase<T3> &abscissa, T &op,
				     bool freeLatents, Eigen::ArrayBase<T1> &latentGradOut);
		template <typename T1, typename T2>
		void addMeanCovLocalToGlobal(Eigen::ArrayBase<T1> &local,
					     Eigen::ArrayBase<T2> &glob);
		void copyStructure(ba81NormalQuad::layer &orig);
		void prepSummary();
		void allocSummary(int numThreads);
		void addSummary(ba81NormalQuad::layer &l1);
		template <typename T1>
		void EAP(const double sampleSize, Eigen::ArrayBase<T1> &scorePad);
		template <typename T1, typename T2>
		void cacheOutcomeProb(const double *ispec, double *iparam,
				      rpf_prob_t prob_fn, int outcomes,
				      Eigen::MatrixBase<T1> &abx,
				      Eigen::MatrixBase<T2> &abscissa);
	};

	double One, ReciprocalOfOne;
	std::vector<layer> layers;
	bool DweightToThread0;

 public:
	struct ifaGroup &ig;            // should be optimized out
	static const double MIN_VARIANCE;

	bool hasBifactorStructure;
	int abilities();                // sum of per-layer abilities
	int abscissaDim() { return std::max(abilities(), 1); };

	ba81NormalQuad(struct ifaGroup *ig);
	ba81NormalQuad(ba81NormalQuad &quad);  // structure only
	void setOne(double one) { One = one; ReciprocalOfOne = 1/one; }
	void setup0();
	template <typename T1, typename T2>
	inline void refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov);
	inline double getReciprocalOfOne() const { return ReciprocalOfOne; };

	template <typename T1>
	void EAP(const double sampleSize, Eigen::ArrayBase<T1> &scorePad);

	void allocBuffers(int numThreads);
	void releaseBuffers();
	inline double computePatternLik(int thrId, int row);
	inline void prepLatentDist(int thrId);  // a.k.a. cai2010part2
	template <typename T1, typename T2, typename T3>
	void setStructure(double Qwidth, int Qpoints,
			  Eigen::ArrayBase<T1> &param,
			  Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
	inline void addToExpected(int thrId, int px);
	bool isAllocated() { return Qpoint.size(); };
	template <typename T>
	void mstepIter(int ix, T &op);
	inline void weightBy(int thrId, double weight);
	inline void weightByAndSummarize(int thrId, double weight);
	template <typename T1, typename T2>
	int cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov);
	template <typename T, typename T1>
	void computeRowDeriv(int thrId, T &op, bool freeLatents, Eigen::ArrayBase<T1> &latentGrad);
	template <typename T>
	void computeRowDeriv(int thrId, T &op);
	void prepSummary();
	void allocSummary(int numThreads);
	void addSummary(ba81NormalQuad &quad);
	void cacheOutcomeProb(double *param, bool wantLog);
	void releaseDerivCoefCache();
	void allocEstep(int numThreads);
	void releaseEstep();
	void prepExpectedTable();
	int getEstepTableSize(int lx) { return layers[lx].expected.rows(); };
	template <typename T1>
	void exportEstepTable(int lx, Eigen::ArrayBase<T1> &out);
	double mstepFit();
};

template <typename T1, typename T2, typename T3, typename T4>
void ba81NormalQuad::layer::calcDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
		   Eigen::MatrixBase<T3> &icov, Eigen::MatrixBase<T4> &where, int qx)
{
	const int pDims = primaryDims;
	const char R='R';
	const char L='L';
	const char U='U';
	const double alpha = 1;
	const double beta = 0;
	const int one = 1;

	std::vector<double> whereDiff(pDims);
	std::vector<double> whereGram(triangleLoc1(pDims));
	for (int d1=0; d1 < pDims; ++d1) {
		whereDiff[d1] = where[d1] - meanVec[d1];
	}
	gramProduct(whereDiff.data(), whereDiff.size(), whereGram.data());

	F77_CALL(dsymv)(&U, &pDims, &alpha, icov.derived().data(), &pDims, whereDiff.data(), &one,
			&beta, &derivCoef.coeffRef(0,qx), &one);

	std::vector<double> covGrad1(pDims * pDims);
	std::vector<double> covGrad2(pDims * pDims);

	int cx=0;
	for (int d1=0; d1 < pDims; ++d1) {
		for (int d2=0; d2 <= d1; ++d2) {
			covGrad1[d2 * pDims + d1] = cov(d2,d1) - whereGram[cx];
			++cx;
		}
	}

	F77_CALL(dsymm)(&R, &L, &pDims, &pDims, &alpha, covGrad1.data(), &pDims, icov.derived().data(),
			&pDims, &beta, covGrad2.data(), &pDims);
	F77_CALL(dsymm)(&R, &L, &pDims, &pDims, &alpha, icov.derived().data(), &pDims, covGrad2.data(),
			&pDims, &beta, covGrad1.data(), &pDims);

	for (int d1=0; d1 < pDims; ++d1) {
		covGrad1[d1 * pDims + d1] *= 0.5;
	}

	cx = pDims;
	for (int d1=0; d1 < pDims; ++d1) {
		int cell = d1 * pDims;
		for (int d2=0; d2 <= d1; ++d2) {
			derivCoef(cx,qx) = -covGrad1[cell + d2];
			++cx;
		}
	}
}

template <typename T1, typename T2, typename T4>
void ba81NormalQuad::layer::calcDerivCoef1(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
					   Eigen::MatrixBase<T4> &where, int qx, int curGroup)
{
	int primaryDeriv = primaryDims + triangleLoc1(primaryDims);
	int base = primaryDeriv + curGroup*2;
	const int specific = maxDims - 1 + curGroup;
	double svar = cov(specific, specific);
	double whereDiff = where[maxDims-1] - meanVec[specific];
	derivCoef(base, qx) = whereDiff / svar;
	derivCoef(base+1, qx) = -(svar - whereDiff * whereDiff) / (2 * svar * svar);
}

template <typename T1>
int ba81quad_InvertSymmetricPosDef(Eigen::MatrixBase<T1> &mat, const char uplo)
{
	if (mat.rows() != mat.cols()) Rf_error("Not square");
	int size = mat.rows();
	int info;
	F77_CALL(dpotrf)(&uplo, &size, mat.derived().data(), &size, &info);
	if (info < 0) Rf_error("Arg %d is invalid", -info);
	if (info > 0) return info;
    
	F77_CALL(dpotri)(&uplo, &size, mat.derived().data(), &size, &info);
	if (info < 0) Rf_error("Arg %d is invalid", -info);
	return info;
}

template <typename T1, typename T2>
int ba81NormalQuad::layer::cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov)
{
	Eigen::MatrixXd priCov = cov.topLeftCorner(primaryDims, primaryDims);
	Eigen::MatrixXd icov = priCov;
	int info = ba81quad_InvertSymmetricPosDef(icov, 'U');
	if (info) return info;
	icov.triangularView<Eigen::Lower>() = icov.transpose().triangularView<Eigen::Lower>();
	
	Eigen::VectorXi abx(numAbil());
	Eigen::VectorXd abscissa(numAbil());
	if (numSpecific == 0) {
		derivCoef.resize(numAbil() + triangleLoc1(numAbil()), totalQuadPoints);
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abx, abscissa);
			calcDerivCoef(meanVec, priCov, icov, abscissa, qx);
		}
	} else {
		derivCoef.resize(primaryDims + triangleLoc1(primaryDims) + 2 * numSpecific, totalQuadPoints);
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abx, abscissa);
			calcDerivCoef(meanVec, priCov, icov, abscissa, qx);
			for (int curGroup=0; curGroup < numSpecific; ++curGroup) {
				calcDerivCoef1(meanVec, cov, abscissa, qx, curGroup);
			}
		}
	}
	return 0;
}

template <typename T1, typename T2>
int ba81NormalQuad::cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov)
{
	int offset=0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		int la = layers[lx].numAbil();
		Eigen::VectorXd meanVec1 = meanVec.segment(offset, la);
		Eigen::MatrixXd cov1 = cov.block(offset, offset, la, la);
		int info = layers[lx].cacheDerivCoef(meanVec1, cov1);
		if (info) return info;
		offset += la;
	}
	return 0;
}

template <typename T2>
void ba81NormalQuad::layer::mapLatentDeriv(double piece, int qx, Eigen::ArrayBase<T2> &derivOut)
{
	int cx = 0;
	for (int d1=0; d1 < primaryDims; ++d1) {
		double amt1 = piece * derivCoef(d1, qx);
		derivOut[d1] += amt1;
		for (int d2=0; d2 <= d1; ++d2) {
			int to = numAbil() + cx;
			double amt2 = piece * derivCoef(primaryDims + cx, qx);
			derivOut[to] += amt2;
			++cx;
		}
	}
}

template <typename T2>
void ba81NormalQuad::layer::mapLatentDerivS(int sgroup, double piece, int qx, int curGroup,
					    Eigen::ArrayBase<T2> &derivOut)
{
	const int priDerivCoef = primaryDims + triangleLoc1(primaryDims);
	int base = priDerivCoef + 2 * curGroup;
	
	int sdim = primaryDims + sgroup;
	double amt3 = piece * derivCoef(base, qx);
	derivOut[sdim] += amt3;

	double amt4 = piece * derivCoef(base+1, qx);
	int to = numAbil() + triangleLoc0(sdim);
	derivOut[to] += amt4;
}

template <typename T, typename T1, typename T2, typename T3>
void ba81NormalQuad::layer::computeRowDeriv(int thrId, Eigen::MatrixBase<T2> &abx,
					    Eigen::MatrixBase<T3> &abscissa, T &op,
					    bool freeLatents, Eigen::ArrayBase<T1> &latentGradOut)
{
	abscissa.setZero();
	const int numLatents = numAbil() + triangleLoc1(numAbil());
	Eigen::ArrayXd latentGrad(numLatents);
	latentGrad.setZero();
	const int specificPoints = quad->gridSize;

	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToGlobalAbscissa(qx, abx, abscissa);
			op.beginQuadPoint(thrId);

			double tmp = Qweight(qx, thrId);
			for (int ix=0; ix < op.getNumItems(); ++ix) {
				op(thrId, abscissa, tmp, ix);
			}
			if (freeLatents) mapLatentDeriv(tmp, qx, latentGrad);

			op.endQuadPoint(thrId);
		}
	} else {
		for (int qloc=0, eisloc=0, qx=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
			for (int sx=0; sx < specificPoints; sx++) {
				pointToGlobalAbscissa(qx, abx, abscissa);
				op.beginQuadPoint(thrId);
				if (freeLatents) mapLatentDeriv(Qweight(qloc, thrId), qx, latentGrad);

				for (int ix=0; ix < op.getNumItems(); ++ix) {
					op(thrId, abscissa, Qweight(qloc + Sgroup[ix], thrId), ix);
				}

				for (int curGroup=0; curGroup < numSpecific; curGroup++) {
					double tmp = Qweight(qloc, thrId);
					if (freeLatents) mapLatentDerivS(curGroup, tmp, qx, curGroup, latentGrad);
					++qloc;
				}
				++qx;
				op.endQuadPoint(thrId);
			}
		}
	}
	if (freeLatents) addMeanCovLocalToGlobal(latentGrad, latentGradOut);
}

template <typename T, typename T1>
void ba81NormalQuad::computeRowDeriv(int thrId, T &op, bool freeLatents, Eigen::ArrayBase<T1> &latentGrad)
{
	Eigen::VectorXi abx;
	Eigen::VectorXd abscissa;
	abx.resize(abscissaDim());
	abscissa.resize(abscissaDim());
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].computeRowDeriv(thrId, abx, abscissa, op, freeLatents, latentGrad);
	}
}

template <typename T>
void ba81NormalQuad::computeRowDeriv(int thrId, T &op)
{
	Eigen::ArrayXd placeholder;
	computeRowDeriv(thrId, op, false, placeholder);
}

namespace ba81quad {

	template<typename T, typename U> struct is_same {
		static const bool value = false;
	};

	template<typename T> struct is_same<T, T> {
		static const bool value = true;
	};

};
using namespace ba81quad;

// dims should correspond to the larger possible qx
template <typename T1>
void ba81NormalQuad::decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out, int dims)
{
	is_same<typename T1::Scalar, int> isInt;
	if (!isInt.value) Rf_error("should be int");

	for (int dx=dims-1; dx >= 0; --dx) {
		out[dx] = qx % base;
		qx = qx / base;
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::pointToGlobalAbscissa(int qx, Eigen::MatrixBase<T1> &abx,
							  Eigen::MatrixBase<T2> &abscissa)
{
	std::vector<double> &Qpoint = quad->Qpoint;
	decodeLocation(qx, quad->gridSize, abx, maxDims);
	for (int dx=0; dx < numAbil(); dx++) {
		abscissa[abilitiesMap[dx]] = Qpoint[abx[std::min(dx, primaryDims)]];
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::pointToLocalAbscissa(int qx, Eigen::MatrixBase<T1> &abx,
						 Eigen::MatrixBase<T2> &abscissa)
{
	std::vector<double> &Qpoint = quad->Qpoint;
	decodeLocation(qx, quad->gridSize, abx, maxDims);
	for (int dx=0; dx < numAbil(); dx++) {
		abscissa[dx] = Qpoint[abx[std::min(dx, primaryDims)]];
	}
}

template <typename T1>
void ba81NormalQuad::layer::mapDenseSpace(double piece, const double *where,
					  const double *whereGram, Eigen::ArrayBase<T1> &latentDist)
{
	int gx = 0;
	int cx = numAbil();
	for (int d1=0; d1 < primaryDims; d1++) {
		double piece_w1 = piece * where[d1];
		latentDist[d1] += piece_w1;
		for (int d2=0; d2 <= d1; d2++) {
			double piece_cov = piece * whereGram[gx];
			latentDist[cx] += piece_cov;
			++cx; ++gx;
		}
	}
}

template <typename T1>
void ba81NormalQuad::layer::mapSpecificSpace(int sgroup, double piece, const double *where,
					     const double *whereGram, Eigen::ArrayBase<T1> &latentDist)
{
	int sdim = primaryDims + sgroup;
	double piece_w1 = piece * where[primaryDims];
	latentDist[sdim] += piece_w1;

	double piece_var = piece * whereGram[triangleLoc0(primaryDims)];
	int to = numAbil() + triangleLoc0(sdim);
	latentDist[to] += piece_var;
}

template <typename T1>
void ba81NormalQuad::layer::finalizeLatentDist(const double sampleSize, Eigen::ArrayBase<T1> &scorePad)
{
	const int padSize = numAbil() + triangleLoc1(numAbil());
	for (int d1=0; d1 < padSize; d1++) {
		scorePad[d1] /= sampleSize;
	}

	int cx = numAbil();
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
		scorePad[numAbil() + triangleLoc0(sdim)] -= ma1 * ma1;
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::addMeanCovLocalToGlobal(Eigen::ArrayBase<T1> &local,
						    Eigen::ArrayBase<T2> &glob)
{
	//mxPrintMat("local", local);
	int totalAbilities = quad->abilities();
	int cx = numAbil();
	for (int a1=0; a1 < numAbil(); ++a1) {   //row
		glob[abilitiesMap[a1]] += local[a1];
		//mxPrintMat("glob", glob);
		for (int a2=0; a2 <= a1; ++a2) { //col
			glob[totalAbilities + triangleLoc1(abilitiesMap[a1]) +
			     abilitiesMap[a2]] += local[cx++];
			//mxPrintMat("glob", glob);
		}
	}
	//mxPrintMat("glob final", glob);
}

template <typename T1>
void ba81NormalQuad::layer::EAP(const double sampleSize,
				Eigen::ArrayBase<T1> &scorePad)
{
	Eigen::ArrayXd layerPad(numAbil() + triangleLoc1(numAbil()));
	layerPad.setZero();
	Eigen::VectorXi abx(numAbil());
	Eigen::VectorXd abscissa(numAbil());
	Eigen::VectorXd whereGram(triangleLoc1(maxDims));
	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			pointToLocalAbscissa(qx, abx, abscissa);
			gramProduct(abscissa.data(), maxDims, whereGram.data());
			mapDenseSpace(Dweight(qx,0), abscissa.data(),
				      whereGram.data(), layerPad);
		}
	} else {
		int qloc=0;
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abx, abscissa);
			gramProduct(abscissa.data(), maxDims, whereGram.data());
			mapDenseSpace(Dweight(qloc,0), abscissa.data(), whereGram.data(), layerPad);
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapSpecificSpace(Sgroup, Dweight(qloc,0), abscissa.data(),
						 whereGram.data(), layerPad);
				++qloc;
			}
		}
	}

	finalizeLatentDist(sampleSize, layerPad);

	addMeanCovLocalToGlobal(layerPad, scorePad);
}

template <typename T1>
void ba81NormalQuad::EAP(const double sampleSize, Eigen::ArrayBase<T1> &scorePad)
{
	scorePad.setZero();
	prepSummary();
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].EAP(sampleSize, scorePad);
	}
}

class ifaGroup {
 private:
	SEXP Rdata;
	void verifyFactorNames(SEXP mat, const char *matName);
 public:
	const int numThreads;

	// item description related
	std::vector<const double*> spec;
	int itemDims;                  // == quad.abilities
	int numItems() const { return (int) spec.size(); }
	int impliedParamRows;          // based on spec set
	int paramRows;
	double *param;  // itemParam->data
	std::vector<const char*> itemNames;
	std::vector<int> itemOutcomes;
	//std::vector<int> cumItemOutcomes;
	int maxOutcomes;
	int totalOutcomes;

	// latent distribution
	double qwidth;
	int qpoints;
	ba81NormalQuad quad;
	ba81NormalQuad &getQuad() { return quad; };
	bool detectIndependence;
	bool twotier;  // rename to detectTwoTier TODO
	double *mean;
	double *cov;
	std::vector<std::string> factorNames;

	// data related
	SEXP dataRowNames;
	std::vector<const int*> dataColumns;
	std::vector<int> rowMap;       // row index into MxData
	int getNumUnique() const { return (int) rowMap.size(); }
	const char *weightColumnName;
	double *rowWeight;
 private:
	int minItemsPerScore;
 public:
	void setMinItemsPerScore(int mips);
	std::vector<bool> rowSkip;     // whether to treat the row as NA

	// workspace
	static const double SmallestPatternLik;
	int excludedPatterns;
	Eigen::ArrayXd patternLik;            // numUnique

	inline static bool validPatternLik(double pl)
	{ return std::isfinite(pl) && pl > SmallestPatternLik; }

	ifaGroup(int cores, bool _twotier);
	~ifaGroup();
	void setGridFineness(double width, int points);
	void import(SEXP Rlist);
	void importSpec(SEXP slotValue);
	void learnMaxAbilities();
	void setLatentDistribution(double *mean, double *cov);
	inline double *getItemParam(int ix) { return param + paramRows * ix; }
	inline const int *dataColumn(int col) { return dataColumns[col]; };
	void buildRowSkip();
	void ba81OutcomeProb(double *param, bool wantLog);
	void setFactorNames(std::vector<const char *> &names);
};

template <typename T1, typename T2>
void ba81NormalQuad::layer::cacheOutcomeProb(const double *ispec, double *iparam,
					     rpf_prob_t prob_fn, int ix,
					     Eigen::MatrixBase<T1> &abx,
					     Eigen::MatrixBase<T2> &abscissa)
{
	ix = glItemsMap[ix];
	if (ix == -1) return;
	abscissa.setZero();
	double *curOutcome = &outcomeProbX.coeffRef(totalQuadPoints * cumItemOutcomes[ix]);
	int outcomes = itemOutcomes[ix];
	for (int qx=0; qx < totalQuadPoints; qx++) {
		pointToGlobalAbscissa(qx, abx, abscissa);
		//mxPrintMat("wh", abscissa);
		(*prob_fn)(ispec, iparam, abscissa.derived().data(), curOutcome);
		//for (int ox=0; ox < outcomes; ++ox) mxLog("%d %d %g", ix, ox, curOutcome[ox]);
		curOutcome += outcomes;
	}
}

double ba81NormalQuad::layer::computePatternLik(int thrId, int row)
{
	double *out = &Qweight.coeffRef(0, thrId);
	double *oProb = outcomeProbX.data();

	double patternLik = 0.0;
	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			out[qx] = priQarea[qx];
		}

		for (int ix=0; ix < numItems(); ix++) {
			int outcomes = itemOutcomes[ix];
			int pick = dataColumns[ix][row];
			if (pick == NA_INTEGER) {
				oProb += outcomes * totalQuadPoints;
				continue;
			}
			pick -= 1;

			for (int qx=0; qx < totalQuadPoints; ++qx) {
				out[qx] *= oProb[pick];
				oProb += outcomes;
			}
		}

		for (int qx=0; qx < totalQuadPoints; ++qx) {
			patternLik += out[qx];
		}
	} else {
		const int specificPoints = quad->gridSize;
		double *Ei = &thrEi.coeffRef(0, thrId);
		double *Eis = &thrEis.coeffRef(0, thrId);

		for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
			for (int sx=0; sx < specificPoints * numSpecific; sx++) {
				out[qloc] = speQarea[sx];
				++qloc;
			}
		}

		for (int ix=0; ix < numItems(); ix++) {
			int outcomes = itemOutcomes[ix];
			int pick = dataColumns[ix][row];
			if (pick == NA_INTEGER) {
				oProb += outcomes * totalQuadPoints;
				continue;
			}
			pick -= 1;
			int Sgroup1 = Sgroup[ix];
			double *out1 = out;
			for (int qx=0; qx < totalQuadPoints; qx++) {
				out1[Sgroup1] *= oProb[pick];
				//mxLog("%d %d %d %g %g", ix, qx, pick, oProb[pick], out1[Sgroup1]);
				oProb += outcomes;
				out1 += numSpecific;
			}
		}

		//mxPrintMat("qweight final", Qweight.col(thrId));

		thrEis.col(thrId).setZero();
		for (int qx=0; qx < totalPrimaryPoints; ++qx) Ei[qx] = priQarea[qx];

		int eisloc = 0;
		for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
			for (int sx=0; sx < specificPoints; sx++) {
				for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
					double piece = out[qloc];
					Eis[eisloc + sgroup] += piece;
					++qloc;
				}
			}
			double roo = quad->getReciprocalOfOne();
			for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
				Ei[qx] *= Eis[eisloc + sgroup] * roo;
			}
			eisloc += numSpecific;
		}
		for (int qx=0; qx < totalPrimaryPoints; ++qx) {
			patternLik += Ei[qx];
		}
		//mxPrintMat("eis", thrEis.col(thrId));
		//mxPrintMat("ei", thrEi.col(thrId));
	}
	return patternLik;
}

double ba81NormalQuad::computePatternLik(int thrId, int row)
{
	double patternLik = 1.0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		patternLik *= l1.computePatternLik(thrId, row);
	}
	return patternLik;
}

void ba81NormalQuad::layer::prepLatentDist(int thrId)
{
	if (0 == numSpecific) return;

	double *Ei = &thrEi.coeffRef(0, thrId);
	double *Eis = &thrEis.coeffRef(0, thrId);

	const int specificPoints = quad->gridSize;

	for (int qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Eis[qloc] = Ei[qx] / Eis[qloc];
			++qloc;
		}
	}

	for (int qloc=0, eisloc=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
		for (int sx=0; sx < specificPoints; sx++) {
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				Qweight(qloc, thrId) *= Eis[eisloc + Sgroup];
				++qloc;
			}
		}
	}
}

void ba81NormalQuad::prepLatentDist(int thrId)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		l1.prepLatentDist(thrId);
	}
}

void ba81NormalQuad::layer::addToExpected(int thrId, int px)
{
	double *out = &expected.coeffRef(0, thrId);

	for (int ix=0; ix < numItems(); ++ix) {
		int outcomes = itemOutcomes[ix];
		int pick = dataColumns[ix][px];
		if (pick == NA_INTEGER) {
			out += outcomes * totalQuadPoints;
			continue;
		}
		pick -= 1;
		double *lastQw = &Qweight.coeffRef(0, thrId) + Qweight.rows();
		if (numSpecific == 0) {
			double *Qw = &Qweight.coeffRef(0, thrId);
			while (Qw < lastQw) {
				out[pick] += *Qw;
				out += outcomes;
				++Qw;
			}
		} else {
			int igroup = Sgroup[ix];
			double *Qw = &Qweight.coeffRef(igroup, thrId);
			while (Qw < lastQw) {
				out[pick] += *Qw;
				out += outcomes;
				Qw += numSpecific;
			}
		}
	}
}

void ba81NormalQuad::addToExpected(int thrId, int px)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].addToExpected(thrId, px);
	}
}

void ba81NormalQuad::layer::weightByAndSummarize(int thrId, double weight)
{
	for (int qx=0; qx < weightTableSize; ++qx) {
		Qweight(qx, thrId) *= weight;
		Dweight(qx, thrId) += Qweight(qx, thrId);
	}
}

void ba81NormalQuad::weightByAndSummarize(int thrId, double weight)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].weightByAndSummarize(thrId, weight);
	}
}

void ba81NormalQuad::layer::weightBy(int thrId, double weight)
{
	for (int qx=0; qx < weightTableSize; ++qx) {
		Qweight(qx, thrId) *= weight;
	}
}

void ba81NormalQuad::weightBy(int thrId, double weight)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].weightBy(thrId, weight);
	}
}

template <typename T>
struct BA81OmitEstep {
	void begin(class ifaGroup *state) {};
	void addRow(class ifaGroup *state, int px, int thrId) {};
	void recordTable(class ifaGroup *state) {};
	bool hasEnd() { return false; }
};

template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename> class EstepPolicy
>
struct BA81Engine : LatentPolicy<T>, EstepPolicy<T> {
	void ba81Estep1(class ifaGroup *state, T extraData);
};

template <
  typename T,
  template <typename> class LatentPolicy,
  template <typename> class EstepPolicy
>
void BA81Engine<T, LatentPolicy, EstepPolicy>::ba81Estep1(class ifaGroup *state, T extraData)
{
	ba81NormalQuad &quad = state->getQuad();
	const int numUnique = state->getNumUnique();
	const int numThreads = state->numThreads;
	state->excludedPatterns = 0;
	state->patternLik.resize(numUnique);
	Eigen::ArrayXd &patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<T>::begin(state);

	quad.allocBuffers(numThreads);
	if (LatentPolicy<T>::wantSummary()) quad.allocSummary(numThreads);

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omp_get_thread_num();
		int mpx = state->rowMap[px];
		double patternLik1 = state->quad.computePatternLik(thrId, mpx);

		if (!ifaGroup::validPatternLik(patternLik1)) {
#pragma omp atomic
			state->excludedPatterns += 1;
			patternLik[px] = 0;
			continue;
		}

		patternLik[px] = patternLik1;

		if (!EstepPolicy<T>::hasEnd() && !LatentPolicy<T>::hasEnd()) continue;

		state->quad.prepLatentDist(thrId);

		LatentPolicy<T>::normalizeWeights(state, extraData, px, patternLik1, thrId);
		EstepPolicy<T>::addRow(state, mpx, thrId);
	}

	if (EstepPolicy<T>::hasEnd() && LatentPolicy<T>::hasEnd()) {
#pragma omp parallel sections
		{
			{ EstepPolicy<T>::recordTable(state); }
#pragma omp section
			{ LatentPolicy<T>::end(state, extraData); }
		}
	} else {
		EstepPolicy<T>::recordTable(state);
		LatentPolicy<T>::end(state, extraData);
	}

	quad.releaseBuffers();
}

template <typename T1, typename T2>
void ba81NormalQuad::refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].refresh(mean, cov);
	}
}

template <typename T1>
void ba81quad_InplaceForcePosDef(Eigen::MatrixBase<T1> &cov)
{
	const double tooSmallEV = 1e-6;

	Eigen::LDLT<T1> ldlt(cov);
	if (ldlt.info() == Eigen::Success && (ldlt.vectorD().array() > tooSmallEV).all()) return;

	// We will rarely need to do the full decomposition
	Eigen::SelfAdjointEigenSolver<T1> esol(cov);
	Eigen::VectorXd ev = esol.eigenvalues().array().max(tooSmallEV).matrix();
	cov.derived() = esol.eigenvectors() * ev * esol.eigenvectors().transpose();
}

template <typename T1, typename T2, typename T3, typename T4>
void ba81NormalQuad::layer::globalToLocalDist(Eigen::MatrixBase<T1> &gmean, Eigen::MatrixBase<T2> &gcov,
					      Eigen::MatrixBase<T3> &mean, Eigen::MatrixBase<T4> &cov)
{
	struct subsetOp {
		std::vector<bool> &abilitiesMask;
		subsetOp(std::vector<bool> &abilitiesMask) : abilitiesMask(abilitiesMask) {};
		bool operator()(int gx) { return abilitiesMask[gx]; };
	} op(abilitiesMask);

	subsetNormalDist(gmean, gcov, op, numAbil(), mean, cov);
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::refresh(Eigen::MatrixBase<T2> &gmeanVec, Eigen::MatrixBase<T1> &gcov)
{
	if (numAbil() == 0) {
		priQarea.clear();
		priQarea.push_back(quad->One);
		return;
	}

	Eigen::VectorXd meanVec;
	Eigen::MatrixXd cov;
	globalToLocalDist(gmeanVec, gcov, meanVec, cov);

	// This is required because EM acceleration can routinely push
	// the covariance matrix to be non-pd.
	if (primaryDims == 1) {
		cov(0, 0) = std::max(cov(0, 0), MIN_VARIANCE);
	} else {
		Eigen::MatrixXd priCov = cov.topLeftCorner(primaryDims, primaryDims);
		ba81quad_InplaceForcePosDef(priCov);
		cov.topLeftCorner(primaryDims, primaryDims) = priCov;
	}

	for (int sx=0; sx < numSpecific; ++sx) {
		int loc = primaryDims + sx;
		cov(loc, loc) = std::max(cov(loc, loc), MIN_VARIANCE);
	}

	std::vector<double> &Qpoint = quad->Qpoint;

	Eigen::VectorXi abscissa(primaryDims);
	Eigen::MatrixXd priCov = cov.topLeftCorner(primaryDims, primaryDims);
	std::vector<double> tmpPriQarea;
	tmpPriQarea.reserve(totalPrimaryPoints);
	{
		Eigen::VectorXd where(primaryDims);
		for (int qx=0; qx < totalPrimaryPoints; qx++) {
			decodeLocation(qx, quad->gridSize, abscissa, primaryDims);
			for (int dx=0; dx < primaryDims; dx++) where[dx] = Qpoint[abscissa[dx]];
			double den = exp(dmvnorm(primaryDims, where.data(),
						 meanVec.derived().data(), priCov.data()));
			tmpPriQarea.push_back(den);
		}
	}

	priQarea.clear();
	priQarea.reserve(totalPrimaryPoints);

	double totalArea = 0;
	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		double den = tmpPriQarea[qx];
		priQarea.push_back(den);
		//double prevTotalArea = totalArea;
		totalArea += den;
		// if (totalArea == prevTotalArea) {
		// 	mxLog("%.4g / %.4g = %.4g", den, totalArea, den / totalArea);
		// }
	}

	for (int qx=0; qx < totalPrimaryPoints; qx++) {
		// must be in correct order to avoid overflow
		priQarea[qx] *= quad->One;
		priQarea[qx] /= totalArea;
		//mxLog("%.5g,", priQarea[qx] / quad->One);
	}
	//pda(priQarea.data(), 1, totalPrimaryPoints);

	if (numSpecific) {
		speQarea.resize(quad->gridSize * numSpecific);
	}

	for (int sgroup=0; sgroup < numSpecific; sgroup++) {
		totalArea = 0;
		double mean = meanVec[primaryDims + sgroup];
		double var = cov.diagonal().coeff(primaryDims + sgroup);
		for (int qx=0; qx < quad->gridSize; qx++) {
			double den = exp(dmvnorm(1, &Qpoint[qx], &mean, &var));
			speQarea[sIndex(sgroup, qx)] = den;
			totalArea += den;
		}
		for (int qx=0; qx < quad->gridSize; qx++) {
			speQarea[sIndex(sgroup, qx)] /= totalArea;
		}
	}
	//pda(speQarea.data(), numSpecific, quad->gridSize);

	for (int sx=0; sx < int(speQarea.size()); ++sx) {
		speQarea[sx] *= quad->One;
	}
	//pda(speQarea.data(), numSpecific, quadGridSize);
}

template <typename T1, typename T2, typename T3>
void ba81NormalQuad::layer::mstepIter(int ix, Eigen::MatrixBase<T1> &abx,
				      Eigen::MatrixBase<T2> &abscissa, T3 &op)
{
	ix = glItemsMap[ix];
	if (ix == -1) return;
	abscissa.setZero();
	int offset = cumItemOutcomes[ix] * totalQuadPoints;
	double *curOutcome = &outcomeProbX.coeffRef(offset);
	double *iexp = &expected.coeffRef(offset, 0);
	int outcomes = itemOutcomes[ix];

	for (int qx=0; qx < totalQuadPoints; ++qx) {
		pointToGlobalAbscissa(qx, abx, abscissa);
		op(abscissa.derived().data(), curOutcome, iexp);
		curOutcome += outcomes;
		iexp += outcomes;
	}
}

template <typename T>
void ba81NormalQuad::mstepIter(int ix, T &op)
{
	Eigen::VectorXi abx;
	Eigen::VectorXd abscissa;
	abx.resize(abscissaDim());
	abscissa.resize(abscissaDim());
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].mstepIter(ix, abx, abscissa, op);
	}
}

template <typename T1>
void ba81NormalQuad::exportEstepTable(int lx, Eigen::ArrayBase<T1> &out)
{
	out.derived() = layers[lx].expected.col(0);
}

#endif
