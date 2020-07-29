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

#ifndef _BA81QUAD_H_
#define _BA81QUAD_H_

#include "omxDefines.h" // need OpenMx's Eigen customizations

#include <Rcpp.h>
using namespace Rcpp;

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include "libifa-rpf.h"
#include "dmvnorm.h"

extern const struct rpf *Glibrpf_model;
extern int Glibrpf_numModels;

namespace ba81quad {

#ifndef _OPENMP
	static inline int omp_get_thread_num() { return 0; }
#endif

	static inline int triangleLoc1(int diag)
	{
		return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
	}

	static inline int triangleLoc0(int diag)
	{
		return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
	}

	static inline bool strEQ(const char *s1, const char *s2) { return strcmp(s1,s2)==0; }

	template<typename _MatrixType, int _UpLo = Eigen::Lower>
	class SimpCholesky : public Eigen::LDLT<_MatrixType, _UpLo> {
	private:
		Eigen::MatrixXd inverse;

	public:
		typedef Eigen::LDLT<_MatrixType, _UpLo> Base;

		SimpCholesky() : Base() {};
		template<typename InputType>
		explicit SimpCholesky(const Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};
		template<typename InputType>
		explicit SimpCholesky(Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};

		double log_determinant() const {
			typename Base::Scalar detL = Base::vectorD().array().log().sum();
			return detL * 0.5;
		}

		void refreshInverse()
		{
			inverse.setIdentity(Base::m_matrix.rows(), Base::m_matrix.rows());
			inverse = Base::solve(inverse);
		};

		const Eigen::MatrixXd &getInverse() const { return inverse; };
	};

	template <typename T1>
	int InvertSymmetricPosDef(Eigen::MatrixBase<T1> &mat)
	{
		SimpCholesky< Eigen::Ref<T1>, Eigen::Upper > sc(mat);
		if (sc.info() != Eigen::Success || !sc.isPositive()) {
			return -1;
		} else {
			sc.refreshInverse();
			mat.derived() = sc.getInverse();
			return 0;
		}
	}

	template <typename T1, typename T2, typename T3, typename T4, typename T5>
	void subsetNormalDist(const Eigen::MatrixBase<T1> &gmean, const Eigen::MatrixBase<T2> &gcov,
												T5 includeTest, int resultSize,
												Eigen::MatrixBase<T3> &mean, Eigen::MatrixBase<T4> &cov)
	{
		mean.derived().resize(resultSize);
		cov.derived().resize(resultSize, resultSize);

		for (int gcx=0, cx=0; gcx < gcov.cols(); gcx++) {
			if (!includeTest(gcx)) continue;
			mean[cx] = gmean[gcx];
			for (int grx=0, rx=0; grx < gcov.rows(); grx++) {
				if (!includeTest(grx)) continue;
				cov(rx,cx) = gcov(grx, gcx);
				rx += 1;
			}
			cx += 1;
		}
	}

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

} // namespace

class ba81NormalQuad {
private:
	template <typename T1>
	static inline void decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out, int dims);

	double width;
	std::vector<double> Qpoint;           // gridSize

public:
	int numThreads;
	void setNumThreads(int nt) { numThreads = nt; }
	int gridSize;

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

		// outcome info
		std::vector<int> itemOutcomes;
		std::vector<int> cumItemOutcomes;
		int totalOutcomes;
		std::vector<const int*> dataColumns;
		std::vector<const double*> spec;
		int paramRows;

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
			//if (sx < 0 || sx >= state->numSpecific) mxThrow("Out of domain");
			//if (qx < 0 || qx >= state->gridSize) mxThrow("Out of domain");
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
											Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov, bool testTwoTier);
		void setupOutcomes(class ifaGroup &ig);
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
		template <typename T1, typename T2>
		void EAP(Eigen::ArrayBase<T1> &wvec, const double sampleSize, Eigen::ArrayBase<T2> &scorePad);
		template <typename T1, typename T2>
		void cacheOutcomeProb(const double *ispec, double *iparam,
													rpf_prob_t prob_fn, int outcomes,
													Eigen::MatrixBase<T1> &abx,
													Eigen::MatrixBase<T2> &abscissa);
	};

private:
	double One, ReciprocalOfOne;
	std::vector<layer> layers;

public:
	static const double MIN_VARIANCE;

	bool hasBifactorStructure;
	int abilities();                // sum of per-layer abilities
	int abscissaDim() { return std::max(abilities(), 1); };

	ba81NormalQuad();
	ba81NormalQuad(ba81NormalQuad &quad);  // structure only
	layer &getLayer() { return layers[0]; }
	void setOne(double one) { One = one; ReciprocalOfOne = 1/one; }
	void setup0();
	template <typename T1, typename T2>
	inline void refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov);
	inline double getReciprocalOfOne() const { return ReciprocalOfOne; };

	template <typename T1, typename T2>
	void EAP(Eigen::ArrayBase<T1> &wvec, const double sampleSize, Eigen::ArrayBase<T2> &scorePad);
	template <typename T2>
	void EAP(const double sampleSize, Eigen::ArrayBase<T2> &scorePad);

	void allocBuffers();
	void releaseBuffers();
	inline double computePatternLik(int thrId, int row);
	inline void prepLatentDist(int thrId);  // a.k.a. cai2010part2
	template <typename T1, typename T2, typename T3>
	void setStructure(double Qwidth, int Qpoints,
										Eigen::ArrayBase<T1> &param,
										Eigen::MatrixBase<T2> &mean,
										Eigen::MatrixBase<T3> &cov, bool testTwoTier);
	void setupOutcomes(class ifaGroup &ig);
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
	void allocSummary();
	void addSummary(ba81NormalQuad &quad);
	void cacheOutcomeProb(double *param, bool wantLog);
	void releaseDerivCoefCache();
	void allocEstep();
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
	using ba81quad::triangleLoc1;
	using ba81quad::gramProduct;
	const int pDims = primaryDims;

	Eigen::VectorXd whereDiff(pDims);
	std::vector<double> whereGram(triangleLoc1(pDims));
	for (int d1=0; d1 < pDims; ++d1) {
		whereDiff[d1] = where[d1] - meanVec[d1];
	}
	gramProduct(whereDiff.data(), whereDiff.size(), whereGram.data());

	Eigen::Map< Eigen::VectorXd > dcol(&derivCoef.coeffRef(0,qx), icov.rows());
	dcol = icov.template selfadjointView<Eigen::Lower>() * whereDiff;

	Eigen::MatrixXd covGrad1(pDims, pDims);

	int cx=0;
	for (int d1=0; d1 < pDims; ++d1) {
		for (int d2=0; d2 <= d1; ++d2) {
			covGrad1(d1, d2) = cov(d2,d1) - whereGram[cx];
			++cx;
		}
	}

	covGrad1 = icov * covGrad1.template selfadjointView<Eigen::Lower>() * icov.template selfadjointView<Eigen::Lower>();
	covGrad1.diagonal() *= 0.5;

	cx = pDims;
	for (int d1=0; d1 < pDims; ++d1) {
		for (int d2=0; d2 <= d1; ++d2) {
			derivCoef(cx,qx) = -covGrad1(d2, d1);
			++cx;
		}
	}
}

template <typename T1, typename T2, typename T4>
void ba81NormalQuad::layer::calcDerivCoef1(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
																					 Eigen::MatrixBase<T4> &where, int qx, int curGroup)
{
	using ba81quad::triangleLoc1;
	int primaryDeriv = primaryDims + triangleLoc1(primaryDims);
	int base = primaryDeriv + curGroup*2;
	const int specific = maxDims - 1 + curGroup;
	double svar = cov(specific, specific);
	double whereDiff = where[maxDims-1] - meanVec[specific];
	derivCoef(base, qx) = whereDiff / svar;
	derivCoef(base+1, qx) = -(svar - whereDiff * whereDiff) / (2 * svar * svar);
}

template <typename T1, typename T2>
int ba81NormalQuad::layer::cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov)
{
	using ba81quad::triangleLoc1;
	Eigen::MatrixXd priCov = cov.topLeftCorner(primaryDims, primaryDims);
	Eigen::MatrixXd icov = priCov;
	int info = ba81quad::InvertSymmetricPosDef(icov);
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
	using ba81quad::triangleLoc1;
	using ba81quad::triangleLoc0;
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
	using ba81quad::triangleLoc1;
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

template<typename T, typename U> struct is_same {
	static const bool value = false;
};

template<typename T> struct is_same<T, T> {
	static const bool value = true;
};

// dims should correspond to the larger possible qx
template <typename T1>
void ba81NormalQuad::decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out, int dims)
{
	is_same<typename T1::Scalar, int> isInt;
	if (!isInt.value) mxThrow("should be int");

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
	using ba81quad::triangleLoc0;
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
	using ba81quad::triangleLoc0;
	scorePad *= (1.0/sampleSize);

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
	using ba81quad::triangleLoc1;
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

template <typename T1, typename T2>
void ba81NormalQuad::layer::EAP(Eigen::ArrayBase<T1> &wvec,
																const double sampleSize,
																Eigen::ArrayBase<T2> &scorePad)
{
	using ba81quad::triangleLoc1;
	using ba81quad::gramProduct;
	Eigen::ArrayXd layerPad(numAbil() + triangleLoc1(numAbil()));
	layerPad.setZero();
	Eigen::VectorXi abx(numAbil());
	Eigen::VectorXd abscissa(numAbil());
	Eigen::VectorXd whereGram(triangleLoc1(maxDims));
	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			pointToLocalAbscissa(qx, abx, abscissa);
			gramProduct(abscissa.data(), maxDims, whereGram.data());
			mapDenseSpace(wvec(qx), abscissa.data(),
										whereGram.data(), layerPad);
		}
	} else {
		int qloc=0;
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abx, abscissa);
			gramProduct(abscissa.data(), maxDims, whereGram.data());
			mapDenseSpace(wvec(qloc), abscissa.data(), whereGram.data(), layerPad);
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapSpecificSpace(Sgroup, wvec(qloc), abscissa.data(),
												 whereGram.data(), layerPad);
				++qloc;
			}
		}
	}

	finalizeLatentDist(sampleSize, layerPad);

	addMeanCovLocalToGlobal(layerPad, scorePad);
}

template <typename T1, typename T2>
void ba81NormalQuad::EAP(Eigen::ArrayBase<T1> &wvec,
												 const double sampleSize, Eigen::ArrayBase<T2> &scorePad)
{
	scorePad.setZero();
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].EAP(wvec, sampleSize, scorePad);
	}
}

template <typename T2>
void ba81NormalQuad::EAP(const double sampleSize, Eigen::ArrayBase<T2> &scorePad)
{
	auto &layer = getLayer();
	Eigen::Map< Eigen::ArrayXd > wvec(layer.Dweight.data(), layer.Dweight.rows());
	scorePad.setZero();
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].EAP(wvec, sampleSize, scorePad);
	}
}

class ifaGroup {
private:
	DataFrame Rdata;
	void verifyFactorNames(const List &dimnames, const char *matName);
public:

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
	bool twotier;  // rename to detectTwoTier TODO
	double *mean;
	double *cov;
	std::vector<std::string> factorNames;

	// data related
	CharacterVector dataRowNames;
	std::vector<const int*> dataColumns;
	std::vector<int> rowMap;       // row index into MxData
	int getNumUnique() const { return (int) rowMap.size(); }
private:
	const char *weightColumnName;
	double *rowWeight;
	const char *freqColumnName;
	int *rowFreq;
	int minItemsPerScore;
	double weightSum;
public:
	void setMinItemsPerScore(int mips);
	std::vector<bool> rowSkip;     // whether to treat the row as NA

	// workspace
	static const double SmallestPatternLik;
	int excludedPatterns;
	double getWeightSum() { return weightSum; }
	Eigen::ArrayXd rowMult;               // numUnique
	Eigen::ArrayXd patternLik;            // numUnique

	inline static bool validPatternLik(double pl)
	{ return std::isfinite(pl) && pl > SmallestPatternLik; }

	ifaGroup(bool _twotier);
	~ifaGroup();
	void setRowWeight(double *_in) { rowWeight = _in; }
	double getRowWeight(int rx) { return rowWeight ? rowWeight[rx] : 1.0; }
	void setRowFreq(int *_in) { rowFreq = _in; }
	void setGridFineness(double width, int points);
	void import(const List &Rlist);
	void importSpec(const List &slotValue);
	void learnMaxAbilities();
	void setLatentDistribution(double *mean, double *cov);
	inline double *getItemParam(int ix) { return param + paramRows * ix; }
	inline const int *dataColumn(int col) { return dataColumns[col]; };
	void buildRowSkip();
	void buildRowMult();
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
	const int numThreads = quad.numThreads;
	state->excludedPatterns = 0;
	state->patternLik.resize(numUnique);
	Eigen::ArrayXd &patternLik = state->patternLik;
	std::vector<bool> &rowSkip = state->rowSkip;

	EstepPolicy<T>::begin(state);

	quad.allocBuffers();
	if (LatentPolicy<T>::wantSummary()) quad.allocSummary();

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

	ba81quad::subsetNormalDist(gmean, gcov, op, numAbil(), mean, cov);
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

template <typename T1, typename T2, typename T3>
void ba81NormalQuad::setStructure(double Qwidth, int Qpoints,
																	Eigen::ArrayBase<T1> &param,
																	Eigen::MatrixBase<T2> &mean,
																	Eigen::MatrixBase<T3> &cov, bool testTwoTier)
{
	hasBifactorStructure = false;
	width = Qwidth;
	gridSize = Qpoints;

	if (int(Qpoint.size()) != gridSize) {
		Qpoint.clear();
		Qpoint.reserve(gridSize);
		double qgs = gridSize-1;
		for (int px=0; px < gridSize; ++px) {
			Qpoint.push_back(px * 2 * width / qgs - width);
		}
	}

	if (!mean.rows()) {
		gridSize = 1;
		layers.resize(1, layer(this));
		layers[0].itemsMask.assign(param.cols(), true);
		layers[0].setStructure(param, mean, cov, false);
		return;
	}

	layers.clear();
	layers.resize(1, layer(this));

	// This overengineering was done in error. We always have 1 layer.
	layers[0].itemsMask.assign(param.cols(), true);
	layers[0].abilitiesMask.assign(mean.rows(), true);
	layers[0].setStructure(param, mean, cov, testTwoTier);
}

template <typename T1, typename T2, typename T3>
void ba81NormalQuad::layer::setStructure(Eigen::ArrayBase<T1> &param,
																				 Eigen::MatrixBase<T2> &gmean, Eigen::MatrixBase<T3> &gcov,
																				 bool testTwoTier)
{
	abilitiesMap.clear();
	for (int ax=0; ax < gmean.rows(); ++ax) {
		if (!abilitiesMask[ax]) continue;
		abilitiesMap.push_back(ax);
	}

	itemsMap.clear();
	glItemsMap.resize(param.cols(), -1);
	for (int ix=0, lx=0; ix < param.cols(); ++ix) {
		if (!itemsMask[ix]) continue;
		itemsMap.push_back(ix);
		glItemsMap[ix] = lx++;
	}

	Eigen::VectorXd mean;
	Eigen::MatrixXd cov;
	globalToLocalDist(gmean, gcov, mean, cov);

	if (mean.size() == 0) {
		numSpecific = 0;
		primaryDims = 0;
		maxDims = 1;
		totalQuadPoints = 1;
		totalPrimaryPoints = 1;
		weightTableSize = 1;
		return;
	}

	numSpecific = 0;

	if (testTwoTier) detectTwoTier(param, mean, cov);
	if (numSpecific) quad->hasBifactorStructure = true;

	primaryDims = cov.cols() - numSpecific;
	maxDims = primaryDims + (numSpecific? 1 : 0);

	totalQuadPoints = 1;
	for (int dx=0; dx < maxDims; dx++) {
		totalQuadPoints *= quad->gridSize;
	}

	totalPrimaryPoints = totalQuadPoints;
	weightTableSize = totalQuadPoints;

	if (numSpecific) {
		totalPrimaryPoints /= quad->gridSize;
		weightTableSize *= numSpecific;
	}
}

template <typename T1, typename T2, typename T3>
void ba81NormalQuad::layer::detectTwoTier(Eigen::ArrayBase<T1> &param,
					  Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov)
{
	if (mean.rows() < 3) return;

	std::vector<int> orthogonal;

	Eigen::Matrix<Eigen::DenseIndex, Eigen::Dynamic, 1>
		numCov((cov.array() != 0.0).matrix().colwise().count());
	std::vector<int> candidate;
	for (int fx=0; fx < numCov.rows(); ++fx) {
		if (numCov(fx) == 1) candidate.push_back(fx);
	}
	if (candidate.size() > 1) {
		std::vector<bool> mask(numItems());
		for (int cx=candidate.size() - 1; cx >= 0; --cx) {
			std::vector<bool> loading(numItems());
			for (int ix=0; ix < numItems(); ++ix) {
				loading[ix] = param(candidate[cx], itemsMap[ix]) != 0;
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

	if (orthogonal.size() == 1) orthogonal.clear();
	if (orthogonal.size() && orthogonal[0] != mean.rows() - int(orthogonal.size())) {
		mxThrow("Independent specific factors must be given after general dense factors");
	}

	numSpecific = orthogonal.size();

	if (numSpecific) {
		Sgroup.assign(numItems(), 0);
		for (int ix=0; ix < numItems(); ix++) {
			for (int dx=orthogonal[0]; dx < mean.rows(); ++dx) {
				if (param(dx, itemsMap[ix]) != 0) {
					Sgroup[ix] = dx - orthogonal[0];
					continue;
				}
			}
		}
		//Eigen::Map< Eigen::ArrayXi > foo(Sgroup.data(), param.cols());
		//mxPrintMat("sgroup", foo);
	}
}

#endif
