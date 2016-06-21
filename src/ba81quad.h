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
	static inline void decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out);

	double width;
	int gridSize;
	std::vector<double> Qpoint;           // gridSize

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

		int abilitiesOffset;                  // index into the global abilities vector
		int abilities;                        // number of logical abilities in this layer
		int maxDims;                          // integration dimension after dimension reduction
		int totalQuadPoints;                  // gridSize ^ maxDims
		int weightTableSize;                  // dense: totalQuadPoints; 2tier: totalQuadPoints * numSpecific
		std::vector<double> priQarea;         // totalPrimaryPoints
		std::vector<double> wherePrep;        // totalQuadPoints * maxDims; better to recompute instead of cache? TODO
		Eigen::MatrixXd whereGram;            // triangleLoc1(maxDims) x totalQuadPoints
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
		// need to dealloc derivCoef TODO

		layer(class ba81NormalQuad *quad)
		: quad(quad), maxDims(-1), numSpecific(-1), primaryDims(-1) {};
		inline int sIndex(int sx, int qx) { // remove? TODO
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

		void allocBuffers(int numThreads, bool wantSummary);
		void releaseBuffers();
		inline double computePatternLik(int thrId, double *oProb, int row);
		inline void prepLatentDist(int thrId);
		template <typename T1, typename T2, typename T3>
		void detectTwoTier(Eigen::ArrayBase<T1> &param,
				   Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
		template <typename T1, typename T2, typename T3>
		void setStructure(Eigen::ArrayBase<T1> &param,
				  Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
		template <typename T1, typename T2>
		inline void refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov);
		inline void addTo(int thrId, int ix, double *out);
		template <typename T1, typename T2>
		inline void foreach(Eigen::MatrixBase<T1> &abscissa, T2 &op);
		inline void weightBy(int thrId, double weight);
		inline void weightByAndSummarize(int thrId, double weight);
		template <typename T1, typename T2>
		void cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov);
		template <typename T1>
		void pointToGlobalAbscissa(int qx, Eigen::MatrixBase<T1> &abscissa);
		template <typename T1>
		void pointToLocalAbscissa(int qx, Eigen::MatrixBase<T1> &abscissa);
		template <typename T, typename T1, typename T2>
		void computeRowDeriv(int thrId, Eigen::MatrixBase<T2> &abscissa, T &op,
				     bool freeLatents, Eigen::ArrayBase<T1> &latentGradOut);
		template <typename T1, typename T2>
		void addMeanCovLocalToGlobal(Eigen::ArrayBase<T1> &local,
					     Eigen::ArrayBase<T2> &glob);
	};

	double One, ReciprocalOfOne;
	std::vector<layer> layers;

 public:
	struct ifaGroup &ig;            // should be optimized out

	// Maybe an abstraction violation to expose weightTableSize, totalQuadPoints
	int totalQuadPoints;            // sum of per-layer totalQuadPoints
	int weightTableSize;            // sum of per-layer weightTableSize

	bool hasBifactorStructure;
	int abilities;                  // sum of per-layer abilities

	ba81NormalQuad(struct ifaGroup *ig);
	void setOne(double one) { One = one; ReciprocalOfOne = 1/one; }
	void setup0();
	template <typename T1, typename T2>
	inline void refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov);
	inline double getReciprocalOfOne() const { return ReciprocalOfOne; };

	template <typename T1>
	void EAP(double *thrDweight, const double sampleSize,  // use Dweight from quad? TODO
		 Eigen::ArrayBase<T1> &scorePad);

	inline void cacheOutcomeProb(const double *ispec, double *iparam,
				     rpf_prob_t prob_fn, double *qProb);
	void allocBuffers(int numThreads, bool wantSummary);
	void releaseBuffers();
	inline double computePatternLik(int thrId, int row);
	inline void prepLatentDist(int thrId);  // a.k.a. cai2010part2
	template <typename T1, typename T2, typename T3>
	void setStructure(double Qwidth, int Qpoints,
			  Eigen::ArrayBase<T1> &param,
			  Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T3> &cov);
	inline void addTo(int thrId, int ix, double *out);
	bool isAllocated() { return Qpoint.size(); };
	template <typename T>
	void foreach(T &op);
	inline void weightBy(int thrId, double weight);
	inline void weightByAndSummarize(int thrId, double weight);
	template <typename T1, typename T2>
	void cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov);
	template <typename T, typename T1>
	void computeRowDeriv(int thrId, T &op, bool freeLatents, Eigen::ArrayBase<T1> &latentGrad);
	template <typename T>
	void computeRowDeriv(int thrId, T &op);
};

template <typename T1, typename T2, typename T3, typename T4>
void ba81NormalQuad::layer::calcDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov,
		   Eigen::MatrixBase<T3> &icov, Eigen::MatrixBase<T4> &where, int qx)
{
	const int pDims = numSpecific? maxDims - 1 : maxDims;
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
void ba81NormalQuad::layer::cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov)
{
	Eigen::MatrixXd icov = cov;
	int info = ba81quad_InvertSymmetricPosDef(icov, 'U');
	if (info) {
		// report error TODO
		//omxRaiseErrorf("%s: latent covariance matrix is not positive definite", oo->name());
		return;
	}
	icov.triangularView<Eigen::Lower>() = icov.transpose().triangularView<Eigen::Lower>();
	
	Eigen::VectorXd abscissa(maxDims);
	if (numSpecific == 0) {
		derivCoef.resize(abilities + triangleLoc1(abilities), totalQuadPoints);
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abscissa);
			calcDerivCoef(meanVec, cov, icov, abscissa, qx);
		}
	} else {
		derivCoef.resize(primaryDims + triangleLoc1(primaryDims) + 2 * numSpecific, totalQuadPoints);
		for (int qx=0; qx < totalQuadPoints; qx++) {
			pointToLocalAbscissa(qx, abscissa);
			calcDerivCoef(meanVec, cov, icov, abscissa, qx);
			for (int curGroup=0; curGroup < numSpecific; ++curGroup) {
				calcDerivCoef1(meanVec, cov, abscissa, qx, curGroup);
			}
		}
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::cacheDerivCoef(Eigen::MatrixBase<T1> &meanVec, Eigen::MatrixBase<T2> &cov)
{
	int offset=0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		int la = layers[lx].abilities;
		Eigen::VectorXd meanVec1 = meanVec.segment(offset, la);
		Eigen::MatrixXd cov1 = cov.block(offset, offset, la, la);
		layers[lx].cacheDerivCoef(meanVec1, cov1);
		offset += la;
	}
}

template <typename T2>
void ba81NormalQuad::layer::mapLatentDeriv(double piece, int qx, Eigen::ArrayBase<T2> &derivOut)
{
	int cx = 0;
	for (int d1=0; d1 < primaryDims; ++d1) {
		double amt1 = piece * derivCoef(d1, qx);
		derivOut[d1] += amt1;
		for (int d2=0; d2 <= d1; ++d2) {
			int to = abilities + cx;
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
	int to = abilities + triangleLoc0(sdim);
	derivOut[to] += amt4;
}

template <typename T, typename T1, typename T2>
void ba81NormalQuad::layer::computeRowDeriv(int thrId, Eigen::MatrixBase<T2> &abscissa, T &op,
					    bool freeLatents, Eigen::ArrayBase<T1> &latentGradOut)
{
	abscissa.setZero();
	const int numLatents = abilities + triangleLoc1(abilities);
	Eigen::ArrayXd latentGrad(numLatents);
	latentGrad.setZero();
	const int specificPoints = quad->gridSize;

	for (int qx=0; qx < totalQuadPoints; qx++) {
		pointToGlobalAbscissa(qx, abscissa);
		op.beginQuadPoint(thrId);
		for (int ix=0; ix < op.getNumItems(); ++ix) {
			op(thrId, abscissa, Qweight(qx, thrId), ix);
		}
		op.endQuadPoint(thrId);
	}
		
	if (!freeLatents) return;

	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; qx++) {
			double tmp = Qweight(qx, thrId);
			mapLatentDeriv(tmp, qx, latentGrad);
		}
	} else {
		for (int qloc=0, eisloc=0, qx=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
			for (int sx=0; sx < specificPoints; sx++) {
				mapLatentDeriv(Qweight(qloc, thrId), qx, latentGrad);

				for (int curGroup=0; curGroup < numSpecific; curGroup++) {
					double tmp = Qweight(qloc, thrId);
					mapLatentDerivS(curGroup, tmp, qx, curGroup, latentGrad);
					++qloc;
				}
				++qx;
			}
		}
	}
	addMeanCovLocalToGlobal(latentGrad, latentGradOut);
}

template <typename T, typename T1>
void ba81NormalQuad::computeRowDeriv(int thrId, T &op, bool freeLatents, Eigen::ArrayBase<T1> &latentGrad)
{
	Eigen::VectorXd abscissa;
	abscissa.resize(abilities);
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].computeRowDeriv(thrId, abscissa, op, freeLatents, latentGrad);
	}
}

template <typename T>
void ba81NormalQuad::computeRowDeriv(int thrId, T &op)
{
	Eigen::ArrayXd placeholder;
	computeRowDeriv(thrId, op, false, placeholder);
}

template <typename T1>
void ba81NormalQuad::decodeLocation(int qx, int base, Eigen::MatrixBase<T1> &out)
{
	for (int dx=out.size()-1; dx >= 0; --dx) {
		out[dx] = qx % base;
		qx = qx / base;
	}
}

template <typename T1>
void ba81NormalQuad::layer::pointToGlobalAbscissa(int qx, Eigen::MatrixBase<T1> &abscissa)
{
	std::vector<double> &Qpoint = quad->Qpoint;
	decodeLocation(qx, quad->gridSize, abscissa);
	for (int dx=0; dx < maxDims; dx++) {
		abscissa[abilitiesOffset + dx] = Qpoint[abscissa[dx]];
	}
}

template <typename T1>
void ba81NormalQuad::layer::pointToLocalAbscissa(int qx, Eigen::MatrixBase<T1> &abscissa)
{
	std::vector<double> &Qpoint = quad->Qpoint;
	decodeLocation(qx, quad->gridSize, abscissa);
	for (int dx=0; dx < maxDims; dx++) {
		abscissa[dx] = Qpoint[abscissa[dx]];
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::foreach(Eigen::MatrixBase<T1> &abscissa, T2 &op)
{
	if (op.wantAbscissa()) abscissa.setZero();

	for (int qx=0; qx < totalQuadPoints; ++qx) {
		if (op.wantAbscissa()) {
			pointToGlobalAbscissa(qx, abscissa);
		}
		op(abscissa.derived().data());
	}
}

template <typename T>
void ba81NormalQuad::foreach(T &op)
{
	Eigen::VectorXd abscissa;
	if (op.wantAbscissa()) {
		abscissa.resize(abilities);
	}
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].foreach(abscissa, op);
	}
}

template <typename T1>
void ba81NormalQuad::layer::mapDenseSpace(double piece, const double *where,
					  const double *whereGram, Eigen::ArrayBase<T1> &latentDist)
{
	int gx = 0;
	int cx = abilities;
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
	int to = abilities + triangleLoc0(sdim);
	latentDist[to] += piece_var;
}

template <typename T1>
void ba81NormalQuad::layer::finalizeLatentDist(const double sampleSize, Eigen::ArrayBase<T1> &scorePad)
{
	const int padSize = abilities + triangleLoc1(abilities);
	for (int d1=0; d1 < padSize; d1++) {
		scorePad[d1] /= sampleSize;
	}

	int cx = abilities;
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
		scorePad[abilities + triangleLoc0(sdim)] -= ma1 * ma1;
	}
}

template <typename T1, typename T2>
void ba81NormalQuad::layer::addMeanCovLocalToGlobal(Eigen::ArrayBase<T1> &local,
						    Eigen::ArrayBase<T2> &glob)
{
	int cx = abilities;
	for (int a1=0; a1 < abilities; ++a1) {
		glob[abilitiesOffset + a1] += local[a1];
		for (int a2=0; a2 < abilities; ++a2) {
			glob[quad->abilities + triangleLoc1(abilitiesOffset+a1) +
			     abilitiesOffset + a2] += local[cx];
		}
	}
}

template <typename T1>
void ba81NormalQuad::EAP(double *thrDweight, const double sampleSize,
			 Eigen::ArrayBase<T1> &scorePad)
{
	int weightTableStart = 0;
	int offset=0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		double *layerDweight = thrDweight + weightTableStart;
		Eigen::ArrayXd layerPad(l1.abilities + triangleLoc1(l1.abilities));
		layerPad.setZero();
		if (l1.numSpecific == 0) {
			for (int qx=0; qx < l1.totalQuadPoints; ++qx) {
				l1.mapDenseSpace(layerDweight[qx], &l1.wherePrep[qx * l1.maxDims],
					      &l1.whereGram.coeffRef(0, qx), layerPad);
			}
		} else {
			int qloc=0;
			for (int qx=0; qx < l1.totalQuadPoints; qx++) {
				const double *whPrep = &l1.wherePrep[qx * l1.maxDims];
				const double *whGram = &l1.whereGram.coeffRef(0, qx);
				l1.mapDenseSpace(layerDweight[qloc], whPrep, whGram, layerPad);
				for (int Sgroup=0; Sgroup < l1.numSpecific; Sgroup++) {
					l1.mapSpecificSpace(Sgroup, layerDweight[qloc], whPrep, whGram, layerPad);
					++qloc;
				}
			}
		}

		l1.finalizeLatentDist(sampleSize, layerPad);

		l1.addMeanCovLocalToGlobal(layerPad, scorePad);

		offset += l1.abilities;
		weightTableStart += l1.weightTableSize;
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
	std::vector<int> cumItemOutcomes;
	int maxOutcomes;
	int totalOutcomes;

	// latent distribution
	double qwidth;
	int qpoints;
	ba81NormalQuad quad;
	ba81NormalQuad &getQuad() { return quad; };
	bool twotier;  // rename to detectTwoTier TODO
	int maxAbilities;
	int numSpecific;
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
	double *outcomeProb;                  // totalOutcomes * totalQuadPoints
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
	void setLatentDistribution(int dims, double *mean, double *cov);
	inline double *getItemParam(int ix) { return param + paramRows * ix; }
	inline const int *dataColumn(int col) { return dataColumns[col]; };
	void buildRowSkip();
	inline void ba81OutcomeProb(double *param, bool wantLog);
};

void ba81NormalQuad::cacheOutcomeProb(const double *ispec, double *iparam,
				      rpf_prob_t prob_fn, double *qProb)
{
	const int outcomes = ispec[RPF_ISpecOutcomes];
	int offset=0;
	Eigen::VectorXd ptheta(abilities);

	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		ptheta.setZero();

		for (int qx=0; qx < l1.totalQuadPoints; qx++) {
			double *where = l1.wherePrep.data() + qx * l1.maxDims;
			for (int dx=0; dx < l1.abilities; dx++) {
				ptheta[offset + dx] = where[std::min(dx, l1.maxDims-1)];
			}
			(*prob_fn)(ispec, iparam, ptheta.data(), qProb);
			qProb += outcomes;
		}
		offset += l1.abilities;
	}
}

// Depends on item parameters, but not latent distribution
void ifaGroup::ba81OutcomeProb(double *param, bool wantLog)
{
	outcomeProb = Realloc(outcomeProb, totalOutcomes * quad.totalQuadPoints, double);

#pragma omp parallel for num_threads(numThreads)
	for (int ix=0; ix < numItems(); ix++) {
		double *qProb = outcomeProb + cumItemOutcomes[ix] * quad.totalQuadPoints;
		const double *ispec = spec[ix];
		int id = ispec[RPF_ISpecID];
		rpf_prob_t prob_fn = wantLog? Glibrpf_model[id].logprob : Glibrpf_model[id].prob;
		quad.cacheOutcomeProb(spec[ix], param + paramRows * ix, prob_fn, qProb);
	}
}

double ba81NormalQuad::layer::computePatternLik(int thrId, double *oProb, int row)
{
	double *out = &Qweight.coeffRef(0, thrId);
	struct ifaGroup &ig = quad->ig;

	double patternLik = 0.0;
	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			out[qx] = priQarea[qx];
		}

		for (int ix=0; ix < ig.numItems(); ix++) {
			int pick = ig.dataColumns[ix][row];
			if (pick == NA_INTEGER) {
				oProb += ig.itemOutcomes[ix] * totalQuadPoints;
				continue;
			}
			pick -= 1;

			for (int qx=0; qx < totalQuadPoints; ++qx) {
				out[qx] *= oProb[pick];
				oProb += ig.itemOutcomes[ix];
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

		for (int ix=0; ix < ig.numItems(); ix++) {
			int pick = ig.dataColumns[ix][row];
			if (pick == NA_INTEGER) {
				oProb += ig.itemOutcomes[ix] * totalQuadPoints;
				continue;
			}
			pick -= 1;
			int Sgroup1 = Sgroup[ix];
			double *out1 = out;
			for (int qx=0; qx < totalQuadPoints; qx++) {
				out1[Sgroup1] *= oProb[pick];
				oProb += ig.itemOutcomes[ix];
				out1 += numSpecific;
			}
		}

		for (int qx=0; qx < totalPrimaryPoints * numSpecific; ++qx) Eis[qx] = 0;
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
	}
	return patternLik;
}

double ba81NormalQuad::computePatternLik(int thrId, int row)
{
	double patternLik = 1.0;
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layer &l1 = layers[lx];
		patternLik *= l1.computePatternLik(thrId, ig.outcomeProb, row);
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

void ba81NormalQuad::layer::addTo(int thrId, int ix, double *out)
{
	std::vector<int> &itemOutcomes = quad->ig.itemOutcomes;

	if (numSpecific == 0) {
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			*out += Qweight(qx, thrId);
			out += itemOutcomes[ix];
		}
	} else {
		int igroup = Sgroup[ix];
		double *Qw = &Qweight.coeffRef(0, thrId);
		for (int qx=0; qx < totalQuadPoints; ++qx) {
			*out += Qw[igroup];
			out += itemOutcomes[ix];
			Qw += numSpecific;
		}
	}
}

void ba81NormalQuad::addTo(int thrId, int ix, double *out)
{
	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].addTo(thrId, ix, out);
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
	void begin(class ifaGroup *state, T extraData) {};
	void addRow(class ifaGroup *state, T extraData, int px, int thrId) {};
	void recordTable(class ifaGroup *state, T extraData) {};
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

	EstepPolicy<T>::begin(state, extraData);

	quad.allocBuffers(numThreads, LatentPolicy<T>::wantSummary());

#pragma omp parallel for num_threads(numThreads)
	for (int px=0; px < numUnique; px++) {
		if (rowSkip[px]) {
			patternLik[px] = 0;
			continue;
		}

		int thrId = omp_get_thread_num();
		double patternLik1 = state->quad.computePatternLik(thrId, state->rowMap[px]);

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
		EstepPolicy<T>::addRow(state, extraData, px, thrId);
	}

	if (EstepPolicy<T>::hasEnd() && LatentPolicy<T>::hasEnd()) {
#pragma omp parallel sections
		{
			{ EstepPolicy<T>::recordTable(state, extraData); }
#pragma omp section
			{ LatentPolicy<T>::end(state, extraData); }
		}
	} else {
		EstepPolicy<T>::recordTable(state, extraData);
		LatentPolicy<T>::end(state, extraData);
	}

	quad.releaseBuffers();
}

template <typename T1, typename T2>
void ba81NormalQuad::refresh(Eigen::MatrixBase<T2> &mean, Eigen::MatrixBase<T1> &cov)
{
	if (abilities == 0) return;

	for (size_t lx=0; lx < layers.size(); ++lx) {
		layers[lx].refresh(mean, cov);
	}
}

/*TODO
	// This is required because the EM acceleration can push the
	// covariance matrix to be slightly non-pd when predictors
	// are highly correlated.
	if (priDims == 1) {
		if (cov(0,0) < BA81_MIN_VARIANCE) cov(0,0) = BA81_MIN_VARIANCE;
	} else {
		Matrix mat(cov.data(), priDims, priDims);
		InplaceForcePosSemiDef(mat, NULL, NULL);
	}

	for (int sx=0; sx < numSpecific; ++sx) {
		int loc = priDims + sx;
		double tmp = fullCov(loc, loc);
		if (tmp < BA81_MIN_VARIANCE) tmp = BA81_MIN_VARIANCE;
		sVar(sx) = tmp;
	}
*/

template <typename T1, typename T2>
void ba81NormalQuad::layer::refresh(Eigen::MatrixBase<T2> &meanVec, Eigen::MatrixBase<T1> &cov)
{
	if (abilities == 0) {
		priQarea.clear();
		priQarea.push_back(quad->One);
		return;
	}

	std::vector<double> &Qpoint = quad->Qpoint;

	Eigen::VectorXi abscissa(primaryDims);
	Eigen::MatrixXd priCov = cov.topLeftCorner(primaryDims, primaryDims);
	std::vector<double> tmpPriQarea;
	tmpPriQarea.reserve(totalPrimaryPoints);
	{
		Eigen::VectorXd where(primaryDims);
		for (int qx=0; qx < totalPrimaryPoints; qx++) {
			decodeLocation(qx, quad->gridSize, abscissa);
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
		//mxLog("%.5g,", priQarea[qx]);
	}

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
	//pda(speQarea.data(), numSpecific, quadGridSize);

	for (int sx=0; sx < int(speQarea.size()); ++sx) {
		speQarea[sx] *= quad->One;
	}
	//pda(speQarea.data(), numSpecific, quadGridSize);
}

#endif
