#ifndef __path_h_
#define __path_h_

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include "polynomial.h"
#include "SelfAdjointEigenSolverNosort.h"

// TODO:
// integrate POV into Rampart
// optimize polynomials
//   try cholesky

// Nice if we could use RowMajor for A matrix, but sparseLU requires ColMajor
// so we do all math transposed instead.

struct PathCalcIO {
	Eigen::SparseMatrix<double> sparse;
	Eigen::MatrixXd full;
	virtual void recompute(FitContext *fc)=0;
	virtual unsigned getVersion(FitContext *fc)=0;
	virtual void refresh(FitContext *fc) {}
	virtual void refreshA(FitContext *fc, double sign) {}
	virtual void refreshSparse1(FitContext *fc, double sign) {}
	virtual PathCalcIO *clone()=0;
	virtual ~PathCalcIO() {};
	void refreshSparse(FitContext *fc, double sign)
	{
		refreshSparse1(fc, sign);
		sparse.makeCompressed();
	}
	void copyLowerToUpper()
	{
		full = full.selfadjointView<Eigen::Lower>(); // TODO optimal?
	}
};

template <typename T>
double polynomialToMoment(Polynomial< double > &polyRep, T& symEv)
{
	double erg = 0;
	for (auto &monom : polyRep.monomials) {
		double zwerg = monom.coeff;
		for (size_t ii=0; ii < monom.exponent.size(); ++ii) {
			if (monom.exponent[ii] % 2 == 1) { zwerg = 0; break; }
			for (int jj=0; jj <= (monom.exponent[ii]/2)-1; ++jj) zwerg *= 2*jj+1;
			zwerg *= pow(symEv[ii], monom.exponent[ii]/2.0);
		}
		erg += zwerg;
	}
	//std::cout << std::string(polyRep) << "\n=" << erg << "\n";
	return erg;
}

class PathCalc {
	std::vector<bool> *latentFilter; // false when latent
	std::vector<bool> *isProductNode; // change to enum?
	bool useSparse;
	unsigned versionM;
	unsigned versionS;
	unsigned versionIA;
	Eigen::MatrixXd IA;
	Eigen::SparseMatrix<double> sparseIA;
	bool sparseLUanal;
	Eigen::SparseLU< Eigen::SparseMatrix<double>,
		Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > sparseLU;
	Eigen::SparseMatrix<double> sparseIdent;
	int numIters;
	bool boker2019;
	int numVars;
	int numObs;
	Eigen::ArrayXi obsMap;
	std::unique_ptr<PathCalcIO> mio, aio, sio;
	bool algoSet;
	Eigen::SelfAdjointEigenSolverNosort< Eigen::MatrixXd > symSolver;
	std::vector< Polynomial< double > > polyRep;
	unsigned versionPoly;
	Eigen::VectorXd meanOut;
	
	void determineShallowDepth(FitContext *fc);
	void evaluate(FitContext *fc, bool filter);
	void prepM(FitContext *fc);
	void prepS(FitContext *fc);
	void prepA(FitContext *fc);
	void buildPolynomial(FitContext *fc);

	void init1();
	void init2();

	void appendPolyRep(int nn, std::vector<int> &status);

	void refreshA(FitContext *fc, double sign)
	{
		if (!useSparse) {
			aio->refreshA(fc, sign);
			if (verbose >= 2) mxPrintMat("fullA", aio->full);
		} else {
			aio->refreshSparse(fc, sign);
			if (verbose >= 2) {
				aio->full = aio->sparse;
				mxPrintMat("fullA", aio->full);
			}
		}
	}

 public:

	const int verbose;
	const bool ignoreVersion;

 PathCalc() :
	 versionM(0), versionS(0), versionIA(0), sparseLUanal(false),
	 numIters(NA_INTEGER),
	 algoSet(false), versionPoly(0), verbose(0), ignoreVersion(false) {}

	void clone(PathCalc &pc)
	{
		if (!pc.algoSet) mxThrow("PathCalc::clone but setAlgo not called yet");
		numVars = pc.numVars;
		numObs = pc.numObs;
		useSparse = pc.useSparse;
		latentFilter = pc.latentFilter;
		isProductNode = pc.isProductNode;
		if (pc.mio) mio = std::unique_ptr<PathCalcIO>(pc.mio->clone());
		aio = std::unique_ptr<PathCalcIO>(pc.aio->clone());
		sio = std::unique_ptr<PathCalcIO>(pc.sio->clone());
		numIters = pc.numIters;
		boker2019 = pc.boker2019;
		init1();
		init2();
	}

	void attach(int _numVars, int _numObs,
							std::vector<bool> &_latentFilter,
							std::vector<bool> &_isProductNode,
							PathCalcIO *_mio,
							PathCalcIO *_aio,
							PathCalcIO *_sio)
	{
		numVars = _numVars;
		numObs = _numObs;
		useSparse = numVars >= 15;
		latentFilter = &_latentFilter;
		isProductNode = &_isProductNode;
		if (_mio) mio = std::unique_ptr<PathCalcIO>(_mio);
		aio = std::unique_ptr<PathCalcIO>(_aio);
		sio = std::unique_ptr<PathCalcIO>(_sio);
	}

	void setAlgo(FitContext *fc, bool _boker2019);

	// called by omxRAMExpectation::populateAttr
	template <typename T>
	void fullCov(FitContext *fc, Eigen::MatrixBase<T> &cov)
	{
		if (!boker2019) {
			evaluate(fc, false);
			prepS(fc);
			if (!useSparse) {
				cov.derived() = IA.transpose() * sio->full.selfadjointView<Eigen::Lower>() * IA;
			} else {
				//sio->copyLowerToUpper();
				cov.derived() = sparseIA.transpose() * sio->sparse.selfadjointView<Eigen::Lower>() * sparseIA;
			}
		} else {
			buildPolynomial(fc);
			auto &symEv = symSolver.eigenvalues();
			auto &symVec = symSolver.eigenvectors();
			for (int ii=0; ii<numVars; ii++) {
				for (int jj=ii; jj<numVars; jj++) {
					auto polyProd = polyRep[ii] * polyRep[jj];
					cov(ii,jj) = polynomialToMoment(polyProd, symEv) - meanOut[ii]*meanOut[jj];
					if (ii != jj) cov(jj,ii) = cov(ii,jj);
				}
			}
		}
	}

	// called by state::computeMean
	template <typename T1>
	Eigen::VectorXd fullMean(FitContext *fc, Eigen::MatrixBase<T1> &meanIn)
	{
		if (!boker2019) {
			evaluate(fc, false);
			if (!useSparse) {
				return IA.transpose() * meanIn; // avoids temporary copy? TODO
			} else {
				return sparseIA.transpose() * meanIn; // avoids temporary copy? TODO
			}
		} else {
			buildPolynomial(fc);
			return meanOut;
		}
	}

	// called by:
	// independentGroup::computeCov
	// omxRAMExpectation::CalculateRAMCovarianceAndMeans
	template <typename T>
	void cov(FitContext *fc, Eigen::MatrixBase<T> &cov)
	{
		if (!boker2019) {
			evaluate(fc, true);
			prepS(fc);
			if (!useSparse) {
				cov.derived() = IA.transpose() * sio->full.selfadjointView<Eigen::Lower>() * IA;
			} else {
				//sio->copyLowerToUpper();
				cov.derived() = sparseIA.transpose() * sio->sparse.selfadjointView<Eigen::Lower>() * sparseIA;
			}
		} else {
			buildPolynomial(fc);
			auto &symEv = symSolver.eigenvalues();
			auto &symVec = symSolver.eigenvectors();
			for (int ii=0; ii<numVars; ii++) {
				for (int jj=ii; jj<numVars; jj++) {
					int oii = obsMap[ii];
					int ojj = obsMap[jj];
					if (oii < 0 || ojj < 0) continue;
					auto polyProd = polyRep[ii] * polyRep[jj];
					cov(oii,ojj) = polynomialToMoment(polyProd, symEv) - meanOut[ii]*meanOut[jj];
					if (oii != ojj) cov(ojj,oii) = cov(oii,ojj);
				}
			}
		}
	}

	// called by omxRAMExpectation::CalculateRAMCovarianceAndMeans
	template <typename T1>
	void mean(FitContext *fc, Eigen::MatrixBase<T1> &copyOut)
	{
		if (!boker2019) {
			evaluate(fc, true);
			prepM(fc);
			if (!useSparse) {
				copyOut.derived() = IA.transpose() * mio->full;
			} else {
				copyOut.derived() = sparseIA.transpose() * mio->full;
			}
		} else {
			buildPolynomial(fc);
			for (int vx=0; vx < numVars; ++vx) {
				int ox = obsMap[vx];
				if (ox < 0) continue;
				copyOut[ox] = meanOut[vx];
			}
		}
	}
};

#endif // __path_h_
