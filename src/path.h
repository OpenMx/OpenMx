#ifndef u__path_h_
#define u__path_h_

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
	int useSparse;
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
	bool algoSet;
	Eigen::SelfAdjointEigenSolverNosort< Eigen::MatrixXd > symSolver;
	std::vector< Polynomial< double > > polyRep;
	unsigned versionPoly;
	Eigen::VectorXd meanOut;
  omxMatrix *fullMeanAccess;
  omxMatrix *fullCovAccess;
  Eigen::MatrixXd tmpFullCov;
  Eigen::VectorXd tmpFullMean;
  omxMatrix *selVecCov;
  omxMatrix *selVecMean;

  struct selPlanRow {
    int step;
    int from;
    int to;
  };
  std::vector<selPlanRow> selPlanCov;
  std::vector<selPlanRow> selPlanMean;

  struct selStep {
    std::vector<bool> selFilter;
    int selDim;
    Eigen::MatrixXd selAdj;
  };
  std::vector<selStep> selStepsCov;
  std::vector<selStep> selStepsMean;

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

  template <typename T1>
  void pearsonSelCov1(Eigen::MatrixBase<T1> &cov);
  void pearsonSelMean1(Eigen::Ref<Eigen::VectorXd> mean);

 public:

	std::unique_ptr<PathCalcIO> mio, aio, sio;
	const int verbose;
	const bool ignoreVersion;

 PathCalc() :
	 useSparse(false), versionM(0), versionS(0), versionIA(0), sparseLUanal(false),
	 numIters(NA_INTEGER),
	 algoSet(false), versionPoly(0), fullMeanAccess(0), fullCovAccess(0),
   selVecCov(0), selVecMean(0), verbose(0), ignoreVersion(false) {}

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
    selVecCov = pc.selVecCov;
	selVecMean = pc.selVecMean;
	
    selPlanCov = pc.selPlanCov;
	selPlanMean = pc.selPlanMean;
	
    selStepsCov = pc.selStepsCov;
	selStepsMean = pc.selStepsMean;
	
    fullMeanAccess = pc.fullMeanAccess;
    fullCovAccess = pc.fullCovAccess;
		init1();
		init2();
	}

	void attach(int u_numVars, int u_numObs,
							std::vector<bool> &u_latentFilter,
							std::vector<bool> &u_isProductNode,
							PathCalcIO *u_mio,
							PathCalcIO *u_aio,
							PathCalcIO *u_sio)
	{
		numVars = u_numVars;
		numObs = u_numObs;
		latentFilter = &u_latentFilter;
		isProductNode = &u_isProductNode;
		if (u_mio) mio = std::unique_ptr<PathCalcIO>(u_mio);
		aio = std::unique_ptr<PathCalcIO>(u_aio);
		sio = std::unique_ptr<PathCalcIO>(u_sio);
	}

  void attachSelectionCov(omxMatrix *u_selVec, DataFrame u_selPlan)
  {
    if (u_selPlan.nrows() == 0) return;
    selPlanCov.resize(u_selPlan.nrows());

    selVecCov = u_selVec;
    IntegerVector step = u_selPlan["step"];
    {
      int selCount = 1;
      int prevStep = step[0];
      for (int sx=1; sx < step.length(); ++sx) {
        if (step[sx] == prevStep) continue;
        prevStep = step[sx];
        ++selCount;
      }
      selStepsCov.resize(selCount);
      for (auto &s1 : selStepsCov) s1.selFilter.assign(numVars, false);
    }

    int curStep = step[0];
    for (int rx=0, sx=0; rx < u_selPlan.nrows(); ++rx) {
      IntegerVector step = u_selPlan["step"];
      IntegerVector from = u_selPlan["from"];
      IntegerVector to = u_selPlan["to"];
      auto &spr = selPlanCov[rx];
      spr.step = step[rx];
      spr.from = from[rx];
      spr.to = to[rx];
      auto &s1 = selStepsCov[sx];
      s1.selFilter[ from[rx] ] = true;
      s1.selFilter[ to[rx] ] = true;
      if (rx == u_selPlan.nrows()-1 || step[rx+1] != curStep) {
        s1.selDim = std::accumulate(s1.selFilter.begin(), s1.selFilter.end(), 0);
        if (rx < u_selPlan.nrows()-1) {
          curStep = step[rx+1];
          ++sx;
        }
      }
    }
  }

  void attachFullMemory(omxMatrix *u_mean, omxMatrix *u_cov)
  {
    fullMeanAccess = u_mean;
    fullCovAccess = u_cov;
  }

	void setAlgo(FitContext *fc, bool u_boker2019, int u_useSparse);

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
      if (selStepsCov.size()) pearsonSelCov1(cov);
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
      Eigen::VectorXd meanOut;
			if (!useSparse) {
				meanOut = IA.transpose() * meanIn; // avoids temporary copy? TODO
			} else {
				meanOut = sparseIA.transpose() * meanIn; // avoids temporary copy? TODO
			}
      //if (selSteps.size()) pearsonSelMean1(meanOut);
      return meanOut;
		} else {
			buildPolynomial(fc);
			return meanOut;
		}
	}

	// called by:
	// independentGroup::computeCov
	// omxRAMExpectation::CalculateRAMCovarianceAndMeans
	template <typename T>
	void cov(FitContext *fc, Eigen::MatrixBase<T> &cov);

	// called by omxRAMExpectation::CalculateRAMCovarianceAndMeans
	template <typename T1>
	void mean(FitContext *fc, Eigen::MatrixBase<T1> &copyOut);

	std::string getPolyRep()
	{
		auto &symEv = symSolver.eigenvalues();
		std::ostringstream temp;
		for (int ii=0; ii < numVars; ++ii) {
			temp << "[" << ii << "] " << symEv[ii] << " : " << std::string(polyRep[ii]) << "\n";
		}
		return temp.str();
	}

};

#include "matrix.h"
#include "Compute.h"

template <typename T1>
void PathCalc::pearsonSelCov1(Eigen::MatrixBase<T1> &cov)
{
  int rx=0;
  for (auto &s1 : selStepsCov) {
  //mxPrintMat("before sel", cov);
    Eigen::MatrixXd v11(s1.selDim, s1.selDim);
    Eigen::MatrixXd v12(s1.selDim, cov.cols() - s1.selDim);
    Eigen::MatrixXd v22(cov.rows() - s1.selDim, cov.cols() - s1.selDim);
    partitionCovariance(cov, [&](int xx){ return s1.selFilter[xx]; }, v11, v12, v22);
    EigenVectorAdaptor EselVec(selVecCov);
    int curStep = selPlanCov[rx].step;
    while (rx < int(selPlanCov.size()) && selPlanCov[rx].step == curStep) {
      int from = selPlanCov[rx].from;
      int to = selPlanCov[rx].to;
      cov(from, to) = EselVec[rx];
      cov(to, from) = EselVec[rx];
      ++rx;
    }
    Eigen::MatrixXd nc(s1.selDim, s1.selDim);
    subsetCovariance(cov, [&](int x)->bool{ return s1.selFilter[x]; }, s1.selDim, nc);
    Eigen::MatrixXd iv11(v11);
    if (InvertSymmetricPosDef(iv11, 'L')) {
      // complain TODO
      return;
    }
    iv11 = iv11.selfadjointView<Eigen::Lower>();
    s1.selAdj = iv11 * v12; // will use for mean too
    Eigen::MatrixXd n12 = nc * s1.selAdj;
    Eigen::MatrixXd n22 = v22 - v12.transpose() * (iv11 - iv11 * nc * iv11) * v12;
    partitionCovarianceSet(cov, [&](int xx){ return s1.selFilter[xx]; }, nc, n12, n22);
    //mxPrintMat("after sel", cov);
  }
}

template <typename T1>
void PathCalc::mean(FitContext *fc, Eigen::MatrixBase<T1> &copyOut)
{
  if (!boker2019) {
    prepM(fc);
    if (selStepsMean.size()) {
      if (!fullMeanAccess) tmpFullMean.resize(numVars);
      omxMatrix *fma = fullMeanAccess;
      if (fc) fma = fc->state->lookupDuplicate(fullMeanAccess);
      Eigen::Map< Eigen::VectorXd > tmpMean(fma? fma->data : tmpFullMean.data(), numVars);
      tmpMean = fullMean(fc, mio->full);
      subsetVector(tmpMean, [&](int cx)->bool{ return (*latentFilter)[cx]; },
                   numObs, copyOut);
    } else {
      evaluate(fc, true);
      if (!useSparse) {
        copyOut.derived() = IA.transpose() * mio->full;
      } else {
        copyOut.derived() = sparseIA.transpose() * mio->full;
      }
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

template <typename T>
void PathCalc::cov(FitContext *fc, Eigen::MatrixBase<T> &cov)
{
  if (!boker2019) {
    if (selStepsCov.size()) {
      if (!fullCovAccess) tmpFullCov.resize(numVars, numVars);
      omxMatrix *fca = fullCovAccess;
      if (fc) fc->state->lookupDuplicate(fullCovAccess);
      Eigen::Map< Eigen::MatrixXd > tmpCov(fca? fca->data : tmpFullCov.data(),
                                           numVars, numVars);
      fullCov(fc, tmpCov);
      subsetCovariance(tmpCov, [&](int cx)->bool{ return (*latentFilter)[cx]; },
                       numObs, cov);
    } else {
      evaluate(fc, true);
      prepS(fc);
      if (!useSparse) {
        cov.derived() = IA.transpose() * sio->full.selfadjointView<Eigen::Lower>() * IA;
      } else {
        //sio->copyLowerToUpper();
        cov.derived() = sparseIA.transpose() * sio->sparse.selfadjointView<Eigen::Lower>() * sparseIA;
      }
    }
  } else {
    cov.derived().resize(numObs, numObs);
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

#endif // u__path_h_
