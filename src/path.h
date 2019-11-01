#ifndef __path_h_
#define __path_h_

struct PathCalcIO {
	virtual unsigned getVersion()=0;
	virtual void refresh(Eigen::Ref<Eigen::MatrixXd> mat, double sign)=0;
};

class PathCalc {
	omxState &state;   // remove? TODO
	std::vector<bool> *latentFilter; // false when latent
	std::vector<bool> *isProductNode; // change to enum?
	Eigen::MatrixXd fullA;
	Eigen::MatrixXd fullS;
	unsigned versionS;
	unsigned versionIA;
	Eigen::MatrixXd IA;  // intermediate result given boker2019=false
	unsigned versionIAF;
	Eigen::MatrixXd IAF;  // intermediate result given boker2019=false
	int numIters; // private TODO
	bool boker2019;
	int numVars;
	int numObs;
	std::unique_ptr<PathCalcIO> aio, sio;
	
	void determineShallowDepth();
	void evaluate();
	void filter();
	void prepS();

 public:

 PathCalc(omxState *st) :
	 state(*st), versionS(0), versionIA(0), versionIAF(0), numIters(NA_INTEGER) {}

	void attach(int _numVars, int _numObs,
							std::vector<bool> &_latentFilter,
							std::vector<bool> &_isProductNode,
							PathCalcIO *_aio,
							PathCalcIO *_sio)
	{
		numVars = _numVars;
		numObs = _numObs;
		latentFilter = &_latentFilter;
		isProductNode = &_isProductNode;
		aio = std::unique_ptr<PathCalcIO>(_aio);
		sio = std::unique_ptr<PathCalcIO>(_sio);
		fullA.resize(numVars, numVars);
		fullA.setZero();
		fullS.resize(numVars, numVars);
		fullS.setZero();
	}

	void setAlgo(bool _boker2019);

	template <typename T>
	void fullCov(Eigen::MatrixBase<T> &cov)
	{
		evaluate();
		if (!boker2019) {
			prepS();
			cov.derived() = IA * fullS.selfadjointView<Eigen::Lower>() * IA.transpose();
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1, typename T2>
	void fullMean(Eigen::MatrixBase<T1> &meanIn,
								Eigen::MatrixBase<T2> &meanOut)
	{
		evaluate();
		if (!boker2019) {
			meanOut.derived() = IA * meanIn;
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1>
	Eigen::VectorXd fullMean(Eigen::MatrixBase<T1> &meanIn)
	{
		evaluate();
		if (!boker2019) {
			return IA * meanIn; // avoids temporary copy? TODO
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T>
	void cov(Eigen::MatrixBase<T> &cov)
	{
		evaluate();
		filter();
		if (!boker2019) {
			prepS();
			cov.derived() = IAF * fullS.selfadjointView<Eigen::Lower>() * IAF.transpose();
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1, typename T2>
	void mean(Eigen::MatrixBase<T1> &meanIn,
						Eigen::MatrixBase<T2> &meanOut)
	{
		evaluate();
		filter();
		if (!boker2019) {
			meanOut.derived() = IAF * meanIn;
		} else {
			mxThrow("Not yet");
		}
	}
};

/*
	NOTES:

void omxRAMExpectation::CalculateRAMCovarianceAndMeans(FitContext *fc) [DONE]
  needs some rewriting because getZ(NULL) then filters Z
  and uses omxDGEMM to do the quadratic product FZSZ'F'
  hence, unfiltered covariance is optimized out
  Maybe need a special API to keep this optimization.

refreshUnfilteredCov(omxExpectation *oo) [DONE]
   reads directly from A, S
    omxMatrix* Z = oro->getZ(NULL);
   writes to Ax

omxMatrix *omxRAMExpectation::getZ(FitContext *fc) --> TO REMOVE
  updates Z using omxShallowInverse conditional on Zversion

	void state::computeMean(FitContext *fc)
   direct input
			omxRecompute(ram->A, fc);
			EigenMatrixAdaptor eZ(ram->getZ(fc));
   mean only output

	void independentGroup::refreshUnitA(FitContext *fc, int px)
    asymT.fullA -- block diagonal aggregation of group units + rampart adjusted
    input entry by entry into sparse matrix

	void independentGroup::determineShallowDepth(FitContext *fc)
    asymT.fullA -- first builds whole A then
		asymT.determineShallowDepth(fc);

	independentGroup::independentGroup(independentGroup *ig)
	& void independentGroup::prep(FitContext *fc)
    asymT initialization

	void independentGroup::computeCov2()
   output dense, filtered covariance
 */

#endif // __path_h_
