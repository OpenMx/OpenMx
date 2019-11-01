#ifndef __path_h_
#define __path_h_

struct PathCalcIO {
	virtual void recompute(FitContext *fc)=0;
	virtual unsigned getVersion(FitContext *fc)=0;
	virtual void refresh(FitContext *fc, Eigen::Ref<Eigen::MatrixXd> mat,
											 double sign)=0;
	virtual ~PathCalcIO() {};
};

class PathCalc {
public:
	std::vector<bool> *latentFilter; // false when latent
	std::vector<bool> *isProductNode; // change to enum?
	Eigen::MatrixXd fullA;
	Eigen::MatrixXd fullS;
	unsigned versionS;
	unsigned versionIA;
	Eigen::MatrixXd IA;  // intermediate result given boker2019=false
	unsigned versionIAF;
	Eigen::MatrixXd IAF;  // intermediate result given boker2019=false
	int numIters;
	bool boker2019;
	int numVars;
	int numObs;
	PathCalcIO *aio, *sio;
	bool algoSet;
	
	void determineShallowDepth(FitContext *fc);
	void evaluate(FitContext *fc);
	void filter();
	void prepS(FitContext *fc);

	void init()
	{
		fullA.resize(numVars, numVars);
		fullA.setZero();
		fullS.resize(numVars, numVars);
		fullS.setZero();
	}

 public:

	int verbose;
	bool ignoreVersion;

 PathCalc() :
	 versionS(0), versionIA(0), versionIAF(0), numIters(NA_INTEGER),
	 algoSet(false), verbose(0), ignoreVersion(false) {}

	void clone(PathCalc &pc)
	{
		if (!pc.algoSet) mxThrow("PathCalc::clone but setAlgo not called yet");
		numVars = pc.numVars;
		numObs = pc.numObs;
		latentFilter = pc.latentFilter;
		isProductNode = pc.isProductNode;
		aio = pc.aio;
		sio = pc.sio;
		numIters = pc.numIters;
		boker2019 = pc.boker2019;
		init();
	}

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
		aio = _aio;
		sio = _sio;
		init();
	}

	void setAlgo(FitContext *fc, bool _boker2019);

	template <typename T>
	void fullCov(FitContext *fc, Eigen::MatrixBase<T> &cov)
	{
		evaluate(fc);
		if (!boker2019) {
			prepS(fc);
			cov.derived() = IA * fullS.selfadjointView<Eigen::Lower>() * IA.transpose();
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1, typename T2>
	void fullMean(FitContext *fc, Eigen::MatrixBase<T1> &meanIn,
								Eigen::MatrixBase<T2> &meanOut)
	{
		evaluate(fc);
		if (!boker2019) {
			meanOut.derived() = IA * meanIn;
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1>
	Eigen::VectorXd fullMean(FitContext *fc, Eigen::MatrixBase<T1> &meanIn)
	{
		evaluate(fc);
		if (!boker2019) {
			return IA * meanIn; // avoids temporary copy? TODO
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T>
	void cov(FitContext *fc, Eigen::MatrixBase<T> &cov)
	{
		evaluate(fc);
		filter();
		if (!boker2019) {
			prepS(fc);
			cov.derived() = IAF * fullS.selfadjointView<Eigen::Lower>() * IAF.transpose();
		} else {
			mxThrow("Not yet");
		}
	}

	template <typename T1, typename T2>
	void mean(FitContext *fc, Eigen::MatrixBase<T1> &meanIn,
						Eigen::MatrixBase<T2> &meanOut)
	{
		evaluate(fc);
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
