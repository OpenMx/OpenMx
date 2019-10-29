#ifndef __path_h_
#define __path_h_

struct coeffLoc {
	double *src;
	int r, c;

coeffLoc(double *_src, int _r, int _c) :
	src(_src), r(_r), c(_c) {}
};

template <typename T>
void loadCoeffVec(std::vector<coeffLoc> &vec, Eigen::MatrixBase<T> &dest)
{
	for (auto &v : vec) dest(v.r, v.c) = *v.src;
}

template <typename T>
void loadCoeffVecNegate(std::vector<coeffLoc> &vec, Eigen::MatrixBase<T> &dest)
{
	for (auto &v : vec) dest(v.r, v.c) = -(*v.src);
}

template <typename avft>
class PathCalc {
	omxState &state;
	std::vector<bool> *latentFilter; // false when latent
	std::vector<bool> *isProductNode; // change to enum?
	Eigen::MatrixXd fullA;
	Eigen::MatrixXd fullS;
	unsigned versionIA;
	Eigen::MatrixXd IA;  // intermediate result given boker2019=false
	unsigned versionIAF;
	Eigen::MatrixXd IAF;  // intermediate result given boker2019=false
	bool boker2019;
	int numVars;
	int numObs;
	std::unique_ptr<avft> avf;
	std::vector<coeffLoc> *aCoeff;
	std::vector<coeffLoc> *sCoeff;
	
	void determineShallowDepth();
	void evaluate();
	void filter();

 public:
	int numIters; // private TODO

 PathCalc(omxState *st) :
	 state(*st), versionIA(0), versionIAF(0), numIters(NA_INTEGER) {}

	void attach(int _numVars, int _numObs,
							std::vector<bool> &_latentFilter,
							std::vector<bool> &_isProductNode,
							avft *_avf,
							std::vector<coeffLoc> &_aCoeff,
							std::vector<coeffLoc> &_sCoeff)
	{
		numVars = _numVars;
		numObs = _numObs;
		latentFilter = &_latentFilter;
		isProductNode = &_isProductNode;
		aCoeff = &_aCoeff;
		sCoeff = &_sCoeff;
		avf = std::unique_ptr<avft>(_avf);
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
			loadCoeffVec(*sCoeff, fullS);  // only need lower triangle TODO
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
			loadCoeffVec(*sCoeff, fullS);  // only need lower triangle TODO
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


template<class avft>
void PathCalc<avft>::setAlgo(bool _boker2019)
{
	if (!_boker2019 && std::any_of(isProductNode->begin(), isProductNode->end(),
																 [](bool x){ return x; })) {
		mxThrow("Must use Boker2019 when product nodes are present");
	}
	boker2019 = _boker2019;

	if (!boker2019) {
		determineShallowDepth();
		//mxLog("PathCalc: depth %d", numIters);
	} else {
		
	}
}

template<class avft>
void PathCalc<avft>::determineShallowDepth()
{
	if (!Global->RAMInverseOpt) return;

	int maxDepth = std::min(numVars, 30);
	if (Global->RAMMaxDepth != NA_INTEGER) {
		maxDepth = std::min(maxDepth, Global->RAMMaxDepth);
	}
	loadCoeffVec(*aCoeff, fullA);
	Eigen::MatrixXd curProd = fullA;
	for (int tx=1; tx <= maxDepth; ++tx) {
		if (false) {
			mxLog("tx=%d", tx);
			mxPrintMat("curProd", curProd);
		}
		curProd = (curProd * fullA.transpose()).eval();
		if ((curProd.array() == 0.0).all()) {
			numIters = tx - 1;
			break;
		}
	}
}

template<class avft>
void PathCalc<avft>::evaluate()
{
	auto v = (*avf)();
	if (versionIA == v) {
		//mxLog("PathCalc<avft>::evaluate() in cache");
		return;
	}
	versionIA = v;

	if (!boker2019) {
		if (numIters >= 0) {
			loadCoeffVec(*aCoeff, fullA);
			// could further optimize using std::swap? (see old code)
			IA = fullA;
			IA.diagonal().array() += 1;
			for (int iter=1; iter <= numIters; ++iter) {
				IA *= fullA;
				IA.diagonal().array() += 1;
				//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
			}
		} else {
			loadCoeffVecNegate(*aCoeff, fullA);
			fullA.diagonal().array() = 1;
			Eigen::FullPivLU< Eigen::MatrixXd > lu(fullA);
			IA.resize(numVars, numVars);
			IA.setIdentity();
			IA = lu.solve(IA);
		}
	} else {
		mxThrow("not impl yet");
	}
}

template<class avft>
void PathCalc<avft>::filter()
{
	if (versionIAF == versionIA) {
		//mxLog("PathCalc<avft>::filter() in cache");
		return;
	}
	versionIAF = versionIA;

	auto &lf = *latentFilter;
	if (!boker2019) {
		IAF.resize(numObs, numVars);
		for (int rx=0, dx=0; rx < IA.rows(); ++rx) {
			if (!lf[rx]) continue;
			// maybe smarter to build transposed so we can filter by column?
			IAF.row(dx) = IA.row(rx);
			dx += 1;
		}
	} else {
		mxThrow("Not yet");
	}
}

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
