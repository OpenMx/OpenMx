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

class PathCalc {
	omxState &state;
	std::vector<bool> *latentFilter; // false when latent
	std::vector<bool> *isProductNode; // change to enum?
	Eigen::MatrixXd fullA;
	Eigen::MatrixXd fullS;
	Eigen::VectorXd fullM;
	Eigen::MatrixXd IA;  // intermediate result given boker2019=false
	Eigen::MatrixXd IAF;  // intermediate result given boker2019=false
	bool boker2019;
	bool algoSelected;
	int numVars;
	int numObs;
	std::vector<coeffLoc> *aCoeff;
	std::vector<coeffLoc> *sCoeff;
	std::vector<coeffLoc> *mCoeff;
	
	void determineShallowDepth();

 public:
	int numIters; // private TODO

 PathCalc(omxState *st) :
	state(*st), algoSelected(false), numIters(NA_INTEGER) {}

	~PathCalc();

	void attach(int _numVars, int _numObs,
							std::vector<bool> &_latentFilter,
							std::vector<bool> &_isProductNode,
							std::vector<coeffLoc> &_aCoeff,
							std::vector<coeffLoc> &_sCoeff,
							std::vector<coeffLoc> &_mCoeff)
	{
		numVars = _numVars;
		numObs = _numObs;
		latentFilter = &_latentFilter;
		isProductNode = &_isProductNode;
		aCoeff = &_aCoeff;
		sCoeff = &_sCoeff;
		// Transpose row vectors
		for (auto &c1 : _mCoeff) {
			if (c1.c) std::swap(c1.r, c1.c);
		}
		mCoeff = &_mCoeff;
		fullA.resize(numVars, numVars);
		fullA.setZero();
		fullS.resize(numVars, numVars);
		fullS.setZero();
		fullM.resize(numVars);
		fullM.setZero();
	}

	void setAlgo(bool _boker2019);

	void evaluate();
	template <typename T> void fullCov(Eigen::MatrixBase<T> &cov);
	template <typename T> void fullMean(Eigen::MatrixBase<T> &mean);
	void filter();
	template <typename T> void cov(Eigen::MatrixBase<T> &cov);
	template <typename T> void mean(Eigen::MatrixBase<T> &mean);
};

template <typename T> void PathCalc::fullCov(Eigen::MatrixBase<T> &cov)
{
	if (!boker2019) {
		loadCoeffVec(*sCoeff, fullS);  // only need lower triangle TODO
		cov.derived() = IA * fullS.selfadjointView<Eigen::Lower>() * IA.transpose();
	} else {
		mxThrow("Not yet");
	}
}

template <typename T> void PathCalc::fullMean(Eigen::MatrixBase<T> &mean)
{
	if (!boker2019) {
		loadCoeffVec(*mCoeff, fullM);
		mean.derived() = IA * fullM.transpose();
	} else {
		mxThrow("Not yet");
	}
}

template <typename T> void PathCalc::cov(Eigen::MatrixBase<T> &cov)
{
	if (!boker2019) {
		loadCoeffVec(*sCoeff, fullS);  // only need lower triangle TODO
		cov.derived() = IAF * fullS.selfadjointView<Eigen::Lower>() * IAF.transpose();
	} else {
		mxThrow("Not yet");
	}
}

template <typename T> void PathCalc::mean(Eigen::MatrixBase<T> &mean)
{
	if (!boker2019) {
		loadCoeffVec(*mCoeff, fullM);
		mean.derived() = IAF * fullM;
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

	void state::analyzeModel1(FitContext *fc)
			propagateDefVar(ram, ram->getZ(fc), ram);
  Uses approximate Z (0/1 only) just to calculate dependencies

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
