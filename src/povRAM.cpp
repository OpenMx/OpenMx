#include <ostream>
#include "omxExpectation.h"
#include "polynomial.h"
#include "path.h"
#include <Eigen/Eigenvalues>
#include "SelfAdjointEigenSolverNosort.h"
#include "EnableWarnings.h"

void PathCalc::prepS(FitContext *fc)
{
	sio->recompute(fc);
	if (ignoreVersion || versionS != sio->getVersion(fc)) {
		sio->refresh(fc, fullS, 0);
		versionS = sio->getVersion(fc);
	}
	if (verbose) mxPrintMat("S", fullS);
}

void PathCalc::setAlgo(FitContext *fc, bool _boker2019)
{
	if (!_boker2019 && std::any_of(isProductNode->begin(), isProductNode->end(),
																 [](bool x){ return x; })) {
		mxThrow("Must use Boker2019 when product nodes are present");
	}
	boker2019 = _boker2019;

	if (!boker2019) {
		determineShallowDepth(fc);
		//mxLog("PathCalc: depth %d", numIters);
	} else {
		
	}
	algoSet = true;
}

void PathCalc::determineShallowDepth(FitContext *fc)
{
	if (!Global->RAMInverseOpt) return;

	int maxDepth = std::min(numVars, 30);
	if (Global->RAMMaxDepth != NA_INTEGER) {
		maxDepth = std::min(maxDepth, Global->RAMMaxDepth);
	}
	aio->refresh(fc, fullA, 1);
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

void PathCalc::evaluate(FitContext *fc)
{
	aio->recompute(fc);
	if (!ignoreVersion && versionIA == aio->getVersion(fc)) {
		//mxLog("PathCalc<avft>::evaluate() in cache");
		return;
	}
	versionIA = aio->getVersion(fc);

	if (!boker2019) {
		if (numIters >= 0) {
			aio->refresh(fc, fullA, 1.0);
			if (verbose) mxPrintMat("A", fullA);
			// could further optimize using std::swap? (see old code)
			IA = fullA;
			IA.diagonal().array() += 1;
			for (int iter=1; iter <= numIters; ++iter) {
				IA *= fullA;
				IA.diagonal().array() += 1;
				//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
			}
		} else {
			aio->refresh(fc, fullA, -1.0);
			if (verbose) mxPrintMat("A", fullA);
			fullA.diagonal().array() = 1;
			Eigen::FullPivLU< Eigen::MatrixXd > lu(fullA);
			IA.resize(numVars, numVars);
			IA.setIdentity();
			IA = lu.solve(IA);
		}
		if (verbose) mxPrintMat("IA", IA);
	} else {
		mxThrow("not impl yet");
	}
}

void PathCalc::filter()
{
	if (!ignoreVersion && versionIAF == versionIA) {
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
		if (verbose) mxPrintMat("new IAF", IAF);
	} else {
		mxThrow("Not yet");
	}
}


// Sparse matrix version below. This is probably
// not worth the trouble because we only do this
// once. Even if it is faster, it might be easier
// to rely on the simpler, dense matrix code.

/*
template <typename T1> void AsymTool<T1>::determineShallowDepth()
{
	int maxDepth = std::min(fullA.cols(), 30);
	if (Global->RAMMaxDepth != NA_INTEGER) maxDepth = std::min(maxDepth, Global->RAMMaxDepth);
	Eigen::SparseMatrix<double> curProd = fullA;
	for (int tx=1; tx <= maxDepth; ++tx) {
		if (false) {
			Eigen::MatrixXd tmp = curProd;
			mxPrintMat("curProd", tmp);
		}
		curProd = (curProd * fullA.transpose()).eval();
		bool allZero = true;
		for (int k=0; k < curProd.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(curProd, k); it; ++it) {
				if (it.value() != 0.0) {
					allZero = false;
					break;
				}
			}
		}
		if (allZero) {
			AshallowDepth = tx - 1;
			break;
		}
	}
	fullA.setZero();

	if (AshallowDepth >= 0) signA = 1.0;
	hasDeterminedDepth = true;
}

	void evaluate()
	{
	if (AshallowDepth >= 0) {
		fullA.makeCompressed();
		IAF = fullA + ident;
		for (int iter=1; iter <= AshallowDepth; ++iter) {
			IAF = (IAF * fullA + ident).eval();
			//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
		}
	} else {
		fullA += ident;
		if (!analyzed) {
			analyzed = true;
			fullA.makeCompressed();
			Asolver.analyzePattern(fullA);
		}
		Asolver.factorize(fullA);
		if (Asolver.info() != Eigen::Success) {
			mxThrow("Failed to invert flattened A matrix; %s",
				 Asolver.lastErrorMessage().c_str());
		}

		IAF = Asolver.solve(ident);
		fullA -= ident;  // leave unchanged
		//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
	}
	filtered = false;
	++invertCount;
}
*/

class povRAMExpectation : public omxExpectation {
	typedef omxExpectation super;
	unsigned Zversion;
	omxMatrix *_Z, *I, *Ax;
	Eigen::VectorXi dataCols;  // composition of F permutation and expectation->dataColumns
	std::vector<const char *> dataColNames;
	std::vector< omxThresholdColumn > thresholds;
	std::vector<bool> isProductNode;

	void appendPolyRep(int nn, std::vector<int> &status,
										 std::vector< Polynomial< double > > &polyRep);

 public:
	std::vector<bool> latentFilter; // false when latent

	povRAMExpectation(omxState *st, int num) : super(st, num), Zversion(0), _Z(0) {};
	virtual ~povRAMExpectation();

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M;

	int verbose;

	void studyF();

	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual omxMatrix *getComponent(const char*);
	virtual const std::vector<const char *> &getDataColumnNames() const { return dataColNames; };
	virtual const Eigen::Map<DataColumnIndexVector> getDataColumns() {
		return Eigen::Map<DataColumnIndexVector>(dataCols.data(), numDataColumns);
	}
	virtual std::vector< omxThresholdColumn > &getThresholdInfo() { return thresholds; }
};

omxExpectation *povRAMExpectationInit(omxState *st, int num)
{ return new povRAMExpectation(st, num); }

povRAMExpectation::~povRAMExpectation()
{
	omxFreeMatrix(I);
	omxFreeMatrix(_Z);
	omxFreeMatrix(Ax);
	omxFreeMatrix(cov);
	omxFreeMatrix(means);
}

void povRAMExpectation::init()
{
	canDuplicate = true;

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose) + OMX_DEBUG;

	M = omxNewMatrixFromSlot(rObj, currentState, "M"); // can be optional? TODO
	A = omxNewMatrixFromSlot(rObj, currentState, "A");
	S = omxNewMatrixFromSlot(rObj, currentState, "S");
	F = omxNewMatrixFromSlot(rObj, currentState, "F");

	int k = A->rows;
	I = omxNewIdentityMatrix(k, currentState);
	_Z = omxInitMatrix(k, k, TRUE, currentState);
	Ax = omxInitMatrix(k, k, TRUE, currentState);

	ProtectedSEXP RprodNode(R_do_slot(rObj, Rf_install("isProductNode")));
	if (Rf_length(RprodNode) != A->cols) {
		mxThrow("isProductNode must be same dimension as A matrix");
	}
	isProductNode.assign(A->cols, false);
	for (int px = 0; px < A->cols; ++px) {
		if (INTEGER(RprodNode)[px]) isProductNode[px] = true;
	}

	int l = F->rows;
	cov = 		omxInitMatrix(l, l, TRUE, currentState);
	means = 	omxInitMatrix(1, l, TRUE, currentState);

	studyF();
}

void povRAMExpectation::studyF()
{
	// Permute the data columns such that the manifest
	// part of F is a diagonal matrix. This permits
	// trivial filtering of the latent variables.
	auto dataColumns = super::getDataColumns();
	auto origDataColumnNames = super::getDataColumnNames();
	auto origThresholdInfo = super::getThresholdInfo();
	EigenMatrixAdaptor eF(F);
	latentFilter.assign(eF.cols(), false);
	dataCols.resize(eF.rows());
	dataColNames.resize(eF.rows(), 0);
	if (!eF.rows()) return;  // no manifests
	for (int cx =0, dx=0; cx < eF.cols(); ++cx) {
		int dest;
		double isManifest = eF.col(cx).maxCoeff(&dest);
		latentFilter[cx] = isManifest;
		if (isManifest) {
			dataColNames[dx] = origDataColumnNames[dest];
			int newDest = dataColumns.size()? dataColumns[dest] : dest;
			dataCols[dx] = newDest;
			if (origThresholdInfo.size()) {
				omxThresholdColumn adj = origThresholdInfo[dest];
				adj.dColumn = dx;
				thresholds.push_back(adj);
			}
			dx += 1;
		}
	}
}

omxMatrix* povRAMExpectation::getComponent(const char* component)
{
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = cov;
	} else if(strEQ("means", component)) {
		retval = means;
	} else if(strEQ("slope", component)) {
		mxThrow("slope not supported");
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}
	
	return retval;
}

void povRAMExpectation::appendPolyRep(int nn, std::vector<int> &status,
																			std::vector< Polynomial< double > > &polyRep)
{
	EigenMatrixAdaptor eA(A);
	if (status[nn] == 2) return;
	if (status[nn] == 1) mxThrow("Asymmetric matrix is cyclic");
	status[nn] = 1;
	
	for (int ii=0; ii < A->rows; ++ii) {
		if (ii == nn || status[ii] == 2 || eA(nn,ii) == 0) continue;
		appendPolyRep(ii, status, polyRep);
	}
	for (int ii=0; ii < A->rows; ++ii) {
		if (ii == nn || eA(nn,ii) == 0) continue;
		Polynomial< double > term(eA(nn,ii));
		//mxLog("A %d %d %f", ii,nn,eA(nn,ii));
		//std::cout << std::string(polyRep[ii]) << "\n";
		term *= polyRep[ii];
		//std::cout << std::string(polyRep[nn]) << "OP " << isProductNode[nn] << " " << std::string(term) << "\n";
		if (isProductNode[nn]) {
			polyRep[nn] *= term;
		} else {
			polyRep[nn] += term;
		}
		//std::cout << "result: " << std::string(polyRep[nn]) << "\n";
	}

	status[nn] = 2;
}

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

void povRAMExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	if (F->rows == 0) return;

	omxRecompute(A, fc);
	omxRecompute(S, fc);
	omxRecompute(F, fc);
	if (M) omxRecompute(M, fc);  // currently required TODO

	EigenMatrixAdaptor eS(S);
	EigenVectorAdaptor eM(M);
	//mxPrintMat("S", eS);

	std::vector< Polynomial< double > > polyRep(S->rows);
	for (int ii=0; ii < S->rows; ++ii) {
		polyRep[ii].addMonomial(Monomial< double >(eM(ii)));
	}

        /* old version with Cholesky, doesn't work with non-positive definite S, but probably easier to extend to derivatives. Also, 
         * polynomials are more sparse with the Cholesky variant. An alternative could be to use a sign switch in the Cholesky decomposition,
         * thereby working with normally distributed variables of variance either 1 or -1. 
        double[][] cholesky = null;
        try {
            cholesky = Statik.choleskyDecompose(symVal,0.001);
        } catch (Exception e) {
            throw new RuntimeException("Polynomial representation is impossible because symmetric matrix is not positive definite: "+e);
        }

        int[] x = new int[anzFac]; int k = 0; for (int i=0; i<anzFac; i++) if (cholesky[i][i]!=0.0) x[i] = k++;
        for (int i=0; i<anzFac; i++) for (int j=0; j<=i; j++) if (cholesky[i][j] != 0.0) polynomialRepresentation[i].addMonomial(cholesky[i][j], x[j]);
        polynomialRepresentationVariances = Statik.ensureSize(polynomialRepresentationVariances, k);
        for (int i=0; i<polynomialRepresentationVariances.length; i++) polynomialRepresentationVariances[i] = 1.0;
        */

	// Add option to use Cholesky with fallback to SelfAdjointEigenSolver TODO
	Eigen::SelfAdjointEigenSolverNosort< Eigen::MatrixXd > sym(eS);
	auto &symEv = sym.eigenvalues();
	auto &symVec = sym.eigenvectors();
	for (int jj=0; jj < S->rows; ++jj) {
		if (symEv(jj) == 0) continue;
		for (int ii=0; ii < S->rows; ++ii) {
			if (symVec(ii,jj) == 0) continue;
			polyRep[ii].addMonomial(symVec(ii,jj), jj);
		}
	}

	// for (int ii=0; ii < S->rows; ++ii) {
	// 	std::cout << ii << ":" << std::string(polyRep[ii]) << "\n";
	// }

	std::vector<int> status(S->rows, 0);
	for (int ii=0; ii<S->rows; ii++) {
		appendPolyRep(ii, status, polyRep);
	}

	// mxPrintMat("vec", symVec);
	 for (int ii=0; ii < S->rows; ++ii) {
	 	std::cout << ii << " " << symEv[ii] << ":" << std::string(polyRep[ii]) << "\n";
	 }

	EigenMatrixAdaptor sigmaBig(Ax);
	Eigen::VectorXd fullMean(S->rows);
	// If don't need covariance then can do mean only for observed variables
	for (int ii=0; ii<S->rows; ++ii) {
		fullMean[ii] = polynomialToMoment(polyRep[ii], symEv);
	}
	// Can do only observed covariance if full is not needed
	for (int ii=0; ii<S->rows; ii++) {
		for (int jj=ii; jj<S->rows; jj++) {
			auto polyProd = polyRep[ii] * polyRep[jj];
			sigmaBig(ii,jj) = polynomialToMoment(polyProd, symEv) - fullMean[ii]*fullMean[jj];
		}
	}
	sigmaBig.derived() = sigmaBig.selfadjointView<Eigen::Upper>();

  //mxPrintMat("full cov", sigmaBig);
	//mxPrintMat("full mean", fullMean);

	//mxThrow("stop");

	EigenMatrixAdaptor eCov(cov);
	EigenVectorAdaptor eMeans(means);
	subsetCovariance(sigmaBig, [&](int x){ return latentFilter[x]; }, F->rows, eCov);
	subsetVector(fullMean, [&](int x){ return latentFilter[x]; }, F->rows, eMeans);
	//mxPrintMat("cov", eCov);
	//mxPrintMat("mean", eMeans);
}
