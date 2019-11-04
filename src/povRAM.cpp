#include <ostream>
#include "omxExpectation.h"
#include "path.h"
#include "EnableWarnings.h"

void PathCalc::prepM(FitContext *fc)
{
	if (!mio) mxThrow("PathCalc::prepM but no PathCalcIO for mean");
	mio->recompute(fc);
	if (ignoreVersion || versionM != mio->getVersion(fc)) {
		mio->refresh(fc, fullM, 0);
		versionM = mio->getVersion(fc);
	}
	if (verbose) mxPrintMat("M", fullM);
}

void PathCalc::prepS(FitContext *fc)
{
	sio->recompute(fc);
	if (ignoreVersion || versionS != sio->getVersion(fc)) {
		sio->refresh(fc, fullS, 0);
		versionS = sio->getVersion(fc);
	}
	if (verbose) mxPrintMat("S", fullS);
}

void PathCalc::prepA(FitContext *fc)
{
	aio->recompute(fc);
	if (ignoreVersion || versionIA != aio->getVersion(fc)) {
		aio->refresh(fc, fullA, 1.0);
		versionIA = aio->getVersion(fc);
	}
	if (verbose) mxPrintMat("A", fullA);
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
	if (boker2019) mxThrow("PathCalc::evaluate but boker2019=TRUE");

	aio->recompute(fc);
	if (!ignoreVersion && versionIA == aio->getVersion(fc)) {
		//mxLog("PathCalc<avft>::evaluate() in cache");
		return;
	}
	versionIA = aio->getVersion(fc);

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
}

void PathCalc::filter()
{
	if (!ignoreVersion && versionIAF == versionIA) {
		//mxLog("PathCalc<avft>::filter() in cache");
		return;
	}
	versionIAF = versionIA;

	auto &lf = *latentFilter;
	IAF.resize(numObs, numVars);
	for (int rx=0, dx=0; rx < IA.rows(); ++rx) {
		if (!lf[rx]) continue;
		// maybe smarter to build transposed so we can filter by column?
		IAF.row(dx) = IA.row(rx);
		dx += 1;
	}
	if (verbose) mxPrintMat("new IAF", IAF);
}

void PathCalc::appendPolyRep(int nn, std::vector<int> &status)
{
	if (status[nn] == 2) return;
	if (status[nn] == 1) mxThrow("Asymmetric matrix is cyclic");
	status[nn] = 1;
	
	for (int ii=0; ii < fullA.rows(); ++ii) {
		if (ii == nn || status[ii] == 2 || fullA(nn,ii) == 0) continue;
		appendPolyRep(ii, status);
	}
	for (int ii=0; ii < fullA.rows(); ++ii) {
		if (ii == nn || fullA(nn,ii) == 0) continue;
		Polynomial< double > term(fullA(nn,ii));
		//mxLog("A %d %d %f", ii,nn,fullA(nn,ii));
		//std::cout << std::string(polyRep[ii]) << "\n";
		term *= polyRep[ii];
		//std::cout << std::string(polyRep[nn]) << "OP " << isProductNode[nn] << " " << std::string(term) << "\n";
		if ((*isProductNode)[nn]) {
			polyRep[nn] *= term;
		} else {
			polyRep[nn] += term;
		}
		//std::cout << "result: " << std::string(polyRep[nn]) << "\n";
	}

	status[nn] = 2;
}

void PathCalc::buildPolynomial(FitContext *fc)
{
	prepS(fc);
	prepA(fc);
	unsigned curV = versionS + versionIA;
	if (mio) {
		prepM(fc);
		curV += versionM;
	}
	if (!ignoreVersion && versionPoly == curV) return;
	versionPoly = curV;

	for (auto &p1 : polyRep) p1.clear();

	if (mio) {
		for (int ii=0; ii < numVars; ++ii) {
			polyRep[ii].addMonomial(Monomial< double >(fullM(ii)));
		}
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
	symSolver.compute(fullS);
	auto &symEv = symSolver.eigenvalues();
	auto &symVec = symSolver.eigenvectors();
	for (int jj=0; jj < numVars; ++jj) {
		if (symEv(jj) == 0) continue;
		for (int ii=0; ii < numVars; ++ii) {
			if (symVec(ii,jj) == 0) continue;
			polyRep[ii].addMonomial(symVec(ii,jj), jj);
		}
	}

	// for (int ii=0; ii < numVars; ++ii) {
	// 	std::cout << ii << ":" << std::string(polyRep[ii]) << "\n";
	// }

	std::vector<int> status(numVars, 0);
	for (int ii=0; ii<numVars; ii++) {
		appendPolyRep(ii, status);
	}

	// mxPrintMat("vec", symVec);
	if (false) {
		for (int ii=0; ii < numVars; ++ii) {
			std::cout << ii << " " << symEv[ii] << ":" << std::string(polyRep[ii]) << "\n";
		}
	}

	// Could be smarter and avoid latents except when needed TODO
	meanOut.resize(numVars);
	for (int ii=0; ii<numVars; ++ii) {
		meanOut[ii] = polynomialToMoment(polyRep[ii], symEv);
	}
}

// Sparse matrix version below. Need to integrate this again.

/*
template <typename T1>
class AsymTool {
	int AshallowDepth;
	double signA;
	std::vector<T1> &latentFilter;
	Eigen::SparseMatrix<double>      ident;
	bool analyzed;
	bool hasDeterminedDepth;
	Eigen::SparseLU< Eigen::SparseMatrix<double>,
		Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > Asolver;
	//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;
	int clumpObs;
	bool filtered;
	int invertCount;
	int filterCount;
 public:
	bool determinedDepth() const { return hasDeterminedDepth; };
	bool isFiltered() const { return filtered; };
	int getDepth() const { return AshallowDepth; };
	void setDepth(int depth) {
		AshallowDepth = depth;
		hasDeterminedDepth = true;
		if (depth >= 0) signA = 1.0;
	};
	double getSign() const { return signA; };
	int getInvertCount() const { return invertCount; };
	int getFilterCount() const { return filterCount; };

	Eigen::SparseMatrix<double> fullA;
	Eigen::SparseMatrix<double> IAF;

	AsymTool(std::vector<T1> &_latentFilter) :
	AshallowDepth(-1), signA(-1.0), latentFilter(_latentFilter), analyzed(false), hasDeterminedDepth(false),
		invertCount(0), filterCount(0) {};
	void resize(int clumpVars, int _clumpObs)
	{
		fullA.resize(clumpVars, clumpVars);
		ident.resize(clumpVars, clumpVars);
		ident.setIdentity();
		clumpObs = _clumpObs;
	};
	void determineShallowDepth();
	void invert();
	void filter();
};

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

template <typename T1> void AsymTool<T1>::invert()
{
	// consider http://users.clas.ufl.edu/hager/papers/Lightning/update.pdf ?
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

template <typename T1> void AsymTool<T1>::filter()
{
	const bool doubleCheck = false;
	Eigen::MatrixXd denseA;
	if (doubleCheck) {
		denseA = IAF;
	}

	// We built A transposed so we can quickly filter columns
	// Switch to filterOuter http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1130 TODO
	IAF.uncompress();
	Eigen::SparseMatrix<double>::Index *op = IAF.outerIndexPtr();
	Eigen::SparseMatrix<double>::Index *nzp = IAF.innerNonZeroPtr();
	int dx = 0;
	for (int cx=0; cx < fullA.cols(); ++cx) {
		if (!latentFilter[cx]) continue;
		op[dx] = op[cx];
		nzp[dx] = nzp[cx];
		++dx;
	}
	op[dx] = op[fullA.cols()];
	IAF.conservativeResize(fullA.rows(), clumpObs);

	// I've screwed this up 3-4 times, better check it!
	if (OMX_DEBUG && dx != clumpObs) mxThrow("latentFilter has wrong count %d != %d",
						  dx, clumpObs);

	if (doubleCheck) {
		Eigen::MatrixXd denseAF;
		denseAF.resize(fullA.rows(), clumpObs);
		int xx=0;
		for (int cx=0; cx < fullA.cols(); ++cx) {
			if (!latentFilter[cx]) continue;
			denseAF.col(xx) = denseA.col(cx);
			++xx;
		}
		if (xx != clumpObs) mxThrow("latentFilter has wrong count %d != %d",
					     xx, clumpObs);

		// ensure inner iterator works
		for (int k=0; k< IAF.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(IAF, k); it; ++it) {
				if (denseAF.coeff(it.row(), it.col()) != it.value()) {
					mxLog("[%d,%d] %f != %f",
					      it.row(), it.col(), denseAF.coeff(it.row(), it.col()), it.value());
				}
			}
		}

		Eigen::MatrixXd denseFilteredA = IAF;
		if ((denseAF.array() != denseFilteredA.array()).any()) {
			for (int rx=0; rx<denseAF.rows(); ++rx) {
				for (int cx=0; cx<denseAF.cols(); ++cx) {
					if (denseAF.coeff(rx,cx) != denseFilteredA.coeff(rx,cx)) {
						mxLog("[%d,%d] %f != %f",
						      rx, cx, denseAF.coeff(rx,cx), denseFilteredA.coeff(rx,cx));
					}
				}
			}
			mxThrow("stop");
		}
	}
	filtered = true;
	++filterCount;
	//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
}

/*
	// Based on lme4CholmodDecomposition.h from lme4
	template<typename _MatrixType, int _UpLo = Eigen::Lower>
	class Cholmod : public Eigen::CholmodDecomposition<_MatrixType, _UpLo> {
	private:
		Eigen::MatrixXd ident;
		cholmod_dense* x_cd;

	protected:
		typedef Eigen::CholmodDecomposition<_MatrixType, _UpLo> Base;
		using Base::m_factorizationIsOk;
		typedef void (*cholmod_error_type)(int status, const char *file, int line, const char *message);

	        cholmod_factor* factor() const { return Base::m_cholmodFactor; }
		cholmod_common& cholmod() const {
			return const_cast<Cholmod<_MatrixType, _UpLo>*>(this)->Base::cholmod();
		}

		cholmod_error_type oldHandler;
		static void cholmod_error(int status, const char *file, int line, const char *message) {
			// cannot throw exception here because we might be in an OpenMP section
			failed = true;
		}

     // * If you are going to factorize hundreds or more matrices with the same
     // * nonzero pattern, you may wish to spend a great deal of time finding a
     // * good permutation.  In this case, try setting Common->nmethods to CHOLMOD_MAXMETHODS
     // * The time spent in cholmod_analysis will be very high, but you need to
     // * call it only once. TODO

	public:
		static int failed;

		Cholmod() : x_cd(NULL) {
			oldHandler = cholmod().error_handler;
			cholmod().error_handler = cholmod_error;
			cholmod().supernodal = CHOLMOD_AUTO;
			cholmod().nmethods = CHOLMOD_MAXMETHODS;
			if (0) {
				cholmod().nmethods = 2;
				cholmod().method[0].ordering = CHOLMOD_NESDIS;
				cholmod().method[1].ordering = CHOLMOD_AMD;
			}
		};
		~Cholmod() {
			if (x_cd) cholmod_free_dense(&x_cd, &cholmod());
			cholmod().error_handler = oldHandler;
		};

		void analyzePattern(const typename Base::MatrixType& matrix)
		{
			if (ident.rows() != matrix.rows()) {
				ident.setIdentity(matrix.rows(), matrix.rows());
			}
			Base::analyzePattern(matrix);
			cholmod_common &cm = cholmod();
			if (OMX_DEBUG) {
				mxLog("Cholmod: selected ordering %d lnz=%f fl=%f super=%d",
				      cm.method[cm.selected].ordering,
				      cm.method[cm.selected].lnz, cm.method[cm.selected].fl, factor()->is_super);
			}
		}

		double log_determinant() const {
			// Based on https://github.com/njsmith/scikits-sparse/blob/master/scikits/sparse/cholmod.pyx
			cholmod_factor *cf = factor();
			if (cf->xtype == CHOLMOD_PATTERN) mxThrow("Cannot extract diagonal from symbolic factor");
			double logDet = 0;
			double *x = (double*) cf->x;
			if (cf->is_super) {
				// This is a supernodal factorization, which is stored as a bunch
				// of dense, lower-triangular, column-major arrays packed into the
				// x vector. This is not documented in the CHOLMOD user-guide, or
				// anywhere else as far as I can tell; I got the details from
				// CVXOPT's C/cholmod.c.

				int *super = (int*) cf->super;
				int *pi = (int*) cf->pi;
				int *px = (int*) cf->px;
				for (size_t sx=0; sx < cf->nsuper; ++sx) {
					int ncols = super[sx + 1] - super[sx];
					int nrows = pi[sx + 1] - pi[sx];

					Eigen::Map<const Eigen::Array<double,1,Eigen::Dynamic>, 0, Eigen::InnerStride<> >
						s1(x + px[sx], ncols, Eigen::InnerStride<>(nrows+1));
					logDet += s1.real().log().sum();
				}
			} else {
				// This is a simplicial factorization, which is simply stored as a
				// sparse CSC matrix in x, p, i. We want the diagonal, which is
				// just the first entry in each column; p gives the offsets in x to
				// the beginning of each column.
				//
				// The ->p array actually has n+1 entries, but only the first n
				// entries actually point to real columns (the last entry is a
				// sentinel)
				int *p = (int*) cf->p;
				for (size_t ex=0; ex < cf->n; ++ex) {
					logDet += log( x[p[ex]] );
				}
			}
			if (cf->is_ll) {
				logDet *= 2.0;
			}
			return logDet;
		};

		double *getInverseData() const
		{
			return (double*) x_cd->x;
		}

		void refreshInverse()
		{
			eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
			cholmod_dense b_cd(viewAsCholmod(ident));
			if (x_cd) cholmod_free_dense(&x_cd, &cholmod());
			x_cd = cholmod_solve(CHOLMOD_A, factor(), &b_cd, &cholmod());
			if(!x_cd) throw std::runtime_error("cholmod_solve failed"); // impossibe?
		};
	};

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
