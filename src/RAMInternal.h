#ifndef _RAMINTERNAL_H_
#define _RAMINTERNAL_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Cholesky>
//#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
#include <Rcpp.h>
//#include <RcppEigenStubs.h>
#include <RcppEigenWrap.h>
//#include <Eigen/UmfPackSupport>
//#include <RcppEigenCholmod.h>

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

	AsymTool(std::vector<T1> &latentFilter) :
	AshallowDepth(-1), signA(-1.0), latentFilter(latentFilter), analyzed(false), hasDeterminedDepth(false),
		invertCount(0), filterCount(0) {};
	void resize(int clumpVars, int _clumpObs)
	{
		fullA.resize(clumpVars, clumpVars);
		ident.resize(clumpVars, clumpVars);
		ident.setIdentity();
		clumpObs = _clumpObs;
	};
	void determineShallowDepth(FitContext *fc);
	void invert();
	void filter();
};

template <typename T1> void AsymTool<T1>::determineShallowDepth(FitContext *fc)
{
	int maxDepth = std::min(fullA.cols(), 30);
	if (Global->RAMMaxDepth != NA_INTEGER) maxDepth = Global->RAMMaxDepth;
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
			Rf_error("Failed to invert flattened A matrix; %s",
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
	if (OMX_DEBUG && dx != clumpObs) Rf_error("latentFilter has wrong count %d != %d",
						  dx, clumpObs);

	if (doubleCheck) {
		Eigen::MatrixXd denseAF;
		denseAF.resize(fullA.rows(), clumpObs);
		int dx=0;
		for (int cx=0; cx < fullA.cols(); ++cx) {
			if (!latentFilter[cx]) continue;
			denseAF.col(dx) = denseA.col(cx);
			++dx;
		}
		if (dx != clumpObs) Rf_error("latentFilter has wrong count %d != %d",
					     dx, clumpObs);

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
			Rf_error("stop");
		}
	}
	filtered = true;
	++filterCount;
	//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
}

template<typename _MatrixType, int _UpLo = Eigen::Lower>
class SimpCholesky : public Eigen::LDLT<_MatrixType, _UpLo> {
 private:
	Eigen::MatrixXd ident;
	Eigen::MatrixXd inverse;

 public:
	typedef Eigen::LDLT<_MatrixType, _UpLo> Base;

	double log_determinant() const {
		typename Base::Scalar detL = Base::vectorD().array().log().sum();
		return detL;
	}

	void refreshInverse()
	{
		if (ident.rows() != Base::m_matrix.rows()) {
			ident.setIdentity(Base::m_matrix.rows(), Base::m_matrix.rows());
		}

		inverse = Base::solve(ident);
	};

	const Eigen::MatrixXd &getInverse() const { return inverse; };
};

/*
template<typename _MatrixType, int _UpLo = Eigen::Lower>
class SimpCholesky : public Eigen::SimplicialLDLT<_MatrixType, _UpLo> {
 private:
	Eigen::MatrixXd ident;
	Eigen::MatrixXd inverse;

 public:
	typedef Eigen::SimplicialLDLT<_MatrixType, _UpLo> Base;
	
	double log_determinant() const {
		typename Base::Scalar detL = Base::vectorD().array().log().sum();
		return detL;
	}

	void refreshInverse()
	{
		if (ident.rows() != Base::m_matrix.rows()) {
			ident.setIdentity(Base::m_matrix.rows(), Base::m_matrix.rows());
		}

		inverse = Base::solve(ident);
	};

 	const Eigen::MatrixXd &getInverse() const { return inverse; };
};
*/

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
			if (cf->xtype == CHOLMOD_PATTERN) Rf_error("Cannot extract diagonal from symbolic factor");
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
*/

class omxRAMExpectation;

namespace RelationalRAMExpectation {

	// addrSetup and addr are conceptual the same object. They are
	// split up because we need the information in addr on an
	// ongoing basis. addrSetup could be discarded after initial
	// analysis of the data.

	struct addrSetup {
		int numKids;
		int numJoins;
		int parent1;  // first parent
		int fk1;      // first foreign key
		bool rotationLeader;

		// clump indexes into the layout for models that
		// are considered a compound component of this model.
		std::vector<int> clump;
		bool clumped;
		int group;
		int copy;
	};

	class addr {
	private:
		omxExpectation *model;  // read-only
	public:
		int row;                 // to load definition variables (never the key)
		struct independentGroup *ig;
		int igIndex;

		int numVars() const;
		int numObsCache;
		int numObs() const { return numObsCache; }
		double rampartScale;

		std::string modelName() {
			std::string tmp = model->data->name;
			tmp = tmp.substr(0, tmp.size() - 5); // remove ".data" suffix
			return tmp;
		};
		void setModel(omxExpectation *ex) { model=ex; };
		omxExpectation *getModel(FitContext *fc) {
			return omxExpectationFromIndex(model->expNum, fc->state);
		};
		int getExpNum() const { return model->expNum; };
		omxData *getData() const { return model->data; };
		std::vector<bool> &getIgnoreDefVar();
		omxRAMExpectation *getRAMExpectation(FitContext *fc) {
			return (omxRAMExpectation*) getModel(fc)->argStruct;
		};
		omxRAMExpectation *getRAMExpectationReadOnly() const {
			// NOTE: not per-thread!
			return (omxRAMExpectation*) model->argStruct;
		};
		std::vector< omxMatrix* > &getBetween() const;
		omxMatrix *getDataColumns() const { return model->dataColumns; };
		void dataRow(omxMatrix *out) const;
	};

	struct RowToLayoutMapCompare {
		bool operator() (const std::pair<omxData*,int> &lhs, const std::pair<omxData*,int> &rhs) const
		{
			if (lhs.first != rhs.first)
				return strcmp(lhs.first->name, rhs.first->name) < 0;
			return lhs.second < rhs.second;
		}
	};

	struct placement {
		int modelStart;  // both latent and obs
		int obsStart;
	};

	struct sufficientSet {
		int                              start;   // index into placements
		int                              length;  // # of clumpSize units
		Eigen::MatrixXd                  dataCov;
		Eigen::VectorXd                  dataMean;
	};

	class independentGroup {
	private:
		class state &st;

		void refreshModel(FitContext *fc);
		void refreshUnitA(FitContext *fc, int px);
		void invertAndFilterA();
	public:
		int arrayIndex;
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToPlacementMapType;
		RowToPlacementMapType            rowToPlacementMap;
		std::vector<int>                 gMap;  // vector of indexes into layout
		std::vector<placement>           placements;
		const int                        clumpSize;
		int                              clumpVars, clumpObs;
		std::vector<sufficientSet>       sufficientSets;
		SEXP                             obsNameVec;
		SEXP                             varNameVec;
		// make dataColumn optional TODO
		Eigen::ArrayXi                   dataColumn; // for OLS profiled constant parameters
		Eigen::VectorXd                  dataVec;
		Eigen::VectorXd                  fullMean;
		Eigen::VectorXd                  rawFullMean;
		Eigen::VectorXd                  expectedVec;
		Eigen::SparseMatrix<double>      fullCov;
		bool                             analyzedCov;
		//Cholmod< Eigen::SparseMatrix<double> > covDecomp;
		SimpCholesky< Eigen::MatrixXd >  covDecomp;
		Eigen::SparseMatrix<double>      fullS;
		std::vector<bool>                latentFilter; // false when latent or missing

		// could store coeff extraction plan in addr TODO
		AsymTool<bool>          asymT;

		independentGroup(class state *st, int size, int clumpSize)
			: st(*st), clumpSize(clumpSize),
			analyzedCov(false), asymT(latentFilter)
		{ placements.reserve(size); };
		independentGroup(independentGroup *ig);
		int numLooseClumps() {
			independentGroup &par = getParent();
			int loose = par.placements.size() / clumpSize;
			if (par.sufficientSets.size()) {
				loose = par.sufficientSets[0].start / clumpSize;
			}
			return loose;
		};
		void place(int ax);
		void prep(FitContext *fc);
		void determineShallowDepth(FitContext *fc);
		int verbose() const;
		void filterFullMean();
		void finalizeData();
		Eigen::SparseMatrix<double> getInputMatrix() const;
		void computeCov1(FitContext *fc);
		void computeCov2();
		void exportInternalState(MxRList &out, MxRList &dbg);
		independentGroup &getParent();
	};

	class state {
	private:
		state *parent;
		std::vector<int>                 rampartUsage;
		std::vector< std::vector<int> >  rotationPlan;

	public:
		typedef std::vector< std::set<int> > ConnectedType;
		struct omxExpectation *homeEx;
		std::set<struct omxExpectation *> allEx;
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToLayoutMapType;
		RowToLayoutMapType               rowToLayoutMap;
		std::vector<addrSetup>		 layoutSetup;
		std::vector<addr>		 layout;

		void clumpWith(int upper, int lower) {
			if (layoutSetup[lower].clumped) Rf_error("%d is already clumped", lower);
			layoutSetup[upper].clump.push_back(lower);
			layoutSetup[lower].clumped = true;
		};

		omxMatrix                       *smallCol;
		std::vector<independentGroup*>   group;
		bool                             doIdentifyZeroVarPred;

	private:
		int flattenOneRow(omxExpectation *expectation, int frow, int &maxSize);
		template <typename T>
		void placeSet(std::set<std::vector<T> > &toPlace, independentGroup *ig, int groupNum, int &copyNum);
		void planModelEval(int maxSize, FitContext *fc);
		void identifyZeroVarPred(FitContext *fc);
		int rampartRotate(int level);
		template <typename T> void oertzenRotate(std::vector<T> &t1);
		template <typename T> void applyRotationPlan(T accessor);
		template <typename T> void appendClump(int ax, std::vector<T> &clump);
		template <typename T> void propagateDefVar(omxRAMExpectation *ram,
							   Eigen::MatrixBase<T> &transition,
							   omxRAMExpectation *ram2, bool within);
		void computeConnected(std::vector<int> &region, ConnectedType &connected);
	public:
		~state();
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		int verbose() const;
		bool hasRotationPlan() const { return rotationPlan.size() != 0; }
		void exportInternalState(MxRList &dbg);
		state &getParent() { return *parent; };
	};
};

typedef std::set< std::pair< omxExpectation*, int> > dvScoreboardSetType;

class omxRAMExpectation {
	bool trivialF;
	int Zversion;
	omxMatrix *_Z;
 public:
	std::vector< dvScoreboardSetType > dvScoreboard;
	Eigen::VectorXd hasVariance;
	std::vector<bool> ignoreDefVar;
	std::vector<bool> latentFilter; // false when latent

	// composition of F permutation and expectation->dataColumns
	Eigen::VectorXi dataCols;

 	omxRAMExpectation(omxMatrix *Z) : trivialF(false), Zversion(0), _Z(Z) {};
	~omxRAMExpectation() {
		omxFreeMatrix(_Z);
	};

	omxMatrix *getZ(FitContext *fc);
	void CalculateRAMCovarianceAndMeans();

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *X, *Y, *Ax;

	int verbose;
	int numIters;
	int rampart;
	bool rampartEnabled() { return rampart == NA_INTEGER || rampart > 0; };
	double logDetObserved;
	double n;
	double *work;
	int lwork;

	std::vector< omxMatrix* > between;
	RelationalRAMExpectation::state *rram;
	bool forceSingleGroup;

	void studyF(omxMatrix *dc);
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx->argStruct)->verbose;
	}

	inline int independentGroup::verbose() const { return st.verbose(); };

};

#endif
