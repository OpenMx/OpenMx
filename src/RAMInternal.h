#ifndef _RAMINTERNAL_H_
#define _RAMINTERNAL_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/CholmodSupport>
#include <Rcpp.h>
#include <RcppEigenStubs.h>
#include <RcppEigenWrap.h>
//#include <Eigen/UmfPackSupport>
#include <RcppEigenCholmod.h>

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
			throw std::runtime_error(message);
		}

     // * If you are going to factorize hundreds or more matrices with the same
     // * nonzero pattern, you may wish to spend a great deal of time finding a
     // * good permutation.  In this case, try setting Common->nmethods to CHOLMOD_MAXMETHODS
     // * The time spent in cholmod_analysis will be very high, but you need to
     // * call it only once. TODO

	public:
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
			m_factorizationIsOk = false; // base class should take care of it TODO
		};
		~Cholmod() {
			if (x_cd) cholmod_free_dense(&x_cd, &cholmod());
			cholmod().error_handler = oldHandler;
		};

		bool analyzedPattern() const { return m_factorizationIsOk; };

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
			if(!x_cd) throw std::runtime_error("cholmod_solve failed");
		};
	};

namespace RelationalRAMExpectation {
	struct addr {
		omxExpectation *model;
		int row;      // to load definition variables (never the key)
		int key;
		int numKids;
		int numJoins;
		int parent1;  // first parent
		int fk1;      // first foreign key

		// clump indexes into the layout for models that
		// are considered a compound component of this model.
		std::vector<int> clump;
		bool clumped;
		int region;
		int group;
		int copy;
		struct independentGroup *ig;
		int igIndex;

		int numVars() const;
		int numObsCache;
		int numObs() const { return numObsCache; }
		double rampartScale;

		std::string modelName() const {
			std::string tmp = model->data->name;
			tmp = tmp.substr(0, tmp.size() - 5); // remove ".data" suffix
			return tmp;
		};
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
		int aIndex;      // index into addr vector
		int modelStart;  // both latent and obs
		int obsStart;
	};

	class independentGroup {
	private:
		class state &st;
		bool analyzed;
		int AshallowDepth;
		double signA;
		Eigen::SparseMatrix<double>      ident;

		void refreshModel(FitContext *fc);
		void refreshUnitA(FitContext *fc, int px);
		void invertAndFilterA();
	public:
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToPlacementMapType;
		RowToPlacementMapType            rowToPlacementMap;
		std::vector<placement>           placements;
		const int                        clumpSize;
		int                              clumpVars, clumpObs;
		SEXP                             obsNameVec;
		SEXP                             varNameVec;
		// make dataColumn optional TODO
		Eigen::ArrayXi                   dataColumn; // for OLS profiled constant parameters
		Eigen::VectorXd                  dataVec;
		Eigen::VectorXd                  fullMean;
		Eigen::VectorXd                  expectedMean;
		Eigen::SparseMatrix<double>      fullCov;
		Cholmod< Eigen::SparseMatrix<double> > covDecomp;
		Eigen::SparseMatrix<double>      fullS;
		std::vector<bool>                latentFilter; // use to reduce the A matrix

		// could store coeff extraction plan in addr TODO
		Eigen::SparseMatrix<double>      fullA;
		Eigen::SparseLU< Eigen::SparseMatrix<double>,
				 Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > Asolver;
		Eigen::SparseMatrix<double>      IAF;
		//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;

		independentGroup(class state *st, int size, int clumpSize)
			: st(*st), analyzed(false), AshallowDepth(-1), signA(-1), clumpSize(clumpSize)
		{ placements.reserve(size); };
		void prep(int maxSize, int totalObserved, FitContext *fc);
		void determineShallowDepth(FitContext *fc);
		int verbose() const;
		void filterFullMean();
		Eigen::SparseMatrix<double> getInputMatrix() const;
		void computeCov(FitContext *fc);
		void exportInternalState(MxRList &out, MxRList &dbg);
	};

	class state {
	private:
		std::vector<int>                 rampartUsage;
		std::vector< std::vector<int> >  rotationPlan;
		Eigen::VectorXd                  expectedMean;  //debug

	public:
		struct omxExpectation *homeEx;
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToLayoutMapType;
		RowToLayoutMapType               rowToLayoutMap;
		std::vector<addr>		 layout;
		omxMatrix                       *smallCol;
		std::vector<independentGroup*>   group;

	private:
		int flattenOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void planModelEval(int maxSize, int totalObserved, FitContext *fc);
		int rampartRotate(int level);
		template <typename T> void oertzenRotate(std::vector<T> &t1);
		template <typename T> void applyRotationPlan(T accessor);
		template <typename T> void appendClump(int ax, std::vector<T> &clump);
	public:
		~state();
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		int verbose() const;
		bool hasRotationPlan() const { return rotationPlan.size() != 0; }
		void exportInternalState(MxRList &dbg);
	};
};

class omxRAMExpectation {
	bool trivialF;
	int Zversion;
	omxMatrix *_Z;
 public:

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

	void ensureTrivialF();
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx->argStruct)->verbose;
	}

	inline int independentGroup::verbose() const { return st.verbose(); };
};

#endif
