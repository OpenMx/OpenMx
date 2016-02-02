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

namespace RelationalRAMExpectation {
	struct addr {
		omxExpectation *model;
		int row;     // to load definition variables (never the key)
		int key;
		int numKids;
		int numJoins;
		int parent1;  // first parent
		int fk1;      // first foreign key

		// clump indexes into the layout for models that
		// are considered a compound component of this model.
		std::vector<int> clump;

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

	class Amatrix {
	private:
		class state &st;
		bool analyzed;
		int AshallowDepth;
		double signA;
		Eigen::SparseMatrix<double>      ident;
	public:
		// could store coeff extraction plan in addr TODO
		Eigen::SparseMatrix<double>      in;
		Eigen::SparseLU< Eigen::SparseMatrix<double>,
				 Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > solver;
		Eigen::SparseMatrix<double>      out;
		//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;

		Amatrix(class state &st) : st(st), analyzed(false), AshallowDepth(-1), signA(-1) {};
		void resize(int dim);
		void determineShallowDepth(FitContext *fc);
		void refreshUnitA(FitContext *fc, struct placement &pl);
		int verbose() const;
		void invertAndFilterA();
		Eigen::SparseMatrix<double> getInputMatrix() const;
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

	struct independentGroup {
		int numRows;
		std::vector<placement> placements;
		omxMatrix *smallCol;
		double signA;
		Eigen::SparseMatrix<double>      fullS;
		std::vector< std::vector<int> >  rotationPlan;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		SEXP obsNameVec;
		SEXP varNameVec;
		Amatrix regularA;
		Amatrix rampartA;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMeans;
		Eigen::SparseMatrix<double>      fullCov;
	};

	class state {
	private:
		Eigen::SparseMatrix<double>      fullS;
		std::vector<int>                 rampartUsage;
		std::vector< std::vector<int> >  rotationPlan;

	public:
		struct omxExpectation *homeEx;
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToLayoutMapType;
		RowToLayoutMapType               rowToLayoutMap;
		std::vector<addr>		 layout;
		omxMatrix *smallCol;
		// below move to independentGroup TODO
		std::vector<placement>           placements;
		std::vector<bool>                latentFilter; // use to reduce the A matrix
		SEXP                             obsNameVec;
		SEXP                             varNameVec;
		Amatrix                          regularA;
		Eigen::ArrayXi                   dataColumn; // for OLS profiled constant parameters
		Eigen::VectorXd                  dataVec;
		Eigen::VectorXd                  fullMean;
		Eigen::VectorXd                  expectedMean;
		Eigen::SparseMatrix<double>      fullCov;

	private:
		void refreshModel(FitContext *fc);
		int flattenOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void planModelEval(int maxSize, int totalObserved);
		int rampartRotate(int level);
		template <typename T> void oertzenRotate(std::vector<T> &t1);
	public:
		state() : regularA(*this) {};
		~state();
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		void exportInternalState(MxRList &dbg);
		int verbose() const;
		template <typename T> void applyRotationPlan(std::vector<placement> &place, Eigen::MatrixBase<T> &resid) const;
		bool hasRotationPlan() const { return rotationPlan.size() != 0; }
	};

	template <typename T>
	void state::applyRotationPlan(std::vector<placement> &place, Eigen::MatrixBase<T> &resid) const
	{
		// maybe faster to do all observations in parallel
		// to allow more possibility of instruction reordering TODO
		//std::string buf;
		for (size_t rx=0; rx < rotationPlan.size(); ++rx) {
			//buf += "rotate";
			const std::vector<int> &units = rotationPlan[rx];

			const addr &specimen = layout[ place[ units[0] ].aIndex ];
			for (int ox=0; ox < specimen.numObs(); ++ox) {
				double partialSum = 0.0;
				for (size_t ux=0; ux < units.size(); ++ux) {
					partialSum += resid[ place[units[ux]].obsStart + ox ];
					//buf += string_snprintf(" %d", 1+ units[ux]);
				}
				double prev = resid[ place[units[0]].obsStart + ox ];
				resid[ place[units[0]].obsStart + ox ] = partialSum / sqrt(units.size());

				for (size_t i=1; i < units.size(); i++) {
					double k=units.size()-i;
					partialSum -= prev;
					double prevContrib = sqrt(k / (k+1)) * prev;
					prev = resid[ place[units[i]].obsStart + ox ];
					resid[ place[units[i]].obsStart + ox ] =
						partialSum * sqrt(1.0 / (k*(k+1))) - prevContrib;
				}
			}
			//buf += "\n";
		}
		//if (buf.size()) mxLogBig(buf);
	}
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

	void ensureTrivialF();
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx->argStruct)->verbose;
	}

	inline int Amatrix::verbose() const { return st.verbose(); };
};

#endif
