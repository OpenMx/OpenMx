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
		int numKids;
		int key;
		int numJoins;
		int parent1;  // first parent
		int fk1;      // first foreign key

		// clump names indexes into the layout for models that
		// are considered a compound component of this model.
		std::vector<int> clump;

		int modelStart, modelEnd;  //both latent and obs
		int numVars() const { return 1 + modelEnd - modelStart; }
		int obsStart, obsEnd;
		int numObsCache;
		int numObs() const { return numObsCache; }
		double rampartScale;

		std::string modelName() const {
			std::string tmp = model->data->name;
			tmp = tmp.substr(0, tmp.size() - 5); // remove ".data" suffix
			return tmp;
		};
		static bool CompareWithModelStart(addr &i, int p1) { return i.modelStart < p1; };
	};

	struct Amatrix {
		bool analyzed;
		Eigen::SparseMatrix<double>      in;
		Eigen::SparseLU< Eigen::SparseMatrix<double>,
				 Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > solver;
		Eigen::SparseMatrix<double>      out;
		//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;

		Amatrix() : analyzed(false) {};
	};

	struct RowToOffsetMapCompare {
		bool operator() (const std::pair<omxData*,int> &lhs, const std::pair<omxData*,int> &rhs) const
		{
			if (lhs.first != rhs.first)
				return strcmp(lhs.first->name, rhs.first->name) < 0;
			return lhs.second < rhs.second;
		}
	};

	struct placement {
		int aIndex;                // index into addr vector
		int modelStart;  //both latent and obs
		int obsStart;
	};

	struct independentGroup {
		int numRows;
		std::vector<placement> placements;
		omxMatrix *smallCol;
		bool AmatDependsOnParameters;
		bool haveFilteredAmat;
		Amatrix testA;
		int AshallowDepth;
		double signA;
		Eigen::SparseMatrix<double>      ident;
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
		struct omxExpectation *homeEx;
		omxMatrix *smallCol;
		bool AmatDependsOnParameters;
		bool haveFilteredAmat;
		Amatrix testA;
		int AshallowDepth;
		double signA;
		Eigen::SparseMatrix<double>      ident;
		Eigen::SparseMatrix<double>      fullS;
		std::vector<int>                 rampartUsage;
		std::vector< std::vector<int> >  rotationPlan;
		typedef std::map< std::pair<omxData*,int>, int, RowToOffsetMapCompare> RowToOffsetMapType;
		RowToOffsetMapType               rowToOffsetMap;

	public:
		std::vector<addr>		 layout;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		SEXP obsNameVec;
		SEXP varNameVec;
		Amatrix regularA;
		Eigen::ArrayXi dataColumn;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMean;
		Eigen::VectorXd expectedMean;
		Eigen::SparseMatrix<double>      fullCov;

	private:
		void refreshLevelTransitions(FitContext *fc, addr &a1, Amatrix &dest, double scale);
		void refreshUnitA(FitContext *fc, addr &a1, Amatrix &dest);
		void invertAndFilterA(Amatrix &Amat);
		void refreshModel(FitContext *fc);
		int placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void examineModel();
		int rampartRotate(int level);
		template <typename T> void oertzenRotate(std::vector<T> &t1);
	public:
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		~state();
		void exportInternalState(MxRList &dbg);
		int verbose() const;
		template <typename T> void applyRotationPlan(Eigen::MatrixBase<T> &resid) const;
		bool hasRotationPlan() const { return rotationPlan.size() != 0; }
	};

	template <typename T>
	void state::applyRotationPlan(Eigen::MatrixBase<T> &resid) const
	{
		//std::string buf;
		for (size_t rx=0; rx < rotationPlan.size(); ++rx) {
			//buf += "rotate";
			const std::vector<int> &om = rotationPlan[rx];
			double partialSum = 0.0;
			for (size_t ox=0; ox < om.size(); ++ox) {
				partialSum += resid[om[ox]];
				//buf += string_snprintf(" %d", 1+ om[ox]);
			}
			double prev = resid[om[0]];
			resid[om[0]] = partialSum / sqrt(om.size());

			for (size_t i=1; i < om.size(); i++) {
				double k=om.size()-i;
				partialSum -= prev;
				double prevContrib = sqrt(k / (k+1)) * prev;
				prev = resid[om[i]];
				resid[om[i]] = partialSum * sqrt(1.0 / (k*(k+1))) - prevContrib;
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
};

#endif
