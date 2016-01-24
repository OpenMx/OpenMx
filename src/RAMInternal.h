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
		bool HEV;
		int row;     // to load definition variables (never the key)
		int numKids;
		int key;
		int numJoins;
		int parent1;  // first parent
		int fk1;      // first foreign key
		std::vector<int> clump;
		int modelStart, modelEnd;  //both latent and obs
		int numVars() const { return 1 + modelEnd - modelStart; }
		int obsStart, obsEnd;
		int numObs() const { return 1 + obsEnd - obsStart; }
		bool rampartUnlinked;
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
		size_t                           maxRotationUnits;

	public:
		std::vector<addr>		 layout;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		SEXP obsNameVec;
		SEXP varNameVec;
		Amatrix regularA;
		Amatrix rampartA;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMeans;
		Eigen::SparseMatrix<double>      fullCov;

	private:
		void refreshLevelTransitions(FitContext *fc, addr &a1, Amatrix &dest, double scale);
		void refreshUnitA(FitContext *fc, addr &a1, Amatrix &dest);
		void invertAndFilterA(Amatrix &Amat);
		void refreshModel(FitContext *fc);
		int placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void examineModel();
		int rampartRotate();
		template <typename T> void oertzenRotateCompound(std::vector<T> &t1);
	public:
		void compute(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		~state();
		void exportInternalState(MxRList &dbg);
		int verbose() const;
		template <typename T> void applyRotationPlan(Eigen::MatrixBase<T> &resid) const;
	};

	template <typename T1, typename T2>
	void oertzenRotate1(const int anzGroups, Eigen::MatrixBase<T1> &vec, Eigen::MatrixBase<T2> &erg)
	{
		double partialSum = vec.sum();
		erg[0] = partialSum / sqrt(anzGroups);

		for (int i=1; i<anzGroups; i++) {
			double k=anzGroups-i;
			partialSum -= vec[i-1];
			erg[i] = partialSum * sqrt(1.0 / (k*(k+1))) - sqrt(k / (k+1)) * vec[i-1];
		}
	}

	template <typename T>
	void state::applyRotationPlan(Eigen::MatrixBase<T> &resid) const
	{
		Eigen::VectorXd obsIn(maxRotationUnits);
		Eigen::VectorXd obsOut(maxRotationUnits);
		for (size_t rx=0; rx < rotationPlan.size(); ++rx) {
			const std::vector<int> &r1 = rotationPlan[rx];
			for (size_t ox=0; ox < r1.size(); ++ox) {
				obsIn[ox] = resid[r1[ox]];
			}
			oertzenRotate1(r1.size(), obsIn, obsOut);
			for (size_t ox=0; ox < r1.size(); ++ox) {
				resid[r1[ox]] = obsOut[ox];
			}
		}
	}
};

class omxRAMExpectation {
	bool determinedHEV;
	bool HEV;    // homogeneous error variance
 public:

	omxRAMExpectation() : determinedHEV(false), HEV(false) {};

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *X, *Y, *Z, *Ax;

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
	bool isHEV();
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx->argStruct)->verbose;
	}
};

#endif
