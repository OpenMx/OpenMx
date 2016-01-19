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
		int chain;    // other addr to consider part of the same unit
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

	class state {
	private:
		struct omxExpectation *homeEx;
		omxMatrix *smallCol;
		bool AmatDependsOnParameters;
		bool haveFilteredAmat;
		Eigen::SparseMatrix<double>      depthTestA;
		int AshallowDepth;
		bool analyzedFullA;
		double signA;
		Eigen::SparseMatrix<double>      fullA;
		Eigen::SparseLU< Eigen::SparseMatrix<double>,
				 Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > Asolver;
		//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;
		Eigen::SparseMatrix<double>      ident;
		Eigen::SparseMatrix<double>      fullS;

	public:
		std::vector<addr>		 layout;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		SEXP obsNameVec;
		SEXP varNameVec;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMeans;
		Eigen::SparseMatrix<double>      filteredA;
		Eigen::SparseMatrix<double>      fullCov;

	private:
		void refreshModel(FitContext *fc);
		int placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void prepOneRow(omxExpectation *expectation, int row_or_key, int &lx, int &dx);
		int rampartRotate();
		template <typename T> void oertzenRotateCompound(std::vector<T> &t1);
	public:
		void compute(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		~state();
		void exportInternalState(MxRList &dbg);
		int verbose() const;
	};
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
