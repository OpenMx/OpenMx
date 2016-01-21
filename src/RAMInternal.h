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
		int numKids; // remove? TODO
		int key;
		int numJoins;
		int fk1;
		int modelStart, modelEnd;  //both latent and obs
		int numVars() const { return modelEnd - modelStart; }
		int obsStart, obsEnd;
		int numObs() const { return obsEnd - obsStart; }

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
		int verbose;
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
		std::vector<addr>		 layout;

	public:
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
	public:
		void compute(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		~state();
		void exportInternalState(MxRList &dbg);
	};
};

struct omxRAMExpectation {

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *X, *Y, *Z, *Ax;

	int numIters;
	double logDetObserved;
	double n;
	double *work;
	int lwork;

	std::vector< omxMatrix* > between;
	RelationalRAMExpectation::state *rram;

	void ensureTrivialF();
};

#endif
