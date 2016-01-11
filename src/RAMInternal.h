#ifndef _RAMINTERNAL_H_
#define _RAMINTERNAL_H_

#include <RcppEigenCholmod.h>
#include <RcppEigenStubs.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
//#include <Eigen/UmfPackSupport>

namespace RelationalRAMExpectation {
	class state {
	private:
		struct omxExpectation *homeEx;
		int verbose;
		omxMatrix *smallCol;
		std::vector<bool> latentFilter; // use to reduce the A matrix
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
		std::vector<const char *> nameVec;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMeans;
		Eigen::SparseMatrix<double>      filteredA;
		Eigen::SparseMatrix<double>      fullCov;

	private:
		void loadOneRow(omxExpectation *expectation, FitContext *fc, int row, int &lx);
		void placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void prepOneRow(omxExpectation *expectation, int row_or_key, int &lx, int &dx);
	public:
		void compute(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		~state();
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
