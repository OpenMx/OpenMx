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
		int numObs() const { return 1 + obsEnd - obsStart; }
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

	public:
		std::vector<addr>		 layout;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		SEXP obsNameVec;
		SEXP varNameVec;
		Amatrix regularA;
		Amatrix rampartA;
		Eigen::VectorXd dataVec;
		Eigen::VectorXd fullMeans;
		Eigen::VectorXd resid;
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
		void compute(FitContext *fc);
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
 public:

	omxRAMExpectation() {};

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
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx->argStruct)->verbose;
	}
};

#endif
