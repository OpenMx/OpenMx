#ifndef u_RAMINTERNAL_H_
#define u_RAMINTERNAL_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Cholesky>
//#include <Eigen/SparseCholesky>
//#include <Eigen/UmfPackSupport>
#include <RcppEigenCholmod.h>
#include <RcppEigenWrap.h>
#include "path.h"
#include "Connectedness.h"

struct coeffLoc {
	int off;
	int r, c;

coeffLoc(int u_off, int u_r, int u_c) :
	off(u_off), r(u_r), c(u_c) {}
};

class omxRAMExpectation;

namespace RelationalRAMExpectation {

	// addrSetup and addr are conceptual the same object. They are
	// split up because we need the information in addr on an
	// ongoing basis. addrSetup could be discarded after initial
	// analysis of the data.

	struct addrSetup {
		int numKids;  // how many lower level units (kids) join with this unit?
		int numJoins; // how many parents does this unit join with?
		int parent1;  // first parent
		int fk1;      // first foreign key

		// clump indexes into the layout for models that
		// are considered a compound component of this model.
		std::vector<int> clump;
		bool clumped;
		int rset; // "rotation set" annotation for debugging
		int skipMean;
		bool heterogenousMean;
	};

	class addr {
	private:
		omxExpectation *model;  // read-only
	public:
		int row;                 // to load definition variables (never the key)
		class independentGroup *ig;
		int igIndex;
		int nextMean;

		int numVars() const;
		int numObsCache;
		int numObs() const { return numObsCache; }
		double rampartScale;
		double quickRotationFactor;

		std::string modelName() {
			std::string tmp = model->data->name;
			tmp = tmp.substr(0, tmp.size() - 5); // remove ".data" suffix
			return tmp;
		};
		void setModel(omxExpectation *ex) { model=ex; };
		omxExpectation *getModel(FitContext *fc);
		int getExpNum() const { return model->expNum; };
		omxData *getData() const { return model->data; };
		const std::vector<bool> &getDefVarInfluenceMean() const;
		const std::vector<bool> &getDefVarInfluenceVar() const;
		omxRAMExpectation *getRAMExpectation(FitContext *fc);
		omxRAMExpectation *getRAMExpectationReadOnly() const {
			// NOTE: not per-thread!
			return (omxRAMExpectation*) model;
		};
		std::vector< omxMatrix* > &getBetween() const;
		const Eigen::Map<DataColumnIndexVector> getDataColumns() const {
			return model->getDataColumns();
		};
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
		// both start & length are in multiples of clumpSize
		int                              start;
		int                              length;
		Eigen::MatrixXd                  dataCov;
		Eigen::VectorXd                  dataMean;
	};

	// This is not really the best organization. It would be better to
	// partition the work into equal size covariance matrices then by
	// identical covariance matrices and then by identical means.

	class independentGroup {
	private:
		class state &st;

		void refreshModel(FitContext *fc);
		void refreshUnitA(FitContext *fc, int px);
		void invertAndFilterA();

		struct MpcIO : PathCalcIO {
			independentGroup &par;
			int clumpSize;
			MpcIO(independentGroup &u_par) : par(u_par), clumpSize(u_par.clumpSize) {}
			virtual void recompute(FitContext *fc) override;
			virtual unsigned getVersion(FitContext *fc) override;
			virtual void refresh(FitContext *fc) override;
			virtual PathCalcIO *clone() override
			{ return new MpcIO(par); }
		};

		struct ApcIO : PathCalcIO {
			independentGroup &par;
			int clumpSize;
			bool useRampart;
			ApcIO(independentGroup &u_par) : par(u_par), clumpSize(u_par.clumpSize), useRampart(true) {}
			virtual void recompute(FitContext *fc) override;
			virtual unsigned getVersion(FitContext *fc) override;
			template <typename T>
			void u_refresh(FitContext *fc, T &mat, double sign);
			virtual void refreshA(FitContext *fc, double sign) override
			{ u_refresh(fc, full, sign); }
			virtual void refreshSparse1(FitContext *fc, double sign) override
			{ u_refresh(fc, sparse, sign); }
			virtual PathCalcIO *clone() override
			{ return new ApcIO(par); }
		};

		struct SpcIO : PathCalcIO {
			independentGroup &par;
			int clumpSize;
			SpcIO(independentGroup &u_par) : par(u_par), clumpSize(u_par.clumpSize) {}
			virtual void recompute(FitContext *fc) override;
			virtual unsigned getVersion(FitContext *fc) override;
			template <typename T>
			void u_refresh(FitContext *fc, T &mat);
			virtual void refresh(FitContext *fc) override
			{ u_refresh(fc, full); }
			virtual void refreshSparse1(FitContext *fc, double sign) override
			{ u_refresh(fc, sparse); }
			virtual PathCalcIO *clone() override
			{ return new SpcIO(par); }
		};

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
		Eigen::VectorXd                  simDataVec;
		Eigen::VectorXd                  fullMean;  // rename, latents are filtered out
		Eigen::VectorXd                  rawFullMean;
		int                              skipMean;
		Eigen::VectorXd                  expectedVec;
		Eigen::MatrixXd                  fullCov;   // rename, latents are filtered out
		std::vector<bool>                latentFilter; // false when latent or missing
    std::vector<bool>                isProductNode;

		PathCalc pcalc;
		double                           fit;  // most recent fit for debugging

		independentGroup(class state *u_st, int size, int u_clumpSize);
		independentGroup(independentGroup *ig);
		int numLooseClumps() {
			independentGroup &par = getParent();
			int loose = par.placements.size() / clumpSize;
			if (par.sufficientSets.size()) {
				loose = par.sufficientSets[0].start;
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
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void simulate();
		void exportInternalState(MxRList &out, MxRList &dbg);
		independentGroup &getParent();
	};

	using namespace UndirectedGraph;

	class state {
	private:
		state *parent;
		std::vector<int>                 rampartUsage;
		std::vector< std::vector<int> >  rotationPlan;
		std::vector< bool >              modelRotationPlanFilter;
		int rotationCount;

	public:
		bool isChild() const { return this != parent; }
		typedef std::vector< std::set<int> > SubgraphType;
		class omxExpectation *homeEx;
		std::set<class omxExpectation *> allEx;
		bool hasProductNodes;
		typedef std::map< std::pair<omxData*,int>, int, RowToLayoutMapCompare> RowToLayoutMapType;
		RowToLayoutMapType               rowToLayoutMap;
		std::vector<addrSetup>		 layoutSetup;
		std::vector<addr>		 layout;

		void clumpWith(int upper, int lower) {
			if (layoutSetup[lower].clumped) mxThrow("%d is already clumped", lower);
			layoutSetup[upper].clump.push_back(lower);
			layoutSetup[lower].clumped = true;
		};

		omxMatrix                       *smallCol;
		std::vector<independentGroup*>   group;
		bool                             doAnalyzeDefVars;

	private:
		int flattenOneRow(omxExpectation *expectation, int frow, int &maxSize);
		template <typename T>
		bool placeSet(std::set<std::vector<T> > &toPlace, independentGroup *ig);
		void planModelEval(int maxSize, FitContext *fc);
		void analyzeModel1(FitContext *fc);
		void analyzeModel2(FitContext *fc);
		int rampartRotate(int level);
		template <typename T> void oertzenRotate(std::vector<T> &t1, bool canOptimize);
		template <typename T> void unapplyRotationPlan(T accessor);
		template <typename T> void applyRotationPlan(T accessor);
		template <typename T> void appendClump(int ax, std::vector<T> &clump);
		template <typename T>
		void propagateDefVar(omxRAMExpectation *to, Eigen::MatrixBase<T> &transition,
												 omxRAMExpectation *from);
		void computeConnected(std::vector<int> &region, SubgraphType &connected);
		void computeMeanByModel(FitContext *fc);
		void computeMeanByGroup(FitContext *fc);
	public:
		~state();
		void computeCov(FitContext *fc);
		void computeMean(FitContext *fc);
		void init(omxExpectation *expectation, FitContext *fc);
		int verbose() const;
		bool hasRotationPlan() const { return rotationPlan.size() != 0; }
		void exportInternalState(MxRList &dbg);
		state &getParent() { return *parent; };
		void simulate(FitContext *fc, MxRList &out);
		int getOptimizeMean();
		void optimizeModelRotation();
    int numObservedStats();
	};
};

class omxRAMExpectation : public MVNExpectation {
	typedef MVNExpectation super;
	Eigen::VectorXi dataCols;  // composition of F permutation and expectation->dataColumns
	std::vector<const char *> dataColNames;
	std::vector< omxThresholdColumn > thresholds;
	std::vector<int> exoDataColumns; // index into omxData
	bool hasProductNodes;
  bool studiedF;
  bool openBox;  // can the user access the expectation during optimization?
	int numExoPred;
	std::vector<int> exoDataColIndex;
  void addSlopeMatrix();

	struct MpcIO : PathCalcIO {
		omxMatrix *M0;
		MpcIO() {}
		virtual void recompute(FitContext *fc) override;
		virtual unsigned getVersion(FitContext *fc) override;
		virtual void refresh(FitContext *fc) override;
		virtual PathCalcIO *clone() override
		{
			auto *mio = new MpcIO;
			mio->M0 = M0;
			return mio;
		}
	};

	struct ApcIO : PathCalcIO {
		omxMatrix *A0;
		std::vector<coeffLoc> &vec;
		ApcIO(std::vector<coeffLoc> &u_vec) : vec(u_vec) {}
		virtual void recompute(FitContext *fc) override;
		virtual unsigned getVersion(FitContext *fc) override;
		template <typename T>
		void u_refresh(FitContext *fc, T &mat, double sign);
		virtual void refreshA(FitContext *fc,double sign) override
		{ u_refresh(fc, full, sign); }
		virtual void refreshSparse1(FitContext *fc, double sign) override
		{ u_refresh(fc, sparse, sign); }
		virtual PathCalcIO *clone() override
		{
			auto *aio = new ApcIO(vec);
			aio->A0 = A0;
			return aio;
		}
	};

	struct SpcIO : PathCalcIO {
		omxMatrix *S0;
		std::vector<coeffLoc> &vec;
		SpcIO(std::vector<coeffLoc> &u_vec) : vec(u_vec) {}
		virtual void recompute(FitContext *fc) override;
		virtual unsigned getVersion(FitContext *fc) override;
		template <typename T>
		void u_refresh(FitContext *fc, T &mat);
		virtual void refresh(FitContext *fc) override
		{ u_refresh(fc, full); }
		virtual void refreshSparse1(FitContext *fc, double sign) override
		{ u_refresh(fc, sparse); }
		virtual PathCalcIO *clone() override
		{
			auto *sio = new SpcIO(vec);
			sio->S0 = S0;
			return sio;
		}
	};

 public:
	typedef std::pair< omxExpectation*, int> dvRefType; // int is offset into data->defVars array
	typedef std::set< dvRefType > dvRefSetType;
	std::vector< dvRefSetType > dvContribution; // per variable dv contribution
	Eigen::VectorXd hasVariance;    // current level, if nonzero then has variance
	Eigen::VectorXd hasMean;        // current level, if nonzero then has means
	// per variable influence of defVars, including lower levels
	std::vector<bool> dvInfluenceMean;
	std::vector<bool> dvInfluenceVar;
	std::vector<bool> latentFilter; // false when latent
	std::vector<bool> isProductNode;
	bool getHasProductNodes() const { return hasProductNodes; }
	std::vector<coeffLoc> ScoeffStorage;
	std::vector<coeffLoc> AcoeffStorage;
	std::vector<coeffLoc> *Scoeff;
	std::vector<coeffLoc> *Acoeff;
	PathCalc pcalc;

	omxRAMExpectation(omxState *st, int num);
	virtual ~omxRAMExpectation();

	void CalculateRAMCovarianceAndMeans(FitContext *fc);
	void analyzeDefVars(FitContext *fc);
	void logDefVarsInfluence();
  bool isOpenBox() const { return openBox; }

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrixPtr covOwner, meanOwner;
	omxMatrix *fullCov, *fullMean;
	omxMatrix *slope;       // exogenous predictor slopes
	omxMatrix *A, *S, *F, *M;
	
	std::vector< Eigen::SparseMatrix<double> > dS_dtheta;
	std::vector< Eigen::SparseMatrix<double> > dA_dtheta;
	std::vector< Eigen::SparseMatrix<double> > dM_dtheta;
	std::vector<bool> alwaysZeroSDeriv;
	std::vector<bool> alwaysZeroADeriv;
	std::vector<bool> alwaysZeroMDeriv;
	
	void provideSufficientDerivs(
			FitContext *fc, std::vector< Eigen::MatrixXd > &u_dSigma_dtheta, std::vector< Eigen::MatrixXd > &u_dNu_dtheta,
			std::vector<bool> &u_alwaysZeroCovDeriv, std::vector<bool> &u_alwaysZeroMeanDeriv, bool wantHess, 
			std::vector< std::vector< Eigen::MatrixXd >> &u_d2Sigma_dtheta1dtheta2, 
			std::vector< std::vector< Eigen::MatrixXd >> &u_d2Mu_dtheta1dtheta2) override;

	int verbose;
	int rampartCycleLimit;
	int rampartUnitLimit;
	int maxDebugGroups;
	bool useSufficientSets;
	int optimizeMean;
	int useSparse;
	bool rampartEnabled() { return (rampartCycleLimit == NA_INTEGER || rampartCycleLimit > 0) && !forceSingleGroup; };
	double logDetObserved;
	double n;
	double *work;
	int lwork;

	std::vector< omxMatrix* > between;
	RelationalRAMExpectation::state *rram; // should use unique_ptr TODO
	bool forceSingleGroup;

	void studyF();
	void studyExoPred();

	virtual void init() override;
	virtual void compute(FitContext *fc, const char *what, const char *how) override;
	virtual omxMatrix *getComponent(const char*) override;
	virtual void populateAttr(SEXP expectation) override;
	virtual const std::vector<const char *> &getDataColumnNames() const override {
    if (studiedF) return dataColNames;
    else return super::getDataColumnNames();
  };
	virtual const Eigen::Map<DataColumnIndexVector> getDataColumns() override {
    if (studiedF)
      return Eigen::Map<DataColumnIndexVector>(dataCols.data(), numDataColumns);
    else
      return super::getDataColumns();
	}
	virtual std::vector< omxThresholdColumn > &getThresholdInfo() override
  {
    if (studiedF) return thresholds;
    else return super::getThresholdInfo();
  }
	virtual void invalidateCache() override;
	virtual void generateData(FitContext *fc, MxRList &out) override;
	void flatten(FitContext *fc);
	virtual void getExogenousPredictors(std::vector<int> &out) override;
  virtual int numLatentVars() const override;
  virtual int numObservedStats() override;
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx)->verbose;
	}

	inline int independentGroup::verbose() const { return st.verbose(); };

};

#endif
