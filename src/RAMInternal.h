#ifndef _RAMINTERNAL_H_
#define _RAMINTERNAL_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Cholesky>
//#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
//#include <RcppEigenStubs.h>
#include <RcppEigenWrap.h>
//#include <Eigen/UmfPackSupport>
//#include <RcppEigenCholmod.h>
#include "path.h"
#include "Connectedness.h"

struct coeffLoc {
	int off;
	int r, c;

coeffLoc(int _off, int _r, int _c) :
	off(_off), r(_r), c(_c) {}
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
		struct independentGroup *ig;
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

		struct ApcIO : PathCalcIO {
			independentGroup &par;
			int clumpSize;
			ApcIO(independentGroup &_par) : par(_par), clumpSize(_par.clumpSize) {}
			virtual void recompute(FitContext *fc);
			virtual unsigned getVersion(FitContext *fc);
			virtual void refresh(FitContext *fc,
													 Eigen::Ref<Eigen::MatrixXd> mat, double sign);
		};

		struct SpcIO : PathCalcIO {
			independentGroup &par;
			int clumpSize;
			SpcIO(independentGroup &_par) : par(_par), clumpSize(_par.clumpSize) {}
			virtual void recompute(FitContext *fc);
			virtual unsigned getVersion(FitContext *fc);
			virtual void refresh(FitContext *fc,
													 Eigen::Ref<Eigen::MatrixXd> mat, double sign);
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
		Eigen::VectorXd                  expectedVec;
		Eigen::MatrixXd                  fullCov;   // rename, latents are filtered out
		bool                             analyzedCov;
		//Cholmod< Eigen::SparseMatrix<double> > covDecomp;
		//SimpCholesky< Eigen::MatrixXd >  covDecomp;
		Eigen::SparseMatrix<double>      fullS;
		std::vector<bool>                latentFilter; // false when latent or missing
    std::vector<bool>                isProductNode;

		ApcIO *aio;
		SpcIO *sio;
		PathCalc pcalc;
		double                           fit;  // most recent fit for debugging

		independentGroup(class state *_st, int size, int _clumpSize);
		independentGroup(independentGroup *ig);
		~independentGroup() {
			if (aio) delete aio;
			if (sio) delete sio;
		}
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
		struct omxExpectation *homeEx;
		std::set<struct omxExpectation *> allEx;
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
	};
};

class omxRAMExpectation : public omxExpectation {
	typedef omxExpectation super;
	Eigen::VectorXi dataCols;  // composition of F permutation and expectation->dataColumns
	std::vector<const char *> dataColNames;
	std::vector< omxThresholdColumn > thresholds;
	std::vector<int> exoDataColumns; // index into omxData
	Eigen::VectorXd exoPredMean;

	struct MpcIO : PathCalcIO {
		omxMatrix *M0;
		MpcIO() {}
		virtual void recompute(FitContext *fc);
		virtual unsigned getVersion(FitContext *fc);
		virtual void refresh(FitContext *fc,
												 Eigen::Ref<Eigen::MatrixXd> mat, double sign);
	};

	struct ApcIO : PathCalcIO {
		omxMatrix *A0;
		std::vector<coeffLoc> &vec;
		ApcIO(std::vector<coeffLoc> &_vec) : vec(_vec) {}
		virtual void recompute(FitContext *fc);
		virtual unsigned getVersion(FitContext *fc);
		virtual void refresh(FitContext *fc,
												 Eigen::Ref<Eigen::MatrixXd> mat, double sign);
	};

	struct SpcIO : PathCalcIO {
		omxMatrix *S0;
		std::vector<coeffLoc> &vec;
		SpcIO(std::vector<coeffLoc> &_vec) : vec(_vec) {}
		virtual void recompute(FitContext *fc);
		virtual unsigned getVersion(FitContext *fc);
		virtual void refresh(FitContext *fc,
												 Eigen::Ref<Eigen::MatrixXd> mat, double sign);
	};

	MpcIO *mio;
	ApcIO *aio;
	SpcIO *sio;

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

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *slope;       // exogenous predictor slopes
	omxMatrix *A, *S, *F, *M;

	int verbose;
	int rampartCycleLimit;
	int rampartUnitLimit;
	int maxDebugGroups;
	bool useSufficientSets;
	int optimizeMean;
	bool rampartEnabled() { return (rampartCycleLimit == NA_INTEGER || rampartCycleLimit > 0) && !forceSingleGroup; };
	double logDetObserved;
	double n;
	double *work;
	int lwork;

	std::vector< omxMatrix* > between;
	RelationalRAMExpectation::state *rram;
	bool forceSingleGroup;

	void studyF();
	void studyExoPred();

	virtual void init();
	virtual void compute(FitContext *fc, const char *what, const char *how);
	virtual omxMatrix *getComponent(const char*);
	virtual void populateAttr(SEXP expectation);
	virtual const std::vector<const char *> &getDataColumnNames() const { return dataColNames; };
	virtual const Eigen::Map<DataColumnIndexVector> getDataColumns() {
		return Eigen::Map<DataColumnIndexVector>(dataCols.data(), numDataColumns);
	}
	virtual std::vector< omxThresholdColumn > &getThresholdInfo() { return thresholds; }
	virtual void invalidateCache();
	virtual void generateData(FitContext *fc, MxRList &out);
	virtual void flatten(FitContext *fc);
	virtual void getExogenousPredictors(std::vector<int> &out);
};

namespace RelationalRAMExpectation {
	inline int state::verbose() const
	{
		return ((omxRAMExpectation*) homeEx)->verbose;
	}

	inline int independentGroup::verbose() const { return st.verbose(); };

};

#endif
