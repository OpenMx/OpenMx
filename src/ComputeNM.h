/*
 *  Copyright 2017-2019 by the individuals mentioned in the source code history
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "Compute.h"
#include <Eigen/Core>
#include <Eigen/Dense>

class omxComputeNM : public omxCompute {
	typedef omxCompute super;
	bool nudge;
	
public:
	omxMatrix *fitMatrix;
	int verbose;
	int maxIter;
	bool defaultMaxIter;
	double alpha, betao, betai, gamma, sigma;
	double bignum;
	int iniSimplexType; //regular=1, right=2, smartRight=3, random=4
	double iniSimplexEdge;
	Eigen::MatrixXd iniSimplexMat;
	std::vector< const char* > iniSimplexColnames;
	bool centerIniSimplex;
	bool greedyMinimize, altContraction;
	double degenLimit;
	Eigen::Vector2i stagnCtrl;
	bool validationRestart;
	double xTolProx, fTolProx;
	bool doPseudoHessian;
	int ineqConstraintMthd, eqConstraintMthd;
	double feasTol;
	double backtrackCtrl1;
	int backtrackCtrl2;

	Eigen::MatrixXd verticesOut;
	Eigen::VectorXd fvalsOut;
	Eigen::VectorXi vertexInfeasOut;
	double fproxOut, xproxOut, bestfitOut;
	//Eigen::MatrixXd equalityOut, inequalityOut;

	Eigen::MatrixXd pseudohess, phpts, phFvals, Xout;
	Eigen::VectorXi phInfeas;

	Eigen::VectorXd simplexGradient;

	omxComputeNM();
	virtual void initFromFrontend(omxState *, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};



class NelderMeadOptimizerContext{
private:
	void copyBounds();
	FitContext *fc;
public:
	NelderMeadOptimizerContext(FitContext* fc, omxComputeNM* nmo);
	void countConstraintsAndSetupBounds();
	void copyParamsFromFitContext(double *ocpars);
	
	template <typename T1>
	void copyParamsFromOptimizer(Eigen::MatrixBase<T1> &x, FitContext* fc2){
		fc2->setEstFromOptimizer(x);
	}
	
	omxComputeNM* NMobj;
	const int numFree;
	int verbose;
	int numIneqC;
	int numEqC;
	int n; //<--number of free parameters minus number of equality constraints
	int backtrackSteps;
	int maxIter;
	int itersElapsed;
	int iniSimplexType; //regular=1, right=2, smartRight=3, random=4
	double iniSimplexEdge;
	bool centerIniSimplex;
	double fr, fe, foc, fic;
	int badr, bade, badoc, badic, badsc;
	int restartsUsed;
	bool failedContraction;
	bool needFullSort;
	//double avgFitValPrev, avgFitValCurr;
	int unchangedx0count;
	double fit2beat, bestfit;
	int estInfeas;
	int statuscode;
	double bignum;
	double rho;
	bool addPenalty;
	int ineqConstraintMthd, eqConstraintMthd;
	bool checkRedundantEqualities;
	
	bool checkBounds(Eigen::VectorXd &x);
	void enforceBounds(Eigen::VectorXd &x);
	void evalIneqC();
	void evalEqC();
	double evalFit(Eigen::VectorXd &x);
	void checkNewPointInfeas(Eigen::VectorXd &x, Eigen::Vector2i &ifcr);
	void evalFirstPoint(Eigen::VectorXd &x, double &fv, int &infeas);
	void evalNewPoint(Eigen::VectorXd &newpt, Eigen::VectorXd oldpt, double &fv, int &newInfeas, int oldInfeas);
	void jiggleCoord(Eigen::VectorXd &xin, Eigen::VectorXd &xout, double scal);
	void invokeNelderMead();
	void initializeSimplex(Eigen::VectorXd startpt, double edgeLength, bool isRestart);
	void fullSort();
	void fastSort();
	void simplexTransformation();
	bool checkConvergence();
	bool checkProgress();
	void printProblemState();
	void printNewPoint(Eigen::VectorXd &x, double fv, int isbad);
	void calculatePseudoHessian();
	void finalize();
	
	std::vector<Eigen::VectorXd> vertices;	
	Eigen::VectorXd est;
	Eigen::VectorXd fvals;
	Eigen::VectorXi vertexInfeas;
	Eigen::VectorXd solLB;
	Eigen::VectorXd solUB;
	Eigen::VectorXd equality;
	Eigen::VectorXd inequality;
	Eigen::Vector2i feasCheckResults;
	Eigen::VectorXd subcentroid;
	Eigen::VectorXd eucentroidPrev, eucentroidCurr;
	Eigen::VectorXd xr, xe, xoc, xic;
	Eigen::VectorXd oldWorstVertex;
	Eigen::MatrixXd iniSimplexMat;
	Eigen::VectorXd tentativpt;
	//Eigen::VectorXd gdpt;

	GradientOptimizerContext subsidiarygoc;
	void *extraData;
	int gdfsIter;
};

double nmgdfso(unsigned n, const double *x, double *grad, void *f_data);
void omxInvokeSLSQPfromNelderMead(NelderMeadOptimizerContext* nmoc, Eigen::VectorXd &gdpt);

struct NldrMd_equality_functional {
	NelderMeadOptimizerContext* nmoc;
	FitContext* fc;
	
	NldrMd_equality_functional(NelderMeadOptimizerContext* _nmoc, FitContext* _fc) : nmoc(_nmoc), fc(_fc) {};
	
	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		nmoc->copyParamsFromOptimizer(x, fc);
		fc->solEqBFun(false, nmoc->verbose);
		result = fc->equality;
	}
	
	template <typename T1, typename T2, typename T3>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result, Eigen::MatrixBase<T3> &jacobian) const {
		nmoc->copyParamsFromOptimizer(x, fc);
		fc->analyticEqJacTmp.resize(jacobian.rows(), jacobian.cols());
		fc->solEqBFun(true, nmoc->verbose);
		result = fc->equality;
		jacobian = fc->analyticEqJacTmp;
	}
};
