/*
 *  Copyright 2017 The OpenMx Project
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
	omxMatrix *fitMatrix;
	bool nudge;
	
public:
	int verbose;
	int maxIter;
	bool defaultMaxIter;
	std::vector<int> excludeVars;
	double alpha, betao, betai, gamma, sigma;
	double bignum;
	int iniSimplexType; //regular=1, right=2, smartRight=3, random=4
	double iniSimplexEdge;
	Eigen::MatrixXd iniSimplexMtx;
	bool greedyMinimize, altContraction;
	double degenLimit;
	Eigen::Vector2i stagnationCtrl;
	bool validationRestart;
	double xTolProx, fTolProx, xTolRelChange, fTolRelChange;
	bool doPseudoHessian;
	int ineqConstraintMthd, eqConstraintMthd;
	double feasTol;
	omxComputeNM();
	virtual void initFromFrontend(omxState *, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	//virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};



class NelderMeadOptimizerContext{
private:
	void copyBounds();
	int countNumFree();
	FitContext *fc;
public:
	NelderMeadOptimizerContext(FitContext* fc, omxComputeNM* nmo);
	void countConstraintsAndSetupBounds();
	void copyParamsFromFitContext(double *ocpars);
	void copyParamsFromOptimizer(Eigen::VectorXd &x, FitContext* fc2);

	
	omxComputeNM* NMobj;
	const int numFree;
	int verbose;
	int numIneqC;
	int numEqC;
	int n; //<--number of free parameters minus number of equality constraints
	int backtrackSteps;
	int itersElapsed;
	double fr, fe, foc, fic;
	int badr, bade, badoc, badic, badsc;
	int restartsUsed;
	bool failedContraction;
	
	bool checkBounds(Eigen::VectorXd &x);
	void enforceBounds(Eigen::VectorXd &x);
	void evalIneqC();
	void evalEqC();
	double evalFit(Eigen::VectorXd &x);
	void checkNewPointInfeas(Eigen::VectorXd &x, Eigen::Vector2i &ifcr);
	void evalFirstPoint(Eigen::VectorXd &x, double fv, int infeas);
	void evalNewPoint(Eigen::VectorXd &newpt, Eigen::VectorXd &oldpt, double fv, int newInfeas, int oldInfeas);
	void jiggleCoord(Eigen::VectorXd &xin, Eigen::VectorXd &xout);
	void invokeNelderMead();
	void initializeSimplex(Eigen::VectorXd startpt, double edgeLength, bool isRestart);
	void fullSort();
	void fastSort();
	void simplexTransformation();
	bool checkConvergence();
	bool checkProgress();
	void printProblemState();
	void printNewPoint(Eigen::VectorXd &x, double fv, int isbad);
	/*void restart();
	void validationRestart();*/
	
	std::vector<Eigen::VectorXd> vertices;	
	Eigen::VectorXd est;
	//Eigen::MatrixXd vertices;
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
	
	bool needFullSort;
	//double avgFitValPrev, avgFitValCurr;
	int restartCount, unchangedx0count;
	
	
	
};
