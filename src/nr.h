/*
 *  Copyright 2013-2019 by the individuals mentioned in the source code history
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

#ifndef _OMX_NR_H_
#define _OMX_NR_H_

struct NewtonRaphsonObjective {
	bool converged;
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;

	virtual void init() {
		converged = false;
	}
	virtual bool isConverged() { return converged; };
	virtual double getFit()=0;
	virtual void resetDerivs() {};
	virtual const char *paramIndexToName(int px)=0;
	virtual void evaluateFit()=0;
	virtual void evaluateDerivs(int want)=0;
	virtual double *getParamVec()=0;
	virtual double *getGrad()=0;
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir)=0;  // ihess * grad
	virtual void reportBadDeriv() {};
	virtual void debugDeriv(const Eigen::Ref<Eigen::VectorXd> searchDir) {};
};

class NewtonRaphsonOptimizer {
	const char *name;
	int maxIter;
	double tolerance;
	int verbose;
	int iter;
	int numParam;
	double refFit; // current fit value
	double priorSpeed;
	double improvement;  // previous fit - refFit
	double maxAdj;
	double maxAdjSigned;
	int maxAdjParam;
	int minorIter;
	Eigen::VectorXd prevEst;
	Eigen::VectorXd searchDir;
	double relImprovement(double im) { return im / (1 + fabs(refFit)); }
	void lineSearch(NewtonRaphsonObjective &nro);
public:
	NewtonRaphsonOptimizer(const char *_name, int _maxIter, double tol, int _verbose) :
		name(_name), maxIter(_maxIter), tolerance(tol), verbose(_verbose) {};
	void operator()(NewtonRaphsonObjective &nro);
	int getIter() { return iter; };
	int getMinorIter() { return minorIter; };
};

#endif
