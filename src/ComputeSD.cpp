/* Steepest Descent optimizer for unconstrained problems*/

#include "ComputeGD.h"
#include "ComputeSD.h"
#include "EnableWarnings.h"

struct fit_functional {
	GradientOptimizerContext &goc;

	fit_functional(GradientOptimizerContext &_goc) : goc(_goc) {};

	double operator()(double *x, int thrId) const {
		int mode = 0;
		return goc.solFun(x, &mode);
	}
};

void omxSD(GradientOptimizerContext &rf)
{
	int maxIter = rf.maxMajorIterations;
	if (maxIter == -1) maxIter = 50000;

	Eigen::VectorXd currEst(rf.numFree);
	rf.copyToOptimizer(currEst.data());

    int iter = 0;
	double priorSpeed = 1.0, shrinkage = 0.7;
    rf.setupSimpleBounds();
    rf.informOut = INFORM_UNINITIALIZED;

    {
	    int mode = 0;
	    rf.solFun(currEst.data(), &mode);
	    if (mode == -1) {
		    rf.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
		    return;
	    }
    }
    double refFit = rf.getFit();

    rf.grad.resize(rf.numFree);

    fit_functional ff(rf);
    Eigen::VectorXd majorEst = currEst;

    while(++iter < maxIter && !isErrorRaised()) {
	    rf.numericalGradientWithRef(majorEst);

	    if (rf.verbose >= 3) mxPrintMat("grad", rf.grad);

        if(rf.grad.norm() == 0)
        {
            rf.informOut = INFORM_CONVERGED_OPTIMUM;
            if(rf.verbose >= 2) mxLog("After %i iterations, gradient achieves zero!", iter);
            break;
        }

        int retries = 300;
        double speed = std::min(priorSpeed, 1.0);
	double bestSpeed = speed;
	bool foundBetter = false;
	Eigen::VectorXd bestEst(majorEst.size());
	Eigen::VectorXd prevEst(majorEst.size());
	Eigen::VectorXd searchDir = rf.grad;
	searchDir /= searchDir.norm();
	prevEst.setConstant(nan("uninit"));
        while (--retries > 0 && !isErrorRaised()){
		Eigen::VectorXd nextEst = majorEst - speed * searchDir;
		nextEst = nextEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);

		if (nextEst == prevEst) break;
		prevEst = nextEst;

		rf.checkActiveBoxConstraints(nextEst);

		int mode = 0;
		double fit = rf.solFun(nextEst.data(), &mode);
		if (fit < refFit) {
			foundBetter = true;
			refFit = rf.getFit();
			bestSpeed = speed;
			bestEst = nextEst;
			break;
		}
		speed *= shrinkage;
        }

	if (false && foundBetter) {
		// In some tests, this did not help so it is not enabled.
		// It might be worth testing more.
		mxLog("trying larger step size");
		retries = 3;
		while (--retries > 0 && !isErrorRaised()){
			speed *= 1.01;
			Eigen::VectorXd nextEst = majorEst - speed * searchDir;
			nextEst = nextEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);
			rf.checkActiveBoxConstraints(nextEst);
			int mode = 0;
			double fit = rf.solFun(nextEst.data(), &mode);
			if (fit < refFit) {
				foundBetter = true;
				refFit = rf.getFit();
				bestSpeed = speed;
				bestEst = nextEst;
			}
		}
	}

        if (!foundBetter) {
            rf.informOut = INFORM_CONVERGED_OPTIMUM;
            if(rf.verbose >= 2) mxLog("After %i iterations, cannot find better estimation along the gradient direction", iter);
            break;
        }

	if (rf.verbose >= 2) mxLog("major fit %f bestSpeed %g", refFit, bestSpeed);
	majorEst = bestEst;
	priorSpeed = bestSpeed * 1.1;
    }
    rf.est = majorEst;
    if ((rf.grad.array().abs() > 0.1).any()) {
	    rf.informOut = INFORM_NOT_AT_OPTIMUM;
    }
    if (iter >= maxIter - 1) {
            rf.informOut = INFORM_ITERATION_LIMIT;
            if(rf.verbose >= 2) mxLog("Maximum iteration achieved!");
    }
    if(rf.verbose >= 1) mxLog("Status code : %i", rf.informOut);
}
