/* Steepest Descent optimizer for unconstrained problems*/

#include "ComputeSD.h"

static void SD_grad(GradientOptimizerContext &rf, double eps)
{
    rf.fc->copyParamToModel();
    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);

    const double refFit = rf.fc->fit;
    Eigen::Map< Eigen::VectorXd > p1(rf.fc->est, rf.fc->numParam);
    Eigen::VectorXd p2 = p1;
    Eigen::VectorXd grad(rf.fc->numParam);

    for (int px = 0; px < int(rf.fc->numParam); px++) {
        p1[px] += eps;
        rf.fc->copyParamToModel();
        ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
        grad[px] = (rf.fc->fit - refFit) / eps;
        p1 = p2;
    }
    rf.fc->copyParamToModel();
    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
    rf.fc->grad = grad;
}

static bool FitCompare(GradientOptimizerContext &rf, double speed)
{
    Eigen::Map< Eigen::VectorXd > currEst(rf.fc->est, rf.fc->numParam);
    Eigen::VectorXd prevEst = currEst;

    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
    if (isnan(rf.fc->fit))
    {
        rf.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
        return FALSE;
    }
    double refFit = rf.fc->fit;

    Eigen::VectorXd searchDir = rf.fc->grad;
    currEst = prevEst - speed * searchDir / searchDir.norm();
    currEst = currEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);
    if(rf.verbose >= 2){
        for(int index = 0; index < int(rf.fc->numParam); index++)
        {
            if(currEst[index] == rf.solLB[index])
                mxLog("paramter %i hit lower bound %f", index, rf.solLB[index]);
            if(currEst[index] == rf.solUB[index])
                mxLog("paramter %i hit upper bound %f", index, rf.solUB[index]);
        }
    }

    rf.fc->copyParamToModel();
    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
    double newFit = rf.fc->fit;

    if(newFit < refFit) return newFit < refFit;
    currEst = prevEst;
    rf.fc->copyParamToModel();
    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
    return newFit < refFit;
}

void omxSD(GradientOptimizerContext &rf, int maxIter)
{
    int iter = 0;
	double priorSpeed = 1.0, shrinkage = 0.7, epsilon = 1e-9;
    rf.setupSimpleBounds();
    rf.informOut = INFORM_UNINITIALIZED;
    rf.fc->copyParamToModel();
    ComputeFit("Steepest Descent", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);  // To check isErrorRaised()

	while(iter < maxIter && !isErrorRaised())
	{
        SD_grad(rf, epsilon);
        if(rf.fc->grad.norm() == 0)
        {
            rf.informOut = INFORM_CONVERGED_OPTIMUM;
            if(rf.verbose >= 2) mxLog("After %i iterations, gradient achieves zero!", iter);
            break;
        }
        bool findit = FitCompare(rf, priorSpeed);

        int retries = 300;
        double speed = priorSpeed;
        while (--retries > 0 && !findit && !isErrorRaised()){
            speed *= shrinkage;
            findit = FitCompare(rf, speed);
        }
        if(findit){
            priorSpeed = speed * 1.1;
            iter++;
            if(iter == maxIter){
                rf.informOut = INFORM_ITERATION_LIMIT;
                if(rf.verbose >= 2) mxLog("Maximum iteration achieved!");
                break;
            }
        }
        else if (iter == 0){
            if(rf.verbose >= 2) mxLog("Infeasbile starting values!");
            break;
        }
        else if (iter == maxIter - 1){
            rf.informOut = INFORM_ITERATION_LIMIT;
            if(rf.verbose >= 2) mxLog("Maximum iteration achieved!");
            break;
        }
        else {
            rf.informOut = INFORM_CONVERGED_OPTIMUM;
            if(rf.verbose >= 2) mxLog("After %i iterations, cannot find better estimation along the gradient direction", iter);
            break;
        }
    }
    if(rf.verbose == 1) mxLog("Status code : %i", rf.informOut);
    return;
}


