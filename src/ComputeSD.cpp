/* Steepest Descent optimizer for unconstrained problems*/

#include "ComputeSD.h"
#include "finiteDifferences.h"

struct fit_functional {
	GradientOptimizerContext &goc;

	fit_functional(GradientOptimizerContext &goc) : goc(goc) {};

	template <typename T1>
	double operator()(Eigen::MatrixBase<T1>& x) const {
		int mode = 0;
		return goc.solFun(x.derived().data(), &mode);
	}
};

void omxSD(GradientOptimizerContext &rf)
{
	int maxIter = rf.maxMajorIterations;
	if (maxIter == -1) maxIter = 50000;

	FitContext *fc = rf.fc;
    int iter = 0;
	double priorSpeed = 1.0, shrinkage = 0.7;
    rf.setupSimpleBounds();
    rf.informOut = INFORM_UNINITIALIZED;

    int mode = 0;
    rf.solFun(fc->est, &mode);
    if (mode == -1) {
	    rf.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
	    return;
    }
    double refFit = fc->fit;

    fc->grad.resize(fc->numParam);

    fit_functional ff(rf);
    Eigen::Map< Eigen::VectorXd > currEst(fc->est, fc->numParam);
    Eigen::VectorXd majorEst = currEst;

    while(++iter < maxIter && !isErrorRaised()) {
	    fc->iterations += 1;
	    gradient_with_ref(rf.gradientAlgo, rf.gradientIterations, rf.gradientStepSize,
			      ff, refFit, majorEst, fc->grad);

	    if (rf.verbose >= 3) mxPrintMat("grad", fc->grad);

        if(fc->grad.norm() == 0)
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
	Eigen::VectorXd searchDir = fc->grad;
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
			refFit = fc->fit;
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
				refFit = fc->fit;
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
    currEst = majorEst;
    if ((fc->grad.array().abs() > 0.1).any()) {
	    rf.informOut = INFORM_NOT_AT_OPTIMUM;
    }
    if (iter >= maxIter - 1) {
            rf.informOut = INFORM_ITERATION_LIMIT;
            if(rf.verbose >= 2) mxLog("Maximum iteration achieved!");
    }
    if(rf.verbose >= 1) mxLog("Status code : %i", rf.informOut);
}

/* Steepest Descent optimizer for constrained problems */

template <typename T1, typename T2>
static void SD_AL(GradientOptimizerContext &rf, double rho,
			Eigen::MatrixBase<T1> &lambda, Eigen::MatrixBase<T2> &mu)
{
    rf.fc->copyParamToModel();
    ComputeFit("SD_AL", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);

    rf.solEqBFun();
    rf.myineqFun();

    for (size_t i = 0; i < unsigned(rf.equality.size()); ++i)
    {
        rf.fc->fit += 0.5 * rho * (rf.equality[i] + lambda[i] / rho) * (rf.equality[i] + lambda[i] / rho);
    }

    for (size_t i = 0; i < unsigned(rf.inequality.size()); ++i)
    {
        rf.fc->fit += 0.5 * rho * std::max(0.0,(rf.inequality[i] + mu[i] / rho)) * std::max(0.0,(rf.inequality[i] + mu[i] / rho));
    }
    return;
}

template <typename T1, typename T2>
static void SD_grad_AL(GradientOptimizerContext &rf, double eps, double rho,
		       Eigen::MatrixBase<T1> &lambda, Eigen::MatrixBase<T2> &mu)
{
    SD_AL(rf, rho, lambda, mu);
    const double refFit = rf.fc->fit;

    Eigen::Map< Eigen::VectorXd > p1(rf.fc->est, rf.fc->numParam);
    Eigen::VectorXd p2 = p1;
    Eigen::VectorXd grad(rf.fc->numParam);

    for (int px = 0; px < int(rf.fc->numParam); px++) {
        p1[px] += eps;
        SD_AL(rf, rho, lambda, mu);
        grad[px] = (rf.fc->fit - refFit) / eps;
        p1 = p2;
    }
    rf.fc->fit = refFit;
    rf.fc->grad = grad;
}

//Armijo's Rule implemented
template <typename T1, typename T2>
static void omxSD_AR(GradientOptimizerContext &rf, int maxIter, double rho,
			Eigen::MatrixBase<T1> &lambda, Eigen::MatrixBase<T2> &mu, double &alpha)
{
	int iter = 0;
	double epsilon = 0.2, eta = 2, eps_grad = 1e-9, grad_tol = 1e-3;
    rf.setupSimpleBounds();
    rf.informOut = INFORM_UNINITIALIZED;
    rf.fc->copyParamToModel();

    SD_grad_AL(rf, eps_grad, rho, lambda, mu);

	while(iter < maxIter && !isErrorRaised())
	{
        alpha = std::min(1.0, alpha);
        iter++;
        SD_grad_AL(rf, eps_grad, rho, lambda, mu);
        Eigen::VectorXd searchDir = rf.fc->grad;

        if (std::isnan(rf.fc->fit))
        {
            rf.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
            break;
        }

        if(rf.fc->grad.norm() < grad_tol)
        {
            rf.informOut = INFORM_CONVERGED_OPTIMUM;
            if(rf.verbose >= 2) mxLog("After %i iterations, gradient achieves zero!", iter);
            break;
        }

        Eigen::Map< Eigen::VectorXd > currEst(rf.fc->est, rf.fc->numParam);
        Eigen::VectorXd prevEst = currEst;

        currEst = prevEst - alpha * searchDir / searchDir.norm();
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

        SD_AL(rf, rho, lambda, mu);
        double phi = rf.fc->fit;

        currEst = prevEst;
        SD_AL(rf, rho, lambda, mu);
        double dashline = rf.fc->fit - epsilon * alpha * searchDir.norm();

        if(phi <= dashline)
        {
            do{
                alpha *= eta;
                currEst = prevEst;
                SD_AL(rf, rho, lambda, mu);
                dashline = rf.fc->fit - epsilon * alpha * searchDir.norm();

                currEst = prevEst - alpha * searchDir / searchDir.norm();
                currEst = currEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);

		rf.checkActiveBoxConstraints(currEst);

                SD_AL(rf, rho, lambda, mu);
                phi = rf.fc->fit;
            }while(phi <= dashline);
            currEst = prevEst - alpha / eta * searchDir / searchDir.norm();
            currEst = currEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);
            SD_AL(rf, rho, lambda, mu);
        }
        else
        {
            do{
                alpha /= eta;
                currEst = prevEst;
                SD_AL(rf, rho, lambda, mu);
                dashline = rf.fc->fit - epsilon * alpha * searchDir.norm();

                currEst = prevEst - alpha * searchDir / searchDir.norm();
                currEst = currEst.cwiseMax(rf.solLB).cwiseMin(rf.solUB);

		rf.checkActiveBoxConstraints(currEst);

                SD_AL(rf, rho, lambda, mu);
                phi = rf.fc->fit;
            }while(phi > dashline);
        }

        if (iter == maxIter - 1)
        {
            rf.informOut = INFORM_ITERATION_LIMIT;
            if(rf.verbose >= 2) mxLog("Maximum iteration achieved!");
            break;
        }

    }
    if(rf.verbose >= 1) mxLog("Status code : %i", rf.informOut);
    return;
}

void omxSD_AL(GradientOptimizerContext &rf)
{
	int maxIter = rf.maxMajorIterations;
	if (maxIter == -1) maxIter = 50000;

    double ICM = HUGE_VAL;
    /* magic parameters from Birgin & Martinez */
    const double tau = 0.5, gam = 10;
    const double lam_min = -1e20, lam_max = 1e20, mu_max = 1e20;
    double alpha = 1;

    rf.fc->copyParamToModel();
    ComputeFit("SD_AL", rf.fitMatrix, FF_COMPUTE_FIT, rf.fc);
    rf.solEqBFun();
    rf.myineqFun();

    size_t ineq_size = rf.inequality.size(), eq_size = rf.equality.size();

    // initialize penalty parameter rho and the Lagrange multipliers lambda and mu
    double eq_norm = 0, ineq_norm = 0;

    for(size_t i = 0; i < eq_size; i++)
    {
      eq_norm += rf.equality[i] * rf.equality[i];
    }

    for(size_t i = 0; i < ineq_size; i++)
    {
      ineq_norm += std::max(0.0, rf.inequality[i]) * std::max(0.0, rf.inequality[i]);
    }

    double rho = std::max(1e-6, std::min(10.0, (2 * std::abs(rf.fc->fit) / (eq_norm + ineq_norm))));

    Eigen::VectorXd lambda(eq_size);
    lambda.setZero();
    Eigen::VectorXd mu(ineq_size);
    mu.setZero();
    int iter = 0, minorIter = 1000;
    Eigen::VectorXd V(ineq_size);

    while (1) {
        iter++;
        double ICM_tol = 1e-4;
        double prev_ICM = ICM;
        ICM = 0;

        omxSD_AR(rf, minorIter, rho, lambda, mu, alpha);

        if(rf.informOut == INFORM_STARTING_VALUES_INFEASIBLE) return;

        rf.fc->copyParamToModel();
        rf.solEqBFun();
        rf.myineqFun();

        for(size_t i = 0; i < eq_size; i++){
            lambda[i] = std::min(std::max(lam_min, (lambda[i] + rho * rf.equality[i])), lam_max);
            ICM = std::max(ICM, std::abs(rf.equality[i]));
        }

        for(size_t i = 0; i < ineq_size; i++){
            mu[i] = std::min(std::max(0.0, (mu[i] + rho * rf.inequality[i])),mu_max);
        }

        for(size_t i = 0; i < ineq_size; i++){
            V[i] = std::max(rf.inequality[i], (-mu[i] / rho));
            ICM = std::max(ICM, std::abs(V[i]));
        }

        if(!(iter == 1 || ICM <= tau * prev_ICM))
        {
            rho *= gam;
        }

        if(ICM < ICM_tol)
        {
            if(rf.verbose >= 1) mxLog("Augmented lagrangian coverges!");
            return;
        }
	if (iter >= maxIter) {
		rf.informOut = INFORM_ITERATION_LIMIT;
		break;
	}
    }
}






