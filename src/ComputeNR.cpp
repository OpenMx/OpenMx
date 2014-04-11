/*
 *  Copyright 2013 The OpenMx Project
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

#include <valarray>

#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "matrix.h"

#include "Eigen/Eigenvalues"

class ComputeNR : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMatrix;

	int maxIter;
	double tolerance;
	int verbose;
	double priorSpeed;

	void lineSearch(FitContext *fc, double *maxAdj, double *maxAdjSigned, int *maxAdjParam, double *improvement);

public:
	virtual void initFromFrontend(SEXP rObj);
	virtual omxFitFunction *getFitFunction();
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class omxCompute *newComputeNewtonRaphson()
{
	return new ComputeNR();
}

void ComputeNR::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	if (!fitMatrix->fitFunction->hessianAvailable ||
	    !fitMatrix->fitFunction->gradientAvailable) {
		Rf_error("Newton-Raphson requires derivatives");
	}

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) Rf_error("tolerance must be positive");

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
}

omxFitFunction *ComputeNR::getFitFunction()
{ return fitMatrix->fitFunction; }

void omxApproxInvertPosDefTriangular(int dim, double *hess, double *ihess, double *stress)
{
	int info;
	int retries = 0;
	const int maxRetries = 31; // assume >=32 bit integers
	double adj = 0;
	do {
		memcpy(ihess, hess, sizeof(double) * dim * dim);

		if (retries >= 1) {
			int th = maxRetries - retries;
			if (th > 0) {
				adj = 1.0/(1 << th);
			} else {
				adj = (1 << -th);
			}
			for (int px=0; px < dim; ++px) {
				ihess[px * dim + px] += adj;
			}
		}

		Matrix ihessMat(ihess, dim, dim);
		info = InvertSymmetricPosDef(ihessMat, 'L');
		if (info == 0) break;
	} while (++retries < maxRetries * 1.5);

	if (info > 0) {
		// or just set stress to something high and return? TODO
		omxRaiseErrorf(globalState, "Hessian is not even close to positive definite (order %d)", info);
		return;
	}

	if (stress) *stress = adj;
}

void omxApproxInvertPackedPosDefTriangular(int dim, int *mask, double *packedHess, double *stress)
{
	int mdim = 0;
	for (int dx=0; dx < dim; ++dx) if (mask[dx]) mdim += 1;
	if (mdim == 0) {
		*stress = 0;
		return;
	}

	std::vector<double> hess(mdim * mdim, 0.0);
	for (int d1=0, px=0, m1=-1; d1 < dim; ++d1) {
		if (mask[d1]) ++m1;
		for (int d2=0, m2=-1; d2 <= d1; ++d2) {
			if (mask[d2]) ++m2;
			if (mask[d1] && mask[d2]) {
				hess[m2 * mdim + m1] = packedHess[px];
			}
			++px;
		}
	}

	std::vector<double> ihess(mdim * mdim);
	omxApproxInvertPosDefTriangular(mdim, hess.data(), ihess.data(), stress);

	for (int d1=0, px=0, m1=-1; d1 < dim; ++d1) {
		if (mask[d1]) ++m1;
		for (int d2=0, m2=-1; d2 <= d1; ++d2) {
			if (mask[d2]) ++m2;
			if (mask[d1] && mask[d2]) {
				packedHess[px] = *stress? 0 : ihess[m2 * mdim + m1];
			}
			++px;
		}
	}
}

void pda(const double *ar, int rows, int cols);

void ComputeNR::lineSearch(FitContext *fc, double *maxAdj, double *maxAdjSigned, int *maxAdjParam, double *improvement)
{
	const size_t numParam = varGroup->vars.size();
	const double epsilon = .5;
	bool steepestDescent = false;

	Eigen::Map<Eigen::VectorXd> prevEst(fc->est, numParam);

	Global->checkpointPrefit(fc, fc->est, false);
	omxFitFunctionCompute(fitMatrix->fitFunction,
			      FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT | FF_COMPUTE_IHESSIAN, fc);
	Global->checkpointPostfit(fc);

	const double refFit = fitMatrix->data[0];
	Eigen::VectorXd searchDir(fc->ihessGradProd());
	double targetImprovement = searchDir.dot(fc->grad);
	if (targetImprovement < tolerance) {
		steepestDescent = true;
		if (0 && fc->grad.norm() > tolerance) {
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
			es.compute(fc->getDenseIHess());
			Eigen::VectorXd ev(es.eigenvalues());
			mxLog("ihess EV");
			pda(ev.data(), 1, numParam);
			mxLog("ihess");
			Eigen::MatrixXd ihess(fc->getDenseIHess());
			ihess.triangularView<Eigen::StrictlyLower>() = ihess.transpose().triangularView<Eigen::StrictlyLower>();
			pda(ihess.data(), numParam, numParam);
		}
		searchDir = fc->grad;
		targetImprovement = searchDir.norm();
		if (targetImprovement < tolerance) return;
	}
	
	// This is based on the Goldstein test. However, we don't enforce
	// a lower bound on the improvement.

	int probeCount = 0;
	double speed = priorSpeed;
	Eigen::VectorXd trial;
	trial.resize(numParam);
	double bestSpeed = 0;
	double bestImproved = 0;
	double goodness = 0;

	while (++probeCount < 16) {
		trial = prevEst - speed * searchDir;
		fc->copyParamToModel(globalState, trial.data());
		omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, fc);
		if (!std::isfinite(fitMatrix->data[0])) {
			speed *= .5;
			continue;
		}
		const double improved = refFit - fitMatrix->data[0];
		if (improved <= 0) {
			speed *= .5;
			continue;
		}
		bestImproved = improved;
		bestSpeed = speed;
		goodness = improved / (speed * targetImprovement);
		if (verbose >= 3) mxLog("initial guess: speed %f for improvement %.3g goodness %f",
					bestSpeed, bestImproved, goodness);
		break;
	}
	if (bestSpeed == 0) return;

	if (speed < 1 && goodness < epsilon) {
		int retries = 15; // search up to 2*speed
		speed *= 1.05;
		while (--retries > 0 && goodness < epsilon) {
			++probeCount;
			trial = prevEst - speed * searchDir;
			fc->copyParamToModel(globalState, trial.data());
			omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_FIT, fc);
			if (!std::isfinite(fitMatrix->data[0])) break;
			const double improved = refFit - fitMatrix->data[0];
			if (bestImproved >= improved) break;
			bestImproved = improved;
			bestSpeed = speed;
			goodness = improved / (speed * targetImprovement);
		}
	}

	if (verbose >= 3) mxLog("Newton-Raphson: steepestDescent %d probes %d speed %f improved %.3g",
				steepestDescent, probeCount, bestSpeed, bestImproved);
	if (!steepestDescent) priorSpeed = bestSpeed;

	trial = prevEst - bestSpeed * searchDir;

	*maxAdj = 0;
	for (size_t px=0; px < numParam; ++px) {
		double oldEst = fc->est[px];
		double badj = fabs(oldEst - trial(px));
		if (*maxAdj < badj) {
			*maxAdj = badj;
			*maxAdjSigned = oldEst - trial(px);
			*maxAdjParam = px;
		}
	}
	memcpy(fc->est, trial.data(), sizeof(double) * numParam);

	*improvement = bestImproved;
}

void ComputeNR::computeImpl(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		Rf_error("Model has no free parameters");
		return;
	}

	fc->flavor.assign(numParam, NULL);

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PARAMFLAVOR, fc);

	// flavor used for debug output only
	for (size_t px=0; px < numParam; ++px) {
		if (!fc->flavor[px]) fc->flavor[px] = "?";
	}

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

	priorSpeed = 1;
	int startIter = fc->iterations;
	bool converged = false;
	double maxAdj = 0;
	double maxAdjSigned = 0;
	int maxAdjParam = -1;
	const char *maxAdjFlavor = "?";

	if (verbose >= 2) {
		mxLog("Welcome to Newton-Raphson (tolerance %.3g, max iter %d)",
		      tolerance, maxIter);
	}
	while (1) {
		fc->iterations += 1;
		int iter = fc->iterations - startIter;
		if (verbose >= 2) {
			if (iter == 1) {
				mxLog("Begin %d/%d iter of Newton-Raphson", iter, maxIter);
			} else {
				const char *pname = "none";
				if (maxAdjParam >= 0) pname = fc->varGroup->vars[maxAdjParam]->name;
				mxLog("Begin %d/%d iter of Newton-Raphson (prev maxAdj %.3g for %s %s)",
				      iter, maxIter, maxAdjSigned, maxAdjFlavor, pname);
			}
		}

		fc->grad = Eigen::VectorXd::Zero(fc->numParam);
		fc->clearHessian();

		maxAdj = 0;
		double improvement = 0;
		lineSearch(fc, &maxAdj, &maxAdjSigned, &maxAdjParam, &improvement);

		converged = improvement < tolerance;
		maxAdjFlavor = fc->flavor[maxAdjParam];

		fc->copyParamToModel(globalState);

		R_CheckUserInterrupt();

		if (converged || iter >= maxIter || isErrorRaised(globalState)) break;
	}

	if (converged) {
		fc->inform = INFORM_CONVERGED_OPTIMUM;
		fc->wanted |= FF_COMPUTE_BESTFIT;
		if (verbose >= 1) {
			int iter = fc->iterations - startIter;
			mxLog("Newton-Raphson converged in %d cycles", iter);
		}
	} else {
		fc->inform = INFORM_ITERATION_LIMIT;
		if (verbose >= 1) {
			int iter = fc->iterations - startIter;
			mxLog("Newton-Raphson failed to converge after %d cycles", iter);
		}
	}
}

void ComputeNR::reportResults(FitContext *fc, MxRList *slots, MxRList *output)
{
	omxPopulateFitFunction(fitMatrix, output);
}
