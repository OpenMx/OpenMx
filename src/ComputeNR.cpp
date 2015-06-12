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

static const char engineName[] = "NnRn";

class ComputeNR : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMatrix;
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;

	int maxIter;
	double tolerance;
	int verbose;
	double priorSpeed;
	int minorIter;
	double refFit;

	void lineSearch(FitContext *fc, int iter, double *maxAdj, double *maxAdjSigned,
			int *maxAdjParam, double *improvement);

public:
	virtual void initFromFrontend(omxState *state, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class omxCompute *newComputeNewtonRaphson()
{
	return new ComputeNR();
}

void ComputeNR::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, state, "fitfunction");
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
		omxRaiseErrorf("Hessian is not even close to positive definite (order %d)", info);
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

void ComputeNR::lineSearch(FitContext *fc, int iter, double *maxAdj, double *maxAdjSigned,
			   int *maxAdjParam, double *improvement)
{
	*maxAdjParam = -1;
	const size_t numParam = varGroup->vars.size();
	bool steepestDescent = false;

	Eigen::VectorXd prevEst(numParam);
	memcpy(prevEst.data(), fc->est, numParam * sizeof(double));

	int want = FF_COMPUTE_GRADIENT | FF_COMPUTE_IHESSIAN;
	if (verbose >= 5) want |= FF_COMPUTE_HESSIAN;
	if (iter == 1) {
		want |= FF_COMPUTE_FIT;
	}

	ComputeFit(engineName, fitMatrix, want, fc);
	if (iter == 1) {
		refFit = fc->fit;
		if (!std::isfinite(refFit)) {
			fc->inform = INFORM_STARTING_VALUES_INFEASIBLE;
			return;
		}
	}

	double speed = std::min(priorSpeed * 1.5, 1.0);
	Eigen::VectorXd searchDir(fc->ihessGradProd());
	double targetImprovement = searchDir.dot(fc->grad);

	if (verbose >= 5) {
		fc->log(FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN);

		std::string buf;
		buf += "searchDir: c(";
		for (int vx=0; vx < searchDir.size(); ++vx) {
			buf += string_snprintf("%.5f", searchDir[vx]);
			if (vx < searchDir.size() - 1) buf += ", ";
		}
		buf += ")\n";
		mxLogBig(buf);
	}

	if (!std::isfinite(targetImprovement)) {
		if (verbose >= 4) mxLog("%s: target improvement %.4g is suspect, using steepest descent",
					name, targetImprovement);
		steepestDescent = true;
		// the steepness really provides no information due to unknown curvature
		searchDir = fc->grad / fc->grad.norm();
		if (!std::isfinite(searchDir.norm())) {
			if (verbose >= 2) {
				for (int px=0; px < fc->grad.size(); ++px) {
					if (std::isfinite(fc->grad[px])) continue;
					omxFreeVar *fv = fc->varGroup->vars[px];
					mxLog("%s=%f is infeasible; try adding bounds", fv->name, fc->est[px]);
				}
			}
			fc->inform = INFORM_BAD_DERIVATIVES;
			return;
		}
		targetImprovement = 1;
		speed = std::max(speed, .1);  // expect steepestDescent
	}
	
	// This is based on the Goldstein test. However, we don't enforce
	// a lower bound on the improvement.

	int probeCount = 0;
	Eigen::Map< Eigen::VectorXd > trial(fc->est, numParam);
	double bestSpeed = 0;
	double bestImproved = 0;
	double goodness = 0;
	double bestFit = 0;

	while (++probeCount < 16) {
		const double scaledTarget = speed * targetImprovement;
		if (!steepestDescent && scaledTarget / fabs(refFit) < tolerance) {
			trial = prevEst;
			return;
		}
		trial = (prevEst - speed * searchDir).cwiseMax(lbound).cwiseMin(ubound);
		++minorIter;
		fc->copyParamToModel();
		ComputeFit(engineName, fitMatrix, FF_COMPUTE_FIT, fc);
		if (verbose >= 4) mxLog("%s: speed %.3g for target %.3g fit %f ref %f",
					name, speed, scaledTarget, fc->fit, refFit);
		if (!std::isfinite(fc->fit)) {
			speed *= .1;
			continue;
		}
		const double improved = refFit - fc->fit;
		if (improved <= 0) {
			const double minSpeedReduction = .1;
			double howBad = refFit / (scaledTarget - improved);
			if (verbose >= 5) {
				mxLog("refFit %.2g / scaledTarget - improved %.2g = %.4g",
				      refFit, scaledTarget - improved, howBad);
			}
			if (howBad < minSpeedReduction && verbose >= 4) {
				mxLog("%s: scaledTarget is %2.g times over optimistic", name, 1/howBad);
			}
			if (howBad > minSpeedReduction) howBad = minSpeedReduction;
			speed *= howBad;
			continue;
		}
		bestImproved = improved;
		bestSpeed = speed;
		bestFit = fc->fit;
		goodness = improved / scaledTarget;
		if (verbose >= 3) mxLog("%s: viable speed %.2g for improvement %.3g goodness %f",
					name, bestSpeed, bestImproved, goodness);
		break;
	}
	if (bestSpeed == 0) {
		trial = prevEst;
		return;
	}

	const double epsilon = .3;
	if (!steepestDescent && speed < 1 && goodness < epsilon) {
		int retries = 8;
		while (--retries > 0) {
			speed *= 1.5;
			++probeCount;
			trial = (prevEst - speed * searchDir).cwiseMax(lbound).cwiseMin(ubound);
			++minorIter;
			fc->copyParamToModel();
			ComputeFit(engineName, fitMatrix, FF_COMPUTE_FIT, fc);
			if (!std::isfinite(fc->fit)) break;
			const double improved = refFit - fc->fit;
			if (bestImproved >= improved) break;
			double improvementOverBest = improved - bestImproved;
			if (verbose >= 4) {
				mxLog("%s: [%d] incr speed for fit improvement of %.2g",
				      name, retries, improvementOverBest);
			}
			bestFit = fc->fit;
			bestImproved = improved;
			bestSpeed = speed;
			if (improvementOverBest / fabs(refFit) < tolerance) break;
		}
	}

	if (verbose >= 3) mxLog("%s: using steepestDescent %d probes %d speed %f improved %.3g",
				name, steepestDescent, probeCount, bestSpeed, bestImproved);
	priorSpeed = bestSpeed;

	trial = (prevEst - bestSpeed * searchDir).cwiseMax(lbound).cwiseMin(ubound);

	*maxAdj = 0;
	for (size_t px=0; px < numParam; ++px) {
		const double oldEst = prevEst[px];
		double badj = fabs(oldEst - trial(px));
		if (*maxAdj < badj) {
			*maxAdj = badj;
			*maxAdjSigned = oldEst - trial(px);
			*maxAdjParam = px;
		}
	}

	*improvement = bestImproved;
	refFit = bestFit;
}

void ComputeNR::computeImpl(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		Rf_error("Model has no free parameters");
		return;
	}

	fc->inform = INFORM_UNINITIALIZED;
	fc->flavor.assign(numParam, NULL);

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PARAMFLAVOR, fc);

	// flavor used for debug output only
	for (size_t px=0; px < numParam; ++px) {
		if (!fc->flavor[px]) fc->flavor[px] = "?";
	}

	omxFitFunctionPreoptimize(fitMatrix->fitFunction, fc);

	// omxFitFunctionPreoptimize can change bounds
	lbound.resize(numParam);
	ubound.resize(numParam);
	for(int px = 0; px < int(numParam); px++) {
		lbound[px] = varGroup->vars[px]->lbound;
		ubound[px] = varGroup->vars[px]->ubound;
	}

	priorSpeed = 1;
	minorIter = 0;
	int startIter = fc->iterations;
	bool converged = false;
	double maxAdj = 0;
	double maxAdjSigned = 0;
	int maxAdjParam = -1;
	const char *maxAdjFlavor = "?";

	if (verbose >= 2) {
		mxLog("Welcome to Newton-Raphson (%d param, tolerance %.3g, max iter %d)",
		      (int)numParam, tolerance, maxIter);
		if (verbose >= 3) {
			mxPrintMat("lbound", lbound);
			mxPrintMat("ubound", ubound);
		}
	}
	while (1) {
		fc->iterations += 1;
		int iter = fc->iterations - startIter;
		if (verbose >= 2) {
			if (iter == 1) {
				mxLog("%s: iter %d/%d", name, iter, maxIter);
			} else {
				const char *pname = "none";
				if (maxAdjParam >= 0) pname = fc->varGroup->vars[maxAdjParam]->name;
				mxLog("%s: iter %d/%d (prev maxAdj %.3g for %s %s)",
				      name, iter, maxIter, maxAdjSigned, maxAdjFlavor, pname);
			}
		}

		fc->grad = Eigen::VectorXd::Zero(fc->numParam);
		fc->clearHessian();

		maxAdj = 0;
		double improvement = 0;
		lineSearch(fc, iter, &maxAdj, &maxAdjSigned, &maxAdjParam, &improvement);
		if (fc->inform != INFORM_UNINITIALIZED) {
			if (verbose >= 2) {
				mxLog("%s: line search failed with code %d", name, fc->inform);
			}
			return;
		}

		converged = improvement / fabs(refFit) < tolerance;
		if (maxAdjParam >= 0) maxAdjFlavor = fc->flavor[maxAdjParam];

		fc->copyParamToModel();

		if (converged || iter >= maxIter || isErrorRaised()) break;
	}

	if (converged) {
		fc->inform = INFORM_CONVERGED_OPTIMUM;
		fc->wanted |= FF_COMPUTE_BESTFIT;
		if (verbose >= 1) {
			int iter = fc->iterations - startIter;
			mxLog("%s: converged in %d cycles (%d minor iterations)", name, iter, minorIter);
		}
	} else {
		fc->inform = INFORM_ITERATION_LIMIT;
		if (verbose >= 1) {
			int iter = fc->iterations - startIter;
			mxLog("%s: failed to converge after %d cycles (%d minor iterations)",
			      name, iter, minorIter);
		}
	}
}

void ComputeNR::reportResults(FitContext *fc, MxRList *slots, MxRList *output)
{
	omxPopulateFitFunction(fitMatrix, output);
}
