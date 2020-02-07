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

#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "matrix.h"
#include "nr.h"
#include "EnableWarnings.h"

void NewtonRaphsonOptimizer::operator()(NewtonRaphsonObjective &nro)
{
	nro.init();
	numParam = nro.lbound.size();
	prevEst.resize(numParam);
	searchDir.resize(numParam);
	maxAdj = 0;
	maxAdjSigned = 0;
	maxAdjParam = -1;
	priorSpeed = 1;
	iter = 0;
	minorIter = 0;

	if (verbose >= 2) {
		mxLog("Welcome to Newton-Raphson (%d param, tolerance %.3g, max iter %d)",
		      (int)numParam, tolerance, maxIter);
		if (verbose >= 3) {
			mxPrintMat("lbound", nro.lbound);
			mxPrintMat("ubound", nro.ubound);
		}
	}

	while (1) {
		iter += 1;
		if (verbose >= 2) {
			if (iter == 1) {
				mxLog("%s: iter %d/%d", name, iter, maxIter);
			} else {
				const char *pname = nro.paramIndexToName(maxAdjParam);
				mxLog("%s: iter %d/%d (prev maxAdj %.3g for %s)",
				      name, iter, maxIter, maxAdjSigned, pname);
			}
		}

		nro.resetDerivs();
		maxAdj = 0;
		maxAdjParam = -1;
		improvement = 0;
		lineSearch(nro);

		if (!std::isfinite(refFit)) return;

		nro.converged =
			relImprovement(improvement) < tolerance || iter >= maxIter;

		if (nro.isConverged()) break;
	}
}

void NewtonRaphsonOptimizer::lineSearch(NewtonRaphsonObjective &nro)
{
	bool steepestDescent = false;

	memcpy(prevEst.data(), nro.getParamVec(), numParam * sizeof(double));

	int want = FF_COMPUTE_GRADIENT | FF_COMPUTE_IHESSIAN;
	if (verbose >= 5) want |= FF_COMPUTE_HESSIAN;
	if (iter == 1) {
		want |= FF_COMPUTE_FIT;
	}

	nro.evaluateDerivs(want);

	double speed = std::min(priorSpeed * 1.5, 1.0);
	nro.setSearchDir(searchDir);
	Eigen::Map<Eigen::VectorXd> grad(nro.getGrad(), numParam);
	double targetImprovement = searchDir.dot(grad);

	if (verbose >= 5) {
		nro.debugDeriv(searchDir);
	}

	if (iter == 1) {
		refFit = nro.getFit();
		if (!std::isfinite(refFit)) return;
	}

	if (!std::isfinite(targetImprovement) || targetImprovement < 0) {
		if (verbose >= 4) mxLog("%s: target improvement %.4g is suspect, using steepest descent",
					name, targetImprovement);
		steepestDescent = true;
		// the steepness really provides no information due to unknown curvature
		searchDir = grad / grad.norm();
		if (!std::isfinite(searchDir.norm())) {
			if (verbose >= 2) {
				for (int px=0; px < numParam; ++px) {
					if (std::isfinite(grad[px])) continue;
					mxLog("%s=%f is infeasible; try adding bounds",
					      nro.paramIndexToName(px), nro.getParamVec()[px]);
				}
			}
			nro.reportBadDeriv();
			return;
		}
		targetImprovement = 1;
		speed = std::max(speed, .1);  // expect steepestDescent
	}
	
	// This is based on the Goldstein test. However, we don't enforce
	// a lower bound on the improvement.

	int probeCount = 0;
	Eigen::Map< Eigen::VectorXd > trial(nro.getParamVec(), numParam);
	double bestSpeed = 0;
	double bestImproved = 0;
	double goodness = 0;
	double bestFit = 0;

	while (++probeCount < 16) {
		const double scaledTarget = speed * targetImprovement;
		if (!steepestDescent && relImprovement(scaledTarget) < tolerance) {
			trial = prevEst;
			return;
		}
		trial = (prevEst - speed * searchDir).
			cwiseMax(nro.lbound).cwiseMin(nro.ubound);
		++minorIter;
		nro.evaluateFit();
		if (verbose >= 4) mxLog("%s: speed %.3g for target %.3g fit %f ref %f",
					name, speed, scaledTarget, nro.getFit(), refFit);
		if (!std::isfinite(nro.getFit())) {
			speed *= .1;
			continue;
		}
		const double improved = refFit - nro.getFit();
		if (improved <= 0) {
			double guess = scaledTarget/(scaledTarget-improved);
			if (verbose >= 4) {
				mxLog("%s: improved %.2g (%.2g), suspect excessive speed",
				      name, improved, guess);
			}
			speed *= std::min(0.1, guess);
			continue;
		}
		bestImproved = improved;
		bestSpeed = speed;
		bestFit = nro.getFit();
		goodness = improved / scaledTarget;
		if (verbose >= 3) mxLog("%s: viable speed %.2g for improvement %.3g goodness %f",
					name, bestSpeed, bestImproved, goodness);
		break;
	}
	if (bestSpeed == 0) {
		trial = prevEst;
		return;
	}

	const double epsilon = 0.5;
	const double accelFactor = 0.5;
	if (!steepestDescent && speed < 1/(1+0.5*accelFactor) && goodness < epsilon) {
		if (verbose >= 3) {
			mxLog("%s: goodness %.2f < %.2f, try to incr speed",
			      name, goodness, epsilon);
		}
		int retries = 8;
		while (--retries > 0) {
			speed = std::min(speed * (1 + accelFactor), 1.0);
			++probeCount;
			trial = (prevEst - speed * searchDir).
				cwiseMax(nro.lbound).cwiseMin(nro.ubound);
			++minorIter;
			nro.evaluateFit();
			if (!std::isfinite(nro.getFit())) break;
			const double improved = refFit - nro.getFit();
			if (bestImproved >= improved) break;
			double improvementOverBest = improved - bestImproved;
			if (verbose >= 4) {
				mxLog("%s: [%d] incr speed for fit incremental improvement of %.2g",
				      name, retries, improvementOverBest);
			}
			bestFit = nro.getFit();
			bestImproved = improved;
			bestSpeed = speed;
			if (speed == 1 || relImprovement(improvementOverBest) < tolerance) break;
		}
	}

	if (verbose >= 3) mxLog("%s: %s, probes %d speed %f improved %.3g",
				name, steepestDescent?"steepestDescent":"normal step",
				probeCount, bestSpeed, bestImproved);
	priorSpeed = bestSpeed;

	trial = (prevEst - bestSpeed * searchDir).
		cwiseMax(nro.lbound).cwiseMin(nro.ubound);

	maxAdj = 0;
	for (int px=0; px < numParam; ++px) {
		const double oldEst = prevEst[px];
		double badj = fabs(oldEst - trial(px));
		if (maxAdj < badj) {
			maxAdj = badj;
			maxAdjSigned = oldEst - trial(px);
			maxAdjParam = px;
		}
	}

	improvement = bestImproved;
	refFit = bestFit;
}

class ComputeNR : public omxCompute {
	typedef omxCompute super;

	int maxIter;
	double tolerance;
	int verbose;

public:
	omxMatrix *fitMatrix;
	int numParam;
	const char *engineName;
	virtual void initFromFrontend(omxState *state, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

struct ComputeNRO : public NewtonRaphsonObjective {
	ComputeNR &nr;
	FitContext *fc;
	ComputeNRO(ComputeNR *_nr, FitContext *_fc) : nr(*_nr), fc(_fc) {};
	virtual bool isConverged() {
		nr.reportProgress(fc);
		return converged || isErrorRaised() ||
			fc->getInform() != INFORM_UNINITIALIZED;
	}
	virtual double getFit() { return fc->fit; };
	virtual void resetDerivs() {
		fc->resetOrdinalRelativeError();
		fc->initGrad(nr.numParam);
		fc->clearHessian();
	};
	virtual const char *paramIndexToName(int px)
	{
		const char *pname = "none";
		if (px >= 0) pname = fc->varGroup->vars[px]->name;
		return pname;
	}
	virtual void evaluateFit() {
		fc->copyParamToModel();
		ComputeFit(nr.engineName, nr.fitMatrix, FF_COMPUTE_FIT, fc);
	}
	virtual void evaluateDerivs(int want) {
		fc->copyParamToModel();
		ComputeFit(nr.engineName, nr.fitMatrix, want, fc);
	}
	virtual double *getParamVec() { return fc->est; };
	virtual double *getGrad() { return fc->gradZ.data(); };
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir) {
		searchDir = fc->ihessGradProd();
	}
	virtual void reportBadDeriv() {
		fc->setInform(INFORM_BAD_DERIVATIVES);
	};
	virtual void debugDeriv(const Eigen::Ref<Eigen::VectorXd> searchDir) {
		fc->log(FF_COMPUTE_FIT | FF_COMPUTE_ESTIMATE | FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN);
		std::string buf;
		buf += "searchDir: c(";
		for (int vx=0; vx < searchDir.size(); ++vx) {
			buf += string_snprintf("%.5f", searchDir[vx]);
			if (vx < searchDir.size() - 1) buf += ", ";
		}
		buf += ")\n";
		mxLogBig(buf);
	}
};

class omxCompute *newComputeNewtonRaphson()
{
	return new ComputeNR();
}

void ComputeNR::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, state, "fitfunction");
	omxCompleteFitFunction(fitMatrix);

	if (!fitMatrix->fitFunction->hessianAvailable ||
	    !fitMatrix->fitFunction->gradientAvailable) {
		mxThrow("Newton-Raphson requires derivatives");
	}

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) mxThrow("tolerance must be positive");

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);

	engineName = "NnRn";
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

		ThinMatrix ihessMat(ihess, dim, dim);
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

void ComputeNR::computeImpl(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	numParam = varGroup->vars.size();
	if (numParam <= 0) { complainNoFreeParam(); return; }

	fc->setInform(INFORM_UNINITIALIZED);

	omxAlgebraPreeval(fitMatrix, fc);

	// bounds might have changed
	ComputeNRO obj(this, fc);
	obj.lbound.resize(numParam);
	obj.ubound.resize(numParam);
	for(int px = 0; px < numParam; px++) {
		obj.lbound[px] = varGroup->vars[px]->lbound;
		obj.ubound[px] = varGroup->vars[px]->ubound;
	}

	NewtonRaphsonOptimizer nro(name, maxIter, tolerance, verbose);
	nro(obj);

	fc->iterations += nro.getIter();

	if (obj.converged) {
		double gradNorm = 0.0;
		double feasibilityTolerance = Global->feasibilityTolerance;
		// factor out simliar code in omxHessianCalculation
		for (int gx=0; gx < numParam; ++gx) {
			if ((fc->gradZ[gx] > 0 && fabs(fc->est[gx] - obj.lbound[gx]) < feasibilityTolerance) ||
			    (fc->gradZ[gx] < 0 && fabs(fc->est[gx] - obj.ubound[gx]) < feasibilityTolerance)) continue;
			double g1 = fc->gradZ[gx];
			gradNorm += g1 * g1;
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm > Global->getGradientThreshold(fc->fit)) {
			if (verbose >= 1) mxLog("gradient norm=%f gradient thresh=%f",
						gradNorm, Global->getGradientThreshold(fc->fit));
			fc->setInform(INFORM_NOT_AT_OPTIMUM);
		} else {
			fc->setInform(INFORM_CONVERGED_OPTIMUM);
			fc->wanted |= FF_COMPUTE_BESTFIT;
		}
		if (verbose >= 1) {
			mxLog("%s: converged in %d cycles (%d minor iterations) inform=%d",
			      name, nro.getIter(), nro.getMinorIter(), fc->getInform());
		}
	} else {
		if (nro.getIter() == 1) {
			fc->setInform(INFORM_STARTING_VALUES_INFEASIBLE);
		} else {
			fc->setInform(INFORM_ITERATION_LIMIT);
			if (verbose >= 1) {
				mxLog("%s: failed to converge after %d cycles (%d minor iterations)",
				      name, nro.getIter(), nro.getMinorIter());
			}
		}
	}
}

void ComputeNR::reportResults(FitContext *fc, MxRList *slots, MxRList *output)
{
	omxPopulateFitFunction(fitMatrix, output);
}
