/*
 *  Copyright 2013-2021 by the individuals mentioned in the source code history
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
		mxLog("Welcome to Newton-Raphson (%d param, fitTol %.3g, gradTol %.3g, max iter %d)",
		      (int)numParam, tolerance, gradTolerance, maxIter);
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
			relImprovement(improvement) <= tolerance || iter >= maxIter;

		if (nro.isConverged()) break;
	}

  Eigen::Map<Eigen::VectorXd> grad(nro.getGrad(), numParam);
  if (grad.norm() > gradTolerance) {
    double bestFit = refFit;
    nro.getParamVec(prevEst);
    Eigen::VectorXd bestEst = prevEst;
    const int want = FF_COMPUTE_FIT | FF_COMPUTE_GRADIENT;
    double speed = priorSpeed;
    nro.adjustSpeed(speed);
    while (++iter < maxIter) {
      Eigen::VectorXd prevGrad = grad;
      double priorGradNorm = grad.norm();
      if (verbose >= 3) mxLog("%s: iter %d/%d polish grad=%.3g", name, iter, maxIter, priorGradNorm);
      int probeCount = 0;
      while (++probeCount < lineSearchMax) {
        Eigen::VectorXd trial = (prevEst - speed * prevGrad).cwiseMax(nro.lbound).cwiseMin(nro.ubound);
        if (!trial.array().isFinite().all()) mxThrow("!trial.array().isFinite().all()");
        nro.setParamVec(trial);
        nro.evaluateDerivs(want); // updates grad
        double curGradNorm = grad.norm();
        if (!std::isfinite(curGradNorm)) {
          double oldSpeed = speed;
          speed *= stepMultiplier;
          if (verbose >= 4) {
            mxLog("%s: grad NaN, suspect excessive speed %.3g->%.3g",
                  name, oldSpeed, speed);
          }
          continue;
        }
        // <= important in next line because nro.getFit() == bestFit is likely!
        if (nro.getFit() <= bestFit) { bestFit = nro.getFit(); bestEst = trial; }
        if (curGradNorm == priorGradNorm || curGradNorm < gradTolerance) {
          const double improved = refFit - bestFit;
          if (improved == 0) {
            if (verbose >= 3) mxLog("%s: grad polish failed", name);
          } else {
            if (verbose >= 3) mxLog("%s: grad < tol at speed %.3g grad %.3g, polish improved fit %.3g",
                                    name, speed, curGradNorm, improved);
          }
          iter = maxIter; break;
        } else if (curGradNorm < priorGradNorm) {
          if (verbose >= 3) mxLog("%s: grad speed %.3g grad %.3g", name, speed, grad.norm());
          //speed = std::min(speed / sqrt(stepMultiplier), 1.0); Seems to not help
          prevEst = trial;
          break;
        } else {
          double oldSpeed = speed;
          speed *= stepMultiplier;
          if (verbose >= 4) {
            mxLog("%s: grad %.3g, suspect excessive speed %.3g->%.3g",
                  name, curGradNorm, oldSpeed, speed);
          }
        }
      }
    }
    nro.setParamVec(bestEst);
  }
}

void NewtonRaphsonOptimizer::setStepMultiplier(double sm)
{
  if (sm >= 1 || sm <= 0) mxThrow("NewtonRaphsonOptimizer::setStepMultiplier must be in (0,1)");
  stepMultiplier = sm;
  lineSearchMax = log10(std::numeric_limits<double>::epsilon()) / log10(sm);
}

void NewtonRaphsonOptimizer::lineSearch(NewtonRaphsonObjective &nro)
{
	bool steepestDescent = false;

  nro.getParamVec(prevEst);

	int want = FF_COMPUTE_GRADIENT | FF_COMPUTE_IHESSIAN;
	if (verbose >= 5) want |= FF_COMPUTE_HESSIAN;
	if (iter == 1) {
		want |= FF_COMPUTE_FIT;
	}

	nro.evaluateDerivs(want);

	double speed = std::min(priorSpeed / sqrt(stepMultiplier), 1.0);
	nro.setSearchDir(searchDir);
  if (priorSpeed == 1) nro.adjustSpeed(speed); // only first line search
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
					      nro.paramIndexToName(px), prevEst[px]);
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
	Eigen::VectorXd trial;
	double bestSpeed = 0;
	double bestImproved = 0;
	double goodness = 0;
	double bestFit = 0;

	while (++probeCount < lineSearchMax) {
		const double scaledTarget = speed * targetImprovement;
		if (!steepestDescent && relImprovement(scaledTarget) <= tolerance) {
			nro.setParamVec(prevEst);
			return;
		}
		trial = (prevEst - speed * searchDir).
			cwiseMax(nro.lbound).cwiseMin(nro.ubound);
    nro.setParamVec(trial);
		++minorIter;
		nro.evaluateFit();
		if (verbose >= 4) mxLog("%s: speed %.3g for target %.3g fit %f ref %f",
					name, speed, scaledTarget, nro.getFit(), refFit);
		if (!std::isfinite(nro.getFit())) {
			speed *= stepMultiplier;
			continue;
		}
		const double improved = refFit - nro.getFit();
    if (improved == 0) {
			if (verbose >= 4) {
				mxLog("%s: fit unchanged at speed %.2g, giving up", name, speed);
			}
      return;
    }
		if (improved < 0) {
      // Can't assume curve is remotely symmetric
			//double guess = scaledTarget/(scaledTarget-improved);
			if (verbose >= 4) {
				mxLog("%s: improved %.2g, suspect excessive speed", name, improved);
			}
			speed *= stepMultiplier;
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
		nro.setParamVec(prevEst);
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
      nro.setParamVec(trial);
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
			if (speed == 1 || relImprovement(improvementOverBest) <= tolerance) break;
		}
	}

	if (verbose >= 3) mxLog("%s: %s, probes %d speed %f improved %.3g",
				name, steepestDescent?"steepestDescent":"normal step",
				probeCount, bestSpeed, bestImproved);
	priorSpeed = bestSpeed;

	trial = (prevEst - bestSpeed * searchDir).
		cwiseMax(nro.lbound).cwiseMin(nro.ubound);
  nro.setParamVec(trial);

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
	virtual void initFromFrontend(omxState *state, SEXP rObj) override;
	virtual void computeImpl(FitContext *fc) override;
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out) override;
};

struct ComputeNRO : public NewtonRaphsonObjective {
	ComputeNR &nr;
	FitContext *fc;
	ComputeNRO(ComputeNR *u_nr, FitContext *u_fc) :
    nr(*u_nr), fc(u_fc) {};
	virtual bool isConverged() override {
		nr.reportProgress(fc);
		return converged || isErrorRaised() ||
			fc->getInform() != INFORM_UNINITIALIZED;
	}
	virtual double getFit() override { return fc->getUnscaledFit(); };
	virtual void resetDerivs() override {
		fc->resetOrdinalRelativeError();
		fc->clearHessian();
	};
	virtual const char *paramIndexToName(int px) override
	{
		const char *pname = "none";
		if (px >= 0) pname = fc->varGroup->vars[ fc->freeToParamMap[px] ]->name;
		return pname;
	}
	virtual void evaluateFit() override {
		ComputeFit(nr.engineName, nr.fitMatrix, FF_COMPUTE_FIT, fc);
	}
	virtual void evaluateDerivs(int want) override {
		ComputeFit(nr.engineName, nr.fitMatrix, want, fc);
	}
	virtual void getParamVec(Eigen::Ref<Eigen::VectorXd> out) override {
    fc->copyEstToOptimizer(out);
  };
	virtual void setParamVec(const Eigen::Ref<const Eigen::VectorXd> in) override {
    fc->setEstFromOptimizer(in);
  };
	virtual double *getGrad() override { return fc->gradZ.data(); };
	virtual void setSearchDir(Eigen::Ref<Eigen::VectorXd> searchDir) override {
		searchDir = fc->ihessGradProd();
	}
	virtual void reportBadDeriv() override {
		fc->setInform(INFORM_BAD_DERIVATIVES);
	};
	virtual void debugDeriv(const Eigen::Ref<Eigen::VectorXd> searchDir) override {
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

	if (!fitMatrix->fitFunction->hessianAvailable) {
		mxThrow("Newton-Raphson requires analytic Hessian");
	}

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance < 0) mxThrow("tolerance must be non-negative");

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
	omxAlgebraPreeval(fitMatrix, fc);

	// complain if there are non-linear constraints TODO

	numParam = fc->getNumFree();
	if (numParam <= 0) { complainNoFreeParam(); return; }

	fc->setInform(INFORM_UNINITIALIZED);

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
		if (fc->isGradientTooLarge()) {
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
