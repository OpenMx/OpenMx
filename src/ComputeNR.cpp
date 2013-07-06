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

#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"

class ComputeNR : public omxComputeOperation {
	typedef omxComputeOperation super;
	omxMatrix *fitMatrix;

	int maxIter;
	double tolerance;
	int inform, iter;

public:
	ComputeNR();
	virtual void initFromFrontend(SEXP rObj);
	virtual void compute(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *out);
	virtual double getOptimizerStatus() { return inform; }  // backward compatibility
};

class omxCompute *newComputeNewtonRaphson()
{
	return new ComputeNR();
}

ComputeNR::ComputeNR()
{
	inform = 0;
	iter = 0;
}

void ComputeNR::initFromFrontend(SEXP rObj)
{
	super::initFromFrontend(rObj);

	fitMatrix = omxNewMatrixFromSlot(rObj, globalState, "fitfunction");
	setFreeVarGroup(fitMatrix->fitFunction, varGroup);
	omxCompleteFitFunction(fitMatrix);

	if (!fitMatrix->fitFunction->hessianAvailable ||
	    !fitMatrix->fitFunction->gradientAvailable) {
		error("Newton-Raphson requires derivatives");
	}

	SEXP slotValue;
	PROTECT(slotValue = GET_SLOT(rObj, install("maxIter")));
	maxIter = INTEGER(slotValue)[0];

	PROTECT(slotValue = GET_SLOT(rObj, install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) error("tolerance must be positive");
}

void ComputeNR::compute(FitContext *fc)
{
	// complain if there are non-linear constraints TODO

	size_t numParam = varGroup->vars.size();
	if (numParam <= 0) {
		error("Model has no free parameters");
		return;
	}

	omxFitFunctionCompute(fitMatrix->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);

	iter = 0;
	double prevLL = nan("unset");
	bool decreasing = TRUE;

	while (1) {
		const int want = FF_COMPUTE_FIT|FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN;

		OMXZERO(fc->grad, numParam);
		OMXZERO(fc->hess, numParam * numParam);

		omxFitFunctionCompute(fitMatrix->fitFunction, want, fc);

		// Only need LL for diagnostics; Can avoid computing it? TODO
		double LL = fitMatrix->data[0];
		if (isfinite(prevLL) && prevLL < LL - tolerance) decreasing = FALSE;
		prevLL = LL;

		fc->fixHessianSymmetry();
		//		fc->log(FF_COMPUTE_ESTIMATE|FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN);

		std::vector<double> ihess(numParam * numParam);
		for (size_t rx=0; rx < numParam; rx++) {
			for (size_t cx=0; cx < numParam; cx++) {
				ihess[rx * numParam + cx] = rx==cx? 1 : 0;
			}
		}
		std::vector<int> ipiv(numParam);
		int info;
		int dim = int(numParam);
		F77_CALL(dgesv)(&dim, &dim, fc->hess, &dim, ipiv.data(),
				ihess.data(), &dim, &info);
		if (info < 0) error("Arg %d is invalid", -info);
		if (info > 0) {
			omxRaiseErrorf(globalState, "Hessian is singular and cannot be inverted");
			break;
		}

		std::vector<double> adj(numParam);
		char trans = 'N';
		double alpha = -1;
		int incx = 1;
		double beta = 0;
		F77_CALL(dgemv)(&trans, &dim, &dim, &alpha, ihess.data(), &dim,
				fc->grad, &incx, &beta, adj.data(), &incx);

		double maxAdj = 0;
		for (size_t px=0; px < numParam; ++px) {
			double param = fc->est[px];
			param += adj[px];
			omxFreeVar *fv = fc->varGroup->vars[px];
			if (param < fv->lbound) param = fv->lbound;
			if (param > fv->ubound) param = fv->ubound;
			double adj = fabs(param - fc->est[px]);
			if (maxAdj < adj)
				maxAdj = adj;
			fc->est[px] = param;
		}
		fc->copyParamToModel(globalState);
		R_CheckUserInterrupt();
		if (maxAdj < tolerance || ++iter > maxIter) break;
	}

	// The check is too dependent on numerical precision to enable by default.
	// Anyway, it's just a tool for developers.
	if (0 && !decreasing) warning("Newton-Raphson iterations did not converge");
}

void ComputeNR::reportResults(FitContext *fc, MxRList *out)
{
	if (Global->numIntervals) {
		warning("Confidence intervals are not implemented for Newton-Raphson");
	}  

	omxPopulateFitFunction(fitMatrix, out);

	size_t numFree = varGroup->vars.size();

	SEXP minimum, estimate;
	PROTECT(minimum = NEW_NUMERIC(1));
	PROTECT(estimate = allocVector(REALSXP, numFree));

	REAL(minimum)[0] = fc->fit;
	memcpy(REAL(estimate), fc->est, sizeof(double) * numFree);

	out->push_back(std::make_pair(mkChar("minimum"), minimum));
	out->push_back(std::make_pair(mkChar("estimate"), estimate));

	SEXP code, iterations;

	PROTECT(code = NEW_NUMERIC(1));
	REAL(code)[0] = inform;
	out->push_back(std::make_pair(mkChar("nr.code"), code));

	PROTECT(iterations = NEW_NUMERIC(1));
	REAL(iterations)[0] = iter;
	out->push_back(std::make_pair(mkChar("nr.iterations"), iterations));
}
