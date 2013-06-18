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

#ifndef _OMX_COMPUTE_H_
#define _OMX_COMPUTE_H_

#include <R.h>
#include <Rinternals.h>

#include "types.h"

class omxCompute {
 public:
        virtual void initFromFrontend(SEXP rObj) = 0;
        virtual void compute(double *startVals) = 0;
        virtual void reportResults(MxRList *out) = 0;
	virtual double getOptimizerStatus() { return NA_REAL; }  // backward compatibility
	virtual double getFit() = 0;
	virtual double *getEstimate() = 0;
        virtual ~omxCompute() {}
};

class omxCompute *omxNewCompute(omxState* os, const char *type);

class omxCompute *newComputeGradientDescent();
class omxCompute *newComputeEstimatedHessian();

// hide impl TODO
class omxComputeGD : public omxCompute {
	omxMatrix *fitMatrix;

	int numFree;
	SEXP minimum, estimate, gradient, hessian;
	SEXP intervals, intervalCodes;
	int inform, iter;

public:
	virtual void initFromFrontend(SEXP rObj);
	virtual void compute(double *startVals);
	virtual void reportResults(MxRList *out);
	virtual double getOptimizerStatus() { return inform; }  // backward compatibility
	virtual double getFit() { return REAL(minimum)[0]; }
	virtual double *getEstimate() { return REAL(estimate); }

	// remove TODO
	omxComputeGD() { init(); }
	void init();
};

#endif
