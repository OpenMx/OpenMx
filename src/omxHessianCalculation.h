/*
 *  Copyright 2007-2013 The OpenMx Project
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

#ifndef _OMXHESSIAN_CALCULATION_H
#define _OMXHESSIAN_CALCULATION_H

#include "types.h"
#include "omxCompute.h"

void omxEstimateHessian(double functionPrecision, int r);

class omxComputeEstimateHessian : public omxCompute {
	const double stepSize;
	const int numIter;

	double *hessian;
	double *gradient;
	double *stdError;

	SEXP calculatedHessian;
	SEXP stdErrors;

	void doHessianCalculation(int numParams, int numChildren, 
				  struct hess_struct *hess_work, omxState* parentState);

 public:
	omxComputeEstimateHessian();
	virtual ~omxComputeEstimateHessian();

        virtual void initFromFrontend(SEXP rObj) {};
        virtual void setStartValues(SEXP startVals) {};
        virtual void compute();
        virtual void saveState() {};
        virtual void reportResults(MxRList *out);
};

class omxCompute *newComputeEstimateHessian();

#endif // _OMXHESSIAN_CALCULATION_H
