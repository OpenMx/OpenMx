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

#ifndef _OMX_NPSOL_SPECIFIC_H
#define _OMX_NPSOL_SPECIFIC_H

#include "omxMatrix.h"

/* NPSOL-specific globals */
extern const double NPSOL_BIGBND, NEG_INF, INF;

void omxInvokeNPSOL(omxMatrix *fitMatrix, double *f, double *x, double *g, double *R, int disableOptimizer);
 
void omxNPSOLConfidenceIntervals(omxMatrix *fitMatrix, double optimum, double *optimalValues,
				 int ciMaxIterations);
 
void omxSetNPSOLOpts(SEXP options, int *numHessians, 
    int *calculateStdErrors, int *ciMaxIterations, int *disableOptimizer,
    int *numThreads, int *analyticGradients, int numFreeParams);

#endif // #define _OMX_NPSOL_SPECIFIC_H
