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

#include "types.h"

/* NPSOL-specific globals */
extern const double NPSOL_BIGBND;

void omxInvokeNPSOL(omxMatrix *fitMatrix, FitContext *fc,
		    int *inform_out, int *iter_out);
 
void omxNPSOLConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc);
 
void omxSetNPSOLOpts(SEXP options, int *ciMaxIterations, int *numThreads,
		     int *analyticGradients);

#endif // #define _OMX_NPSOL_SPECIFIC_H
