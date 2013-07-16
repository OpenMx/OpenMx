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

#ifndef _OPENMX_TYPES_H_
#define _OPENMX_TYPES_H_

// Put forward type declarations here

#include <vector>

enum omxFFCompute {
	FF_COMPUTE_FIT          = 1<<0,
	FF_COMPUTE_PREOPTIMIZE  = 1<<1,
	FF_COMPUTE_ESTIMATE     = 1<<2,
	FF_COMPUTE_GRADIENT     = 1<<3,
	FF_COMPUTE_HESSIAN      = 1<<4,
	FF_COMPUTE_POSTOPTIMIZE = 1<<5
};

typedef struct omxMatrix omxMatrix;
typedef struct omxState omxState;
class FitContext;
struct FreeVarGroup;
typedef struct omxFitFunction omxFitFunction;
typedef struct omxExpectation omxExpectation;
typedef struct omxDefinitionVar omxDefinitionVar;
typedef struct omxRFitFunction omxRFitFunction;
typedef std::vector< std::pair<SEXP, SEXP> > MxRList;
class omxCompute;

#endif
