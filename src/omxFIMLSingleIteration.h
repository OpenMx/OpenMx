/*
 *  Copyright 2007-2014 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#ifndef _OMX_FIML_SINGLE_ITERATION
#define _OMX_FIML_SINGLE_ITERATION

#include "omxFIMLFitFunction.h"

bool omxFIMLSingleIterationJoint(omxFitFunction *localobj,
	omxFitFunction *sharedobj, int rowbegin, int rowcount);

bool omxFIMLSingleIterationOrdinal(omxFitFunction *localobj,
	omxFitFunction *sharedobj, int rowbegin, int rowcount);

bool omxFIMLSingleIteration(omxFitFunction *localobj,
	omxFitFunction *sharedobj, int rowbegin, int rowcount);



#endif // _OMX_FIML_SINGLE_ITERATION
