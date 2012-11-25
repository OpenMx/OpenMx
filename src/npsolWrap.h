/*
 *  Copyright 2007-2012 The OpenMx Project
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

#ifndef _NPSOLWRAP_H
#define _NPSOLWRAP_H

#include "omxState.h"

/* Functions for Export */
SEXP callNPSOL(SEXP fitfunction, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList, SEXP expList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state);  // Calls NPSOL.  Duh.

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options);
SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
	SEXP skipMissingness, SEXP skipDefvars);

#endif // #define _NPSOLWRAP_H
