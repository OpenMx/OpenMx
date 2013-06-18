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

#ifndef _NPSOLWRAP_H
#define _NPSOLWRAP_H

#include <exception>
#include <string>

#include "omxState.h"
#include "omxOpenmpWrap.h"

/* Functions for Export */
SEXP omxBackend(SEXP fitfunction, SEXP startVals, SEXP constraints,
		SEXP matList, SEXP varList, SEXP algList, SEXP expList, SEXP computeList,
		SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options);

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options);
SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
	SEXP skipMissingness, SEXP skipDefvars);

/** R is not thread-safe. This lock must be held while manipulating R
 * data structures from a context where there may be multiple threads
 * active.
 */
extern omp_lock_t GlobalRLock;

class omxManageProtectInsanity {
	PROTECT_INDEX initialpix;
 public:
	omxManageProtectInsanity() {
		PROTECT_WITH_INDEX(R_NilValue, &initialpix);
		UNPROTECT(1);
	}
	~omxManageProtectInsanity() {
		PROTECT_INDEX pix;
		PROTECT_WITH_INDEX(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		UNPROTECT(1 + diff);
	}
};

void string_to_try_error( const std::string& str) __attribute__ ((noreturn));

void exception_to_try_error( const std::exception& ex ) __attribute__ ((noreturn));

#endif // #define _NPSOLWRAP_H
