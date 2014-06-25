/*
 *  Copyright 2007-2014 The OpenMx Project
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
SEXP omxBackend(SEXP constraints, SEXP matList,
		SEXP varList, SEXP algList, SEXP expectList, SEXP computeList,
		SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options);

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options);
SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
	SEXP skipMissingness, SEXP skipDefvars);

class omxManageProtectInsanity {
	PROTECT_INDEX initialpix;
 public:
	omxManageProtectInsanity() {
		PROTECT_WITH_INDEX(R_NilValue, &initialpix);
		UNPROTECT(1);
	}
	PROTECT_INDEX getDepth() {
		PROTECT_INDEX pix;
		PROTECT_WITH_INDEX(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		UNPROTECT(1);
		return diff;
	}
	~omxManageProtectInsanity() {
		UNPROTECT(getDepth());
	}
};

typedef std::vector< std::pair<const char *, SEXP> > MxRListBase;
class MxRList : private MxRListBase {
 public:
	size_t size() const { return MxRListBase::size(); }
	SEXP asR();
	void add(const char *key, SEXP val) {
		Rf_protect(val);
		push_back(std::make_pair(key, val));
	};
};

void string_to_try_Rf_error( const std::string& str) __attribute__ ((noreturn));

void exception_to_try_Rf_error( const std::exception& ex ) __attribute__ ((noreturn));

void getMatrixDims(SEXP r_theta, int *rows, int *cols);

#endif // #define _NPSOLWRAP_H
