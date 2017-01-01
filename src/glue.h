/*
 *  Copyright 2007-2017 The OpenMx Project
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

#ifndef _GLUE_H
#define _GLUE_H

#include <exception>
#include <string>

#include "omxState.h"

// Can uncomment these to ensure that they are not used in the OpenMx
// source code.
//
//#undef PROTECT_WITH_INDEX
//#undef UNPROTECT

class omxManageProtectInsanity {
	PROTECT_INDEX initialpix;
 public:
	omxManageProtectInsanity() {
		R_ProtectWithIndex(R_NilValue, &initialpix);
		Rf_unprotect(1);
	}
	PROTECT_INDEX getDepth() {
		PROTECT_INDEX pix;
		R_ProtectWithIndex(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		Rf_unprotect(1);
		return diff;
	}
	~omxManageProtectInsanity() {
		Rf_unprotect(getDepth());
	}
};

typedef std::vector< std::pair<SEXP, SEXP> > MxRListBase;
class MxRList : private MxRListBase {
 public:
	size_t size() const { return MxRListBase::size(); }
	SEXP asR();
	void add(const char *key, SEXP val) {
		SEXP rkey = Rf_mkChar(key);
		Rf_protect(rkey);
		Rf_protect(val);
		push_back(std::make_pair(rkey, val));
	};
};

class ScopedProtect { // DEPRECATED, use ProtectedSEXP
	PROTECT_INDEX initialpix;
 public:
	ScopedProtect(SEXP &var, SEXP src) {
		R_ProtectWithIndex(R_NilValue, &initialpix);
		Rf_unprotect(1);
		Rf_protect(src);
		var = src;
	}
	~ScopedProtect() {
		PROTECT_INDEX pix;
		R_ProtectWithIndex(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		if (diff != 1) Rf_error("Depth %d != 1, ScopedProtect was nested", diff);
		Rf_unprotect(2);
	}
};

class ProtectedSEXP {
	PROTECT_INDEX initialpix;
	SEXP var;
 public:
	ProtectedSEXP(SEXP src) {
		R_ProtectWithIndex(R_NilValue, &initialpix);
		Rf_unprotect(1);
		Rf_protect(src);
		var = src;
	}
	~ProtectedSEXP() {
		PROTECT_INDEX pix;
		R_ProtectWithIndex(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		if (diff != 1) Rf_error("Depth %d != 1, ProtectedSEXP was nested", diff);
		Rf_unprotect(2);
	}
        operator SEXP() const { return var; }
 private:
        ProtectedSEXP( const ProtectedSEXP& );
        ProtectedSEXP& operator=( const ProtectedSEXP& );
};

void string_to_try_Rf_error( const std::string& str) __attribute__ ((noreturn));

void exception_to_try_Rf_error( const std::exception& ex ) __attribute__ ((noreturn));

void getMatrixDims(SEXP r_theta, int *rows, int *cols);

SEXP makeFactor(SEXP vec, int levels, const char **labels);
void markAsDataFrame(SEXP list, int rows);
inline void markAsDataFrame(SEXP list) { markAsDataFrame(list, -1); }

SEXP dtmvnorm_marginal(SEXP xn, SEXP n, SEXP sigma, SEXP lower, SEXP upper);
SEXP dtmvnorm_marginal2(SEXP Rxq, SEXP Rxr, SEXP Rq, SEXP Rr,
			SEXP Rsigma, SEXP Rlower, SEXP Rupper);
SEXP mtmvnorm(SEXP sigma, SEXP lower, SEXP upper);

#ifndef M_LN_2PI
#define M_LN_2PI        1.837877066409345483560659472811        /* log(2*pi) */
#endif

#endif // #define _GLUE_H
