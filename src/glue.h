/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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

void string_to_Rf_error( const std::string& str) __attribute__ ((noreturn));

void exception_to_Rf_error( const std::exception& ex ) __attribute__ ((noreturn));

void getMatrixDims(SEXP r_theta, int *rows, int *cols);

SEXP makeFactor(SEXP vec, int levels, const char **labels);
void markAsDataFrame(SEXP list, int rows);
inline void markAsDataFrame(SEXP list) { markAsDataFrame(list, -1); }

SEXP dtmvnorm_marginal(SEXP xn, SEXP n, SEXP sigma, SEXP lower, SEXP upper);
SEXP dtmvnorm_marginal2(SEXP Rxq, SEXP Rxr, SEXP Rq, SEXP Rr,
			SEXP Rsigma, SEXP Rlower, SEXP Rupper);
SEXP mtmvnorm(SEXP sigma, SEXP lower, SEXP upper);

void friendlyStringToLogical(const char *key, SEXP rawValue, int *out);

void loadCharVecFromR(const char *context, SEXP names, std::vector<const char *> &dest);

#ifndef M_LN_2PI
#define M_LN_2PI        1.837877066409345483560659472811        /* log(2*pi) */
#endif

#endif // #define _GLUE_H
