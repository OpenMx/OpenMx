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
 *
 */

/***********************************************************
 * 
 *  omxDefines.h
 *
 *  Created: Timothy R. Brick 	Date: 2009-09-23
 *
 *	Contains #define information for debugging purposes.
 *
 **********************************************************/
#ifndef _OMXDEFINES_H_
#define _OMXDEFINES_H_

#include <string.h>

#define MIN_ROWS_PER_THREAD 8
#define EPSILON 1e-16
#define TRUE 1
#define FALSE 0

#define OMXINLINE inline

#define OMXZERO(ptr, num)	memset(ptr, 0, sizeof(*ptr)*num)

static inline bool strEQ(const char *s1, const char *s2) { return strcmp(s1,s2)==0; }

#define OMX_STATIC_ARRAY_SIZE(ar) (sizeof(ar)/sizeof(ar[0]))

#ifdef EIGEN_WORLD_VERSION
#error "omxDefines.h must be included before Eigen"
#endif

#define EIGEN_NO_DEBUG 1
#define EIGEN_DONT_PARALLELIZE

#ifdef DEBUGMX
//#ifdef NDEBUG
//#error "Undefine NDEBUG when debugging"
//#endif // NDEBUG
#define OMX_DEBUG 1
#define OMX_BOUNDS_CHECK 1
#else
#define OMX_DEBUG 0
#endif /* DEBUGMX */

#ifdef OMX_BOUNDS_CHECK
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#undef EIGEN_NO_DEBUG
//#define _GLIBCXX_DEBUG  // but gives link errors without -D_GLIBCXX_DEBUG on command line
#endif // OMX_BOUNDS_CHECK

#ifdef DEBUGMX_ROWS
#define OMX_DEBUG_ROWS(row) ((row < 10) || (row % 50 == 0))
#else
#define OMX_DEBUG_ROWS(row) 0
#endif /* DEBUGMX_ROWS */

#ifdef DEBUGMX_MATRIX
#define OMX_DEBUG_MATRIX 1
#else
#define OMX_DEBUG_MATRIX 0
#endif /* DEBUGMX_MATRIX */

#ifdef DEBUGMX_ALGEBRA
#define OMX_DEBUG_ALGEBRA 1
#else
#define OMX_DEBUG_ALGEBRA 0
#endif /* DEBUGMX_ALGEBRA */

// Put forward type declarations here

#include <vector>

enum omxFFCompute {
	FF_COMPUTE_PARAMFLAVOR  = 1<<0,
	FF_COMPUTE_PREOPTIMIZE  = 1<<1,
	FF_COMPUTE_MAXABSCHANGE = 1<<2,
	FF_COMPUTE_FIT          = 1<<3,
	FF_COMPUTE_ESTIMATE     = 1<<4,
	FF_COMPUTE_GRADIENT     = 1<<5,
	FF_COMPUTE_HESSIAN      = 1<<6,
	FF_COMPUTE_IHESSIAN     = 1<<7,

	// Use this to obtain a Hessian or Inverse Hessian evaluated at the MLE.
	// Check FitContext::wanted to see which one you got. It may be
	// more efficient to compute one or the other depending on the
	// estimation method. The information matrix is -1 * Hessian.

	FF_COMPUTE_INFO         = 1<<8,   // Fisher information
	FF_COMPUTE_BESTFIT      = 1<<9,
	FF_COMPUTE_STARTING     = 1<<10,   // for special hacks, not for routine use
	FF_COMPUTE_INITIAL_FIT  = 1<<11,   // omxInitialMatrixAlgebraCompute
	FF_COMPUTE_DIMS         = 1<<12
};

typedef struct omxMatrix omxMatrix;
typedef struct omxState omxState;
class FitContext;
struct FreeVarGroup;
typedef struct omxFitFunction omxFitFunction;
typedef struct omxExpectation omxExpectation;
typedef struct omxDefinitionVar omxDefinitionVar;
typedef struct omxRFitFunction omxRFitFunction;
typedef struct SEXPREC *SEXP;
class MxRList;
class omxCompute;
struct Matrix;
struct Param_Obj;

// debug tools
void pda(const double *ar, int rows, int cols);
void pia(const int *ar, int rows, int cols);

static inline int triangleLoc1(int diag)
{
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

static inline int triangleLoc0(int diag)
{
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

#endif /* _OMXDEFINES_H_ */
