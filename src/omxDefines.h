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
#include <string>

#define R_NO_REMAP
#include <Rcpp.h>
#include <Rmath.h>
#include <Rinternals.h>

typedef uint64_t nanotime_t;
nanotime_t get_nanotime(void);

#define MIN_ROWS_PER_THREAD 8

#define OMXINLINE inline

#define OMXZERO(ptr, num)	memset(ptr, 0, sizeof(*ptr)*num)

static inline bool strEQ(const char *s1, const char *s2) { return strcmp(s1,s2)==0; }

#define OMX_STATIC_ARRAY_SIZE(ar) (sizeof(ar)/sizeof(ar[0]))

#ifdef EIGEN_WORLD_VERSION
#error "omxDefines.h must be included before Eigen"
#endif

#define EIGEN_NO_DEBUG 1
#define EIGEN_DONT_PARALLELIZE
#define  _OpenMx_Compilation_ 1  // work around bug in Eigen 3.2.5
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int   // default is 8 but 4 bytes is plenty for us

#ifdef OMX_BOUNDS_CHECK
#ifdef NDEBUG
#error "Undefine NDEBUG when using OMX_BOUNDS_CHECK"
#endif
#endif

#ifdef DEBUGMX
#ifdef NDEBUG
#error "Undefine NDEBUG when debugging"
#endif // NDEBUG
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
	// In preoptimize, we determine which fit functions are
	// actually in use so we can made a decision about the best
	// multithreading strategy.
	FF_COMPUTE_PREOPTIMIZE  = 1<<1,
	FF_COMPUTE_MAXABSCHANGE = 1<<2,
	FF_COMPUTE_FIT          = 1<<3,
	FF_COMPUTE_ESTIMATE     = 1<<4,
	FF_COMPUTE_GRADIENT     = 1<<5,
	FF_COMPUTE_HESSIAN      = 1<<6,
	FF_COMPUTE_IHESSIAN     = 1<<7,
	FF_COMPUTE_DERIV        = FF_COMPUTE_GRADIENT | FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN,

	// Use this to obtain a Hessian or Inverse Hessian evaluated at the MLE.
	// Check FitContext::wanted to see which one you got. It may be
	// more efficient to compute one or the other depending on the
	// estimation method. The information matrix is -1 * Hessian.

	FF_COMPUTE_INFO         = 1<<8,   // Fisher information
	FF_COMPUTE_BESTFIT      = 1<<9,
	FF_COMPUTE_STARTING     = 1<<10,   // for special hacks, not for routine use
	FF_COMPUTE_INITIAL_FIT  = 1<<11,    // omxInitialMatrixAlgebraCompute
	FF_COMPUTE_FINAL_FIT    = 1<<12
};

enum FitStatisticUnits {
	FIT_UNITS_UNINITIALIZED=0,
	FIT_UNITS_UNKNOWN,
	FIT_UNITS_PROBABILITY,
	FIT_UNITS_MINUS2LL,
	FIT_UNITS_SQUARED_RESIDUAL  // OK?
};

#define GRADIENT_FUDGE_FACTOR(x) (pow(10.0,x))

typedef struct omxMatrix omxMatrix;
typedef struct omxState omxState;
typedef struct omxAlgebra omxAlgebra;
class FitContext;
struct FreeVarGroup;
typedef struct omxFitFunction omxFitFunction;
typedef struct omxExpectation omxExpectation;
typedef struct omxDefinitionVar omxDefinitionVar;
typedef struct omxRFitFunction omxRFitFunction;
typedef struct SEXPREC *SEXP;
class MxRList;
class omxCompute;
class FitContext;
class GradientOptimizerContext;
struct Matrix;
struct Param_Obj;

// debug tools
void pda(const double *ar, int rows, int cols);
void pia(const int *ar, int rows, int cols);
std::string string_snprintf(const char *fmt, ...) __attribute__((format (printf, 1, 2)));
void mxLog(const char* msg, ...) __attribute__((format (printf, 1, 2)));   // thread-safe
void mxLogSetCurrentRow(int row);
void mxLogBig(const std::string &str);

static inline int triangleLoc1(int diag)
{
	return (diag) * (diag+1) / 2;   // 0 1 3 6 10 15 ..
}

static inline int triangleLoc0(int diag)
{
	return triangleLoc1(diag+1) - 1;  // 0 2 5 9 14 ..
}

static inline bool doubleEQ(double d1, double d2)
{
	// memcmp is required here because NaN != NaN always
	// but careful because there is more than 1 bit pattern for NaN
	return memcmp(&d1, &d2, sizeof(double)) == 0;
}

#ifdef _OPENMP

#include <omp.h>

#if _OPENMP <= 200505

static OMXINLINE int omx_absolute_thread_num(void) {
   return(omp_get_thread_num());
}

#else  // _OPENMP <= 200505

static OMXINLINE int omx_absolute_thread_num(void) {
   int retval = 0;
   int level = omp_get_level();
   int scale = 1;
   for(int i = level; i > 0; i--) {
       retval += scale * omp_get_ancestor_thread_num(i);
       scale *= omp_get_team_size(i);
   }
   return(retval);
}

#endif // _OPENMP <= 200505

static OMXINLINE void omx_omp_init_lock(omp_lock_t* lock) {
   omp_init_lock(lock);
}

static OMXINLINE void omx_omp_destroy_lock(omp_lock_t* lock) {
   omp_destroy_lock(lock);
}

static OMXINLINE void omx_omp_set_lock(omp_lock_t* lock) {
   omp_set_lock(lock);
}

static OMXINLINE void omx_omp_unset_lock(omp_lock_t* lock) {
   omp_unset_lock(lock);
}

#else  // _OPENMP

typedef int omp_lock_t;

static OMXINLINE int omx_absolute_thread_num(void) {
   return(0);
}

static OMXINLINE void omx_omp_init_lock(omp_lock_t* __attribute__((unused)) lock) {}

static OMXINLINE void omx_omp_destroy_lock(omp_lock_t* __attribute__((unused)) lock) {}

static OMXINLINE void omx_omp_set_lock(omp_lock_t* __attribute__((unused)) lock) {}

static OMXINLINE void omx_omp_unset_lock(omp_lock_t* __attribute__((unused)) lock) {}


#endif // #ifdef _OPENMP


#ifndef _OPENMP
static inline int omp_get_thread_num() { return 0; }
#endif

#include <Eigen/Core>

// Refactor as a single split function that pulls out all 3 parts
// of the covariance matrix in one iteration through the elements?
template <typename T1, typename T2, typename T3, typename T4, typename T5>
void subsetNormalDist(const Eigen::MatrixBase<T1> &gmean, const Eigen::MatrixBase<T2> &gcov,
		      T5 includeTest, int resultSize,
		      Eigen::MatrixBase<T3> &mean, Eigen::MatrixBase<T4> &cov)
{
	mean.derived().resize(resultSize);
	cov.derived().resize(resultSize, resultSize);

	for (int gcx=0, cx=0; gcx < gcov.cols(); gcx++) {
		if (!includeTest(gcx)) continue;
		mean[cx] = gmean[gcx];
		for (int grx=0, rx=0; grx < gcov.rows(); grx++) {
			if (!includeTest(grx)) continue;
			cov(rx,cx) = gcov(grx, gcx);
			rx += 1;
		}
		cx += 1;
	}
}

// Refactor as a single split function that pulls out all 3 parts
// of the covariance matrix in one iteration through the elements?
template <typename T2, typename T4, typename T5>
void subsetCovariance(const Eigen::MatrixBase<T2> &gcov,
		      T5 includeTest, int resultSize,
		      Eigen::MatrixBase<T4> &cov)
{
	cov.derived().resize(resultSize, resultSize); // can avoid reallocation? TODO

	for (int gcx=0, cx=0; gcx < gcov.cols(); gcx++) {
		if (!includeTest(gcx)) continue;
		for (int grx=0, rx=0; grx < gcov.rows(); grx++) {
			if (!includeTest(grx)) continue;
			cov(rx,cx) = gcov(grx, gcx);
			rx += 1;
		}
		cx += 1;
	}
}

template <typename T1, typename T2, typename T3>
void subsetVector(const Eigen::MatrixBase<T1> &gmean, T2 includeTest,
		  int resultSize, Eigen::MatrixBase<T3> &out)
{
	out.derived().resize(resultSize); // can avoid reallocation? TODO

	for (int gcx=0, cx=0; gcx < gmean.size(); gcx++) {
		if (!includeTest(gcx)) continue;
		out[cx] = gmean[gcx];
		cx += 1;
	}
}

template <typename T2, typename T4, typename T5>
void subsetCovarianceStore(Eigen::MatrixBase<T2> &gcov,
		      T5 includeTest, const Eigen::MatrixBase<T4> &cov)
{
	for (int gcx=0, cx=0; gcx < gcov.cols(); gcx++) {
		if (!includeTest(gcx)) continue;
		for (int grx=0, rx=0; grx < gcov.rows(); grx++) {
			if (!includeTest(grx)) continue;
			gcov(grx, gcx) = cov(rx,cx);
			rx += 1;
		}
		cx += 1;
	}
}

template<typename _MatrixType, int _UpLo = Eigen::Lower>
class SimpCholesky : public Eigen::LDLT<_MatrixType, _UpLo> {
 private:
	Eigen::MatrixXd inverse;

 public:
	typedef Eigen::LDLT<_MatrixType, _UpLo> Base;

	double log_determinant() const {
		typename Base::Scalar detL = Base::vectorD().array().log().sum();
		return detL;
	}

	void refreshInverse()
	{
		inverse.setIdentity(Base::m_matrix.rows(), Base::m_matrix.rows());
		inverse = Base::solve(inverse);
	};

	const Eigen::MatrixXd &getInverse() const { return inverse; };
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

#endif /* _OMXDEFINES_H_ */
