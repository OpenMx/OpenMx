/*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
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

#ifndef u_OMXDEFINES_H_
#define u_OMXDEFINES_H_

#include <Rcpp.h>
using namespace Rcpp;

template <typename... Args>
inline void NORET mxThrow(const char* fmt, Args&&... args) {
    throw std::runtime_error( tfm::format(fmt, std::forward<Args>(args)... ).c_str() );
}

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
//#define EIGEN_DONT_PARALLELIZE
#define  u_OpenMx_Compilation_ 1  // work around bug in Eigen 3.2.5
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
//#define u_GLIBCXX_DEBUG  // but gives link errors without -D_GLIBCXX_DEBUG on command line
#endif // OMX_BOUNDS_CHECK

#ifdef DEBUGMX_ROWS
#define OMX_DEBUG_ROWS(row) ((row < 10) || (row % 50 == 0))
#else
#define OMX_DEBUG_ROWS(row) 0
#endif /* DEBUGMX_ROWS */

#ifdef DEBUGMX_ALGEBRA
#define OMX_DEBUG_ALGEBRA 1
#else
#define OMX_DEBUG_ALGEBRA 0
#endif /* DEBUGMX_ALGEBRA */

#ifdef DEBUGMX_NEWSTUFF
#define OMX_DEBUG_NEWSTUFF 1
#else
#define OMX_DEBUG_NEWSTUFF 0
#endif /* DEBUGMX_NEWSTUFF */

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
	FIT_UNITS_SQUARED_RESIDUAL,
	FIT_UNITS_SQUARED_RESIDUAL_CHISQ,  // full weight matrix
	FIT_UNITS_ANY,
};

#define GRADIENT_FUDGE_FACTOR(x) (pow(10.0,x))

class omxMatrix;
class omxState;
typedef struct omxAlgebra omxAlgebra;
class FitContext;
struct FreeVarGroup;
typedef struct omxFitFunction omxFitFunction;
class omxExpectation;
typedef struct omxDefinitionVar omxDefinitionVar;
typedef struct omxRFitFunction omxRFitFunction;
typedef struct SEXPREC *SEXP;
class MxRList;
class omxCompute;
class FitContext;
class GradientOptimizerContext;
struct ThinMatrix;
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
static inline int omp_get_num_threads(void) { return 1; }
#endif

#include <stan/math/version.hpp>
#if STAN_MATH_MAJOR >= 4
#include <stan/math/prim/fun/Eigen.hpp>
#else
#include <stan/math/prim/mat/fun/Eigen.hpp>
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
//^^^This is what filters out rows and columns of the mean vector and covariance matrix due to missingness, when the
//sufficient-statistics FIML speedup is in use.

// Refactor as a single split function that pulls out all 3 parts
// of the covariance matrix in one iteration through the elements?
// Maybe rename to subsetSymmetric?
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

template <typename T1, typename T2, typename T3>
void subsetVector(const Eigen::MatrixBase<T1> &in, T2 filter, Eigen::MatrixBase<T3> &out)
{
	int ox = 0;
	for (int ix=0; ix < in.size(); ++ix) {
		if (!filter(ix)) continue;
		out[ox++] = in[ix];
	}
}

template <typename T1, typename T2, typename T3>
void subsetVector(const Eigen::ArrayBase<T1> &in, T2 filter, Eigen::ArrayBase<T3> &out)
{
	int ox = 0;
	for (int ix=0; ix < in.size(); ++ix) {
		if (!filter(ix)) continue;
		out[ox++] = in[ix];
	}
}

template <typename T1, typename T3>
void subsetVector(const Eigen::ArrayBase<T1> &in, const std::vector<int> &ind, Eigen::ArrayBase<T3> &out)
{
	for (int ix=0; ix < int(ind.size()); ++ix) {
		out[ix] = in[ ind[ix] ];
	}
}

template <typename T1, typename T3>
void subsetVector(const Eigen::MatrixBase<T1> &in, const std::vector<int> &ind, Eigen::MatrixBase<T3> &out)
{
	for (int ix=0; ix < int(ind.size()); ++ix) {
		out[ix] = in[ ind[ix] ];
	}
}

template <typename T1, typename T3>
void subsetMatrix(const Eigen::ArrayBase<T1> &in, const std::vector<int> &ind, Eigen::ArrayBase<T3> &out)
{
	for (int ix=0; ix < int(ind.size()); ++ix) {
		out.row(ix) = in.row(ind[ix]);
	}
}

template <typename T1, typename T3>
void subsetMatrix(const Eigen::MatrixBase<T1> &in, const std::vector<int> &ind, Eigen::MatrixBase<T3> &out)
{
	for (int ix=0; ix < int(ind.size()); ++ix) {
		out.row(ix) = in.row(ind[ix]);
	}
}

template <typename T1, typename T2>
void subsetVectorStore(Eigen::MatrixBase<T1> &in, T2 filter, double val)
{
	for (int ix=0; ix < in.size(); ++ix) {
		if (!filter(ix)) continue;
		in[ix] = val;
	}
}

template <typename T1, typename T2>
void subsetVectorStore(Eigen::ArrayBase<T1> &in, T2 filter, double val)
{
	for (int ix=0; ix < in.size(); ++ix) {
		if (!filter(ix)) continue;
		in[ix] = val;
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

template <typename T1, typename T2, typename T3>
void subsetColumns(const Eigen::MatrixBase<T1> &in, T2 includeTest,
		   Eigen::MatrixBase<T3> &out)
{
	for (int gcx=0, cx=0; gcx < in.cols(); gcx++) {
		if (!includeTest(gcx)) continue;
		out.col(cx) = in.col(gcx);
		cx += 1;
	}
}

template <typename T1, typename T2, typename T3>
void subsetRows(const Eigen::MatrixBase<T1> &in, T2 includeTest,
		Eigen::MatrixBase<T3> &out)
{
	for (int gcx=0, cx=0; gcx < in.rows(); gcx++) {
		if (!includeTest(gcx)) continue;
		out.row(cx) = in.row(gcx);
		cx += 1;
	}
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void partitionCovariance(const Eigen::MatrixBase<T1> &gcov,
		     T2 filterTest,
		     Eigen::MatrixBase<T3> &v11,
		     Eigen::MatrixBase<T4> &v12,
		     Eigen::MatrixBase<T5> &v22)
{
	for (int gcx=0, c1=0,c2=0,c3=0; gcx < gcov.cols(); gcx++) {
		for (int grx=0, r1=0,r2=0,r3=0; grx < gcov.rows(); grx++) {
			if (filterTest(grx)) {
				if (filterTest(gcx)) {
					v11(r1++, c1) = gcov(grx, gcx);
				} else {
					v12(r2++, c2) = gcov(grx, gcx);
				}
			} else {
				if (!filterTest(gcx)) {
					v22(r3++, c3) = gcov(grx, gcx);
				}
			}
		}
		if (filterTest(gcx)) {
			c1 += 1;
		} else {
			c2 += 1;
			c3 += 1;
		}
	}
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void partitionCovarianceSet(Eigen::MatrixBase<T1> &gcov,
			    T2 filterTest,
			    const Eigen::MatrixBase<T3> &v11,
			    const Eigen::MatrixBase<T4> &v12,
			    const Eigen::MatrixBase<T5> &v22)
{
	for (int gcx=0, c1=0,c2=0,c3=0,c4=0; gcx < gcov.cols(); gcx++) {
		for (int grx=0, r1=0,r2=0,r3=0,r4=0; grx < gcov.rows(); grx++) {
			if (filterTest(grx)) {
				if (filterTest(gcx)) {
					gcov(grx, gcx) = v11(r1++, c1);
				} else {
					gcov(grx, gcx) = v12(r2++, c2);
				}
			} else {
				if (!filterTest(gcx)) {
					gcov(grx, gcx) = v22(r3++, c3);
				} else {
					gcov(grx, gcx) = v12(c4, r4++);
				}
			}
		}
		if (filterTest(gcx)) {
			c1 += 1;
			c4 += 1;
		} else {
			c2 += 1;
			c3 += 1;
		}
	}
}

template<typename u_MatrixType, int u_UpLo = Eigen::Lower>
class SimpCholesky : public Eigen::LDLT<u_MatrixType, u_UpLo> {
 private:
	Eigen::MatrixXd inverse;

 public:
	typedef Eigen::LDLT<u_MatrixType, u_UpLo> Base;

	SimpCholesky() : Base() {};
	template<typename InputType>
	explicit SimpCholesky(const Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};
	template<typename InputType>
	explicit SimpCholesky(Eigen::EigenBase<InputType>& matrix) : Base(matrix) {};

	double log_determinant() const {
		typename Base::Scalar detL = Base::vectorD().array().log().sum();
		return detL * 0.5;
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
		Rf_protect(val);
		SEXP rkey = Rf_mkChar(key);
		Rf_protect(rkey);
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
		if (diff != 1) mxThrow("Depth %d != 1, ScopedProtect was nested", diff);
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
		if (diff != 1) mxThrow("Depth %d != 1, ProtectedSEXP was nested", diff);
		Rf_unprotect(2);
	}
        operator SEXP() const { return var; }
 private:
        ProtectedSEXP( const ProtectedSEXP& );
        ProtectedSEXP& operator=( const ProtectedSEXP& );
};

class ProtectAutoBalanceDoodad {
	PROTECT_INDEX initialpix;
 public:
	ProtectAutoBalanceDoodad() {
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
	~ProtectAutoBalanceDoodad() {
		Rf_unprotect(getDepth());
	}
};

class AssertProtectStackBalanced {
	const char *context;
	int preDepth;
	PROTECT_INDEX initialpix;

	PROTECT_INDEX getDepth() {
		PROTECT_INDEX pix;
		R_ProtectWithIndex(R_NilValue, &pix);
		PROTECT_INDEX diff = pix - initialpix;
		Rf_unprotect(1);
		return diff;
	}
 public:
 	AssertProtectStackBalanced(const char *u_context) : context(u_context) {
		R_ProtectWithIndex(R_NilValue, &initialpix);
		Rf_unprotect(1);
		preDepth = getDepth();
	};
	~AssertProtectStackBalanced() {
		int postDepth = getDepth();
		if (preDepth != postDepth) {
			Rf_warning("%s: "
				   "protect stack usage %d > 0, PLEASE REPORT TO OPENMX DEVELOPERS",
				   context, postDepth - preDepth);
		}
	}
};

#define OOPS mxThrow("%s at %d: oops", __FILE__, __LINE__)

inline int cast_with_NA(double val)
{
  if (std::isfinite(val)) return int(val);
  return NA_INTEGER;
}

inline double cast_with_NA(int val)
{
  if (val == NA_INTEGER) return NA_REAL;
  return double(val);
}

struct cstrCmp {
	bool operator() (const char *s1, const char *s2) const
	{ return strcmp(s1,s2) < 0; }
};

#endif /* u_OMXDEFINES_H_ */
