#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
  double *half = NULL;
  double *cloc = Realloc(NULL, dim, double);
  for (int dx=0; dx < dim; dx++) {
    cloc[dx] = loc[dx] - center[dx];
  }

  double *cov = Realloc(NULL, dim * dim, double);
  memcpy(cov, origCov, sizeof(double) * dim * dim);

  double *icov = Realloc(NULL, dim * dim, double);
  for (int rx=0; rx < dim; rx++) {
    for (int cx=0; cx < dim; cx++) {
      icov[rx * dim + cx] = rx==cx? 1 : 0;
    }
  }
  
  int *ipiv = Realloc(NULL, dim, int);
  int info;
  F77_CALL(dgesv)(&dim, &dim, cov, &dim, ipiv, icov, &dim, &info);
  if (info < 0) {
    snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
    goto errout;
  }
  if (info > 0) {
    snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
    goto errout;
  }

  half = Realloc(NULL, dim, double);
  char trans='n';
  double alpha=1;
  double beta=0;
  int inc=1;
  F77_CALL(dgemv)(&trans, &dim, &dim, &alpha, icov, &dim, cloc, &inc, &beta, half, &inc);

  double got=0;
  for (int dx=0; dx < dim; dx++) got += half[dx] * cloc[dx];

 errout:
  Free(half);
  Free(cloc);
  Free(cov);
  Free(icov);
  Free(ipiv);
  return got;
}

static double
mahalanobis(int dim, double *loc, double *center, double *origCov)
{
  char err[ERROR_LEN];
  err[0] = 0;
  double ret = _mahalanobis(err, dim, loc, center, origCov);
  if (err[0]) error("%s", err);
  return ret;
}

static double
_dmvnorm(char *err, int dim, double *loc, double *mean, double *origSigma)
{
  double dist = mahalanobis(dim, loc, mean, origSigma);

  double *sigma = Realloc(NULL, dim * dim, double);
  memcpy(sigma, origSigma, sizeof(double) * dim * dim);
  double *work = NULL;
  int *iwork = NULL;

  char jobz = 'N';
  char range = 'A';
  char uplo = 'U';
  double vunused;
  int iunused;
  double abstol = 0;
  int m;
  double w[dim];
  double Z[dim];
  int ldz=1;
  int isuppz[2*dim];
  int lwork = -1;
  double optlWork;
  int optliWork;
  int liwork = -1;
  int info;

  F77_CALL(dsyevr)(&jobz, &range, &uplo,
		   &dim, sigma, &dim,
		   &vunused, &vunused,
		   &iunused, &iunused,
		   &abstol, &m, w,
		   Z, &ldz, isuppz,
		   &optlWork, &lwork,
		   &optliWork, &liwork, &info);
  if (info != 0) {
    snprintf(err, ERROR_LEN, "dsyevr failed when requesting work space size");
    goto errout;
  }

  lwork = optlWork;
  work = Realloc(NULL, lwork, double);
  liwork = optliWork;
  iwork = Realloc(NULL, liwork, int);

  F77_CALL(dsyevr)(&jobz, &range, &uplo, &dim, sigma, &dim,
		   &vunused, &vunused, &iunused, &iunused, &abstol, &m, w, Z, &ldz, isuppz,
		   work, &lwork, iwork, &liwork, &info);
  if (info < 0) {
    snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
    goto errout;
  }
  if (info > 0) {
    snprintf(err, ERROR_LEN, "dsyevr: internal error");
    goto errout;
  }
  if (m < dim) {
    snprintf(err, ERROR_LEN, "Sigma not of full rank");
    goto errout;
  }

  for (int dx=0; dx < dim; dx++) dist += log(w[dx]);
  double got = -(dim * M_LN_SQRT_2PI*2 + dist)/2;

 errout:
  Free(sigma);
  Free(work);
  Free(iwork);
  return got;
}

double
dmvnorm(int dim, double *loc, double *mean, double *sigma)
{
  char err[ERROR_LEN];
  err[0] = 0;
  double ret = _dmvnorm(err, dim, loc, mean, sigma);
  if (err[0]) error("%s", err);
  return ret;
}
