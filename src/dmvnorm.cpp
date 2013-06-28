#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <vector>

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
	std::vector<double> cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	std::vector<double> cov(dim * dim);
	memcpy(cov.data(), origCov, sizeof(double) * dim * dim);

	std::vector<double> icov(dim * dim);
	for (int rx=0; rx < dim; rx++) {
		for (int cx=0; cx < dim; cx++) {
			icov[rx * dim + cx] = rx==cx? 1 : 0;
		}
	}
  
	std::vector<int> ipiv(dim);
	int info;
	F77_CALL(dgesv)(&dim, &dim, cov.data(), &dim, ipiv.data(), icov.data(), &dim, &info);
	if (info < 0) {
		snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
		return nan("error");
	}
	if (info > 0) {
		snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
		return nan("error");
	}

	std::vector<double> half(dim);
	char trans='n';
	double alpha=1;
	double beta=0;
	int inc=1;
	F77_CALL(dgemv)(&trans, &dim, &dim, &alpha, icov.data(), &dim, cloc.data(), &inc, &beta, half.data(), &inc);

	double got=0;
	for (int dx=0; dx < dim; dx++) got += half[dx] * cloc[dx];
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

	std::vector<double> sigma(dim * dim);
	memcpy(sigma.data(), origSigma, sizeof(double) * dim * dim);

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
			 &dim, sigma.data(), &dim,
			 &vunused, &vunused,
			 &iunused, &iunused,
			 &abstol, &m, w,
			 Z, &ldz, isuppz,
			 &optlWork, &lwork,
			 &optliWork, &liwork, &info);
	if (info != 0) {
		snprintf(err, ERROR_LEN, "dsyevr failed when requesting work space size");
		return nan("error");
	}

	lwork = optlWork;
	std::vector<double> work(lwork);
	liwork = optliWork;
	std::vector<int> iwork(liwork);

	F77_CALL(dsyevr)(&jobz, &range, &uplo, &dim, sigma.data(), &dim,
			 &vunused, &vunused, &iunused, &iunused, &abstol, &m, w, Z, &ldz, isuppz,
			 work.data(), &lwork, iwork.data(), &liwork, &info);
	if (info < 0) {
		snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
		return nan("error");
	}
	if (info > 0) {
		snprintf(err, ERROR_LEN, "dsyevr: internal error");
		return nan("error");
	}
	if (m < dim) {
		snprintf(err, ERROR_LEN, "Sigma not of full rank");
		return nan("error");
	}

	for (int dx=0; dx < dim; dx++) dist += log(w[dx]);
	double got = -(dim * M_LN_SQRT_2PI*2 + dist)/2;
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
