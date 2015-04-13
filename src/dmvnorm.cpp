#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "omxDefines.h"
#include <Eigen/Core>
#include "omxBuffer.h"
#include "matrix.h"
#include "glue.h"

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
	std::vector<double> cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	Matrix covMat(origCov, dim, dim);
	omxBuffer<double> icov(dim * dim);
	Matrix icovMat(icov.data(), dim, dim);
	int info = MatrixSolve(covMat, icovMat, true); // can optimize for symmetry TODO
	if (info) {
		snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
		return nan("Rf_error");
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
	if (err[0]) Rf_error("%s", err);
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
	Eigen::VectorXd w(dim);
	Eigen::VectorXd Z(dim);
	int ldz=1;
	Eigen::VectorXi isuppz(2*dim);
	int lwork = -1;
	double optlWork;
	int optliWork;
	int liwork = -1;
	int info;

	F77_CALL(dsyevr)(&jobz, &range, &uplo,
			 &dim, sigma.data(), &dim,
			 &vunused, &vunused,
			 &iunused, &iunused,
			 &abstol, &m, w.data(),
			 Z.data(), &ldz, isuppz.data(),
			 &optlWork, &lwork,
			 &optliWork, &liwork, &info);
	if (info != 0) {
		snprintf(err, ERROR_LEN, "dsyevr failed when requesting work space size");
		return nan("Rf_error");
	}

	lwork = optlWork;
	std::vector<double> work(lwork);
	liwork = optliWork;
	std::vector<int> iwork(liwork);

	F77_CALL(dsyevr)(&jobz, &range, &uplo, &dim, sigma.data(), &dim,
			 &vunused, &vunused, &iunused, &iunused, &abstol, &m, w.data(), Z.data(), &ldz, isuppz.data(),
			 work.data(), &lwork, iwork.data(), &liwork, &info);
	if (info < 0) {
		snprintf(err, ERROR_LEN, "Arg %d is invalid", -info);
		return nan("Rf_error");
	}
	if (info > 0) {
		snprintf(err, ERROR_LEN, "dsyevr: internal Rf_error");
		return nan("Rf_error");
	}
	if (m < dim) {
		snprintf(err, ERROR_LEN, "Sigma not of full rank");
		return nan("Rf_error");
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
	if (err[0]) Rf_error("%s", err);
	return ret;
}

SEXP dmvnorm_wrapper(SEXP Rloc, SEXP Rmean, SEXP Rsigma)
{
	SEXP ret;
	ScopedProtect p1(ret, Rf_allocVector(REALSXP, 1));
	REAL(ret)[0] = dmvnorm(Rf_length(Rloc), REAL(Rloc), REAL(Rmean), REAL(Rsigma));
	return ret;
}
