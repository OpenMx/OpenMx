#include "omxDefines.h"
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Eigen/Core>
#include "omxBuffer.h"
#include "matrix.h"
#include "glue.h"
#include "EnableWarnings.h"

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
		return nan("mxThrow");
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
	if (err[0]) mxThrow("%s", err);
	return ret;
}

static double
_dmvnorm(char *err, int dim, double *loc, double *mean, double *origSigma)
{
	Eigen::Map< Eigen::MatrixXd > Esigma(origSigma, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Upper > sc(Esigma);

	double dist = mahalanobis(dim, loc, mean, origSigma);

	double got = -(dim * M_LN_SQRT_2PI + 0.5*dist + sc.log_determinant());
	return got;
}

double
dmvnorm(int dim, double *loc, double *mean, double *sigma)
{
	char err[ERROR_LEN];
	err[0] = 0;
	double ret = _dmvnorm(err, dim, loc, mean, sigma);
	if (err[0]) mxThrow("%s", err);
	return ret;
}

SEXP dmvnorm_wrapper(SEXP Rloc, SEXP Rmean, SEXP Rsigma)
{
	SEXP ret;
	ScopedProtect p1(ret, Rf_allocVector(REALSXP, 1));
	REAL(ret)[0] = dmvnorm(Rf_length(Rloc), REAL(Rloc), REAL(Rmean), REAL(Rsigma));
	return ret;
}
