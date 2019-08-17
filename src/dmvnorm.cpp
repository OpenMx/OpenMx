#include "omxDefines.h"
#include <Eigen/Core>
#include "matrix.h"
#include "glue.h"
#include "EnableWarnings.h"

static const int ERROR_LEN = 80;

static double
_mahalanobis(char *err, int dim, double *loc, double *center, double *origCov)
{
	Eigen::VectorXd cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	Eigen::Map<Eigen::MatrixXd> covMat(origCov, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Lower > sc(covMat);
	if (sc.info() != Eigen::Success || !sc.isPositive()) {
		snprintf(err, ERROR_LEN, "Sigma is singular and cannot be inverted");
		return nan("mxThrow");
	}

	sc.refreshInverse();
	return cloc.transpose() * sc.getInverse() * cloc;
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
