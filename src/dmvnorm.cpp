#include "omxDefines.h"
#include <Eigen/Core>
#include "matrix.h"
#include "glue.h"
#include "EnableWarnings.h"

static double
mahalanobis(int dim, double *loc, double *center, double *origCov)
{
	Eigen::VectorXd cloc(dim);
	for (int dx=0; dx < dim; dx++) {
		cloc[dx] = loc[dx] - center[dx];
	}

	Eigen::Map<Eigen::MatrixXd> covMat(origCov, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Lower > sc(covMat);
	if (sc.info() != Eigen::Success || !sc.isPositive()) {
		mxThrow("mahalanobis: sigma is singular and cannot be inverted");
	}

	sc.refreshInverse();
	return cloc.transpose() * sc.getInverse() * cloc;
}

double dmvnorm(int dim, double *loc, double *mean, double *origSigma)
{
	Eigen::Map< Eigen::MatrixXd > Esigma(origSigma, dim, dim);
	SimpCholesky< Eigen::MatrixXd, Eigen::Upper > sc(Esigma);

	double dist = mahalanobis(dim, loc, mean, origSigma);

	double got = -(dim * M_LN_SQRT_2PI + 0.5*dist + sc.log_determinant());
	return got;
}

SEXP dmvnorm_wrapper(SEXP Rloc, SEXP Rmean, SEXP Rsigma)
{
	SEXP ret;
	ScopedProtect p1(ret, Rf_allocVector(REALSXP, 1));
	REAL(ret)[0] = dmvnorm(Rf_length(Rloc), REAL(Rloc), REAL(Rmean), REAL(Rsigma));
	return ret;
}
