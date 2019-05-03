#undef OMX_BOUNDS_CHECK  // Very, very expensive here

#include "omxDefines.h"
#include "unsupported/Eigen/MatrixFunctions"

// For background, see
// http://epubs.siam.org/doi/abs/10.1137/090768539

void logm_eigen(int n, double *rz, double *out)
{
    Eigen::Map< Eigen::MatrixXd > inMat(rz, n, n);
    Eigen::Map< Eigen::MatrixXd > outMat(out, n, n);
    outMat = inMat.log();
}

void expm_eigen(int n, double *rz, double *out)
{
    Eigen::Map< Eigen::MatrixXd > inMat(rz, n, n);
    Eigen::Map< Eigen::MatrixXd > outMat(out, n, n);
    outMat = inMat.exp();
}
