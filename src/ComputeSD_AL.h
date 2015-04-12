#include "ComputeSD.h"

template <typename T1, typename T2>
static void SD_AL(GradientOptimizerContext &, double rho, Eigen::MatrixBase<T1> &, Eigen::MatrixBase<T2> &);

template <typename T1, typename T2>
static void SD_grad_AL(GradientOptimizerContext &, double, double, Eigen::MatrixBase<T1> &, Eigen::MatrixBase<T2> &);

template <typename T1, typename T2>
static void omxSD_AR(GradientOptimizerContext &, int, double, Eigen::MatrixBase<T1> &, Eigen::MatrixBase<T2> &, double &);

void omxSD_AL(GradientOptimizerContext &);

