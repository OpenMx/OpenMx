//
//  finiteDifferences_csolnp.h
//  
//
//  Created by Mahsa Zahery on 9/1/16.
//
//

#ifndef _finiteDifferences_csolnp_h
#define _finiteDifferences_csolnp_h
// See http://en.wikipedia.org/wiki/Finite_difference

struct forward_difference_grad_c {
    template <typename T1, typename T2>
    void operator()(T1 ff, double refFit, int thrId, double *point,
                    double offset, int px, int numIter, double *Gaprox, Eigen::MatrixBase<T2>& vscale_e, bool flag)
    {
        double orig = point[px];
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            Gaprox[k] = (ff(point, thrId, vscale_e, flag) - refFit) / offset;
        }
        point[px] = orig;
    };
};

struct central_difference_grad_c {
    template <typename T1, typename T2>
    void operator()(T1 ff, double refFit, int thrId, double *point,
                    double offset, int px, int numIter, double *Gaprox, Eigen::MatrixBase<T2>& vscale_e, bool flag)
    {
        double orig = point[px];
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            double f1 = ff(point, thrId, vscale_e, flag);
            point[px] = orig - offset;
            double f2 = ff(point, thrId, vscale_e, flag);
            Gaprox[k] = (f1 - f2) / (2.0 * offset);
        }
        point[px] = orig;
    };
};

struct forward_difference_grad_csolnp {
    template <typename T1, typename T2, typename T3, typename T4>
    void operator()(T3 ff, double refFit, int thrId, double *point, int nineq, int np, double offset, int px, int numIter, double *Gaprox, Eigen::MatrixBase<T4> const & fullPointV, Eigen::MatrixBase<T2>& yy_e, double rho, Eigen::MatrixBase<T1>& vscale_e, Eigen::MatrixBase<T2>& a_e, Eigen::MatrixBase<T2>& b_e, bool flag)
    {
        double orig = point[px];
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            const_cast< Eigen::MatrixBase<T4>& >(fullPointV)[nineq + px] = orig + offset;
            Gaprox[k] = (ff(point, thrId, px, np, fullPointV, yy_e, rho, vscale_e, a_e, b_e, flag) - refFit) / offset;
            //offset *= .5;
        }
        point[px] = orig;
        const_cast< Eigen::MatrixBase<T4>& >(fullPointV)[nineq + px] = orig;
    };
};

struct central_difference_grad_csolnp {
    template <typename T1, typename T2, typename T3, typename T4>
    void operator()(T3 ff, double refFit, int thrId, double *point, int nineq, int np, double offset, int px, int numIter, double *Gaprox, Eigen::MatrixBase<T4> const & fullPointV, Eigen::MatrixBase<T2>& yy_e, double rho, Eigen::MatrixBase<T1>& vscale_e, Eigen::MatrixBase<T2>& a_e, Eigen::MatrixBase<T2>& b_e, bool flag)
    {
        double orig = point[px];
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            const_cast< Eigen::MatrixBase<T4>& >(fullPointV)[nineq + px] = orig + offset;
            double f1 = ff(point, thrId, px, np, fullPointV, yy_e, rho, vscale_e, a_e, b_e, flag);
            point[px] = orig - offset;
            const_cast< Eigen::MatrixBase<T4>& >(fullPointV)[nineq + px] = orig - offset;
            double f2 = ff(point, thrId, px, np, fullPointV, yy_e, rho, vscale_e, a_e, b_e, flag);
            Gaprox[k] = (f1 - f2) / (2.0 * offset);
            //offset *= .5;
        }
        point[px] = orig;
        const_cast< Eigen::MatrixBase<T4>& >(fullPointV)[nineq + px] = orig;
    };
};

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void gradientImpl_c(T1 ff, int numThreads, double refFit, Eigen::MatrixBase<T2> &point, int numIter,
                  const double eps, T4 dfn, Eigen::MatrixBase<T3> &gradOut, Eigen::MatrixBase<T5>& vscale_e, bool flag)
{
    Eigen::MatrixXd grid(numIter, point.size());
    numThreads = std::min(numThreads, point.size()); // could break work into smaller pieces TODO
    Eigen::MatrixXd thrPoint(point.size(), numThreads);
    thrPoint.colwise() = point;
    double offset = eps;
#pragma omp parallel for num_threads(numThreads)
    for (int px=0; px < int(point.size()); ++px) {
        int thrId = omp_get_thread_num();
        int thrSelect = numThreads==1? -1 : thrId;
        dfn(ff, refFit, thrSelect, &thrPoint.coeffRef(0, thrId), offset, px, numIter, &grid.coeffRef(0,px), vscale_e, flag);
        for(int m = 1; m < numIter; m++) {						// Richardson Step
            for(int k = 0; k < (numIter - m); k++) {
                // NumDeriv Hard-wires 4s for r here. Why?
                grid(k,px) = (grid(k+1,px) * pow(4.0, m) - grid(k,px))/(pow(4.0, m)-1);
            }
        }
        gradOut[px] = grid(0,px);
    }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void gradientImpl_csolnp(T1 ff, int numThreads, double refFit, Eigen::MatrixBase<T2> &fullPoint, int numIter,
                         const double eps, T4 dfn, Eigen::MatrixBase<T3> &gradOut, int nineq, int nc, int np, Eigen::MatrixBase<T5>& yy_e, double rho, Eigen::MatrixBase<T6>& vscale_e, Eigen::MatrixBase<T5>& a_e, Eigen::MatrixBase<T5>& b_e, bool flag)
{
    Eigen::MatrixXd grid(numIter, np);
    numThreads = std::min(numThreads, np); // could break work into smaller pieces TODO
    Eigen::MatrixXd thrFullPointV(fullPoint.size(), numThreads);
    thrFullPointV.colwise() = fullPoint;
    Eigen::MatrixXd thrPoint(np, numThreads);
    Eigen::VectorXd point = (fullPoint.transpose().block(0, nineq, 1, np).array() * vscale_e.block(0, nc+1, 1, np).array()).transpose();
    thrPoint.colwise() = point;
    double offset = eps;
#pragma omp parallel for num_threads(numThreads)
    for (int px=0; px < int(point.size()); ++px) {
        int thrId = omp_get_thread_num();
        int thrSelect = numThreads==1? -1 : thrId;
        dfn(ff, refFit, thrSelect, &thrPoint.coeffRef(0, thrId), nineq, np, offset, px, numIter, &grid.coeffRef(0,px), thrFullPointV.col(thrId).eval(), yy_e, rho, vscale_e, a_e, b_e, flag);
        for(int m = 1; m < numIter; m++) {                                              // Richardson Step
            for(int k = 0; k < (numIter - m); k++) {
                // NumDeriv Hard-wires 4s for r here. Why?
                grid(k,px) = (grid(k+1,px) * pow(4.0, m) - grid(k,px))/(pow(4.0, m)-1);
            }
        }
        gradOut[px] = grid(0,px);
    }
}

template <typename T1, typename T2, typename T3, typename T4>
void gradient_with_ref_c(GradientAlgorithm algo, int numThreads, int order, double eps, T1 ff, double refFit,
                       Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &gradOut, Eigen::MatrixBase<T4>& vscale_e, bool flag)
{
    switch (algo) {
        case GradientAlgorithm_Forward:{
            forward_difference_grad_c dfn;
            gradientImpl_c(ff, numThreads, refFit, point, order, eps, dfn, gradOut, vscale_e, flag);
            break;}
        case GradientAlgorithm_Central:{
            central_difference_grad_c dfn;
            gradientImpl_c(ff, numThreads, refFit, point, order, eps, dfn, gradOut, vscale_e, flag);
            break;}
        default: Rf_error("Unknown gradient algorithm %d", algo);
    }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void gradient_with_ref_csolnp(GradientAlgorithm algo, int numThreads, int order, double eps, T1 ff, double refFit, Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &gradOut, int nineq, int nc, int np, Eigen::MatrixBase<T4>& yy_e, double rho, Eigen::MatrixBase<T5>& vscale_e, Eigen::MatrixBase<T4>& a_e, Eigen::MatrixBase<T4>& b_e, bool flag)
{
    switch (algo) {
        case GradientAlgorithm_Forward:{
            forward_difference_grad_csolnp dfn;
            gradientImpl_csolnp(ff, numThreads, refFit, point, order, eps, dfn, gradOut, nineq, nc, np, yy_e, rho, vscale_e, a_e, b_e, flag);
            break;}
        case GradientAlgorithm_Central:{
            central_difference_grad_csolnp dfn;
            gradientImpl_csolnp(ff, numThreads, refFit, point, order, eps, dfn, gradOut, nineq, nc, np, yy_e, rho, vscale_e, a_e, b_e, flag);
            break;}
        default: Rf_error("Unknown gradient algorithm %d", algo);
    }
}

struct forward_difference_jacobi_c {
    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    void operator()(T1 ff, Eigen::MatrixBase<T4> &refFit, Eigen::MatrixBase<T2> &point,
                    double offset, int px, int numIter, Eigen::MatrixBase<T3> &Gaprox, Eigen::MatrixBase<T5> &vscale_e, bool flag)
    {
        double orig = point[px];
        Eigen::VectorXd result(refFit.size());
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            ff(point, result, vscale_e, flag);
            Gaprox.col(k) = (result - refFit) / offset;
        }
        point[px] = orig;
    };
};

struct central_difference_jacobi_c {
    template <typename T1, typename T2, typename T3, typename T4, typename T5>
    void operator()(T1 ff, Eigen::MatrixBase<T4> &refFit, Eigen::MatrixBase<T2> &point,
                    double offset, int px, int numIter, Eigen::MatrixBase<T3> &Gaprox, Eigen::MatrixBase<T5> &vscale_e, bool flag)
    {
        double orig = point[px];
        Eigen::VectorXd result1(refFit.size());
        Eigen::VectorXd result2(refFit.size());
        for(int k = 0; k < numIter; k++) {
            point[px] = orig + offset;
            ff(point, result1, vscale_e, flag);
            point[px] = orig - offset;
            ff(point, result2, vscale_e, flag);
            Gaprox.col(k) = (result1 - result2) / (2.0 * offset);
        }
        point[px] = orig;
    };
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void jacobianImpl_c(T1 ff,  Eigen::MatrixBase<T2> &ref, Eigen::MatrixBase<T3> &point,
                  int numIter, const double eps, T4 dfn, Eigen::MatrixBase<T5> &jacobiOut, Eigen::MatrixBase<T6> &vscale_e, bool flag)
{
    double offset = eps;
    // TODO evaluate jacobian in parallel
    for (int px=0; px < int(point.size()); ++px) {
        Eigen::MatrixXd Gaprox(ref.size(), numIter);
        dfn(ff, ref, point, offset, px, numIter, Gaprox, vscale_e, flag);
        for(int m = 1; m < numIter; m++) {						// Richardson Step
            for(int k = 0; k < (numIter - m); k++) {
                // NumDeriv Hard-wires 4s for r here. Why?
                Gaprox.col(k) = (Gaprox.col(k+1) * pow(4.0, m) - Gaprox.col(k))/(pow(4.0, m)-1);
            }
        }
        jacobiOut.row(px) = Gaprox.col(0).transpose();
    }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
void fd_jacobian_c(GradientAlgorithm algo, int numIter, double eps, T1 ff, Eigen::MatrixBase<T2> &ref,
                 Eigen::MatrixBase<T3> &point, Eigen::MatrixBase<T4> &jacobiOut, Eigen::MatrixBase<T5> &vscale_e, bool flag)
{
    switch (algo) {
        case GradientAlgorithm_Forward:{
            forward_difference_jacobi_c dfn;
            jacobianImpl_c(ff, ref, point, numIter, eps, dfn, jacobiOut, vscale_e, flag);
            break;}
        case GradientAlgorithm_Central:{
            central_difference_jacobi_c dfn;
            jacobianImpl_c(ff, ref, point, numIter, eps, dfn, jacobiOut, vscale_e, flag);
            break;}
        default: Rf_error("Unknown gradient algorithm %d", algo);
    }
}

#endif
