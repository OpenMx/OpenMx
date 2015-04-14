#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include <Eigen/Core>

typedef struct Matrix Matrix;

struct Matrix {
    int rows;
    int cols;
    double *t;
    
    Matrix() : rows(0), cols(0), t(NULL) {};
    Matrix(double *_t, int _r, int _c) : rows(_r), cols(_c), t(_t) {}
    Matrix(omxMatrix *mat);

	template <typename T1> void copyDimsFromEigen(Eigen::MatrixBase<T1> &mb) {
		if (mb.rows() == 1 || mb.cols() == 1) {
			// transpose column vectors (Eigen) to row vectors (CSOLNP)
			cols = mb.rows();
			rows = mb.cols();
		} else {
			rows = mb.rows();
			cols = mb.cols();
		}
	}
	template <typename T1> Matrix(Eigen::MatrixBase<T1> &mb) : t(mb.derived().data()) {
		copyDimsFromEigen(mb);
	};
	template <typename T1> void operator=(Eigen::MatrixBase<T1> &mb) {
		t = mb.derived().data();
		copyDimsFromEigen(mb);
	};
};

template <typename T1>
double solvecond(Eigen::MatrixBase<T1> & inMat)
{
    Eigen::MatrixXd result = inMat;
    int l;
    char JOBZ = 'S';
    int lwork = -1;
    double wkopt;
    int dim_s = std::max(result.cols(), result.rows()); // maybe min is sufficient
    Eigen::ArrayXi iwork(8 * dim_s);
    Eigen::ArrayXd sv(dim_s);
    Eigen::ArrayXd u(dim_s * result.rows());
    Eigen::ArrayXd vt(dim_s * result.cols());
    int result_row = result.rows();
    int result_col = result.cols();
    F77_CALL(dgesdd)(&JOBZ, &result_row, &result_col, result.data(), &result_row, sv.data(), u.data(), &result_row, vt.data(), &result_col, &wkopt, &lwork, iwork.data(), &l);
    lwork = (int)wkopt;
    Eigen::ArrayXd work(lwork);
    F77_CALL(dgesdd)(&JOBZ, &result_row, &result_col, result.data(), &result_row, sv.data(), u.data(), &result_row, vt.data(), &result_col, work.data(), &lwork, iwork.data(), &l);
    
    if (l < 0) Rf_error("the i-th argument had an illegal value");
    else if (l > 0) Rf_error("DBDSDC did not converge, updating process failed.");
    else
    {
        if ((sv == 0).count()) return std::numeric_limits<double>::infinity();
        else return sv.maxCoeff() / sv.minCoeff();
    }
}

template <typename T1>
Eigen::MatrixXd QRdsolve(Eigen::MatrixBase<T1> &mainMat, Eigen::MatrixBase<T1> &RHSMat)

{
    int lwork = 4 * mainMat.rows() * mainMat.cols();
    int l;
    char TRANS = 'N';
    int LDB = std::max(mainMat.cols(), mainMat.rows());
    Eigen::MatrixXd input = mainMat;
    Eigen::MatrixXd result(LDB, RHSMat.cols());
    result.setZero();
    for (int i = 0; i < RHSMat.rows(); i++)
        for (int j = 0; j < RHSMat.cols(); j++)
            result(i, j) = RHSMat(i, j);
    int input_row = input.rows();
    int input_col = input.cols();
    int result_col = result.cols();
    Eigen::ArrayXd work(lwork);
    F77_CALL(dgels)(&TRANS, &input_row, &input_col, &result_col, input.data(), &input_row, result.data(), &LDB, work.data(), &lwork, &l);
    Eigen::MatrixXd Final_result(mainMat.cols(), RHSMat.cols());
    for (int i = 0; i < mainMat.cols(); i++)
    {
        for(int j = 0; j < RHSMat.cols(); j++)
        {
            Final_result(i, j) = result(i, j);
        }
    }
    
    return Final_result;
}

void InplaceForcePosSemiDef(Matrix mat, double *origEv, double *condnum);
Matrix MatrixInvert(Matrix inMat);
int MatrixInvert1(Matrix result);
int InvertSymmetricPosDef(Matrix matd, const char uplo);
int InvertSymmetricIndef(Matrix mat, const char uplo);
void SymMatrixMultiply(char side, char uplo, double alpha, double beta,
                       Matrix amat, Matrix bmat, Matrix cmat);
void MeanSymmetric(Matrix mat);
int MatrixSolve(Matrix mat1, Matrix mat2, bool identity);

template <typename T2>
void rowSort_e(Eigen::MatrixBase<T2>& mat)
{
    int r, i, j;
    for ( r = 0; r < mat.rows(); r++ )
    {
        for ( i = 0; i < mat.cols(); i++ )
        {
            for ( j = 0; j < mat.cols(); j++ )
            {
                if (mat(r, i) < mat(r, j)){
                    double a = mat(r, i);
                    mat(r, i) = mat(r, j);
                    mat(r, j) = a;
                }
            }
        }
    }
}

template <typename T1>
void minMaxAbs(Eigen::MatrixBase<T1> &t, double tol){
    int c;
    for ( c = 0; c < t.size(); c++ )
    {
        double buf = fabs(t[c]);
        buf = std::max(buf, tol);
        buf = std::min(buf, 1.0/(tol));
        t[c] = buf;
    }
}

#endif
