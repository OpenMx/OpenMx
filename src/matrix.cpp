#include "omxDefines.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <iostream>
using std::cout;
using std::endl;
#include <list>
#include <algorithm>
#include <iterator>
#include <limits>
#include "omxBuffer.h"
#include "matrix.h"
#include "omxMatrix.h"
#include <Eigen/LU>
#include "EnableWarnings.h"

Matrix::Matrix(omxMatrix *mat)
: rows(mat->rows), cols(mat->cols), t(mat->data) {}

int InvertSymmetricIndef(Matrix mat, const char uplo)
{
    if (mat.rows != mat.cols) mxThrow("Not square");
    int info;
    omxBuffer<int> ipiv(mat.rows);
    double temp;
    int lwork = -1;
    F77_CALL(dsytrf)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), &temp, &lwork, &info);
    if (info < 0) mxThrow("Arg %d is invalid", -info);
    if (info > 0) return info;
    
    if (lwork < mat.rows) lwork = mat.rows; // for dsytri
    omxBuffer<double> work(lwork);
    F77_CALL(dsytrf)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), work.data(), &lwork, &info);
    if (info < 0) mxThrow("Arg %d is invalid", -info);
    if (info > 0) return info;
    
    F77_CALL(dsytri)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), work.data(), &info);
    if (info < 0) mxThrow("Arg %d is invalid", -info);
    return info;
}

void MeanSymmetric(Matrix mat)
{
    if (mat.rows != mat.cols) mxThrow("Not conformable");
    const int len = mat.rows;
    
    for (int v1=1; v1 < len; ++v1) {
        for (int v2=0; v2 < v1; ++v2) {
            int c1 = v1 * len + v2;
            int c2 = v2 * len + v1;
            double mean = (mat.t[c1] + mat.t[c2])/2;
            mat.t[c1] = mean;
            mat.t[c2] = mean;
        }
    }
}

void SymMatrixMultiply(char side, Matrix amat, Matrix bmat, Matrix cmat)
{
	using Eigen::Map;
	using Eigen::MatrixXd;
	Map< MatrixXd > Ea(amat.t, amat.rows, amat.cols);
	Map< MatrixXd > Eb(bmat.t, bmat.rows, bmat.cols);
	Map< MatrixXd > Ec(cmat.t, cmat.rows, cmat.cols);

    if (side == 'R') {
	Ec.derived() = Eb * Ea.selfadjointView<Eigen::Upper>();
    } else if (side == 'L') {
	Ec.derived() = Ea.selfadjointView<Eigen::Upper>() * Eb;
    } else {
        mxThrow("Side of %c is invalid", side);
    }
}

int MatrixSolve(Matrix mat1, Matrix mat2, bool identity)
{
    if (mat1.rows != mat1.cols ||
        mat2.rows != mat2.cols ||
        mat1.rows != mat2.rows) mxThrow("Not conformable");
    const int dim = mat1.rows;
    
    omxBuffer<double> pad(dim * dim);
    memcpy(pad.data(), mat1.t, sizeof(double) * dim * dim);
    
    if (identity) {
        for (int rx=0; rx < dim; rx++) {
            for (int cx=0; cx < dim; cx++) {
                mat2.t[rx * dim + cx] = rx==cx? 1 : 0;
            }
        }
    }
    
    std::vector<int> ipiv(dim);
    int info;
    F77_CALL(dgesv)(&dim, &dim, pad.data(), &dim, ipiv.data(), mat2.t, &dim, &info);
    if (info < 0) {
        mxThrow("Arg %d is invalid", -info);
    }
    return info;
}

int InvertSymmetricPosDef(Matrix mat, char uplo)
{
	Eigen::Map<Eigen::MatrixXd> Emat(mat.t, mat.rows, mat.cols);
	if (uplo == 'L') {
		SimpCholesky< Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower > sc(Emat);
		if (sc.info() != Eigen::Success || !sc.isPositive()) {
			return -1;
		} else {
			sc.refreshInverse();
			Emat.derived() = sc.getInverse();
			return 0;
		}
	} else if (uplo == 'U') {
		SimpCholesky< Eigen::Ref<Eigen::MatrixXd>, Eigen::Upper > sc(Emat);
		if (sc.info() != Eigen::Success || !sc.isPositive()) {
			return -1;
		} else {
			sc.refreshInverse();
			Emat.derived() = sc.getInverse();
			return 0;
		}
	} else {
		mxThrow("uplo invalid");
	}
}

