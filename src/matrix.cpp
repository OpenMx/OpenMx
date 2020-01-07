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
#include "matrix.h"
#include "omxMatrix.h"
#include <Eigen/LU>
#include "EnableWarnings.h"

ThinMatrix::ThinMatrix(omxMatrix *mat)
: rows(mat->rows), cols(mat->cols), t(mat->data) {}

int InvertSymmetricIndef(ThinMatrix mat, const char uplo)
{
	// Not as efficient as dsytrf/dsytri, but we generally
	// use this only when InvertSymmetricPosDef fails or
	// in non-performance critical paths.
	Eigen::Map< Eigen::MatrixXd > Emat(mat.t, mat.rows, mat.cols);
	if (uplo == 'L') {
		Emat.derived() = Emat.selfadjointView<Eigen::Lower>();
	} else if (uplo == 'U') {
		Emat.derived() = Emat.selfadjointView<Eigen::Upper>();
	} else {
		mxThrow("uplo='%c'", uplo);
	}
	Eigen::FullPivLU< Eigen::MatrixXd > lu(Emat);
	if (lu.rank() < mat.rows) return -1;
	Emat.derived() = lu.inverse();
	return 0;
}

void MeanSymmetric(ThinMatrix mat)
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

void SymMatrixMultiply(char side, ThinMatrix amat, ThinMatrix bmat, ThinMatrix cmat)
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

int MatrixSolve(ThinMatrix mat1, ThinMatrix mat2, bool identity)
{
	Eigen::Map< Eigen::MatrixXd > Emat1(mat1.t, mat1.rows, mat1.cols);
	Eigen::Map< Eigen::MatrixXd > Emat2(mat2.t, mat2.rows, mat2.cols);
	Eigen::FullPivLU< Eigen::MatrixXd > lu(Emat1);
	if (lu.rank() < mat1.rows) return -1;

	if (identity) Emat2.setIdentity();
	Emat2 = lu.solve(Emat2);
    
	return 0;
}

int InvertSymmetricPosDef(ThinMatrix mat, char uplo)
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

