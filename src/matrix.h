#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "omxDefines.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

struct ThinMatrix { // deprecated
    int rows;
    int cols;
    double *t;
    
    ThinMatrix() : rows(0), cols(0), t(NULL) {};
    ThinMatrix(double *_t, int _r, int _c) : rows(_r), cols(_c), t(_t) {}
    ThinMatrix(omxMatrix *mat);

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
	template <typename T1> ThinMatrix(Eigen::MatrixBase<T1> &mb) : t(mb.derived().data()) {
		copyDimsFromEigen(mb);
	}
	template <typename T1> void operator=(Eigen::MatrixBase<T1> &mb) {
		t = mb.derived().data();
		copyDimsFromEigen(mb);
	}
};

template <typename T1>
void ForceInvertSymmetricPosDef(Eigen::MatrixBase<T1> &mat, double *origEv, double *condnum)
{
	// lower triangle is referenced
	Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > es(mat);
	if (origEv) memcpy(origEv, es.eigenvalues().data(), sizeof(double) * mat.rows());

	const double tooSmallEV = 1e-6;
	// The abs() preserves the measure of curvature. This is conservative
	// compared to setting the curvature to some value near zero.
	Eigen::VectorXd ev = 1.0 / es.eigenvalues().array().abs().max(tooSmallEV);
	if (condnum) *condnum = ev.maxCoeff() / ev.minCoeff();
	mat.derived() = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().inverse();
}

template <typename T1>
void ForceInvertSymmetricPosDef(Eigen::MatrixBase<T1> &mat)
{ ForceInvertSymmetricPosDef(mat, 0, 0); }

ThinMatrix MatrixInvert(ThinMatrix inMat);
int InvertSymmetricIndef(ThinMatrix mat, const char uplo);
void SymMatrixMultiply(char side, ThinMatrix amat, ThinMatrix bmat, ThinMatrix cmat);
void MeanSymmetric(ThinMatrix mat);
int MatrixSolve(ThinMatrix mat1, ThinMatrix mat2, bool identity);

int InvertSymmetricPosDef(ThinMatrix mat, char uplo);

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
