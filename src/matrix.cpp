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

Matrix::Matrix(omxMatrix *mat)
: rows(mat->rows), cols(mat->cols), t(mat->data) {}

void InplaceForcePosSemiDef(Matrix mat, double *origEv, double *condnum)
{
    {
        // Variances must be positive so the diagonal must be
        // non-zero. If there are no covariances then dsyevr
        // triggers a valgrind error. It's probably harmless
        // but annoying. This check avoids it.
        Eigen::Map< Eigen::ArrayXXd > tmp(mat.t, mat.rows, mat.cols);
        if ((tmp != 0).count() == mat.rows) return;
    }
    
    const double tooSmallEV = 1e-6;
    double *target = mat.t;
    int numParams = mat.rows;
    if (mat.rows != mat.cols) Rf_error("InplaceForcePosDef must be square");
    
    omxBuffer<double> hessWork(numParams * numParams);
    memcpy(hessWork.data(), target, sizeof(double) * numParams * numParams);
    
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    double abstol = 0;
    int m;
    omxBuffer<double> w(numParams);
    omxBuffer<double> z(numParams * numParams);
    double optWork;
    int optIwork;
    int lwork = -1;
    int liwork = -1;
    int info;
    double realIgn = 0;
    int intIgn = 0;
    omxBuffer<int> isuppz(numParams * 2);
    
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &numParams, hessWork.data(),
                     &numParams, &realIgn, &realIgn, &intIgn, &intIgn, &abstol, &m, w.data(),
                     z.data(), &numParams, isuppz.data(), &optWork, &lwork, &optIwork, &liwork, &info);
    
    lwork = optWork;
    omxBuffer<double> work(lwork);
    liwork = optIwork;
    omxBuffer<int> iwork(liwork);
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &numParams, hessWork.data(),
                     &numParams, &realIgn, &realIgn, &intIgn, &intIgn, &abstol, &m, w.data(),
                     z.data(), &numParams, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info < 0) {
        Rf_error("dsyevr %d", info);
    } else if (info) {
        return;
    }
    
    std::vector<double> evalDiag(numParams * numParams);
    double minEV = 0;
    double maxEV = 0;
    if (origEv) memcpy(origEv, w.data(), sizeof(double) * numParams);
    for (int px=0; px < numParams; ++px) {
        // record how many eigenvalues are zeroed TODO
        if (w[px] < tooSmallEV) {
            evalDiag[px * numParams + px] = tooSmallEV; // exactly zero can still fail
            continue;
        }
        evalDiag[px * numParams + px] = w[px];
        if (w[px] > 0) {
            if (minEV == 0) minEV = w[px];
            else minEV = std::min(minEV, w[px]);
            maxEV = std::max(maxEV, w[px]);
        }
    }
    
    //fc->infoDefinite = true;  actually we don't know!
    if (condnum) *condnum = maxEV/minEV;
    
    Matrix evMat(z.data(), numParams, numParams);
    Matrix edMat(evalDiag.data(), numParams, numParams);
    omxBuffer<double> prod1(numParams * numParams);
    Matrix p1Mat(prod1.data(), numParams, numParams);
    SymMatrixMultiply('R', 'U', 1.0, 0, edMat, evMat, p1Mat);
    char transa = 'N';
    char transb = 'T';
    double alpha = 1.0;
    double beta = 0;
    F77_CALL(dgemm)(&transa, &transb, &numParams, &numParams, &numParams, &alpha,
                    prod1.data(), &numParams, z.data(), &numParams, &beta, target, &numParams);
}

int InvertSymmetricPosDef(Matrix mat, const char uplo)
{
    if (mat.rows != mat.cols) Rf_error("Not square");
    int info;
    F77_CALL(dpotrf)(&uplo, &mat.rows, mat.t, &mat.rows, &info);
    if (info < 0) Rf_error("Arg %d is invalid", -info);
    if (info > 0) return info;
    
    F77_CALL(dpotri)(&uplo, &mat.rows, mat.t, &mat.rows, &info);
    if (info < 0) Rf_error("Arg %d is invalid", -info);
    return info;
}

int InvertSymmetricIndef(Matrix mat, const char uplo)
{
    if (mat.rows != mat.cols) Rf_error("Not square");
    int info;
    omxBuffer<int> ipiv(mat.rows);
    double temp;
    int lwork = -1;
    F77_CALL(dsytrf)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), &temp, &lwork, &info);
    if (info < 0) Rf_error("Arg %d is invalid", -info);
    if (info > 0) return info;
    
    if (lwork < mat.rows) lwork = mat.rows; // for dsytri
    omxBuffer<double> work(lwork);
    F77_CALL(dsytrf)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), work.data(), &lwork, &info);
    if (info < 0) Rf_error("Arg %d is invalid", -info);
    if (info > 0) return info;
    
    F77_CALL(dsytri)(&uplo, &mat.rows, mat.t, &mat.rows, ipiv.data(), work.data(), &info);
    if (info < 0) Rf_error("Arg %d is invalid", -info);
    return info;
}

void MeanSymmetric(Matrix mat)
{
    if (mat.rows != mat.cols) Rf_error("Not conformable");
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

void SymMatrixMultiply(char side, char uplo, double alpha, double beta,
                       Matrix amat, Matrix bmat, Matrix cmat)
{
    if (amat.rows != amat.cols) Rf_error("Not conformable");
    if (bmat.rows != cmat.rows || bmat.cols != cmat.cols) Rf_error("Not conformable");
    int lda;
    if (side == 'R') {
        if (amat.cols != cmat.rows) Rf_error("Not conformable");
        lda = cmat.cols;
    } else if (side == 'L') {
        if (amat.cols != cmat.cols) Rf_error("Not conformable");
        lda = cmat.rows;
    } else {
        Rf_error("Side of %c is invalid", side);
    }
    F77_CALL(dsymm)(&side, &uplo, &cmat.rows, &cmat.cols,
                    &alpha, amat.t, &lda, bmat.t, &bmat.rows,
                    &beta, cmat.t, &cmat.rows);
}

int MatrixSolve(Matrix mat1, Matrix mat2, bool identity)
{
    if (mat1.rows != mat1.cols ||
        mat2.rows != mat2.cols ||
        mat1.rows != mat2.rows) Rf_error("Not conformable");
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
        Rf_error("Arg %d is invalid", -info);
    }
    return info;
}

int MatrixInvert1(Matrix result)
{
    omxBuffer<int> ipiv(result.rows);
    int info;
    F77_CALL(dgetrf)(&(result.cols), &(result.rows), result.t, &(result.rows), ipiv.data(), &info);
    if (info < 0) Rf_error("dgetrf info %d", info);
    if (info > 0) return info;
    
    int opt_lwork = -1;
    double opt_work;
    F77_CALL(dgetri)(&(result.cols), result.t, &(result.rows), ipiv.data(), &opt_work, &opt_lwork, &info);
    if (info != 0) Rf_error("dgetri workspace query failed %d", info);
    
    opt_lwork = opt_work;
    omxBuffer<double> work(opt_lwork);
    F77_CALL(dgetri)(&(result.cols), result.t, &(result.rows), ipiv.data(), work.data(), &opt_lwork, &info);
    if (info < 0) Rf_error("dgetri info %d", info);
    if (info > 0) return info;   // probably would fail at dgetrf already
    
    return 0;
}
