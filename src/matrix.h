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

double rnd_double();

typedef struct Matrix Matrix;

struct Matrix {
    int rows;
    int cols;
    double *t;
    
    Matrix(): t(NULL) {}
    Matrix(double *_t, int _r, int _c) : rows(_r), cols(_c), t(_t) {}
    Matrix(omxMatrix *mat);
};

typedef struct Param_Obj Param_Obj;

struct Param_Obj {
    Matrix parameter;
    double objValue;
};


#define M(m,x,y) m.t[y+x*m.rows]

bool all(Matrix x);

Matrix MatrixToVector(Matrix mat);

Matrix fillMatrix(int cols, int rows, double* array);

void fillMatrix_t(Matrix t, int cols, int rows, double* array);

double solvecond(Matrix inMat);

Matrix diag2(Matrix A);

void chol_lpk(Matrix mainMat);

Matrix QRdsolve(Matrix mainMat, Matrix RHSMat);

void QRdsolve_t(Matrix Final_result, Matrix mainMat, Matrix RHSMat);

void solveinv(Matrix inMat);

void printMatrices();

void freeMatrices();

Matrix QRd(Matrix mainMat, Matrix RHSMat);

void InplaceForcePosSemiDef(Matrix mat, double *origEv, double *condnum);
Matrix MatrixInvert(Matrix inMat);
int MatrixInvert1(Matrix result);
int InvertSymmetricPosDef(Matrix matd, const char uplo);
int InvertSymmetricIndef(Matrix mat, const char uplo);
void SymMatrixMultiply(char side, char uplo, double alpha, double beta,
                       Matrix amat, Matrix bmat, Matrix cmat);
void MeanSymmetric(Matrix mat);
int MatrixSolve(Matrix mat1, Matrix mat2, bool identity);

Matrix condNumPurpose(Matrix inMat);

struct Matrix new_matrix(int cols,int rows);

struct Matrix transposeDotProduct(Matrix t);

struct Matrix fill(int cols, int rows, double value);

void fill_t(Matrix t, int cols, int rows, double value);

struct Matrix getRow(struct Matrix t, int row);

void getRow_t (Matrix toReturn, Matrix t, int row);

void setRow(struct Matrix x, int row, struct Matrix y);

struct Matrix getColumn(struct Matrix t, int colNum);

void getColumn_t (Matrix toReturn, Matrix t, int colNum);

void setColumn(struct Matrix x, struct Matrix y, int colNum);

void print(struct Matrix t);

struct Matrix matrix_mult(struct Matrix a,struct Matrix b);

struct Matrix solve(struct Matrix A, struct Matrix b);

struct Matrix cholesky(struct Matrix A);

struct Matrix diag(struct Matrix A);

void diag_t(Matrix result, Matrix A);

bool allGreaterThan(struct Matrix t, double value);

double vnorm(struct Matrix t);

double findMin(struct Matrix t);

double findMax(struct Matrix t);

double min(double a, double b);

double ourAbs(double a);

double max(double a, double b);

double dotProduct(struct Matrix a, struct Matrix b);

Matrix matrixDotProduct(struct Matrix a, struct Matrix b);

void minMaxAbs(struct Matrix t, double tol);

void add(struct Matrix x, struct Matrix y);

void subtract(struct Matrix x, struct Matrix y);

void multiply(struct Matrix x, struct Matrix y);

void divide(struct Matrix x, struct Matrix y);

void copyInto(struct Matrix x, struct Matrix y, int rowNum, int colStart, int colStop);

struct Matrix rowWiseMin(struct Matrix t);

struct Matrix transposeDP(struct Matrix t);

void transposeDP_t(Matrix toReturn, Matrix t);

struct Matrix transpose2D(struct Matrix t);

struct Matrix transpose(struct Matrix t);

void transpose_t(Matrix t_t, Matrix t);

void negate(struct Matrix t);

struct Matrix duplicateIt(struct Matrix t);

void duplicateIt_t(Matrix result, Matrix t);

void matrixAbs(struct Matrix t);

void multiplyByScalar2D(struct Matrix t, double multiplier);

void divideByScalar2D(struct Matrix t, double divisor);

struct Matrix checkControlList(struct Matrix t);

struct Matrix subset(struct Matrix t, int row, int colStart, int colStop);

void subset_t(Matrix result, Matrix t, int row, int colStart, int colStop);

Matrix copy(Matrix x,  Matrix y);

void copy_t(Matrix result, Matrix x,  Matrix y);

struct Matrix rbind(struct Matrix x, struct Matrix y);

//struct Matrix copyThree(struct Matrix a, struct Matrix b, struct Matrix c);

//struct Matrix copyFive(struct Matrix a, struct Matrix b, struct Matrix c, struct Matrix d, struct Matrix e);

struct Matrix timess(struct Matrix a, struct Matrix b);

void timess_t(Matrix result, Matrix a,  Matrix b);

struct Matrix luSolve(struct Matrix a, struct Matrix b);

struct Matrix qrSolve(struct Matrix a, struct Matrix b);

struct Matrix qrDecomposition(struct Matrix t, bool rDecomp);

void rowSort(struct Matrix t);

void logm_eigen(int n, double *rz, double *out);
void expm_eigen(int n, double *rz, double *out);

#endif
