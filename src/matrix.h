#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "types.h"

double rnd_double();

typedef struct Matrix Matrix;

struct Matrix {
	int rows;
	int cols;
	double *t;

	Matrix() {}
	Matrix(double *_t, int _r, int _c) : rows(_r), cols(_c), t(_t) {}
	Matrix(omxMatrix *mat);
};

typedef struct Param_Obj Param_Obj;

struct Param_Obj {
    Matrix parameter;
    double objValue;
};


#define M(m,x,y) m.t[x+y*m.cols]

void printMatrices();

void freeMatrices();

Matrix QRd(Matrix mainMat, Matrix RHSMat);

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

struct Matrix getRow(struct Matrix t, int row);

struct Matrix setRow(struct Matrix x, int row, struct Matrix y);

struct Matrix getColumn(struct Matrix t, int colNum);

struct Matrix setColumn(struct Matrix x, struct Matrix y, int colNum);

void print(struct Matrix t);

struct Matrix matrix_mult(struct Matrix a,struct Matrix b);

struct Matrix solve(struct Matrix A, struct Matrix b);

struct Matrix cholesky(struct Matrix A);

struct Matrix diag(struct Matrix A);

bool allGreaterThan(struct Matrix t, double value);

double vnorm(struct Matrix t);

double findMin(struct Matrix t);

double findMax(struct Matrix t);

double min(double a, double b);

double ourAbs(double a);

double max(double a, double b);

double dotProduct(struct Matrix a, struct Matrix b);

struct Matrix matrixDotProduct(struct Matrix a, struct Matrix b);

struct Matrix minMaxAbs(struct Matrix t, double tol);

struct Matrix add(struct Matrix x, struct Matrix y);

struct Matrix subtract(struct Matrix x, struct Matrix y);

struct Matrix multiply(struct Matrix x, struct Matrix y);

struct Matrix divide(struct Matrix x, struct Matrix y);

struct Matrix copyInto(struct Matrix x, struct Matrix y, int rowNum, int colStart, int colStop);

struct Matrix rowWiseMin(struct Matrix t);

struct Matrix transposeDP(struct Matrix t);

struct Matrix transpose2D(struct Matrix t);

struct Matrix transpose(struct Matrix t);

struct Matrix negate(struct Matrix t);

struct Matrix duplicateIt(struct Matrix t);

struct Matrix matrixAbs(struct Matrix t);

struct Matrix multiplyByScalar2D(struct Matrix t, double multiplier);

struct Matrix divideByScalar2D(struct Matrix t, double divisor);

struct Matrix checkControlList(struct Matrix t);

struct Matrix subset(struct Matrix t, int row, int colStart, int colStop);

struct Matrix copy(struct Matrix x, struct Matrix y);

struct Matrix rbind(struct Matrix x, struct Matrix y);

//struct Matrix copyThree(struct Matrix a, struct Matrix b, struct Matrix c);

//struct Matrix copyFive(struct Matrix a, struct Matrix b, struct Matrix c, struct Matrix d, struct Matrix e);

struct Matrix timess(struct Matrix a, struct Matrix b);

struct Matrix luSolve(struct Matrix a, struct Matrix b);

struct Matrix qrSolve(struct Matrix a, struct Matrix b);

struct Matrix qrDecomposition(struct Matrix t, bool rDecomp);

struct Matrix rowSort(struct Matrix t);

#endif
