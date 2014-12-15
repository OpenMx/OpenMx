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
#include "unsupported/Eigen/MatrixFunctions"

template <typename T> void printList( const std::list< T > &listRef);

static std::list< double* > matrices;

void freeMatrices(){
    while (!matrices.empty()){
        free(matrices.front());
        matrices.pop_front();
    }
}

Matrix::Matrix(omxMatrix *mat)
: rows(mat->rows), cols(mat->cols), t(mat->data) {}

Matrix new_matrix(int cols,int rows)
{
    if (rows < 0 || cols < 0) Rf_error("Cannot create matrix smaller than 0,0");
    Matrix t;
    t.rows=rows;
    t.cols=cols;
    if (rows == 0 || cols == 0) {
        t.t = NULL;
    } else {
        t.t=(double *)malloc(sizeof(double)*cols*rows);
        matrices.push_front(t.t);
    }
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++) {
            M(t,j,i)=nan("uninit");
        }
    }
    return t;
}

Matrix fill(int cols, int rows, double value){
    Matrix t = new_matrix(cols, rows);
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++) {
            M(t,j,i)=value;
        }
    }
    return t;
}

void fill_t(Matrix t, int cols, int rows, double value){
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++) {
            M(t,j,i)=value;
        }
    }
}

Matrix getRow (Matrix t, int row){
    Matrix toReturn = fill(t.cols, 1, (double)0.0);
    int i;
    for (i=0;i < t.cols; i++){
        M(toReturn,i,0) = M(t,i,row);
    }
    return toReturn;
}

void getRow_t (Matrix toReturn, Matrix t, int row){
    int i;
    for (i=0;i < t.cols; i++){
        M(toReturn,i,0) = M(t,i,row);
    }
}

void setRow (Matrix x, int row,  Matrix y){
    
    int i;
    for (i=0;i < x.cols; i++){
        M(x,i,row) = M(y,i,0);
    }
}

void setRowInplace( Matrix x, int cc,  Matrix y)
{
    Eigen::Map< Eigen::MatrixXd > xx(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::VectorXd > yy(y.t, y.rows * y.cols);
    xx.row(cc) = yy;
}

Matrix getColumn (Matrix t, int colNum)
{
    int r, c;
    int index = 0;
    Matrix toReturn = fill(t.rows, 1, (double)0.0);
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            if (c==colNum){
                M(toReturn, index, 0) = M(t, c, r);
                index++;
            }
        }
    }
    return toReturn;
}

void getColumn_t (Matrix toReturn, Matrix t, int colNum)
{
    int r, c;
    int index = 0;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            if (c==colNum){
                M(toReturn, index, 0) = M(t, c, r);
                index++;
            }
        }
    }
}

void setColumnInplace( Matrix x, Matrix y, int cc)
{
    Eigen::Map< Eigen::MatrixXd > xx(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::VectorXd > yy(y.t, y.rows * y.cols);
    xx.col(cc) = yy;
}

void print(Matrix t) {
    int i,j;
    for(i=0;i<t.rows;i++) {
        printf("| ");
        for(j=0;j<t.cols;j++)
            printf("%.20f ",M(t,j,i));
        printf("|\n");
    }
    printf("\n");
}

Matrix matrix_mult(Matrix a, Matrix b) {
    int x,y,z;
    Matrix r;
    r=new_matrix(a.rows,a.rows);
    
    for(x=0;x<a.rows;x++)
        for(y=0;y<b.cols;y++) {
            M(r,y,x)=0.0;
            for(z=0;z<a.cols;z++) {
                M(r,y,x)+=M(a,z,x)*M(b,y,z);
            }
        }
    return r;
}

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

Matrix solve(Matrix X,  Matrix y){
    Matrix A = duplicateIt(X);
    Matrix b = duplicateIt(y);
    int N = b.cols;
    /** b is a 1 row matrix **/
    int p, i, j;
    for (p=0; p<N; p++){
        // find pivot row and swap
        int max = p;
        for (i = p + 1; i<N; i++) {
            if ( fabs(M(A,p,i)) > fabs(M(A,p,max)) ) {
                max = i;
            }
        }
        
        Matrix temp = fill(A.cols, 1, (double)0.0);
        
        temp = getRow(A,p); setRow(A,p, getRow(A,max));
        setRow(A,max, temp);
        double t = M(b,p,0); M(b,p,0) = M(b,max,0); M(b,max,0)= t;
        
        
        // pivot within A and b
        for (i = p + 1; i < N; i++) {
            double alpha = (M(A,p,i)/M(A,p,p));
            M(b,i,0) -= alpha * M(b,p,0);
            for (j = p; j < N; j++) {
                M(A,j,i) -= alpha * M(A,j,p);
            }
        }
    }
    
    // back substitution
    Matrix x = fill(N, 1, (double)0.0);
    for (i = (N - 1); i >= 0; i--) {
        double sum = 0.0;
        for (j = i + 1; j < N; j++) {
            sum += M(A,j,i) * M(x,j,0);
        }
        M(x, i, 0) = (M(b, i, 0) - sum) / M(A, i, i);
    }
    return duplicateIt(x);
}

// return Cholesky factor L of psd matrix A = L L^T
Matrix cholesky(Matrix A){
    Matrix L = fill(A.rows, A.cols, (double)0.0);
    int i,j, k;
    for (i=0; i<A.cols; i++){
        for (j=0; j<= i; j++){
            double sum = 0.0;
            for (k=0; k<j; k++){
                sum += (M(L,k,i) * M(L,k,j));
            }
            if (i == j) {
                M(L,i,i) = sqrt(M(A,i,i) - sum);
            }
            else{
                M(L,j,i) = 1.0 / M(L,j,j) * (M(A,j,i) - sum);
            }
        }
    }
    return L;
}

Matrix diag(Matrix A){
    int i,j;
    double len;
    Matrix result;
    if (A.cols > A.rows)
    {
        len = A.cols;
        result = fill(len, len, (double)0.0);
        for (i=0; i<len; i++){
            for (j=0; j<=len; j++){
                if (i==j)
                    M(result, j, i) = M(A, j, 0);
            }
        }
    }
    
    else
    {
        len = A.rows;
        result = fill(len, len, (double)0.0);
        for (i=0; i<len; i++){
            for (j=0; j<=len; j++){
                if (i==j)
                    M(result, j, i) = M(A, 0, j);
            }
        }
    }
    return result;
}

void diag_t(Matrix result, Matrix A){
    int i,j;
    double len;
    if (A.cols > A.rows)
    {
        len = A.cols;
        for (i=0; i<len; i++){
            for (j=0; j<=len; j++){
                if (i==j)
                    M(result, j, i) = M(A, j, 0);
            }
        }
    }
    
    else
    {
        len = A.rows;
        for (i=0; i<len; i++){
            for (j=0; j<=len; j++){
                if (i==j)
                    M(result, j, i) = M(A, 0, j);
            }
        }
    }
}


Matrix diag2(Matrix A){
    int i,j;
    Matrix result = fill(A.cols, A.cols, (double)0.0);
    for (i=0; i<A.cols; i++){
        for (j=0; j<=A.cols; j++){
            if (i==j)
                M(result, j, 0) = M(A, j, i);
        }
    }
    return result;
}

bool allGreaterThan(Matrix t, double value){
    int i,j;
    for (i=0; i<t.rows; i++){
        for (j=0; j<t.cols; j++){
            if (value >= M(t,j,i)){
                return false;
            }
        }
    }
    return true;
}

double vnorm(Matrix t){
    double result = 0;
    int i,j;
    for ( i = 0; i < t.rows; i++ )
    {
        for ( j = 0; j < t.cols; j++ )
        {
            result = result + (M(t,j,i) * M(t,j,i));
        }
    }
    
    return sqrt(result);
}

double findMin(Matrix t){
    int i,j;
    double min = DBL_MAX;
    for ( i = 0; i < t.rows; i++ )
    {
        for ( j = 0; j < t.cols; j++ )
        {
            if ( min > M(t, j, i) ){
                min = M(t, j, i);
            }
        }
    }
    return min;
}

double findMax(Matrix t){
    int i,j;
    double max = -999999.0;
    
    for ( i = 0; i < t.rows; i++ )
    {
        for ( j = 0; j < t.cols; j++ )
        {
            if (max < M(t,j,i) ){
                max = M(t,j,i);
            }
        }
    }
    return max;
}

double min(double a, double b){
    if (a < b){
        return a;
    }
    else {
        return b;
    }
}

double ourAbs(double a){
    if (a < 0.0){
        return (a * -1.0);
    }
    else {
        return a;
    }
}

double max(double a, double b){
    if (a > b){
        return a;
    }
    else {
        return b;
    }
}


double dotProduct(Matrix a, Matrix b) { //Mahsa: a.cols element of matrix a is multiplied by a.cols first element of matrix b. e.g. if matrix a is 2 by 4, and matrix b is the transpose of matrix a (hence 4 by 2), then dotProduct multiplied first row of matrix a by first row and second row of matrix b (2+2 = 4 first elements of matrix b), and adds them all together.
    double runningSum = 0;
    int index;
    for (index = 0; index < a.cols; index++)
    {
        runningSum += M(a, index, 0) * M(b, index, 0);
    }
    return runningSum;
}

Matrix matrixDotProduct(Matrix a,  Matrix b){
    int r;
    Matrix result = fill(a.cols, 1, (double)0.0);
    for ( r = 0; r < a.rows; r++ )
    {
        
        M(result, r, 0) = dotProduct(getRow(a, r), getRow(b, 0));
    }
    return result;
}

void minMaxAbs(Matrix t, double tol)
{
    // for each x[i]: min( max( abs(x[i]), tol ), 1/tol )
    int r,c;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            double buf = ourAbs(M(t,c,r));
            buf = max(buf, tol);
            buf = min(buf, 1.0/(tol));
            M(t,c,r) = buf;
        }
    }
}

void addEigen(Matrix x,  Matrix y)
{
    Eigen::Map< Eigen::MatrixXd > firstM(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::MatrixXd > secondM(y.t, y.rows, y.cols);
    
    if (x.cols != y.cols || x.rows != y.rows) {
        if (x.cols == y.rows) {
            firstM += secondM.transpose();
        }
        else Rf_error("CSOLNP BUG: noncomformant matrices are added");
    } else {
	    firstM += secondM;
    }
}

void subtractEigen(Matrix x,  Matrix y)
{
    Eigen::Map< Eigen::MatrixXd > firstM(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::MatrixXd > secondM(y.t, y.rows, y.cols);
    firstM -= secondM;
}

void multiplyEigen(Matrix x,  Matrix y)
{
    
    Eigen::Map< Eigen::MatrixXd > firstM(x.t, x.rows, x.cols);
    if (x.cols == y.cols && x.rows == y.rows){
        Eigen::Map< Eigen::MatrixXd > secondM(y.t, y.rows, y.cols);
        firstM = firstM.cwiseProduct(secondM);
    }
    else if (y.cols > 1 && y.rows > 1)
        Rf_error("Only a vector is acceptable");
    else if ((x.cols * x.rows) % (y.cols * y.rows) != 0)
        Rf_error("longer object length is not a multiple of shorter object length");
    else{
        Eigen::Map< Eigen::VectorXd > secondM(y.t, y.cols);
        firstM = secondM.asDiagonal() * firstM;
    }
}

void divideEigen(Matrix x,  Matrix y)
{
    if (x.cols != y.cols) {
        if (x.cols > y.cols) y = copy(y, fill(x.cols - y.cols, 1, (double)1.0));
    }
    Eigen::Map< Eigen::MatrixXd > firstM(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::VectorXd > secondM(y.t, y.cols);
    firstM = firstM * secondM.asDiagonal().inverse();
}

void copyInto(Matrix x,  Matrix y, int rowNum, int colStart, int colStop){
    
    int r, c, index;
    index = 0;
    for ( r = 0; r < x.rows; r++ )
    {
        for ( c = 0; c < x.cols; c++ )
        {
            if (r == rowNum){
                if (c >= colStart && c <= colStop){
                    M(x, c, r) = M(y, index, 0);
                    index++;
                }
            }
        }
    }
}

void copyIntoInplace(Matrix x,  Matrix y, int rowNum, int colStart, int colStop)
{
    Eigen::Map< Eigen::MatrixXd > xx(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::MatrixXd > yy(y.t, y.rows, y.cols);
    int len = 1 + colStop - colStart;
    xx.block(rowNum, colStart, 1, len) = yy.block(rowNum, 0, 1, len);
}

Matrix rowWiseMin(Matrix t)
{
    Matrix mins = fill(t.rows, 1, (double)0.0);
    int r,c;
    for ( r = 0; r < t.rows; r++ )
    {
        double min = DBL_MAX;
        for ( c = 0; c < t.cols; c++ )
        {
            if ( M(t,c,r) < min){
                min = M(t,c,r);
            }
        }
        M(mins, r, 0) = min;
    }
    return mins;
}

Matrix transpose2D(Matrix t){ //Mahsa: 1) first row is copied into a matrix with the size equal to the argument matrix. then each element of this row is multiplied by the corresponding row in the copied matrix. i.e first element is multiplied by the first row. second element is multiplied by the second row, and so on.
    int r, c;
    Matrix result = fill(t.cols,t.cols,(double)0.0);
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(result,c,r) = M(t,c,0);
        }
    }
    
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(result,c,r) = M(result,c,r) * M(t,r,0);
        }
    }
    
    return result;
    
}

Matrix transposeDotProduct(Matrix t){ //Mahsa: 1) takes the minimum dimension. 2) multiplies the first row by itself and the other rows(as far as the minimum dimension allows). multiplies the second row by itself and the rest, and so on. This gives the same answer as "a%*%t(a)" in R.
    int r,c;
    double minDim = min(t.cols, t.rows);
    Matrix result = fill(minDim,minDim,(double)0.0);
    
    for ( r = 0; r < minDim; r++ )
    {
        for ( c = 0; c < minDim; c++ )
        {
            M(result,c,r) = dotProduct(getRow(t,r), getRow(t,c));
        }
    }
    return result;
}

Matrix transposeDP(Matrix t){ //Mahsa: same as transpose2D function
    Matrix toReturn = fill(t.cols, t.cols, (double)0.0);
    int r, c;
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(toReturn, c, r) = M(t, r, 0) * M(t, c, 0);
        }
    }
    return toReturn;
}

void transposeDP_t(Matrix toReturn, Matrix t){ //Mahsa: same as transpose2D function
    int r, c;
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(toReturn, c, r) = M(t, r, 0) * M(t, c, 0);
        }
    }
}


Matrix transpose(Matrix t){ // Mahsa: simply transposes the matrix
    Matrix result = fill(t.rows, t.cols, (double)0.0);
    int r,c;
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.rows; c++ )
        {
            M(result, c, r) = M(t, r, c);
        }
    }
    return result;
}

void negate(Matrix t)
{
    int r, c;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(t, c, r) = (M(t, c, r) * -1.0);
        }
    }
}

Matrix duplicateIt(Matrix t)
{
    int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(result, c, r) = M(t, c, r);
        }
    }
    return result;
}

void duplicateIt_t(Matrix result, Matrix t)
{
    int r, c;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(result, c, r) = M(t, c, r);
        }
    }
}

double matrixMaxAbs(Matrix t)
{
    Eigen::Map< Eigen::ArrayXXd > tt(t.t, t.rows, t.cols);
    return tt.abs().maxCoeff();
}

void multiplyByScalar2D(Matrix t, double multiplier)
{
    int r, c;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(t, c, r) = M(t, c, r)*multiplier;
        }
    }
}

void divideByScalar2D(Matrix t, double divisor)
{
    int r, c;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(t, c, r) = M(t, c, r)/divisor;
        }
    }
}

Matrix multiplyByScalar2D_Eigen(Matrix t, double multiplier)
{
    Eigen::Map< Eigen::ArrayXXd > tt(t.t, t.rows, t.cols);
    Matrix result = new_matrix(t.cols, t.rows);
    Eigen::Map< Eigen::ArrayXXd > dest(result.t, result.rows, result.cols);
    dest = tt * multiplier;
    return result;
}

Matrix divideByScalar2D_Eigen(Matrix t, double divisor)
{
    return multiplyByScalar2D_Eigen(t, 1.0/divisor);
}


Matrix checkControlList(Matrix t){
    Matrix result = fill(t.cols, t.rows, (double)0.0);
    int Rf_length = t.cols;
    int c;
    
    M(result, 0, 0)=1;
    M(result, 1, 0)=400;
    M(result, 2, 0)=800;
    M(result, 3, 0)=1.0e-7;
    M(result, 4, 0)=1.0e-8;
    M(result, 5, 0)=1;
    
    for ( c = 0; c < Rf_length; c++ )
    {
        M(result, c, 0) = M(t, c, 0);
    }
    return result;
}

Matrix subset(Matrix t, int row, int colStart, int colStop)
{//subset(p0, 0, 0, npic-1)
    
    int r, c;
    int count = (colStop - colStart)+1;
    int index = 0;
    
    Matrix result = fill(count, 1, (double)0.0);
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            if (r == row){
                if ( (c >= colStart) && (c <= colStop) ){
                    M(result, index, 0) = M(t, c, r);
                    index++;
                }
            }
        }
    }
    return result;
}

void subset_t(Matrix result, Matrix t, int row, int colStart, int colStop)
{
    int r, c;
    int index = 0;
    
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            if (r == row){
                if ( (c >= colStart) && (c <= colStop) ){
                    M(result, index, 0) = M(t, c, r);
                    index++;
                }
            }
        }
    }
}

void subsetEigen(Matrix result, Matrix x, int rowNum, int colStart, int colStop){
    Eigen::Map < Eigen::MatrixXd> src(x.t, x.rows, x.cols);
    Eigen::Map < Eigen::MatrixXd> dest(result.t, result.rows, result.cols);
    dest.row(rowNum) = src.block(rowNum, colStart, 1, colStop - colStart + 1);
}

void copy_t(Matrix result, Matrix x,  Matrix y){
    
    int totalRows = x.rows;
    int totalCols = x.cols + y.cols;
    
    int r, c;
    for ( r = 0; r < totalRows; r++ )
    {
        for ( c = 0; c < totalCols; c++ )
        {
            if (c < x.cols){
                M(result, c, r) = M(x, c, r);
            }
            else {
                M(result, c, r) = M(y, (c-x.cols), r);
            }
        }
    }
}

void copyEigen(Matrix result, Matrix x,  Matrix y)
{
    Eigen::Map< Eigen::MatrixXd > resultEigen(result.t, result.rows, result.cols);
    Eigen::Map< Eigen::MatrixXd > firstM(x.t, x.rows, x.cols);
    Eigen::Map< Eigen::MatrixXd > secondM(y.t, y.rows, y.cols);
    resultEigen << firstM, secondM;
}

Matrix copy(Matrix x,  Matrix y){
    
    int totalRows = x.rows;
    int totalCols = x.cols + y.cols;
    Matrix result = fill(totalCols, totalRows, (double)0.0);
    
    int r, c;
    for ( r = 0; r < totalRows; r++ )
    {
        for ( c = 0; c < totalCols; c++ )
        {
            if (c < x.cols){
                M(result, c, r) = M(x, c, r);
            }
            else {
                M(result, c, r) = M(y, (c-x.cols), r);
            }
        }
    }
    
    return result;
}

Matrix rbind(Matrix x,  Matrix y){
    Matrix result = copy(x, y);
    return result;
}

/*Matrix copyThree(Matrix a,  Matrix b,  Matrix c){
 Matrix result = copy(a, b);
 Matrix result_three = copy(result, c);
 }*/

/*Matrix copyFive(Matrix a,  Matrix b,  Matrix c,  Matrix d,  Matrix e){
 copy(copy(copy(a,b), copy(c,d)), e);
 }*/

Matrix timess(Matrix a,  Matrix b){
    int i, j, k;
    if (a.cols != b.rows)
    {
        Rf_error("CSOLNP BUG: noncomformant matrices");
    }
    Matrix result = fill(b.cols, a.rows, (double)0.0);
    
    Eigen::ArrayXd Bcolj;
    Bcolj.resize(a.cols);
    
    for (j=0; j<b.cols; j++){
        for (k=0; k<a.cols; k++){
            Bcolj[k] = M(b, j, k);
        }
        for (i=0; i<a.rows; i++){
            double s = 0;
            for (k=0; k<a.cols; k++){
                
                s+= M(a, k, i) * Bcolj[k];
            }
            M(result, j, i) = s;
        }
    }
    return result;
}

void timess_t(Matrix result, Matrix a,  Matrix b){
    int i, j, k;
    if (a.cols != b.rows)
    {
        Rf_error("CSOLNP BUG: noncomformant matrices");
    }
    Matrix Bcolj = fill(a.cols, 1, (double)0.0);
    for (j=0; j<b.cols; j++){
        for (k=0; k<a.cols; k++){
            M(Bcolj, k, 0) = M(b, j, k);
        }
        for (i=0; i<a.rows; i++){
            double s = 0;
            for (k=0; k<a.cols; k++){
                
                s+= M(a, k, i) * M(Bcolj, k, 0);
            }
            M(result, j, i) = s;
        }
    }
}

void timessEigen(Matrix result, Matrix a,  Matrix b){
    Eigen::Map< Eigen::MatrixXd > firstM(a.t, a.rows, a.cols);
    Eigen::Map< Eigen::MatrixXd > secondM(b.t, b.rows, b.cols);
    Eigen::Map< Eigen::MatrixXd > resultM(result.t, result.rows, result.cols);
    resultM = firstM * secondM;
    Eigen::Map< Eigen::MatrixXd >(result.t, resultM.rows(), resultM.cols()) = resultM;
}

Matrix luSolve(Matrix a,  Matrix b){
    int m, n, pivsign;
    double EMPTY = -999999.0;
    m = a.rows;
    n = a.cols;
    
    Matrix LU = duplicateIt(a);
    Matrix piv = fill(m, 1, EMPTY);
    
    int i, j, k;
    for (i = 0; i < m; i++) {
        M(piv, i, 0) = i;
    }
    
    pivsign = 1;
    
    // Outer loop
    
    for (j = 0; j < n; j++) {
        
        // Apply previous transformations.
        for (i = 0; i < m; i++) {
            
            // Most of the time is spent in the following dot product.
            int kmax = min(i,j);
            double s = 0.0;
            for (k = 0; k < kmax; k++) {
                s += M(LU, k, i) * M(LU,j, k);
            }
            M(LU,j,i) -= s;
        }
        
        // Find pivot and exchange if necessary.
        int p = j;
        for (i = j+1; i < m; i++) {
            if ( fabs(M(LU,j, i)) > fabs(M(LU,j, p)) ){
                p = i;
            }
        }
        
        if (p != j) {
            for (k = 0; k < n; k++) {
                double t = M(LU, k, p); M(LU, k, p) = M(LU, k, j); M(LU, k, j) = t;
            }
            int k = M(piv, p, 0); M(piv, p, 0) = M(piv, j, 0); M(piv, j, 0) = k;
            pivsign = -1 * pivsign;
        }
        
        
        // Compute multipliers.
        
        if (j < m && M(LU, j, j) != 0.0) {
            for (i = j+1; i < m; i++) {
                M(LU, j, i) /= M(LU, j, j);
            }
        }
    }
    
    
    int nx = b.rows;
    // array of row indices
    int pivLength = m;
    
    for (i=0; i<m; i++){
        if (M(piv, i, 0) == EMPTY){
            pivLength = i;
        }
    }
    
    // Copy right hand side with pivoting
    Matrix X = fill(nx, pivLength, (double)0.0);
    
    
    for (i=0; i<pivLength; i++){
        for (j=0; j<nx; j++){
            int index = M(piv, i, 0);
            M(X, j, i) = M(b, j, index);
        }
    }
    
    
    // Solve L*Y = B(piv,:)
    for (k = 0; k < n; k++) {
        for (i = k+1; i < n; i++) {
            for (j = 0; j < nx; j++) {
                M(X, j, i) -= M(X, j, k)*M(LU, k, i);
            }
        }
    }
    
    // Solve U*X = Y;
    for (k = n-1; k >= 0; k--) {
        for (j = 0; j < nx; j++) {
            M(X, j, k) /= M(LU, k, k);
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < nx; j++) {
                M(X, j, i) -= M(X, j, k)*M(LU, k, i);
            }
        }
    }
    return duplicateIt(X);
}

Matrix qrSolve(Matrix a,  Matrix b){
    // Initialize.
    int m = a.rows;
    int n = a.cols;
    
    Matrix QR = duplicateIt(a);
    Matrix Rdiag = fill(n, 1, (double)0.0);
    
    // Main loop.
    int k, i, j;
    for (k = 0; k < n; k++) {
        // Compute 2-norm of k-th column without under/overflow.
        double nrm = 0;
        for (i = k; i < m; i++) {
            nrm = hypot(nrm, M(QR, k, i));
        }
        if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (M(QR, k, k) < 0) {
                nrm = -1*nrm;
            }
            for (i = k; i < m; i++) {
                M(QR, k, i) /= nrm;
            }
            M(QR,k,k) += 1.0;
            
            // Apply transformation to remaining columns.
            for (j = k+1; j < n; j++) {
                double s = 0.0;
                for (i = k; i < m; i++) {
                    s += M(QR, k, i) * M(QR, j, i);
                }
                s = -s/M(QR, k, k);
                for (i = k; i < m; i++) {
                    M(QR, j, i) += s * M(QR, k, i);
                }
            }
        }
        M(Rdiag, k, 0) = -1 * nrm;
    }
    int nx = b.cols;
    Matrix X = duplicateIt(b);
    // Compute Y = transpose(Q)*B
    for (k = 0; k < n; k++) {
        for (j = 0; j < nx; j++) {
            double s = 0.0;
            for (i = k; i < m; i++) {
                s += M(QR, k, i) * M(X, j, i);
            }
            s = (-1*s)/M(QR, k, k);
            for (i = k; i < m; i++) {
                M(X, j, i) += s * M(QR, k, i);
            }
        }
    }
    // Solve R*X = Y;
    for (k = n-1; k >= 0; k--) {
        for (j = 0; j < nx; j++) {
            M(X, j, k) /= M(Rdiag, k, 0);
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < nx; j++) {
                M(X, j, i) -= M(X, j, k) * M(QR, k, i);
            }
        }
    }
    return duplicateIt(X);
}

//
// Creates an orthogonal matrix Q such that QA = R, i.e. A = Q^T R. A may be rectangular.
// @param matrix    Matrix of size n x m
// @param q         Result Matrix of size n x n (orthognal)
// @param r         Result Matrix of size n x m (upper right)
// @param work      Work Vector of size n

Matrix qrDecomposition(Matrix t, bool rDecomp){
    
    int m = t.rows;
    int n = t.cols;
    
    Matrix r = fill(m, n, (double)0.0);
    Matrix q = fill(m, n, (double)0.0);
    
    Matrix work = fill(n, 1, (double)0.0);
    
    int i, j, k;
    
    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            M(r, j, i) = M(t, j, i);
        }
    }
    
    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            if (i==j){
                M(q, j, i) = 1;
            }
        }
    }
    
    // Looping through columns for Householder Transformations
    for (i=0; i<min(m,n); i++) {
        // Preparing lambda fro Householder
        double lambda = 0; for (j=i; j<n; j++) lambda += M(r, i, j) * M(r, i, j);
        lambda = sqrt(lambda);
        if (M(r, i, i) < 0) lambda = -lambda;
        
        // compute Householder
        M(work, i, 0) = M(r, i, i) + lambda; for (j=i+1; j<n; j++) M(work, j, 0) = M(r, i, j);
        double denom = 0; for (j=i; j<n; j++) denom += M(work, j, 0) * M(work, j, 0);
        
        // Apply Householder to result triangular matrix
        for (j=i; j<m; j++) {
            double p = 0; for (k=i; k<n; k++) p+= M(work, k, 0) * M(r, j, k);
            for (k=i; k<n; k++) M(r, j, k) -= 2*M(work, k, 0) * p/denom;
        }
        
        // Apply Householder to result orthogoal matrix
        for (j=0; j<n; j++) {
            double p = 0; for (k=i; k<n; k++) p += M(work, k, 0) * M(q, j, k);
            for (k=i; k<n; k++) M(q, j, k) -= 2*M(work, k, 0) * p/denom;
        }
    }
    
    if (rDecomp){
        return r;
    }
    else{
        return q;
    }
}

void rowSort(Matrix t)
{
    int r, c, i, j;
    for ( r = 0; r < t.rows; r++ )
    {
        for ( c = 0; c < t.cols; c++ )
        {
            M(t, c, r) = M(t, c, r);
        }
        
        for ( i = 0; i < t.cols; i++ )
        {
            for ( j = 0; j < t.cols; j++ )
            {
                if (M(t, i, r) < M(t, j, r)){
                    double a = M(t, i, r);
                    M(t, i, r) = M(t, j, r);
                    M(t, j, r) = a;
                }
            }
        }
    }
}

Matrix QRd(Matrix mainMat, Matrix RHSMat)
{
    int lwork = 4 * mainMat.rows * mainMat.cols;
    int l = 0;
    char TRANS = 'N';
    Matrix result;
    result = duplicateIt(RHSMat);
    double* work = (double*) malloc(lwork * sizeof(double));
    
    F77_CALL(dgels)(&TRANS, &(mainMat.rows), &(mainMat.cols), &(result.cols), mainMat.t, &(mainMat.rows), result.t, &(result.rows), work, &lwork, &l);
    return result;
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

Matrix MatrixInvert(Matrix inMat)
{
    Matrix result = duplicateIt(inMat);
    int info = MatrixInvert1(result);
    if (info) Rf_error("MatrixInvert: attempt to invert singular matrix (info=%d)", info);
    return result;
}



Matrix condNumPurpose(Matrix inMat)
{
    //printf("inMat is: \n");
    //print(inMat); putchar('\n');
    Matrix result = duplicateIt(inMat);
    char jobz = 'N';
    char range = 'A';
    char uplo = 'U';
    double vunused;
    int iunused;
    double abstol = 0;
    int m = inMat.cols;
    Eigen::VectorXd w(inMat.cols);
    Eigen::VectorXd Z(inMat.cols);
    int ldz=1;
    Eigen::VectorXi isuppz(2*inMat.cols);
    int lwork = -1;
    double optlWork;
    int optliWork;
    int liwork = -1;
    int info;
    
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &(result.cols), result.t, &(result.cols),
                     &vunused, &vunused, &iunused, &iunused, &abstol, &m, w.data(), Z.data(), &ldz, isuppz.data(),
                     &optlWork, &lwork, &optliWork, &liwork, &info);
    lwork = optlWork;
    std::vector<double> work(lwork);
    liwork = optliWork;
    std::vector<int> iwork(liwork);
    
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &(result.cols), result.t, &(result.cols),
                     &vunused, &vunused, &iunused, &iunused, &abstol, &m, w.data(), Z.data(), &ldz, isuppz.data(),
                     work.data(), &lwork, iwork.data(), &liwork, &info);
    
    Matrix eigenVals = fill(result.cols, 1, (double)0.0);
    for (int i = 0; i <result.cols; i++)
    {
        M(eigenVals, i, 0) = w[i];
        //printf("w[i] is: \n");
        //printf("%2f", w[i]);
    }
    return eigenVals;
}

void solveinv(Matrix inMat)
{
    //Matrix result;
    //result = duplicateIt(inMat);
    MatrixSolve(inMat, inMat, true);
    //return result;
}

/*int this_main(int argc,char *argv[]) {
 int i,j;
 int size=atoi(argv[1]);
 Matrix a,b,c;
 
 a=new_matrix(size,size);
 b=new_matrix(size,size);
 for(i=0;i<size;i++){
 for(j=0;j<size;j++) {
 M(a,i,j)=rnd_double();
 M(b,i,j)=rnd_double();
 }
 }
 c=matrix_mult(a,b);
 print(c);
 }*/
void QRdsolve_t(Matrix Final_result, Matrix mainMat, Matrix RHSMat)

{
    int lwork = 4 * mainMat.rows * mainMat.cols;
    int l;
    char TRANS = 'N';
    int LDB = max(mainMat.cols, mainMat.rows);
    
    Eigen::ArrayXd work(lwork);
    F77_CALL(dgels)(&TRANS, &(mainMat.rows), &(mainMat.cols), &(RHSMat.cols), mainMat.t, &(mainMat.rows), RHSMat.t, &LDB, work.data(), &lwork, &l);
    
    for (int i = 0; i < mainMat.cols; i++)
    {
        for(int j = 0; j < RHSMat.cols; j++)
        {
            M(Final_result, j, i) = M(RHSMat, j, i);
        }
    }
}



Matrix QRdsolve(Matrix mainMat, Matrix RHSMat)

{
    int lwork = 4 * mainMat.rows * mainMat.cols;
    int l;
    char TRANS = 'N';
    int LDB = max(mainMat.cols, mainMat.rows);
    Matrix result, input;
    input = duplicateIt(mainMat);
    result = fill(RHSMat.cols, LDB, (double)0.0);
    for (int i = 0; i < RHSMat.rows; i++)
        for (int j = 0; j < RHSMat.cols; j++)
            M(result, j, i) = M(RHSMat, j, i);
    
    Eigen::ArrayXd work(lwork);
    F77_CALL(dgels)(&TRANS, &(input.rows), &(input.cols), &(result.cols), input.t, &(input.rows), result.t, &LDB, work.data(), &lwork, &l);
    Matrix Final_result = new_matrix(RHSMat.cols, mainMat.cols);
    for (int i = 0; i < mainMat.cols; i++)
    {
        for(int j = 0; j < RHSMat.cols; j++)
        {
            M(Final_result, j, i) = M(result, j, i);
        }
    }
    
    return Final_result;
}

void chol_lpk(Matrix mainMat)
{
    int l;
    char UPLO = 'U';
    
    F77_CALL(dpotrf)(&UPLO, &(mainMat.cols), mainMat.t, &(mainMat.rows), &l);
    
    for(int i = 0; i < mainMat.rows; i++) {
        for(int j = i+1; j < mainMat.cols; j++) {
            M(mainMat, i, j) = 0;
        }
    }
    
}

double solvecond(Matrix inMat)
{
    Matrix result = duplicateIt(inMat);
    int l;
    char JOBZ = 'S';
    int lwork = -1;
    double wkopt;
    int dim_s = std::max(result.cols, result.rows); // maybe min is sufficient
    Eigen::ArrayXi iwork(8 * dim_s);
    Eigen::ArrayXd sv(dim_s);
    Eigen::ArrayXd u(dim_s * result.rows);
    Eigen::ArrayXd vt(dim_s * result.cols);
    F77_CALL(dgesdd)(&JOBZ, &(result.rows), &(result.cols), result.t, &(result.rows), sv.data(), u.data(), &(result.rows), vt.data(), &(result.cols), &wkopt, &lwork, iwork.data(), &l);
    lwork = (int)wkopt;
    Eigen::ArrayXd work(lwork);
    F77_CALL(dgesdd)(&JOBZ, &(result.rows), &(result.cols), result.t, &(result.rows), sv.data(), u.data(), &(result.rows), vt.data(), &(result.cols), work.data(), &lwork, iwork.data(), &l);
    
    if (l < 0) Rf_error("the i-th argument had an illegal value");
    else if (l > 0) Rf_error("DBDSDC did not converge, updating process failed.");
    else
    {
        if ((sv == 0).count()) return std::numeric_limits<double>::infinity();
        else return sv.maxCoeff() / sv.minCoeff();
    }
}

Matrix fillMatrix(int cols, int rows, double* array)
{
    Matrix t = new_matrix(cols, rows);
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++) {
            M(t,j,i)=array[j];
        }
    }
    return t;
}

Matrix MatrixToVector(Matrix mat)
{
    Matrix result = new_matrix(mat.rows*mat.cols, 1);
    
    int ind_hess = 0;
    for (int i = 0; i < mat.cols; i++)
    {
        for (int j = 0; j < mat.rows; j++)
        {
            M(result, ind_hess, 0) = M(mat, i, j);
            ind_hess = ind_hess + 1;
        }
    }
    return result;
}

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

bool all(Matrix x)
{
    for(int i = 0; i < x.cols * x.rows; i++)
    {    if (x.t[i] != 0) return false;}
    return true;
}

void fillMatrix_t(Matrix t, int cols, int rows, double* array)
{
    int i,j;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++) {
            M(t,j,i)=array[j];
        }
    }
}

void transpose_t(Matrix t_t, Matrix t){
    int r,c;
    for ( r = 0; r < t.cols; r++ )
    {
        for ( c = 0; c < t.rows; c++ )
        {
            M(t_t, c, r) = M(t, r, c);
        }
    }
}
