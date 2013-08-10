#include <math.h>
#include <float.h>               
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "matrix.h"

double rnd_double() { return (double)1.0; }

Matrix new_matrix(int cols,int rows)
{
	Matrix t;
	t.sentinel1 = 0xdeadbeef;
	t.rows=rows;
	t.cols=cols;
	t.isColMajor = 1;
	t.t=(double *)malloc(sizeof(double)*cols*rows);
	t.sentinel2 = 0xdeadbeef;
	
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++) {
			M(t,j,i)=rnd_double();
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

Matrix getRow( Matrix t, int row){
	printf("t.cols is: \n");
    printf("%d", t.cols); putchar('\n');
    printf("t,0,0 is: \n");
    printf("%2f", M(t,0,0)); putchar('\n');
    

	Matrix toReturn = fill(t.cols, 1, (double)0.0);
    printf("toReturn is: \n");
    printf("%2f", M(toReturn, 0, 0)); putchar('\n');

	int i;
	for (i=0;i < t.cols; i++){
		M(toReturn,i,0) = M(t,i,row);
       // printf("toReturn is: \n");
       // print(toReturn); putchar('\n');
	}
	
	return toReturn;
}

Matrix setRow( Matrix x, int row,  Matrix y){
	
	Matrix toReturn = duplicateIt(x);
	
	int i,j;
	
	for (i=0;i < x.cols; i++){
		M(toReturn,i,row) = M(y,i,0);
	}
	
	return toReturn;
}

Matrix getColumn( Matrix t, int colNum)
{	
   int r, c;
   int index = 0;
    Matrix result = fill(t.rows, 1, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
       for ( c = 0; c < t.cols; ++c )
       {
	        if (c==colNum){
				M(result, index, 0) = M(t, c, r);
				index++; 
			}
       }
   }
   return result;
}

Matrix setColumn( Matrix x,  Matrix y, int colNum)
{
   int r, c;
   int index = 0;
    Matrix result = fill(x.cols, x.rows, (double)0.0);
   for ( r = 0; r < x.rows; ++r )
   {
      for ( c = 0; c < x.cols; ++c )
      {
		if (c==colNum){
			M(result, c, r) = M(y, index, 0);
			index++; 
		}
		else{
			M(result, c, r) = M(x, c, r);
		}
      }
   }
   return result;
}

void print(Matrix t) {
	int i,j;
	for(i=0;i<t.rows;i++) {
		printf("| ");
		for(j=0;j<t.cols;j++) 
			printf("%.16f ",M(t,j,i));
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
		
		temp = getRow(A,p); A = setRow(A,p, getRow(A,max)); 
		A = setRow(A,max, temp);
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
	    Matrix result = fill(A.cols, A.cols, (double)0.0);
	   for (i=0; i<A.cols; i++){
			for (j=0; j<=A.cols; j++){
				 if (i==j)
					M(result, j, i) = M(A, j, 0);
			}
	   }		
	   return result;
}

//Mahsa: this was wrong. it was "value <= M(t,j,i)", and "j<=t.cols". both are wrong.
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

   for ( i = 0; i < t.rows; ++i )
   {
      for ( j = 0; j < t.cols; ++j )
      {
		result += (M(t,j,i) * M(t,j,i));
      }
	}
	
	return sqrt(result);
}

double findMin(Matrix t){
   int i,j;
   double min = DBL_MAX;
   for ( i = 0; i < t.rows; ++i )
   {
      for ( j = 0; j < t.cols; ++j )
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

   for ( i = 0; i < t.rows; ++i )
   {
      for ( j = 0; j < t.cols; ++j )
      {
         if (max < M(t,j,i) ){
			max = M(t,j,i);
             //printf("max is:"); putchar('\n');
             //printf("%2f", max); putchar('\n');
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
   for ( r = 0; r < a.rows; ++r )
   {
       
	    M(result, r, 0) = dotProduct(getRow(a, r), getRow(b, 0));
   }
   return result;	
}

Matrix minMaxAbs(Matrix t, double tol)
{
   // for each x[i]: min( max( abs(x[i]), tol ), 1/tol )	
   int r,c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
      {
		double buf = ourAbs(M(t,c,r));
		buf = max(buf, tol);
		buf = min(buf, 1.0/(tol));
		M(result,c,r) = buf;
      }
   }
   return result;
}

Matrix add(Matrix x,  Matrix y)
{
    Matrix result = fill(x.cols, x.rows, (double)0.0);
   int r,c;
   for ( r = 0; r < x.rows; ++r )
   {
      for ( c = 0; c < x.cols; ++c )
      {
		  M(result, c, r) = M(x, c, r) + M(y, c, r);
      }
   }
   return result;
}

Matrix subtract(Matrix x,  Matrix y)
{
    Matrix result = fill(x.cols, x.rows, (double)0.0);
   int r,c;
   for ( r = 0; r < x.rows; ++r )
   {
      for ( c = 0; c < x.cols; ++c )
      {
		  M(result, c, r) = M(x, c, r) - M(y, c, r);
      }
   }
   return result;
}

Matrix multiply(Matrix x,  Matrix y)
{
    Matrix result = fill(x.cols, x.rows, (double)0.0);
   int r,c;
   for ( r = 0; r < x.rows; ++r )
   {
      for ( c = 0; c < x.cols; ++c )
      {
		  M(result, c, r) = M(x, c, r) * M(y, c, r);
      }
   }
   return result;
}

Matrix divide(Matrix x,  Matrix y)
{
    Matrix result = fill(x.cols, x.rows, (double)0.0);
   int r,c;
   for ( r = 0; r < x.rows; ++r )
   {
      for ( c = 0; c < x.cols; ++c )
      {
		  M(result, c, r) = M(x, c, r) / M(y, c, r);
      }
   }
   return result;
}

Matrix copyInto(Matrix x,  Matrix y, int rowNum, int colStart, int colStop){
	
	   int r, c, index;
	   index = 0;
	    Matrix result = fill(x.cols, x.rows, (double)0.0);
	   for ( r = 0; r < x.rows; ++r )
	   {
	      for ( c = 0; c < x.cols; ++c )
	      {
			if (r == rowNum){
				if (c >= colStart && c <= colStop){
					M(result, c, r) = M(y, index, 0);
					index++;
				}
				else{
					M(result, c, r) = M(x, c, r);
				}
			}
			else{
					M(result, c, r) = M(x, c, r);
			}  
	      }
	   }
	   return result;
}

Matrix rowWiseMin(Matrix t)
{
    Matrix mins = fill(t.rows, 1, (double)0.0);
   int r,c;
   for ( r = 0; r < t.rows; ++r )
   {
      double min = DBL_MAX;
      for ( c = 0; c < t.cols; ++c )
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
   for ( r = 0; r < t.cols; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
      {
			M(result,c,r) = M(t,c,0);
      }
   }
    //printf("result: transpose2D: \n");
    //print(result); putchar('\n');
   for ( r = 0; r < t.cols; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
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
	
    for ( r = 0; r < minDim; ++r )
    {
        for ( c = 0; c < minDim; ++c )
        {
			M(result,c,r) = dotProduct(getRow(t,r), getRow(t,c));
		}
	}
  	return result;
}

Matrix transposeDP(Matrix t){ //Mahsa: same as transpose2D function
    Matrix toReturn = fill(t.cols, t.cols, (double)0.0);
   int r, c;
   for ( r = 0; r < t.cols; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
	  {
		  M(toReturn, c, r) = M(t, r, 0) * M(t, c, 0);
	  }
   }
   return toReturn;
}

Matrix transpose(Matrix t){ // Mahsa: simply transposes the matrix
	 Matrix result = fill(t.rows, t.cols, (double)0.0);
	int r,c;
	for ( r = 0; r < t.cols; ++r )
    {
		for ( c = 0; c < t.rows; ++c )
         {
			M(result, c, r) = M(t, r, c);
         }
    }
	return result;
}

Matrix negate(Matrix t)
{
   int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
        for ( c = 0; c < t.cols; ++c )
        {
			M(result, c, r) = (M(t, c, r) * -1.0);
        }
   }
   return result;
}

Matrix duplicateIt(Matrix t)
{
   int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
        for ( c = 0; c < t.cols; ++c )
        {
			M(result, c, r) = M(t, c, r);
        }
   }
   return result;
}

Matrix matrixAbs(Matrix t)
{
   int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
        for ( c = 0; c < t.cols; ++c )
        {
			M(result, c, r) = ourAbs(M(t, c, r));
        }
   }
   return result;
}

Matrix multiplyByScalar2D(Matrix t, double multiplier)
{
   int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
        for ( c = 0; c < t.cols; ++c )
        {
			M(result, c, r) = M(t, c, r)*multiplier;
        }
   }
   return result;
}

Matrix divideByScalar2D(Matrix t, double divisor)
{
   int r, c;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
        for ( c = 0; c < t.cols; ++c )
        {
			M(result, c, r) = M(t, c, r)/divisor;
        }
   }
   return result;
}

Matrix checkControlList(Matrix t){
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   int length = t.cols;
   int c;   

   M(result, 0, 0)=1;
   M(result, 1, 0)=400;
   M(result, 2, 0)=800;
   M(result, 3, 0)=1.0e-7;
   M(result, 4, 0)=1.0e-8;
   M(result, 5, 0)=1;
   
	for ( c = 0; c < length; ++c )
    {
		M(result, c, 0) = M(t, c, 0);
    }
	return result;
}

Matrix subset(Matrix t, int row, int colStart, int colStop)
{

   int r, c;
   int count = (colStop - colStart)+1;
   int index = 0;
   
    Matrix result = fill(count, 1, (double)0.0);	
   for ( r = 0; r < t.rows; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
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


Matrix copy(Matrix x,  Matrix y){
	
	int totalRows = x.rows;
    int totalCols = x.cols + y.cols;
    Matrix result = fill(totalCols, totalRows, (double)0.0);
	
	int r, c;
	for ( r = 0; r < totalRows; ++r )
     {
        for ( c = 0; c < totalCols; ++c )
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

Matrix copyThree(Matrix a,  Matrix b,  Matrix c){
    Matrix result = copy(a, b);
    Matrix result_three = copy(result, c);
}

Matrix copyFive(Matrix a,  Matrix b,  Matrix c,  Matrix d,  Matrix e){
	copy(copy(copy(a,b), copy(c,d)), e);
}

Matrix timess(Matrix a,  Matrix b){
	printf("times is: \n");
	int i, j, k;
	 Matrix result = fill(b.cols, a.rows, (double)0.0);
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
	printf("result is: \n");
	print(result);
	return result;
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

Matrix rowSort(Matrix t)
{
   int r, c, i, j;
    Matrix result = fill(t.cols, t.rows, (double)0.0);
   for ( r = 0; r < t.rows; ++r )
   {
      for ( c = 0; c < t.cols; ++c )
      {
		M(result, c, r) = M(t, c, r);
      }
      	//for ( r = 0; r < t.rows; ++r )
      	//{
         	for ( i = 0; i < t.cols; ++i )
         	{
			 	for ( j = 0; j < t.cols; ++j )
	         	{
					if (M(result, i, r) < M(result, j, r)){
						double a = M(result, i, r);
						M(result, i, r) = M(result, j, r);
						M(result, j, r) = a;
					}
		     	}
         	//}
      	}
   }
   return result;
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
    //result = subset(result, 0, 0, 1);
    return result;
}


Matrix MatrixInvert(Matrix inMat)
{
	Matrix result;
	
	int lwork = 4 * inMat.rows * inMat.cols;
	int l = 0;

	int*    ipiv = (int*) malloc(inMat.rows * sizeof(int));
	double* work = (double*) malloc(lwork * sizeof(double));

	result = duplicateIt(inMat);
	F77_CALL(dgetrf)(&(result.cols), &(result.rows), result.t, &(result.rows), ipiv, &l);
	if(l != 0) {
		error("Attempted to invert non-invertable matrix.");
	} else {
		F77_CALL(dgetri)(&(result.cols), result.t, &(result.rows), ipiv, work, &lwork, &l);
	}

	free(ipiv);
	free(work);
	return result;
}

int this_main(int argc,char *argv[]) {
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
}


