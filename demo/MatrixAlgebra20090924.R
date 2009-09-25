require(OpenMx)

# -----------------------------------------------------------------------
# Program: MatrixAlgebra20090924.R  
#  Author: Hermine Maes
#    Date: Wed Sep 24 11:45:52 EDT 2009
#
# Revision History
#   Hermine Maes -- Wed Sep 24 11:45:52 EDT 2009 MatrixAlgebra20090924.R
# -----------------------------------------------------------------------

# Specify Model Matrices and Algebra
# -----------------------------------------------------------------------
algebraExercises <- mxModel(
	mxMatrix(
    	type="Full", 
    	nrow=3, 
    	ncol=1, 
    	values=c(1,2,3), 
    	name='A'
	),
	mxMatrix(
    	type="Full", 
    	nrow=3, 
    	ncol=1, 
    	values=c(1,2,3), 
    	name='B'
	),
	mxAlgebra(
    	A + B, 
    	name='q1'		# addition
	),
	mxAlgebra(
    	A * A, 
    	name='q2'		# dot multiplication
	),
	mxAlgebra(
	    t(A), 
	    name='q3'		# transpose
	),
    mxAlgebra(
        A %*% t(A), 
        name='q4'		# inner product
    ),
    mxAlgebra(
        t(A) %*% A, 
        name='q5'		# outer product
    )
)

# Run Model and Generate Output
# -----------------------------------------------------------------------
answers <- mxRun(algebraExercises)
answers@algebras
result <- mxEval(list(q1,q2,q3,q4,q5),answers)	

