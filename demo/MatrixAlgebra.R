#
#   Copyright 2007-2010 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# -----------------------------------------------------------------------
# Program: MatrixAlgebra.R  
#  Author: Hermine Maes
#    Date: 08 01 2009 
#
# Matrix Algebra in OpenMx: Basic matrix algebra operations
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

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

