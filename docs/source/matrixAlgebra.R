
require(OpenMx)

algebraExercises <- mxModel(
	mxMatrix(type="Full", values=c(1,2,3), nrow=3, ncol=1, name='A'),
	mxMatrix(type="Full", values=c(1,2,3), nrow=3, ncol=1, name='B'),
	mxAlgebra(A%*%t(A), name='q1'),
	mxAlgebra(t(A)%*%A, name='q2'),
	mxAlgebra(A*A, name='q3'),
	mxAlgebra(A+B, name='q4'))

answers <- mxRun(algebraExercises)
answers@algebras
result <- MxEvaluate(list(q1,q2,q3,q4),answers)
