library(OpenMx)
library(testthat)

if(mxOption(key="Default optimizer") == 'CSOLNP') stop("SKIP")

mxOption(key="feasibility tolerance", value = .001)

#mxOption(NULL,"Print level",20)
#mxOption(NULL,"Print file",1)
#mxOption(NULL,"Verify level",3)

powellmod1 <- mxModel(
	"PowellBenchmarkNoJacobians",
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),
           labels=paste("x",1:5,sep=""),name="X"),
	mxAlgebra( exp(prod(X)), name="powellfunc"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,values=NA_real_,name="nana"),
	mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5]),
                   #powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]),
						 name="objgrad",
						 dimnames=list(NULL,paste("x",1:4,sep="")) ),
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1"),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2"),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3"),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad"),
  mxComputeGradientDescent()
)

if (0) {
  powellmod1 <- mxOption(powellmod1,"Always Checkpoint","Yes")
  powellmod1 <- mxOption(powellmod1,"Checkpoint Units","evaluations")
  powellmod1 <- mxOption(powellmod1,"Checkpoint Count",1)
  powellmod1 <- mxOption(powellmod1, "Checkpoint Fullpath", "/dev/fd/2")
}

powellrun1 <- mxRun(powellmod1)

powellmod1 <- mxOption(powellmod1,"Analytic gradients","No")
powellrun2 <- mxRun(powellmod1)

expect_equal(coef(powellrun1), coef(powellrun2), 1e-6)

expect_equal(powellrun2$output$evaluations - powellrun1$output$evaluations,
             148, 10)

#cat(deparse(round(coef(powellrun1), 3)))
expect_equal(coef(powellrun1),
             c(x1 = -1.717, x2 = 1.596, x3 = 1.827, x4 = -0.764, x5 = -0.764 ),
             1e-3)
