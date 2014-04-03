#options(error = browser)
require(OpenMx)
require(numDeriv)

set.seed(1)
if (0) {
  
mat1 <- mxMatrix("Full", rnorm(1), free=TRUE, nrow=1, ncol=1, labels="m1", name="mat1")

mu <- 0
sigma <- 2
Scale <- -2
obj <- mxAlgebra(Scale * -.5 * (log(2*pi) + log(sigma) + (mat1[1,1] - mu)^2/sigma), name = "obj")
grad <- mxAlgebra(Scale * -(mat1[1,1] - mu)/sigma, name = "grad", dimnames=list("m1"))
hess <- mxAlgebra(Scale * -1/sigma, name = "hess", dimnames=list("m1", "m1"))

#----------------------------------------------------------- test failures

model2 <- mxModel("model2", mat1, obj, grad, hess,
                  mxFitFunctionAlgebra("obj"))

omxCheckError(mxRun(mxModel(model2, mxComputeOnce('fitfunction', 'gradient')), silent=TRUE),
              "Gradient requested but not available")

omxCheckError(mxRun(mxModel(model2, mxComputeOnce('fitfunction', 'hessian')), silent=TRUE),
              "Hessian requested but not available")

#----------------------------------------------------------- single variable

model1 <- mxModel("model1", mat1, obj, grad, hess,
                  mxFitFunctionAlgebra("obj", gradient="grad", hessian="hess"),
                  mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian', 'ihessian')))
got <- mxRun(model1, silent=TRUE)
omxCheckCloseEnough(got$output$fit, -2 * log(dnorm(got$output$estimate, sd=sqrt(sigma))), 1e-4)
omxCheckCloseEnough(got$output$gradient, 2*(mat1$values[1]-mu)/sigma, 1e-4)
omxCheckCloseEnough(got$output$hessian, 2/sigma)
omxCheckCloseEnough(got$output$ihessian, sigma/2)

numer <- mxModel(model1, mxComputeNumericDeriv())
got <- mxRun(numer, silent=TRUE)
omxCheckCloseEnough(got$output$hessian, 1, 1e-3)

model3 <- mxModel(model1, mxComputeNewtonRaphson(sparseProduct=FALSE))
got <- mxRun(model3, silent=TRUE)
omxCheckCloseEnough(got$output$estimate, 0, 1e-4)
omxCheckCloseEnough(got$compute$output$inform, 0)
omxCheckCloseEnough(got$compute$output$iterations, 2L)
}

#----------------------------------------------------------- multivariate

mat2 <- mxMatrix("Full", c(50,50), free=TRUE, nrow=2, ncol=1, labels=paste("x",1:2,sep=""), name="x")
#x = matrix(c(50,50), ncol=1)

obj <- mxAlgebra(x[1,1]^2 + x[2,1]^2 + sin(x[1,1]+x[2,1]) + x[1,1] - x[2,1], name = "obj")
grad <- mxAlgebra(rbind(2*cos(2*x[1,1]+x[2,1])+2*x[1,1]+1,
                        cos(2*x[1,1]+x[2,1])+4*x[2,1]-1),
                  dimnames=list(paste("x",1:2,sep=""),c()), name = "grad")
hess <- mxAlgebra(rbind(cbind(2-4*sin(2*x[1,1]+x[2,1]), -2*sin(2*x[1,1]+x[2,1])),
                        cbind(-2*sin(2*x[1,1]+x[2,1]), 4-sin(2*x[1,1]+x[2,1]))), name = "hess",
                  dimnames=list(paste("x",1:2,sep=""), paste("x",1:2,sep="")))
mv1 <- mxModel("mv1", mat2, obj, grad, hess,
                  mxFitFunctionAlgebra("obj", gradient="grad", hessian="hess"),
               mxComputeSequence(list(
                 mxComputeOnce('fitfunction', c('gradient', 'hessian', 'ihessian')),
                 mxComputeReportDeriv())))
mv1.fit <- mxRun(mv1, silent=TRUE)
omxCheckCloseEnough(mv1.fit$output$gradient, c(102.4, 199.7), .1)
omxCheckCloseEnough(c(mv1.fit$output$hessian), c(4.86, 1.43, 1.43, 4.71), .01)
omxCheckCloseEnough(mv1.fit$output$hessian, solve(mv1.fit$output$ihessian), 1e-2)

mv2 <- mxModel(mv1, mxComputeNewtonRaphson())
mv2.fit <- mxRun(mv2, silent=TRUE)
omxCheckCloseEnough(mv2.fit$output$estimate, rep(0, 2), 1)  # probably in a local minimum

mat3 <- mat2
mat3$free[1,1] <- FALSE
mv3 <- mxRun(mxModel(mv1, mat3), silent=TRUE)
omxCheckCloseEnough(c(mv3$output$gradient), mv1.fit$output$gradient["x2"])
omxCheckCloseEnough(c(mv3$output$hessian), mv1.fit$output$hessian[2,2])
