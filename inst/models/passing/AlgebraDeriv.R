#
#   Copyright 2007-2018 The OpenMx Project
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


#options(error = browser)
require(OpenMx)
require(numDeriv)

set.seed(1)

mat1 <- mxMatrix("Full", rnorm(1), free=TRUE, nrow=1, ncol=1, labels="m1", name="mat1")

mu <- 0
sigma <- 2
Scale <- -2
obj <- mxAlgebra(Scale * -.5 * (log(2*pi) + log(sigma) + (mat1[1,1] - mu)^2/sigma), name = "obj")
grad <- mxAlgebra(Scale * -(mat1[1,1] - mu)/sigma, name = "grad", dimnames=list("m1", NULL))
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
                  mxComputeSequence(list(
                    mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian', 'ihessian')),
                    mxComputeReportDeriv()
                  )))
got <- omxCheckWarning(mxRun(model1, silent=TRUE, useOptimizer=FALSE),
                       "mxRun(..., useOptimizer=FALSE) ignored due to custom compute plan")
omxCheckCloseEnough(got$output$fit, -2 * log(dnorm(got$output$estimate, sd=sqrt(sigma))), 1e-4)
omxCheckCloseEnough(got$output$gradient, 2*(mat1$values[1]-mu)/sigma, 1e-4)
omxCheckCloseEnough(got$output$hessian, 2/sigma)
omxCheckCloseEnough(got$output$ihessian, sigma/2)

numer <- mxModel(model1, mxComputeSequence(list(
  mxComputeNumericDeriv(checkGradient = FALSE),
  mxComputeReportDeriv())))
got <- mxRun(numer, silent=TRUE)
omxCheckCloseEnough(got$output$hessian, 1, 1e-3)

model3 <- mxModel(model1, mxComputeNewtonRaphson())
model3 <- mxRun(model3, silent=TRUE)
omxCheckCloseEnough(model3$output$estimate, 0, 1e-4)
omxCheckCloseEnough(model3$output$status$code, 0)
omxCheckCloseEnough(model3$output$iterations, 2L)

model3 <- mxModel(model3, mxComputeSequence(list(
    mxComputeOnce('fitfunction', 'information', 'hessian'),
    mxComputeReportDeriv())))
model3 <- mxRun(model3)
omxCheckCloseEnough(model3$output$hessian, 1, 1e-3)

#----------------------------------------------------------- code 6

mat1 <- mxMatrix("Full", rnorm(1), free=TRUE, nrow=1, ncol=1, labels="m1", name="mat1")
obj <- mxAlgebra(abs(mat1) + mat1^2, name = "obj")
grad <- mxAlgebra(mat1/abs(mat1) + 2 * mat1, name = "grad", dimnames=list("m1", NULL))
hess <- mxAlgebra(2, name = "hess", dimnames=list("m1", "m1"))
code6 <- mxModel("code6", mat1, obj, grad, hess,
                 mxFitFunctionAlgebra("obj", gradient="grad", hessian="hess"))
m1 <- mxRun(mxModel(code6, mxComputeSequence(list(
  mxComputeNumericDeriv(checkGradient=FALSE),
  mxComputeReportDeriv()
))))
m2 <- mxRun(mxModel(code6, mxComputeSequence(list(
  mxComputeOnce('fitfunction', c('fit','gradient','hessian')),
  mxComputeReportDeriv()
))))
omxCheckCloseEnough(c(m1$output$gradient), c(m2$output$gradient), 1e-6)
omxCheckCloseEnough(c(m1$output$hessian), c(m2$output$hessian), 1e-5)

code6 <- mxRun(mxModel(code6, mxComputeSequence(list(
  mxComputeNewtonRaphson(),
  mxComputeReportDeriv()
))), suppressWarnings = TRUE)

omxCheckEquals(code6$output$status$code, 6)
omxCheckCloseEnough(code6$output$estimate, 0, 1e-6)
omxCheckCloseEnough(abs(code6$output$gradient), 1, 1e-6)

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
mv2.fit <- mxRun(mv2, silent=TRUE, suppressWarnings = TRUE)
omxCheckEquals(mv2.fit$output$status$code, 6)
omxCheckCloseEnough(mv2.fit$output$estimate, rep(0, 2), 1)  # probably in a local minimum
