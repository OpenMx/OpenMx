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


library(OpenMx)
#mxOption(NULL, "Default optimizer", "NPSOL")

covariance <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow=2, dimnames=list(c("a", "b"),
                                                                  c("a", "b")))
means <- c(-1,.5)
names(means) <- c('a','b')

model <- mxModel("CIExample",
                 mxMatrix(name="expectedCov", "Symm", 2, 2, free=T, values = c(1, 0, 1),
                          labels = c("var1", "cov12", "var2")),
                 mxMatrix(name="expectedMean", "Full", 1, 2, free=T, labels=c('m1','m2')),
                 mxExpectationNormal("expectedCov", "expectedMean", dimnames=c("a", "b")),
                 mxFitFunctionML(),
                 mxData(covariance, "cov", means, numObs=10000)
)
diag(model$expectedCov$lbound) <- .1

model <- mxOption(model,"Checkpoint Units",'iterations')
model <- mxOption(model,"Checkpoint Count",1)

fit1 <- mxRun(model, silent=TRUE)

if (mxOption(NULL, 'Default optimizer') != "SLSQP") {ctype = 'none'} else {ctype = 'ineq'}

cimodel <- mxModel(fit1,
                   mxCI("var1", type="lower", boundAdj=FALSE),
                   mxCI("cov12", type="upper"),
                   mxCI("m1", type="both"),
                   mxComputeConfidenceInterval(verbose=0,plan=mxComputeGradientDescent(verbose=0), constraintType = ctype))

fit2 <- mxRun(cimodel,
              intervals = TRUE, silent=TRUE, checkpoint=FALSE)
print(fit2$compute$output$detail)

fit3 <- mxRun(cimodel, intervals = FALSE)
omxCheckTrue(is.null(fit3$output$confidenceIntervals))

# For multivariate normal means, SEs match likelihood-based CIs
omxCheckCloseEnough(fit2$output$estimate['m1'] + fit1$output$standardErrors['m1',] * qnorm(.025),
                    fit2$output$confidenceIntervals['m1', 'lbound'], .0001)
omxCheckCloseEnough(fit2$output$estimate['m1'] - fit1$output$standardErrors['m1',] * qnorm(.025),
                    fit2$output$confidenceIntervals['m1', 'ubound'], .0001)

# cat(deparse(round(model$output$confidenceIntervals, 3)))
omxCheckCloseEnough(fit2$output$confidenceIntervals['var1','lbound'], c(0.9727), .001)
omxCheckCloseEnough(fit2$output$confidenceIntervals['cov12','ubound'], c(0.522), .001)

omxCheckCloseEnough(fit1$output$fit, fit2$output$fit, 1e-6)

fit4 <- omxCheckWarning(mxRun(mxModel(cimodel, mxCI('expectedMean[1,1]', interval=runif(1,.9,.95))),
			      intervals = TRUE, silent=TRUE, checkpoint=FALSE),
			"Different confidence intervals 'CIExample.expectedMean[1,1]' and 'm1' refer to the same thing")

# ensure the [1,] syntax is supported
data(demoOneFactor)
factorModel <- mxModel("One Factor",
      mxMatrix("Full", 5, 1, values=0.2, lbound=0, ubound=5,
           free=TRUE, name="A"),
      mxMatrix("Symm", 1, 1, values=1,
           free=FALSE, name="L"),
      mxMatrix("Diag", 5, 5, values=1,
           free=TRUE, name="U"),
      mxAlgebra(A %*% L %*% t(A) + U, name="R"),  mxCI("A[1,]"),
      mxExpectationNormal("R", dimnames = names(demoOneFactor)),
      mxFitFunctionML(),
      mxData(cov(demoOneFactor), type="cov", numObs=500))
factorModel <- mxRun(factorModel, intervals=TRUE)
if (0) {
  # NPSOL can't find 'em
  factorModel <- mxRun(mxModel(factorModel,
                               mxComputeConfidenceInterval(verbose=3, constraintType = "none",
                                                           plan=mxComputeGradientDescent())))
  print(factorModel$compute$output)
}
ci <- factorModel$output$confidenceIntervals
omxCheckEquals(nrow(ci), 1)
omxCheckCloseEnough(c(0.397), ci[1,'estimate'], .001)
if (mxOption(NULL, "Default optimizer") != "NPSOL") {
  omxCheckCloseEnough(c(0.3679), ci[1,'lbound'], .02)
  omxCheckCloseEnough(c(0.429), ci[1,'ubound'], .02)
}
