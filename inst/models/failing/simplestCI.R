library(OpenMx)
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

model <- mxOption(model,"Checkpoint Units",'iterations')
model <- mxOption(model,"Checkpoint Count",1)

fit1 <- mxRun(model, silent=TRUE)

mle <- fit1$output$fit

var1L <- mxModel("var1L", model,
                 mxAlgebra(CIExample.expectedCov[1,1], "param"),
                 mxConstraint(0 > CIExample.fitfunction - (mle + 3.84)),
                 mxFitFunctionAlgebra("param"))
var1L <- mxRun(var1L)
omxCheckCloseEnough(var1L$output$estimate['var1'], c(0.973), .01)
#mxOption(NULL, "Default optimizer", "CSOLNP")

cimodel <- mxModel(model,
                   mxCI("var1", type="lower"),
                   mxCI("cov12", type="upper"),
                   mxCI("m1", type="both"))
fit2 <- mxRun(cimodel,
              intervals = TRUE, silent=TRUE, checkpoint=FALSE)

# For multivariate normal means, SEs match likelihood-based CIs
omxCheckCloseEnough(fit2$output$estimate['m1'] + fit2$output$standardErrors['m1',] * qnorm(.025),
                    fit2$output$confidenceIntervals['m1', 'lbound'], .0001)
omxCheckCloseEnough(fit2$output$estimate['m1'] - fit2$output$standardErrors['m1',] * qnorm(.025),
                    fit2$output$confidenceIntervals['m1', 'ubound'], .0001)

# cat(deparse(round(model$output$confidenceIntervals, 3)))
omxCheckCloseEnough(fit2$output$confidenceIntervals['var1','lbound'], c(0.973), .01)
omxCheckCloseEnough(fit2$output$confidenceIntervals['cov12','ubound'], c(0.522), .01)

omxCheckCloseEnough(fit1$output$fit, fit2$output$fit, 1e-6)
omxCheckCloseEnough(fit1$output$standardErrors, fit2$output$standardErrors, 1e-6)

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
factorModel <- mxRun(factorModel, intervals=T)
omxCheckEquals(nrow(factorModel$output$confidenceIntervals), 1)
omxCheckCloseEnough(c(0.368, 0.397, 0.429),
                    factorModel$output$confidenceIntervals[1,], .01)
