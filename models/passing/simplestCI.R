library(OpenMx)
covariance <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow=2, dimnames=list(c("a", "b"),
                                                                  c("a", "b")))
model <- mxModel("CI Example",
                 mxMatrix(name="expectedCov", "Symm", 2, 2, free=T, values = c(1, .5, 1),
                          labels = c("var1", "cov12", "var2")),
                 mxExpectationNormal("expectedCov", dimnames=c("a", "b")),
                 mxFitFunctionML(),
                 mxData(covariance, "cov", numObs=10000)
)

fit1 <- mxRun(model, silent=TRUE)

fit2 <- mxRun(mxModel(model, mxCI(c("var1", "cov12"))), intervals = TRUE, silent=TRUE)
# cat(deparse(round(model@output$confidenceIntervals, 3)))
omxCheckCloseEnough(c(fit2@output$confidenceIntervals),
                    c(0.973, 0.478, 1.028, 0.522), .01)

omxCheckCloseEnough(fit1@output$fit, fit2@output$fit, 1e-6)
omxCheckCloseEnough(fit1@output$standardErrors, fit2@output$standardErrors, 1e-6)
