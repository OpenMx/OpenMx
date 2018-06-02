library(OpenMx)

data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("G")

m1 <- mxModel(
    "m1",
    mxMatrix("Symm",5,5, name="S"),
    mxMatrix("Full",1,5, free=TRUE, name="M"),
    mxExpectationNormal("S","M", dimnames = paste0('x',1:5)),
    mxFitFunctionML(),
    mxData(observed=cov(demoOneFactor), means=colMeans(demoOneFactor),
	    type="cov", numObs=nrow(demoOneFactor)))
diag(m1$S$free) <- TRUE
diag(m1$S$values) <- 1

m1f <- mxRun(m1)

m2 <- m1
m2$expectation$covariance <- "Sprime"
m2 <- mxModel(m2, mxAlgebra(S, name="Sprime"))
m2f <- mxRun(m2)

m3 <- m1
m3$expectation$covariance <- "Sprime"
m3 <- mxModel(m3, mxMatrix(type = "Symm", nrow = 5, ncol=5,
                          name="Sprime"))
for (rx in 1:5) {
    m3$Sprime$labels[rx,rx] <- paste0('S[',rx,',',rx,']')
}
m3f <- mxRun(m3)

omxCheckCloseEnough(m1f$output$fit, m2f$output$fit, 1e-6)
omxCheckCloseEnough(m1f$output$fit, m3f$output$fit, 1e-6)
