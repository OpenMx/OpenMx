# Demonstrate how to convert the derivative from log likelihood
# units to probability units.

library(rpf)
library(OpenMx)
library(numDeriv)

numItems <- 5
spec <- list()
spec[1:numItems] <- rpf.grm()
param <- sapply(spec, rpf.rparam)

pat <- data.frame(t(round(runif(numItems))))
for (c in colnames(pat)) {
  pat[[c]] <- ordered(pat[[c]], levels=0:1)
}

if (0) {
  data.fm <- sapply(rpf.sample(1000, spec, param), unclass)-1
  write.table(data.fm, file="info.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ip.mat <- mxMatrix(name="ItemParam", nrow=dim(param)[1], ncol=numItems,
                   values=param, free=FALSE)

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

objective <- function(x) {
  ip.mat@values[,] <- x
  m1 <- mxModel(model="latent",
                ip.mat, m.mat, cov.mat,
                mxData(observed=pat, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=spec,
                                  ItemParam="ItemParam"),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', fit=TRUE))))
  m1.fit <- mxRun(m1, silent=TRUE)
  lik <- exp(m1.fit@output$minimum/-2)
  lik
}

deriv <- grad(objective, c(param), method="simple")

m1 <- mxModel(model="latent",
              ip.mat, m.mat, cov.mat,
              mxData(observed=pat, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=spec,
                                ItemParam="ItemParam"),
              mxFitFunctionML(),
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', context="EM"),
                mxComputeOnce('fitfunction', gradient=TRUE))))
m1@matrices$ItemParam@free[,] <- TRUE
m1.fit <- mxRun(m1, silent=TRUE)

ref <- m1.fit@expectation@patternLikelihood
grad <- exp(ref) - exp(ref + m1.fit@output$gradient)
omxCheckCloseEnough(grad, deriv, .01)   # poor accuracy
if (0) {
  max(abs(grad- deriv))
}
