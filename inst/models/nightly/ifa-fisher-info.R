# Here is some code to compute the expected Fisher information matrix,
# Equation 9 in Bock & Aitken (1981).
#
# The Appendix of Tian, Cai, Thissen, & Xin (2012) has a good discussion.

library(rpf)
library(OpenMx)
library(numDeriv)

mxOption(NULL, 'loglikelihoodScale', -1)

sample.size <- 1000
numItems <- 5
spec <- list()
spec[1:numItems] <- list(rpf.grm())
param <- structure(c(1.76438, -1.02925, 0.676778, -0.639092, 0.832068,  1.44757, 2.01329, -1.52509, 1.00949, 0.10616),
                   .Dim = c(2L, 5L ), .Dimnames = list(NULL, c("i1", "i2", "i3", "i4", "i5")))

if (0) {
  data.fm <- sapply(rpf.sample(sample.size, spec, param), unclass)-1
  write.table(data.fm, file="info.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ip.mat <- mxMatrix(name="item", nrow=dim(param)[1], ncol=numItems,
                   values=param, free=TRUE)
rownames(ip.mat) <- c('f1', 'b')

compute.gradient <- function(v) {
  pat <- data.frame(t(v))
  for (c in colnames(pat)) {
    pat[[c]] <- ordered(pat[[c]], levels=0:1)
  }
  
  colnames(ip.mat) <- colnames(pat)

  m1 <- mxModel(model="latent", ip.mat,
                mxData(observed=pat, type="raw"),
                mxExpectationBA81(ItemSpec=spec, debugInternal=TRUE),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation', 'scores'),
                  mxComputeOnce('fitfunction', 'gradient'),
				      mxComputeReportDeriv(),
				      mxComputeReportExpectation()
		)))
  m1.fit <- mxRun(m1, silent=TRUE)
  
  ref <- m1.fit$expectation$debug$patternLikelihood
  
  # The gradient is computed in log likelihood units
  # but we need probability units. In log likelihood units,
  # the gradient matches numDeriv to around 1e-8. After
  # converting units, the precision goes way down to
  # around 0.01.
  grad <- exp(ref) - exp(ref + m1.fit$output$gradient)
  
  list(l=exp(ref), g=grad)
}

fisher.info <- function() {
  numPat <- prod(vapply(spec, function(s) s$outcomes, 1))
  patLik <- c()
  grad <- matrix(NA, nrow=numPat, sum(ip.mat$free))
  for (h in 1:numPat) {
    index <- h - 1
    v <- rep(NA, numItems)
    for (i in 1:numItems) {
      v[i] <- index %% spec[[i]]$outcomes
      index <- index %/% spec[[i]]$outcomes
    }
    stuff <- compute.gradient(v)
    patLik[h] <- stuff$l
    grad[h,] <- stuff$g
  }
  list(l=patLik, g=grad)
}

fi <- fisher.info()
omxCheckCloseEnough(sum(fi$l), 1, .01)

se <- sqrt(diag(solve(sample.size * t(fi$g) %*% diag(1/fi$l) %*% fi$g)))
fm.se <- c(.25, .13, .11, .07, .14, .10, .32, .18, .13, .08)
omxCheckCloseEnough(se, fm.se, .04)
