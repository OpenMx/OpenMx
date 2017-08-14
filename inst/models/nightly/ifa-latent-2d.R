# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

compute.factored.ll <- function(m, obs) {
  # Factored form of the complete data likelihood
  # Equation 8 from Cai (2010, p. 37)
  
	grp <- as.IFAgroup(m, obs)
	spec <- grp$spec
	ip <- grp$param
	sc <- EAPscores(grp)
	scores <- as.matrix(sc[,1:2])
  
  ll <- 0
  for (ii in 1:nrow(obs)) {
    for (xx in 1:ncol(obs)) {
      ll <- ll + rpf.logprob(spec[[xx]], ip[,xx], scores[ii,])[obs[ii,xx]]
    }
  }
  ld1 <- apply(scores, 1,
               function (sc) dmvnorm(sc, c(m$mean$values),
                                     m$cov$values, log=TRUE))
  ll <- ll + sum(ld1)
  ll
}

set.seed(1)
numItems <- 128  # need lots of items to get precise factor scores
numPeople <- 100

items <- list()
items[1:numItems] <- list(rpf.grm(factors=2, outcomes=5))
correct.mat <- sapply(items, rpf.rparam, version=1)

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))

ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                   values=correct.mat, free=FALSE)
true.mean <- c(.2,-.1)
true.cov <- matrix(c(.5, -.2, -.2, .5), nrow=2)

data <- rpf.sample(numPeople, items, correct.mat, mean=true.mean, cov=true.cov)
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c(paste('f', 1:2, sep=""), rep('p', nrow(ip.mat) - 2))
  
m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=c(0,0), free=TRUE)
colnames(m.mat) <- paste('f', 1:2, sep="")
cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2),
                    free=TRUE, labels=c("v1","c12","c12","v2"))
dimnames(cov.mat) <- list(paste('f', 1:2, sep=""), paste('f', 1:2, sep=""))

m2 <- mxModel(model="latent",
              ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=items),
              mxFitFunctionML(),
	      mxComputeIterate(steps=list(
				   mxComputeOnce('fitfunction', "fit"),
				   mxComputeOnce('expectation', "latentDistribution", "copy"))))
m2 <- mxRun(m2, silent=TRUE)

if (1) {
  omxCheckCloseEnough(-2 * compute.factored.ll(m2, data) /  m2$output$fit, 1, .02)
}

tri1 <- function(k) k*(k+1)/2

mkfilter <- function(factors, r,c) {
  f <- matrix(0, factors, factors)
  f[r,c] <- 1
  f[c,r] <- 1
  f
}

sigma.coef <- function(Icov) {
  factors <- nrow(Icov)
  f1 <- 1
  coef <- matrix(NA, tri1(factors), tri1(factors))
  for (r1 in 1:factors) {
    for (c1 in 1:r1) {
      f2 <- f1
      for (r2 in r1:factors) {
        cstart <- ifelse(r1==r2, c1, 1)
        for (c2 in cstart:r2) {
          coef[f1,f2] <- sum(diag(Icov %*% mkfilter(factors, r1,c1) %*% Icov %*% mkfilter(factors, r2,c2)))
          coef[f2,f1] <- coef[f1,f2]
          f2 <- f2+1
        }
      }
      f1 <- f1+1
    }
  }
  coef
}

Icov <- solve(m2$cov$values)

dimnames(m2$mean) <- list(NULL,paste("f",1:2,sep=""))
dimnames(m2$cov) <- list(paste("f",1:2,sep=""),paste("f",1:2,sep=""))
mean <- c(m2$mean$values)
names(mean) <- paste("f",1:2,sep="")
cov <- m2$cov$values

m3 <- mxModel("en", m2$mean, m2$cov,
              mxExpectationNormal("cov", "mean"),
              mxData(observed=cov, type="cov", means=mean, numObs=numPeople),
              mxFitFunctionML(),
              mxComputeSequence(list(
                mxComputeOnce('fitfunction', 'information', "hessian"),
                mxComputeReportDeriv())))
m3Fit <- mxRun(m3, silent=TRUE)
omxCheckCloseEnough(m3Fit$output$hessian[1:2,1:2], -2 * -numPeople * Icov, .001)
# not sure if numPeople-2 is correct
omxCheckCloseEnough(m3Fit$output$hessian[3:5,3:5], -2 * (numPeople-2) * -.5 * sigma.coef(Icov), .001)
