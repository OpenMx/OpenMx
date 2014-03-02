library(OpenMx)
library(rpf)
library(mvtnorm)

q()

pcov <- function(v) {
  m <- apply(v, 2, mean)
  v %*% t(v)
  (t(v) %*% v) / nrow(v) - m %*% t(m)
}

perDim <- 10
discr <- expand.grid(a1=seq(0,8,length.out=perDim), a2=seq(0,8,length.out=perDim))
discr <- discr[apply(discr, 1, max) > 4,]
if (nrow(discr) %% 2) {
  discr <- discr[-1,]
}

outcomes <- 10
spec <- list()
spec[1:nrow(discr)] <- rpf.nrm(outcomes=outcomes, factors=2)

ip <- matrix(0, 2 + (outcomes-1)*2, nrow(discr))
rownames(ip) <- names(rpf.rparam(spec[[1]]))
ip[1:2,] <- t(as.matrix(discr))
ip["alf1",] <- 1.5
ip["gam1",] <- seq(-3,3, length.out=nrow(discr)/2)
ip['gam2',] <- 20

plot.icc <- function(item, param, basis, width=3) {
  require(ggplot2)
  require(reshape2)
  pm <- t(rpf.prob(item, param, kronecker(as.matrix(basis), t(seq(-width, width, .1)))))
  icc <- as.data.frame(melt(pm, varnames=c("theta",'category')))
  icc$theta <- seq(-width, width, .1)
  icc$category <- as.factor(icc$category - 1)
  ggplot(icc, aes(theta, value)) +
    geom_line(aes(color=category, linetype=category)) +
    ylim(0,1) + xlim(-width,width) + labs(y="Probability", x="Theta")
}
#plot.icc(spec[[1]], ip[,5], c(1,0))

numPeople <- 35
g1.mean <- c(-.5,.8)
g1.cov <- matrix(c(1.8,.5,.5,2),nrow=2)
g2.mean <- c(-.1,-.8)
g2.cov <- matrix(c(.9,-.5,-.5,1.2),nrow=2)

if (0) {
  
  nn <- 2000
  svar1 <- matrix(0, 2, 2)
  pcov1 <- matrix(0, 2, 2)
  for (i in 1:nn) {
    scores <- rmvnorm(numPeople, g1.mean, g1.cov)
    svar1 <- svar1 + var(scores)
    pcov1 <- pcov1 + pcov(scores)
  }
  svar1/nn - g1.cov
  pcov1/nn - g1.cov
}

mkmodel <- function(name, theta, sampleCov) {
  dat <- rpf.sample(theta, spec, ip)
  
  ip.mat <- mxMatrix(name="ItemParam", values=ip, free=FALSE)
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=TRUE)
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=TRUE)
  cov.mat@labels[2,1] <- cov.mat@labels[1,2] <- paste(name,"cov",sep="_")
  
  m1 <- mxModel(name, ip.mat, m.mat, cov.mat,
                mxData(observed=dat, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov", sampleCov=sampleCov,
                                  ItemSpec=spec, ItemParam="ItemParam",
                                  qpoints=31, qwidth=5),
                mxFitFunctionML())
  m1
}

mkgroup <- function(seed, sampleCov) {
  set.seed(seed)
  g1 <- rmvnorm(numPeople, mean=g1.mean, sigma=g1.cov)
  m1 <- mkmodel("g1", t(g1), sampleCov)
  g2 <- rmvnorm(numPeople, mean=g2.mean, sigma=g2.cov)
  m2 <- mkmodel("g2", t(g2), sampleCov)
  groups <- paste("g", 1:2, sep="")
  mxModel("grp", m1, m2,
          mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
          mxComputeSequence(list(
            mxComputeOnce(paste(groups, 'expectation', sep='.')),
            mxComputeIterate(list(
              mxComputeOnce('fitfunction', 'fit',
                            free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.'))),
              tolerance=1e-6))))
}

fit1 <- function(seed) {
  set.seed(seed)
  g1 <- rmvnorm(numPeople, mean=g1.mean, sigma=g1.cov)
  m1 <- mkmodel("g1", t(g1), FALSE)
  g2 <- rmvnorm(numPeople, mean=g2.mean, sigma=g2.cov)
  m2 <- mkmodel("g2", t(g2), FALSE)
  groups <- paste("g", 1:2, sep="")
  gm <- mxModel("grp", m1, m2,
          mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
          mxComputeSequence(list(
            mxComputeOnce(paste(groups, 'expectation', sep='.')),
            mxComputeIterate(list(
              mxComputeOnce('fitfunction', 'fit',
                            free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.'))),
              tolerance=1e-6))))
  gm.fit <- mxRun(gm, silent=TRUE)
  
  g1.cov <- gm.fit@submodels$g1@matrices$cov@values
  g2.cov <- gm.fit@submodels$g2@matrices$cov@values
  
  list(got=cbind(g1.samp=c(g1.cov - cov(g1)),
                 g1.pop=c(g1.cov - pcov(g1)),
                 g2.samp=c(g2.cov - cov(g2)),
                 g2.pop=c(g2.cov - pcov(g2))))
}

correct <- c(g1.mean, g1.cov[1:2,1], g1.cov[2,2], g2.mean, g2.cov[1:2,1], g2.cov[2,2])
isCov <- c(rep(F,2), rep(T,3))

bank <- list()
for (rep in 1:500) {
  bank[[rep]] <- fit1(rep)

  report <- rep
  for (type in c("g1.samp","g1.pop","g2.samp","g2.pop")) {
    est <- sapply(bank, function (t) t$got[,type])
    report <- c(report, norm(apply(est, 1, mean), "2"))
  }
  names(report) <- c("rep","g1.samp","g1.pop","g2.samp","g2.pop")
  print(report)
}

