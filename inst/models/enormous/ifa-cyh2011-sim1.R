# This is a replication of Cai, Yang, & Hansen (2011) simulation study #1.

#options(error = browser)
library(OpenMx)
library(rpf)
library(mvtnorm)

mk.model <- function(name, numItems, latent.free) {
  spec <- list()
  dims <- (1 + numItems/4)
  spec[1:numItems] <- rpf.grm(factors = dims)
  
  imat <- mxMatrix(name="item", nrow=rpf.numParam(spec[[1]]), ncol=numItems,
                   values=0, free=FALSE, dimnames=list(names(rpf.rparam(spec[[1]])),
                                                       paste("i", 1:numItems, sep="")))
  
  for (ix in 1:numItems) {
    sg <- (ix-1) %/% 4
    for (px in c(1,2+sg,nrow(imat))) {
      if (px < nrow(imat)) imat$values[px,ix] <- 1
      imat$free[px,ix] <- TRUE
      imat$labels[px,ix] <- paste(rownames(imat)[px], colnames(imat)[ix], sep=",")
    }
  }

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  colnames(m.mat) <- paste("a", 1:dims, sep="")
  cov.mat.free <- FALSE
  if (latent.free) {
    cov.mat.free <- diag(dims)==1
  }
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=cov.mat.free, lbound=1e-2,
                      dimnames=list(colnames(m.mat), colnames(m.mat)))
  
  lname <- paste(name, "latent", sep="")
  latent <- mxModel(lname, m.mat, cov.mat,
		    mxDataDynamic("cov", expectation=paste(name, "expectation", sep=".")),
		    mxExpectationNormal(covariance="cov", means="mean"),
		    mxFitFunctionML())

  m1 <- mxModel(model=name, imat,
                mxExpectationBA81(spec,
                  mean=paste(lname, "mean",sep="."),
                  cov=paste(lname, "cov", sep="."),
                  qpoints=21, qwidth=5),
                mxFitFunctionML())
  
  list(ifa=m1, latent=latent)
}

g2.mean <- c(1, -.5, 0, .5)
g2.cov <- diag(c(.8, 1.2, 1.5, 1))
correct <- matrix(0, 6, 16)
correct[1,] <- c(1, 1.4, 1.7, 2, 1.4, 1.7, 2, 1, 1.7, 2, 1, 1.4, 2, 1, 1.4, 1.7)
correct[2,1:4] <- c(.8, 1.5, 1.2, 1)
correct[3,5:8] <- c(1, .8, 1.5, 1.2)
correct[4,9:12] <- c(1.2, 1, .8, 1.5)
correct[5,13:16] <- c(1.5, 1.2, 1, .8)
correct[6,] <- rep(c(1, .25, -.25, -1), 4)
groups <- paste("g", 1:2, sep="")

modelGen <- function() {
  g1 <- mk.model("g1", 16, FALSE)
  g2 <- mk.model("g2", 12, TRUE)
  
  data.g1 <- rpf.sample(1000, g1$ifa$expectation$ItemSpec, correct)
  data.g2 <- rpf.sample(1000, g2$ifa$expectation$ItemSpec, correct[-5,1:12],
                        mean=g2.mean, cov=g2.cov)
  
  g1$ifa <- mxModel(g1$ifa, mxData(observed=data.g1, type="raw"))
  g2$ifa <- mxModel(g2$ifa, mxData(observed=data.g2, type="raw"))
  
  latent <- mxModel("latent", mxFitFunctionMultigroup("g2latent.fitfunction"))
  
  if (0) {
    # groups are swapped in my flexmirt specs
    write.table(sapply(data.g1, unclass)-1, file=paste("cyh",seed,"g2.csv",sep="-"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(sapply(data.g2, unclass)-1, file=paste("cyh",seed,"g1.csv",sep="-"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  
  omxIFAComputePlan <- function(groups) {
    latent.plan <- NULL
    latentFG <- apply(expand.grid(paste(groups,"latent",sep=""), c('mean','cov')), 1, paste, collapse='.')
    if (0) {
      latent.plan <- mxComputeSequence(list(mxComputeOnce(paste(groups, 'expectation', sep='.')),
                                            mxComputeOnce(paste(groups, 'expectation', sep='.'),
                                                          "latentDistribution", "copy"),  # c('mean','covariance')
                                            mxComputeOnce('fitfunction', "set-starting-values")),
                                       freeSet=latentFG)
    } else {
      latent.plan <- mxComputeNewtonRaphson(latentFG, fitfunction="latent.fitfunction")
    }
    
      mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                  mxComputeSequence(list(
                    mxComputeNewtonRaphson(freeSet=paste(groups, 'item', sep=".")),
                    latent.plan)))
  }
  
  grpModel <- mxModel(model="groupModel", g1, g2, latent,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))
}

# ----------------------------------------------------------------------------

if (file.exists("models/enormous/lib/stderrlib.R")) {
  source("models/enormous/lib/stderrlib.R")
} else if (file.exists("lib/stderrlib.R")) {
  source("lib/stderrlib.R")
} else {
  stop("Cannot find stderrlib.R")
}

#got <- MCphase(modelGen, reps=5, verbose=TRUE, checkCondNum=FALSE)

name <- "cyh2011-sim1"
getMCdata(name, modelGen, c(correct[correct != 0], g2.mean, diag(g2.cov)),
          maxCondNum=5000, recompute=FALSE)

omxCheckCloseEnough(norm(mcBias, "2"), .14194, .001)
omxCheckCloseEnough(max(abs(mcBias)), .0552, .001)
omxCheckCloseEnough(log(det(mcHessian)), 281.37, .1)

if (0) {
  set.seed(500-484)
  model <- modelGen()
  em <- model$compute
  fitfun <- c()
  if (is(em$mstep, "MxComputeSequence")) {
    fitfun <- sapply(em$mstep$steps, function(step) step$fitfunction)
  } else {
    fitfun <- em$mstep$fitfunction
  }
  em$accel <- ""
  em$tolerance <- 1e-11
  em$maxIter <- 750L
#  em$verbose <- 2L
  em$information <- "mr1991"
  em$infoArgs <- list(fitfunction=fitfun, semMethod='mr', semTolerance=sqrt(1e-6))
  plan <- mxComputeSequence(list(
    em,
    mxComputeStandardError(),
    mxComputeReportDeriv()
  ))
  model$compute <- plan
  fit <- try(mxRun(model, silent=TRUE, suppressWarnings=TRUE), silent=TRUE)
  
}

rda <- paste(name, "-result.rda", sep="")
load(rda)

detail <- testPhase(modelGen, 500, methods=c('re', 'estepH', 'mr', 'tian', 'agile', 'meat', 'oakes'))

if (0) {
  asem <- studyASEM(modelGen)
  smooth <- checkSmoothness(modelGen)
}

save(detail, asem, smooth, file=rda)

if (0) {
  condnum <- detail['condnum','agile',]
  head(order(condnum))
  seed <- 101
  set.seed(seed)
  m1 <- modelGen()
  f1 <- paste("cyh-", seed, "-g1.csv", sep="")
  f2 <- paste("cyh-", seed, "-g2.csv", sep="")
  write.table(sapply(m1$g1$data$observed, unclass)-1, f1, quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(m1$g2$data$observed, unclass)-1, f2, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  covname <- paste("~/ifa/cai2009/cyh-", seed, "-cov.txt", sep="")
  H <- read.csv(covname, header = FALSE, row.names=NULL)
  H <- H[,-ncol(H)]
  H <- solve(H)
  labels <- m1$g1$item$labels
  fmOrder <- c(labels[6,13:16], labels[1,13:16], labels[5,13:16],
               labels[1,1:12], labels[2,1:4], labels[3,5:8], labels[4,9:12],
               labels[6,1:12],
               "g2latent.mean[1,1]",  "g2latent.mean[1,2]", "g2latent.mean[1,3]", "g2latent.mean[1,4]",
               "g2latent.cov[1,1]", "g2latent.cov[2,2]", "g2latent.cov[3,3]",  "g2latent.cov[4,4]")
  fmPerm <- match(colnames(mcHessian), fmOrder)
  H <- H[fmPerm, fmPerm]
  H <- 2 * H
  mvt_KL_D(H, solve(H))
  norm((sqrt(2 * diag(solve(H))) - mcSE) / mcSE, "2")
}

if (0) {
  # maybe these are too precise to reproduce?
  EMcycles <- fivenum(sapply(bank[estmask], function (t) t$em$EMcycles))
  # cat(deparse(EMcycles))
  #omxCheckCloseEnough(EMcycles, c(104, 144, 177, 222, 476), 8)  # maybe check against sum TODO
  
  semProbeCount <- fivenum(sapply(bank[estmask], function (t) t$em$semProbeCount))
  # cat(deparse(semProbeCount))
  #omxCheckCloseEnough(semProbeCount, c(126, 137, 140.5, 144, 157), 2)
  omxCheckTrue(all(semProbeCount == 168))
  
  omxCheckCloseEnough(norm(emp$bias, "2"), 0.14211, .01)
  omxCheckCloseEnough(norm(emp$se, "2"), 1.26, .015)
  omxCheckCloseEnough(norm(apply(to.rd(bank[estmask], "meat", emp), 1, mean), "2"), .3157, .001)
  omxCheckCloseEnough(norm(apply(to.rd(bank[estmask], "sem", emp), 1, mean), "2"), .2484, .001)
  omxCheckCloseEnough(cor(log(meat.condnum[estmask]), log(sem.condnum[estmask])), .873, .1)
  
  hist(cputime)
  which(!estmask)
  hist(emp$bias)
  round(emp$bias[order(emp$bias)],3)
  hist(apply(to.rd(bank[estmask], "sem", emp), 1, mean))
  hist(apply(to.rd(bank[estmask], "meat", emp), 1, mean))
  hist(apply(se.bias(bank[estmask], "sem", emp), 1, mean))
  hist(apply(se.bias(bank[estmask], "meat", emp), 1, mean))
  sum(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  which(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  plot(log(meat.condnum[estmask]), log(sem.condnum[estmask]))
  head(order(-apply(to.rd(bank[estmask], "sem", emp), 2, norm, type="2"))) #worst
}
