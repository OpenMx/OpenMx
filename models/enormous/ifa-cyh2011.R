# This is a replication of Cai, Yang, & Hansen (2011) simulation study #1.

#options(error = utils::recover)
library(OpenMx)
library(rpf)
library(mvtnorm)

mk.model <- function(model.name, numItems, latent.free) {
  spec <- list()
  spec[1:numItems] <- rpf.grm(factors = 2)
  
  dims <- (1 + numItems/4)
  design <- matrix(as.integer(c(rep(1L,numItems),
                     kronecker(2L:dims,rep(1L,4)))), byrow=TRUE, ncol=numItems)
  
  ip.mat <- mxMatrix(name="ItemParam", nrow=3, ncol=numItems,
                     values=c(1.4,1,0),
                     free=c(TRUE,TRUE,TRUE))
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  cov.mat.free <- FALSE
  if (latent.free) {
    cov.mat.free <- diag(dims)==1
  }
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=cov.mat.free)
  
  m1 <- mxModel(model=model.name, ip.mat, m.mat, cov.mat,
                mxExpectationBA81(
                  ItemSpec=spec,
                  design=design, ItemParam="ItemParam",
                  mean="mean", cov="cov",
                  qpoints=21, qwidth=5),
                mxFitFunctionML())
  m1
}

g2.mean <- c(1, -.5, 0, .5)
g2.cov <- diag(c(.8, 1.2, 1.5, 1))
correct <- matrix(c(1, 1.4, 1.7, 2, 1.4, 1.7, 2, 1, 1.7, 2, 1, 1.4, 2, 1, 1.4, 1.7,
		    .8,1.5, 1.2, 1, 1, .8, 1.5, 1.2, 1.2, 1, .8, 1.5, 1.5, 1.2, 1, .8,
		    rep(c(1, .25, -.25, -1), 4)), byrow=TRUE, ncol=16)

groups <- paste("g", 1:2, sep="")

fit1 <- function(seed, ntarget) {
  set.seed(seed)

  g1 <- mk.model("g1", 16, FALSE)
  g2 <- mk.model("g2", 12, TRUE)

  data.g1 <- rpf.sample(1000, g1@expectation@ItemSpec, correct, g1@expectation@design)
  data.g2 <- rpf.sample(1000, g2@expectation@ItemSpec, correct[,1:12], g2@expectation@design,
                        mean=g2.mean, cov=g2.cov)
  
  if (0) {
    # groups are swapped in my flexmirt specs
    write.table(sapply(data.g1, unclass)-1, file=paste("cyh",seed,"g2.csv",sep="-"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(sapply(data.g2, unclass)-1, file=paste("cyh",seed,"g1.csv",sep="-"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  
  g1 <- mxModel(g1, mxData(observed=data.g1, type="raw"))
  g2 <- mxModel(g2, mxData(observed=data.g2, type="raw"))
  
  omxIFAComputePlan <- function(groups) {
    mxComputeSequence(steps=list(
      mxComputeEM(paste(groups, 'expectation', sep='.'),
                  mxComputeNewtonRaphson(free.set=paste(groups, 'ItemParam', sep=".")),
                  mxComputeOnce('fitfunction', free.set=apply(expand.grid(groups, c('mean','cov')),
                                                              1, paste, collapse='.'), fit=TRUE),
                  information=TRUE, info.method="hessian", noiseTarget=ntarget, agileMaxIter=3L),
      mxComputeStandardError(forcePositiveDefinite=TRUE),
      mxComputeHessianQuality()))
  }

  grpModel <- mxModel(model="groupModel", g1, g2,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))

  grpModel <- try(mxRun(grpModel, silent=TRUE), silent=TRUE)
  if (inherits(grpModel, "try-error")) return(NULL)

  if (0) {
    i1 <- mxModel(grpModel,
                  mxComputeSequence(steps=list(
                    mxComputeOnce(paste(groups, 'expectation', sep='.')),
                    mxComputeOnce('fitfunction', information=TRUE, info.method="meat"),
                    mxComputeStandardError(),
                    mxComputeHessianQuality())))
    i1 <- mxRun(i1, silent=TRUE)
  }
  
  got <- cbind(grpModel@output$estimate,
               grpModel@output$standardErrors)
  colnames(got) <- c("est", "sem")
  
  sem.condnum <- grpModel@output$conditionNumber
  
  list(condnum=c(sem=ifelse(is.null(sem.condnum), NA, sem.condnum)),
       got=got,
       em=grpModel@compute@steps[[1]]@output,
       cpuTime=grpModel@output$cpuTime,
       ev=grpModel@compute@steps[[2]]@output$eigenvalues)
}

estmask <- rep(TRUE, 500)
emp <- list()

if (0) {
  meat.condnum <- sapply(bank, function (t) t$condnum['meat'])
  estmask <- !is.na(meat.condnum) & meat.condnum < 2300
  emp <- list(bias=apply(diff[,estmask], 1, mean), se=apply(diff[,estmask], 1, sd))
} else {
  emp <- list(se=c(c(0.08944, 0.11908, 0.08761, 0.14751, 0.2364, 0.11452,  0.14415,
                     0.16342, 0.10032, 0.16025, 0.13619, 0.11682, 0.11997,  0.15394,
                     0.09972, 0.12668, 0.1275, 0.09362, 0.197, 0.25594, 0.117,  0.10577,
                     0.18637, 0.10332, 0.14281, 0.20174, 0.10412, 0.16481,  0.16631, 0.10784,
                     0.08349, 0.13497, 0.07527, 0.15189, 0.3054,  0.13913, 0.26249, 0.33094,
                     0.15823, 0.1349, 0.27334, 0.09186,  0.15155, 0.20331, 0.09334, 0.17177,
                     0.18151, 0.1149, 0.1202,  0.17129, 0.15999, 0.16465, 0.09639, 0.26893, 0.31799, 0.26872 )))
  estmask[c(62,  133,  140, 184,  356,  380,  408 )] <- FALSE
}

to.rd <- function(bank, type) {
  se <- sapply(bank[estmask[1:length(bank)]], function (t) t$got[,type])
  mask <- apply(se, 2, function(c) all(!is.na(c)))
  rd <- apply(se[,mask], 2, function (c) (c - emp$se)/emp$se)
  rd
}

if (0) {
  setwd("/opt/OpenMx")
  rda <- "cyh2011.rda"
  save(mbank, file=rda)
  if (file.exists(rda)) load(rda)
}

mbank <- list()
Targets <- seq(-4.7, -1.5, .1)
for (seed in 249:500) {
  if (!estmask[seed]) next
  rds <- matrix(NA, length(Targets), 6)
  rds[,1] <- Targets
  for (tx in 1:length(Targets)) {
    if (length(mbank) < tx) mbank[[tx]] <- list()
    mbank[[tx]][[seed]] <- fit1(seed, exp(Targets[tx]))
    if (seed >= 3) {
      rd <- to.rd(mbank[[tx]], "sem")
      bank <- mbank[[tx]][estmask[1:seed]]
      rds[tx,2:6] <- c(norm(apply(rd, 1, mean),"2"),
                       norm(apply(rd, 1, sd), "2"),
                       mean(sapply(bank, function (t) t$em$semProbeCount / nrow(rd))),
                       mean(sapply(bank, function (t) t$condnum['sem'])),
                       mean(sapply(bank, function (t) ifelse(t$ev[1] < 0,
                                                                    which(order(abs(t$ev))==1) / length(t$ev), NA)), na.rm=TRUE))
    }
  }
  print(seed)
  print(rds)
}

est <- sapply(bank, function (t) t$got[,'est'])
diff <- apply(est, 2, function(c) (c - c(correct, g2.mean, diag(g2.cov))))

if (0) {
  which(!estmask) #   62  133  184  356
  hist(emp$bias)
  hist(check.se("sem"))
  hist(check.se("pdse"))
  hist(check.rd("meat"))
  hist(check.rd("sem"))
  hist(check.rd("pdse"))
  norm(check.rd("meat"), "2")  # 0.313
  norm(check.rd("sem"), "2")
  #exp(-.5) was 0.3
  #exp(-1) was 1.02
  #exp(-2) was 0.83
  #exp(-2.5) was .607
  #exp(-2.75) was .479
  #exp(-3) was 0.37
  #exp(-3.25) was 0.394
  #exp(-3.5) was 0.44
  norm(check.rd("pdse"), "2")
  fivenum(sapply(bank, function (t) t$em$semProbeCount) / length(c(correct, g2.mean, diag(g2.cov))))
  #exp(-2) 237 were non-pd
  #exp(-3) 296 non-pd
  sum(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  which(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  sum(sapply(bank, function (t) any(is.na(t$got[,'pdse']))))
  # target=exp(-4), 1 10 23 25 27 28 31 38 61 64 65 68 72 73 78 90
  # target=exp(-4.5), 10 38 61 68 90
}
