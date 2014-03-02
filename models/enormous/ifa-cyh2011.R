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

fit1 <- function(seed) {
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
      mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                  mxComputeNewtonRaphson(free.set=paste(groups, 'ItemParam', sep=".")),
                  mxComputeOnce('fitfunction', 'fit', free.set=apply(expand.grid(groups, c('mean','cov')),
                                                              1, paste, collapse='.')),
                  tolerance=1e-5, information=TRUE),
      mxComputeStandardError(),
      mxComputeHessianQuality()))
  }

  grpModel <- mxModel(model="groupModel", g1, g2,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))

  grpModel <- try(mxRun(grpModel, silent=TRUE), silent=TRUE)
  if (inherits(grpModel, "try-error")) return(NULL)

    i1 <- mxModel(grpModel,
                  mxComputeSequence(steps=list(
                    mxComputeOnce(paste(groups, 'expectation', sep='.')),
                    mxComputeOnce('fitfunction', 'information', "meat"),
                    mxComputeStandardError(),
                    mxComputeHessianQuality())))
    i1 <- mxRun(i1, silent=TRUE)
  
  got <- cbind(grpModel@output$estimate,
               grpModel@output$standardErrors,
               i1@output$standardErrors)
  colnames(got) <- c("est", "sem", "meat")
  
  sem.condnum <- grpModel@output$conditionNumber
  
  list(condnum=c(sem=ifelse(is.null(sem.condnum), NA, sem.condnum), meat=i1@output$conditionNumber),
       got=got,
       em=grpModel@compute@steps[[1]]@output,
       cpuTime=grpModel@output$cpuTime)
}

to.rd <- function(bank, type, emp) {
  se <- sapply(bank, function (t) t$got[,type])
  mask <- apply(se, 2, function(c) all(!is.na(c)))
  rd <- apply(se[,mask], 2, function (c) (c - emp$se)/emp$se)
  rd
}

se.bias <- function(bank, type, emp) {
  se <- sapply(bank, function (t) t$got[,type])
  mask <- apply(se, 2, function(c) all(!is.na(c)))
  bias <- apply(se[,mask], 2, function (c) (c - emp$se))
  bias
}

bank <- list()

for (seed in 1:500) {
  bank[[seed]] <- fit1(seed)
  print(seed)
}

if (0) {
  setwd("/opt/OpenMx")
  rda <- "cyh2011.rda"
  save(bank, file=rda)
  if (file.exists(rda)) load(rda)
}

cputime <- sapply(bank, function (t) t$cpuTime)
est <- sapply(bank, function (t) t$got[,'est'])
diff <- apply(est, 2, function(c) (c - c(correct, g2.mean, diag(g2.cov))))
meat.condnum <- sapply(bank, function (t) t$condnum['meat'])
sem.condnum <- sapply(bank, function (t) t$condnum['sem'])
estmask <- !is.na(sem.condnum) & !is.na(meat.condnum) & meat.condnum < 2300
emp <- list(bias=apply(diff[,estmask], 1, mean), se=apply(diff[,estmask], 1, sd))

if (0) {
  hist(cputime)
  which(!estmask) #   62  133  184  356
  hist(emp$bias)
  round(emp$bias[order(emp$bias)],3)
  hist(apply(to.rd(bank[estmask], "sem", emp), 1, mean))
  hist(apply(to.rd(bank[estmask], "meat", emp), 1, mean))
  hist(apply(se.bias(bank[estmask], "sem", emp), 1, mean))
  hist(apply(se.bias(bank[estmask], "meat", emp), 1, mean))
  fivenum(sapply(bank, function (t) t$em$semProbeCount) / length(c(correct, g2.mean, diag(g2.cov))))
  sum(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  which(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  plot(log(meat.condnum[estmask]), log(sem.condnum[estmask]))
}

omxCheckCloseEnough(norm(emp$bias, "2"), 0.14211, .001)
omxCheckCloseEnough(norm(emp$se, "2"), 1.2709, .001)
omxCheckCloseEnough(norm(apply(to.rd(bank[estmask], "meat", emp), 1, mean), "2"), .3141, .001)
omxCheckCloseEnough(norm(apply(to.rd(bank[estmask], "sem", emp), 1, mean), "2"), .377039, .001)
omxCheckCloseEnough(cor(log(meat.condnum[estmask]), log(sem.condnum[estmask])), .36, .02)
