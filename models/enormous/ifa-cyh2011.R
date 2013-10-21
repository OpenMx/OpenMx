# This is a replication of Cai, Yang, & Hansen (2011) simulation study #1.

#options(error = utils::recover)
library(OpenMx)
library(rpf)
library(mvtnorm)

mk.model <- function(model.name, numItems, latent.free) {
  spec <- list()
  spec[1:numItems] <- rpf.grm(factors = 2)
  
  dims <- (1 + numItems/4)
  design <- matrix(c(rep(1,numItems),
                     kronecker(2:dims,rep(1,4))), byrow=TRUE, ncol=numItems)
  
  ip.mat <- mxMatrix(name="ItemParam", nrow=3, ncol=numItems,
                     values=c(1.4,1,0),
                     free=c(TRUE,TRUE,TRUE))
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }
  eip.mat <- mxAlgebra(ItemParam, name="EItemParam")

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  cov.mat.free <- FALSE
  if (latent.free) {
    cov.mat.free <- diag(dims)==1
  }
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=cov.mat.free)
  
  m1 <- mxModel(model=model.name, ip.mat, eip.mat, m.mat, cov.mat,
                mxExpectationBA81(
                  ItemSpec=spec,
                  design=design,
                  EItemParam="EItemParam", ItemParam="ItemParam",
                  mean="mean", cov="cov",
                  qpoints=21, qwidth=5),
                mxFitFunctionML())
  m1
}

omxIFAComputePlan <- function(groups) {
  mxComputeIterate(steps=list(
    mxComputeOnce(paste(groups, 'expectation', sep='.'), context='EM'),
    mxComputeNewtonRaphson(free.set=paste(groups, 'ItemParam', sep=".")),
    mxComputeOnce(paste(groups, 'expectation', sep=".")),
    mxComputeOnce(adjustStart=TRUE, 'fitfunction')
  ))
}

g2.mean <- c(1, -.5, 0, .5)
g2.cov <- diag(c(.8, 1.2, 1.5, 1))
correct <- matrix(c(1, 1.4, 1.7, 2, 1.4, 1.7, 2, 1, 1.7, 2, 1, 1.4, 2, 1, 1.4, 1.7,
		    .8,1.5, 1.2, 1, 1, .8, 1.5, 1.2, 1.2, 1, .8, 1.5, 1.5, 1.2, 1, .8,
		    rep(c(1, .25, -.25, -1), 4)), byrow=TRUE, ncol=16)

groups <- paste("g", 1:2, sep="")

fit1 <- function(seed) {
  result <- list(seed=seed)
  
  set.seed(seed)

  g1 <- mk.model("g1", 16, FALSE)
  g2 <- mk.model("g2", 12, TRUE)

  data.g1 <- rpf.sample(1000, g1@expectation@ItemSpec, correct, g1@expectation@design)
  data.g2 <- rpf.sample(1000, g2@expectation@ItemSpec, correct[,1:12], g2@expectation@design,
                        mean=g2.mean, cov=g2.cov)
  
  g1 <- mxModel(g1, mxData(observed=data.g1, type="raw"))
  g2 <- mxModel(g2, mxData(observed=data.g2, type="raw"))
  
  grpModel <- mxModel(model="groupModel", g1, g2,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))

  grpModel <- mxRun(grpModel)
  
  result$cpuTime <- grpModel@output$cpuTime
  result$LL <- grpModel@output$Minus2LogLikelihood
  result$param <- grpModel@submodels$g1@matrices$ItemParam@values
  result$mean <- grpModel@submodels$g2@matrices$mean@values
  result$cov <- grpModel@submodels$g2@matrices$cov@values
  result
}

mc.estimate <- function (bank, sl) {
  example <- bank[[1]][[sl]]
  bias <- matrix(0, nrow=dim(example)[1], ncol=dim(example)[2])
  for (sx in 1:length(bank)) {
    bias <- bias + bank[[sx]][[sl]]
  }
  bias / length(bank)
}

bank <- list()
#setwd("/opt/OpenMx")
rda <- "ifa-cyh2011.rda"
if (file.exists(rda)) load(rda)
for (seed in 1:500) {
  if (length(bank)) {
    if (any(seed == sapply(bank, function (b) b$seed))) next;
  }
  bi <- length(bank) + 1
  bank[[bi]] <- fit1(seed)
  save(bank, file=rda)

  if (0) {
    cur <- bank[[bi]]
    print(cor(c(cur$param, cur$mean, diag(cur$cov)), c(correct, g2.mean, diag(g2.cov))))
  }
}

if (1) {
  bias <- c(mc.estimate(bank, 'param') - correct,
            mc.estimate(bank, 'mean') - g2.mean,
            mc.estimate(bank, 'cov') - g2.cov)
  omxCheckTrue(abs(bias) < .061)
}

if (0) {
  require(ggplot2)
  sbank <- bank[300:800]
  df <- rbind(
    data.frame(true=c(correct), bias=c(mc.estimate(sbank, 'param') - correct),
               type=c("primary", "specific", "diff")),
    data.frame(true=c(g2.mean), bias=c(mc.estimate(sbank, 'mean') - g2.mean), type="mean"),
    data.frame(true=diag(g2.cov), bias=diag(mc.estimate(sbank, 'cov') - g2.cov), type="var"))
  
  df$type <- factor(df$type)
  ggplot(df, aes(true, bias, color=type)) + geom_point(size=3) + xlab("true parameter value") +
     ylab(paste("bias (", length(sbank), "replications)"))
}
