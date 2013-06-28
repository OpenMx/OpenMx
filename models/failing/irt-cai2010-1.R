# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

options(error = utils::recover)
library(OpenMx)
library(rpf)
library(mvtnorm)

numPersons <- 500

items.3d <- list()
items.3d[1:20] <- rpf.drm(dimensions=3)
#lapply(items, function (i) i@spec[4] <- 100)  # flatten prior on slope
numItems <- length(items.3d)
maxParam.3d <- max(vapply(items.3d, function(i) i@numParam, 0))

if(0) {
  correct <- sapply(items.3d, rpf.rparam)
  correct[5,] <- 0
  correct[1,11:20] <- 0
  correct[2,1:10] <- 0
  correct[3,c(3:5,8:10,13:15,18:20)] <- 0
  for (d in c(1,6,11,16)) correct[2,d+1] <- correct[2,d]
} else {
  # Cai (2010, p. 592)
  correct <- matrix(c(1.7,1.5,2.06,2.56,1.96,2.82,1,2.6,1.15,2.76,rep(0,10),
                      rep(0,10),2.62,2.3,1.94,1.43,1.95,1.67,2.67,1.98,2.22,2.19,
                      1.63,1.64,rep(0,3),1.56,1.56,rep(0,3),1.89,1.89,rep(0,3),2,2,rep(0,3),
                      .91,.61,.02,.04,.02,.21,-.4,-.24,-.74,.3,1.85,
                      1.12,.58,-1.26,.07,1.05,2.05,.16,-.44,.3,
                      rep(0,20)), byrow=TRUE, ncol=numItems)
}

design.3d <- matrix(c(rep(1,10),rep(NA,10),
                   rep(NA,10),rep(2,10),
                   3,3,rep(NA,3),4,4,rep(NA,3),5,5,rep(NA,3),6,6,rep(NA,3)), byrow=TRUE, nrow=3)

cov <- matrix(c(1, .68, .68, 1), nrow=2)

fit1 <- function(seed=5, ghp=11) {
  result <- list(seed=seed, ghp=ghp)
  
  set.seed(seed)
  ability <- cbind(rmvnorm(numPersons, sigma=cov),
                   matrix(rnorm(numPersons * 4), ncol=4))
  design.full <- design.3d
  design.full[is.na(design.3d)] <- 1  # ignored because of item parameters
  data.3d <- rpf.sample(ability, items.3d, correct, design.full)
  
  ip.labels <- matrix(NA, nrow=maxParam.3d, ncol=numItems)
  ip.labels[3,1:2] <- "s1"
  ip.labels[3,6:7] <- "s2"
  ip.labels[3,11:12] <- "s3"
  ip.labels[3,16:17] <- "s4"
  
  ip.free <- correct != 0
  ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam.3d, ncol=numItems,
                     values=c(1.414, 1.414, 1, 0, 0),
                     free=ip.free, labels=ip.labels,
                     lbound=c(1*10^-9, 1*10^-9, 1*10^-9, -10^9,0))
  ip.mat@values[!ip.free] <- 0
  
  m1 <- mxModel(model="cai1",
                mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
                         values=sapply(items.3d, function(m) slot(m,'spec')),
                         free=FALSE, byrow=TRUE),
                mxMatrix(name="Design", nrow=dim(design.3d)[1], ncol=numItems, values=design.3d),
                ip.mat,
                mxData(observed=data.3d, type="raw"),
                mxExpectationBA81(ItemSpec="ItemSpec", Design="Design", ItemParam="ItemParam",
                                  GHpoints=ghp),
                mxFitFunctionBA81()
  )
  
  m1 <- mxOption(m1, "Analytic Gradients", 'no')
  if (1) {
    m1 <- mxOption(m1, "Analytic Gradients", 'yes')
    m1 <- mxOption(m1, "Verify level", '-1')
  }
  m1 <- mxOption(m1, "Function precision", '1.0E-5')
  m1 <- mxOption(m1, "Calculate Hessian", "No")
  m1 <- mxOption(m1, "Standard Errors", "No")
  
  m1 <- mxRun(m1, silent=TRUE)
  
  result$cpuTime <- m1@output$cpuTime
  result$LL <- m1@output$Minus2LogLikelihood
  result$param <- m1@matrices$ItemParam@values
  result
}

calc.bias <- function (bank) {
  bias <- matrix(0, nrow=dim(correct)[1], ncol=dim(correct)[2])
  for (sx in 1:length(bank)) {
    bias <- bias + bank[[sx]]$param
  }
  bias <- (bias / length(bank)) - correct
  bias
}

if (0) {
  bank <- list()
}
setwd("/opt/OpenMx")
rda <- "irt-cai2010.rda"
load(rda)
for (seed in 1:300) {
  if (any(seed == sapply(bank, function (b) b$seed))) next;
  bank[[length(bank) + 1]] <- fit1(ghp=13, seed=seed)
  save(bank, file=rda)
}

calc.bias(bank[sapply(bank, function (b) b$ghp == 13)])
