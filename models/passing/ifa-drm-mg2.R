library(OpenMx)
library(rpf)

mxOption(NULL, 'loglikelihoodScale', -2)  # default
set.seed(9)

numItems <- 20
i1 <- rpf.grm(outcomes=2)
items <- list()
items[1:numItems] <- i1
correct <- matrix(NA, 4, numItems)
for (x in 1:numItems) correct[,x] <- rpf.rparam(i1)

data.g1 <- rpf.sample(500, items, correct)
data.g2 <- rpf.sample(500, items, correct, mean=-1, cov=matrix(5,1,1))
data.g3 <- rpf.sample(500, items, correct, mean=1, cov=matrix(.5,1,1))

if (0) {
  # for flexMIRT, write CSV
  write.table(sapply(data.g1, unclass)-1, file="drm-mg2-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g2, unclass)-1, file="drm-mg2-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g3, unclass)-1, file="drm-mg2-g3.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

mkgroup <- function(model.name, data, latent.free) {
  ip.mat <- mxMatrix(name="ItemParam", nrow=2, ncol=numItems,
                     values=c(1,0), free=TRUE)
  
  for (ix in 1:numItems) {
    for (px in 1:2) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }
  
  dims <- 1
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=latent.free)
  
  m1 <- mxModel(model=model.name, ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(
                  ItemSpec=items,
                  ItemParam="ItemParam",
                  mean="mean", cov="cov"),
                mxFitFunctionML())
  m1
}

g1 <- mkgroup("g1", data.g1, FALSE)
g2 <- mkgroup("g2", data.g2, TRUE)
g3 <- mkgroup("g3", data.g3, TRUE)

groups <- paste("g", 1:3, sep="")

# load flexmirt fit and compare TODO

if (1) {
  # Cannot test derivatives at starting values because Hessian starts very close to singular.
  
  grpModel <- mxModel(model="groupModel", g1, g2, g3,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      mxComputeSequence(steps=list(
		      mxComputeEM(paste(groups, 'expectation', sep='.'),
		                  mxComputeNewtonRaphson(free.set=paste(groups,'ItemParam',sep=".")),
		                  mxComputeOnce('fitfunction', fit=TRUE,
		                                free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.'))))))
  
  #grpModel <- mxOption(grpModel, "Number of Threads", 1)
  
  grpModel <- mxRun(grpModel)
  #grpModel@output$minimum TODO
  omxCheckCloseEnough(grpModel@submodels$g2@matrices$mean@values, -.834, .01)
  omxCheckCloseEnough(grpModel@submodels$g2@matrices$cov@values, 3.93, .01)
  omxCheckCloseEnough(grpModel@submodels$g3@matrices$mean@values, .933, .01)
  omxCheckCloseEnough(grpModel@submodels$g3@matrices$cov@values, .444, .01)
  
  mxOption(NULL, 'loglikelihoodScale', 1)
  i1 <- mxModel(grpModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce(paste(groups, 'expectation', sep='.'), context="EM"),
                  mxComputeOnce('fitfunction', information=TRUE, info.method="meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i1 <- mxRun(i1)
  
  #cat(deparse(round(i1@output$standardErrors,3)))
  se <- c(0.071, 0.078, 0.076, 0.079, 0.097, 0.099, 0.132,  0.117, 0.075,
          0.077, 0.135, 0.121, 0.081, 0.083, 0.215, 0.169,  0.111, 0.141,
          0.121, 0.113, 0.213, 0.159, 0.074, 0.082, 0.077,  0.139, 0.084,
          0.087, 0.095, 0.09, 0.064, 0.07, 0.135, 0.115,  0.091, 0.095, 0.097,
          0.098, 0.096, 0.093, 0.12, 0.512, 0.072,  0.057)
  omxCheckCloseEnough(c(i1@output$standardErrors), se, .01)
  omxCheckCloseEnough(i1@output$conditionNumber, 281, 1)

  mxOption(NULL, 'loglikelihoodScale', -2)
  i1 <- mxModel(grpModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce(paste(groups, 'expectation', sep='.'), context="EM"),
                  mxComputeOnce('fitfunction', information=TRUE, info.method="meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i1 <- mxRun(i1)
  omxCheckCloseEnough(c(i1@output$standardErrors), se, .01)
  omxCheckCloseEnough(i1@output$conditionNumber, 281, 1)
}

if (0) {
  library(mirt)
  rdata <- sapply(data, unclass)-1
  # for flexMIRT, write CSV
  #write.table(rdata, file="ifa-drm-mg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  pars <- mirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars='values')
  pars[pars$name=="a1",'value'] <- 1
  pars[pars$name=="a1",'est'] <- FALSE
  pars[pars$name=="COV_11",'est'] <- TRUE
  fit <- mirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars=pars)
  # LL -7064.519 * -2 = 14129.04
  got <- coef(fit)
  print(got$GroupPars)
  # COV 4.551
  got$GroupPars <- NULL
  round(m2@matrices$itemParam@values - simplify2array(got), 2)
}
