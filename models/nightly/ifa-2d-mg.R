# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
#library(mvtnorm)

set.seed(7)
correct.LL <- 48990.18

numItems <- 30
numPeople <- 500
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.grm(factors=2)
	correct[[ix]] <- rpf.rparam(items[[ix]])
}
correct.mat <- simplify2array(correct)
correct.mat[1,1:10] <- 0
correct.mat[2,20:30] <- 0

if (0) {
  startpar <- sapply(items, rpf.rparam)
} else {
  startpar <- matrix(c(1.4, 1, 0), nrow=3, ncol=numItems)
}
startpar[correct.mat==0] <- 0

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))
maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))

data.g1 <- rpf.sample(numPeople, items, correct.mat)
data.g2 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.5,.8), cov=matrix(c(2,.5,.5,2),nrow=2))
data.g3 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.1,-.8), cov=matrix(c(.9,-.5,-.5,.9),nrow=2))

if (0) {
  # for flexMIRT, write CSV
  write.table(sapply(data.g1, unclass)-1, file="2d-mg-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g2, unclass)-1, file="2d-mg-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g3, unclass)-1, file="2d-mg-g3.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  fm <- read.flexmirt("/home/joshua/irt/ifa-2d-mg/2d-mg-prm.txt")
}

mkgroup <- function(model.name, data, latent.free) {  
  ip.mat <- mxMatrix(name="ItemParam", values=startpar, free=TRUE)
  ip.mat$free[correct.mat==0] <- FALSE
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat$labels[px,ix] <- name
    }
  }
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0)
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2))
  mean <- "mean"
  cov <- "cov"
  if (latent.free) {
    lm <- paste(model.name, "latent", sep="")
    mean <- paste(lm, "mean", sep=".")
    cov <- paste(lm, "cov", sep=".")
  }
  
  m1 <- mxModel(model=model.name,
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean=mean, cov=cov,
                                  ItemSpec=items, ItemParam="ItemParam",
                                  qpoints=21, qwidth=5, verbose=0L),
                mxFitFunctionML())
  m1
}

mklatent <- function(name) {
  cov <- mxMatrix(type="Symm", nrow=2, ncol=2, free=T, values=diag(2), name="cov")
  mask <- c(FALSE,TRUE,TRUE,FALSE)
  cov$labels[mask] <- paste(name, "cov",sep="")
  mxModel(paste(name, "latent", sep=""),
          mxDataDynamic('cov', expectation=paste(name, "expectation", sep=".")),
          mxMatrix(nrow=1, ncol=2, free=T, values=0, name="mean"),
          cov,
          mxExpectationNormal(covariance="cov", means="mean"),
          mxFitFunctionML()
  )
}

groups <- paste("g", 1:3, sep="")

latent <- mxModel("latent",
                  mxFitFunctionMultigroup(paste(paste(groups[-1],"latent",sep=""), "fitfunction", sep=".")))

g2.latent <- mklatent("g2")
g3.latent <- mklatent("g3")

latent.vargroup <- apply(expand.grid(paste(groups[-1], "latent", sep=""), c('mean','cov')),
                         1, paste, collapse='.')

latent.plan <- NULL  # need a plan for latent distribution parameters

if (0) {
  # Copy latent distribution parameters from current estimates without transformation.
  latent.plan <- mxComputeSequence(list(mxComputeOnce(paste(groups, 'expectation', sep='.'),
                                                      "latentDistribution", "copy"),  # c('mean','covariance')
                                        mxComputeOnce('fitfunction', "starting")),
                                   free.set=latent.vargroup)
} else {
  # Obtain latent distribution parameters via mxExpectationNormal.
  # This permits equality constraints (and potentially more complex latent structure).
  latent.plan <- mxComputeGradientDescent(latent.vargroup, fitfunction="latent.fitfunction")

  mask <- c(FALSE, TRUE, TRUE, FALSE)
  g2.latent$cov$labels[mask] <- 'eq1'
  g3.latent$cov$labels[mask] <- 'eq1'
}

g1 <- mkgroup("g1", data.g1, FALSE)
  g2 <- mkgroup("g2", data.g2, TRUE)
  g3 <- mkgroup("g3", data.g3, TRUE)

  grpModel <- mxModel(model="groupModel", g1, g2, g3, g2.latent, g3.latent, latent,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                                  mxComputeNewtonRaphson(paste(groups, 'ItemParam', sep=".")),
                                  latent.plan,
                                  mxComputeOnce('fitfunction', 'fit')))
  
grpModel <- mxRun(grpModel, silent=TRUE)
  
omxCheckCloseEnough(grpModel$output$minimum, correct.LL, .01)
omxCheckCloseEnough(c(grpModel$submodels$g2latent$mean$values), c(-.516, .708), .01)
omxCheckCloseEnough(c(grpModel$submodels$g2latent$cov$values), c(2.114, -.279,-.279, 2.259), .01)
omxCheckCloseEnough(c(grpModel$submodels$g3latent$mean$values), c(-.027, -.823), .01)
omxCheckCloseEnough(c(grpModel$submodels$g3latent$cov$values), c(.779, -.279, -.279, .738), .01)

emstat <- grpModel$compute$output
if (0) {
# TODO too inconsistent
omxCheckCloseEnough(emstat$EMcycles, 70, 2)
omxCheckCloseEnough(emstat$totalMstep, 222, 10)
}
