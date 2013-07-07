# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
#library(mvtnorm)

set.seed(7)

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

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

data.g1 <- rpf.sample(numPeople, items, correct.mat)
data.g2 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.5,.8), cov=matrix(c(2,.5,.5,2),nrow=2))
data.g3 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.1,-.8), cov=matrix(c(.9,-.5,-.5,.9),nrow=2))

if (0) {
  # for flexMIRT, write CSV
  write.table(sapply(data.g1, unclass)-1, file="2d-mg-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g2, unclass)-1, file="2d-mg-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g3, unclass)-1, file="2d-mg-g3.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

mkgroup <- function(model.name, data, latent.free) {  
  spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
                   values=sapply(items, function(m) slot(m,'spec')),
                   free=FALSE, byrow=TRUE)
  
  ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                     values=c(1, 1.4, 0),
                     lbound=c(1e-6, 1e-6, NA), free=TRUE)
  ip.mat@free.group <- "param"
  ip.mat@values[correct.mat==0] <- 0
  ip.mat@free[correct.mat==0] <- FALSE
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }
  
  eip.mat <- mxAlgebra(ItemParam, name="EItemParam")
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=latent.free)
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=latent.free)
  
  m1 <- mxModel(model=model.name,
                spec, ip.mat, m.mat, cov.mat, eip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec="ItemSpec",
                                  EItemParam="EItemParam",
                                  qpoints=21, qwidth=5),
                mxFitFunctionBA81(ItemParam="ItemParam"),
                mxComputeIterate(steps=list(
                  mxComputeOnce("EItemParam"),
                  mxComputeOnce('expectation', context='E'),
                  mxComputeNewtonRaphson(free.group='param'),
                  #				 mxComputeGradientDescent(free.group='param'),
                  mxComputeOnce('expectation', context='M'),
                  mxComputeOnce('fitfunction'))))
  m1
}

g1 <- mkgroup("g1", data.g1, FALSE)
g2 <- mkgroup("g2", data.g2, TRUE)
g3 <- mkgroup("g3", data.g3, TRUE)

groups <- paste("g", 1:3, sep="")

if (1) {
  grpModel <- mxModel(model="groupModel", g1, g2, g3,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      mxComputeIterate(steps=list(
                        mxComputeOnce(paste(groups, "EItemParam", sep=".")),
                        mxComputeOnce(paste(groups, 'expectation', sep='.'), context='E'),
                        #                      mxComputeGradientDescent(free.group='param'),
                        mxComputeNewtonRaphson(free.group='param'),
                        mxComputeOnce(paste(groups, 'expectation', sep="."), context='M'),
                        mxComputeOnce('fitfunction')
                      )))

  grpModel <- mxOption(grpModel, "Analytic Gradients", 'Yes')
	grpModel <- mxOption(grpModel, "Verify level", '-1')
  grpModel <- mxOption(grpModel, "Function precision", '1.0E-5')

  grpModel <- mxRun(grpModel, silent=TRUE)
#  grpModel@output$minimum   # 48939.35 TODO
  omxCheckCloseEnough(c(grpModel@submodels$g2@matrices$mean@values), c(-.507, .703), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g2@matrices$cov@values), c(1.877, .491, .491, 2.05), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g3@matrices$mean@values), c(-.028, -.827), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g3@matrices$cov@values), c(.892, -.422, -.422, .836), .01)
}
