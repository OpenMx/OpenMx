#options(error = browser)
library(OpenMx)
library(rpf)

numItems <- 12
spec <- list()
spec[1:numItems] <- rpf.drm()

g341 <- suppressWarnings(try(read.table("models/passing/data/g341-19.dat"), silent=TRUE))
if (is(g341, "try-error")) g341 <- read.table("data/g341-19.dat")

if (0) {
  # no longer implemented in MIRT
  require(mirt)
  pars <- mirt(data=g341, 1, itemtype='3PL', pars='values', parprior=list(c(seq(3,47,4), 'beta', 2, 5)))
  fit <- mirt(data=g341, 1, pars=pars)
}

g341 <- as.data.frame(lapply(g341, mxFactor, levels=0:1))

ip.mat <- mxMatrix(name="ItemParam", nrow=4, ncol=numItems,
                   values=c(1,0, NA, logit(1)), free=c(TRUE, TRUE, TRUE, FALSE))
plabel <- paste("g", 1:numItems, sep="")
ip.mat$labels[3,] <- plabel
ip.mat$values[3,] <- logit(seq(.1, .2, length.out=numItems))
ip.mat$ubound[3,] <- logit(.9)
# try random starting parameters? TODO

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

g.mat <- mxMatrix(name="gparam", nrow=1, ncol=numItems, free=TRUE, labels=plabel,
                  values=ip.mat$values[3,], ubound=ip.mat$ubound[3,])

# plot(function (x) dbeta(1/(1+exp(-x)),2,5), -10,10)
#beta.const <- 1/beta(2,5)  #=30

prior <- mxAlgebra(-2* sum(log(30 * exp(gparam) / (1 + exp(gparam))^5 )), name="prior")
prior.grad <- mxAlgebra(2 * (4 * exp(gparam) - 1) / (exp(gparam)+1), name="pgrad",
                        dimnames=list(c(),plabel))
prior.hess <- mxAlgebra(vec2diag(10*exp(gparam)/(exp(gparam)+1)^2),
                        name="phess", dimnames=list(plabel, plabel))

pm <- mxModel(model="pmodel", g.mat, prior, prior.grad, prior.hess,
              mxFitFunctionAlgebra("prior", gradient="pgrad", hessian="phess"))

m1 <- mxModel(model="item", ip.mat, m.mat, cov.mat,
              mxData(observed=g341, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam",
                                mean="mean", cov="cov"),
              mxFitFunctionML())

gm <- mxModel(model="gm", pm,
  mxComputeSequence(list(
    mxComputeOnce('pmodel.fitfunction', c('fit', 'gradient', 'hessian', 'ihessian')),
    mxComputeReportDeriv())))
testDeriv <- mxRun(gm)
omxCheckCloseEnough(testDeriv$output$fit, -2 * sum(dbeta(1/(1+exp(-ip.mat$values[3,])), 2, 5, log=TRUE)), .01)
#cat(deparse(round(fivenum(testDeriv$output$gradient), 3)))
omxCheckCloseEnough(fivenum(testDeriv$output$gradient),
                    c(-1, -0.773, -0.5, -0.227, 0), .01)
#cat(deparse(round(fivenum(diag(testDeriv$output$hessian)), 3)))
omxCheckCloseEnough(fivenum(diag(testDeriv$output$hessian)),
                    c(0.9, 1.076, 1.275, 1.458, 1.6), .01)
omxCheckCloseEnough(testDeriv$output$hessian %*% testDeriv$output$ihessian, diag(12), 1e-2)

if (0) {
  require("numDeriv")
  got <- genD(function(x) {
    pm$matrices$gparam$values[,seq(1,12,3)] <- x
    gm <- mxRun(mxModel(model="gm", pm,
                        mxComputeOnce('pmodel.fitfunction', 'fit')), silent=TRUE)
    gm$output$fit
  }, pm$matrices$gparam$values[,seq(1,12,4)], method.args=list(r=2))
}

if (0) {
  pm <- mxModel(model="pmodel", g.mat, prior, prior.grad, prior.hess,
                mxFitFunctionAlgebra(NULL, gradient="pgrad", hessian="phess"))
}

ponly <- mxModel("ex", pm,
                 mxComputeNewtonRaphson(fitfunction="pmodel.fitfunction"))
ponly.fit <- mxRun(ponly)
omxCheckCloseEnough(ponly.fit$submodels$pmodel$matrices$gparam$values[,],
                    rep(logit(.2), numItems), .001)

if (1) {
  fm1 <- mxModel("fm", m1, pm,
                 mxFitFunctionMultigroup(groups=c('pmodel.fitfunction', 'item.fitfunction')),
                 mxComputeSequence(list(
                   mxComputeOnce('fitfunction', c('fit','gradient')),
                   mxComputeReportDeriv())))
  citem <- fm1$submodels$item
  citem$ItemParam$values[1:3,] <-
    c(1.8391, -0.5642, -1.0869, 1.2723, 2.7656, -2.0005,  1.504, -0.9809, -1.193, 1.3107,
      0.4609, -1.1923, 1.348, -0.422,  -1.076, 1.3897, 1.9059, -1.9338, 1.9577, 4.2709,
      -2.2242, 1.0483,  1.2432, -1.5608, 1.9411, 3.3146, -2.5864, 1.3474, 2.2773, -2.2306,
      1.7741, -1.398, -1.254, 1.5506, -0.5158, -1.3686)
  cpmodel <- fm1$submodels$pmodel
  cpmodel$gparam$values[,] <- citem$ItemParam$values[3,]
  fm1 <- mxModel(fm1, citem, cpmodel)
  fm1 <- mxRun(fm1, silent=TRUE)
  omxCheckCloseEnough(max(abs(fm1$output$gradient)), 1.29, .1)
  omxCheckCloseEnough(fm1$output$fit - fm1$submodels$pmodel$fitfunction$result, 33335.75, .01)
  
  if (0) {
    dm <- fm1
    dm$submodels$item$expectation$EItemParam <- dm$submodels$item$matrices$ItemParam$values[,]
    dm$submodels$item$matrices$ItemParam$free[,3:12] <- FALSE
    dm$submodels$item$matrices$ItemParam$values[1:3,1:2] <- sapply(spec[1:2], rpf.rparam)[1:3,]
    dm$submodels$pmodel$matrices$gparam$free[1,3:12] <- FALSE
    dm$submodels$pmodel$matrices$gparam$values[,] <-
      dm$submodels$item$matrices$ItemParam$values[3,]
    dm$compute <- mxComputeSequence(list(
      mxComputeOnce('fitfunction', c('fit','gradient', 'hessian')),
      mxComputeReportDeriv()))
    dm.fit <- mxRun(dm, silent=TRUE)
    
    dm$compute <- mxComputeSequence(list(
      mxComputeOnce('fitfunction', c('fit'))))
    require("numDeriv")
    got <- genD(function(x) {
      dm$submodels$item$matrices$ItemParam$values[1:3,1:2] <- x
      dm$submodels$pmodel$matrices$gparam$values[,] <-
        dm$submodels$item$matrices$ItemParam$values[3,]
      fit <- mxRun(dm, silent=TRUE)
      fit$output$fit
    }, dm$submodels$item$matrices$ItemParam$values[1:3,1:2], method.args=list(r=2))
    
    unpackHessian <- function(deriv, np) {
      hess <- matrix(NA, nrow=np, ncol=np)
      dx <- np+1
      for (hr in 1:np) {
        hess[1:hr,hr] <- hess[hr,1:hr] <- deriv[dx:(dx+hr-1)]
        dx <- dx + hr
      }
      hess
    }
    omxCheckCloseEnough(got$D[1:6], dm.fit$output$gradient, 1e-6)
    # max(abs(unpackHessian(got$D, 6) - dm.fit$output$hessian))
    omxCheckCloseEnough(unpackHessian(got$D, 6), dm.fit$output$hessian, .2)
  }
}

m2 <- mxRun(mxModel("ex", m1, pm,
                    mxFitFunctionMultigroup(groups=c('pmodel.fitfunction', 'item.fitfunction'),
                                            verbose=0L),
                    mxCI(c('g1')),
                    mxComputeSequence(list(
                      mxComputeEM('item.expectation', 'scores',
                                  mxComputeNewtonRaphson(verbose=0L, maxIter=50L),
                                tolerance=1e-5, verbose=0L),
                      mxComputeConfidenceInterval(),
                      mxComputeOnce('fitfunction', c('fit','gradient')),  # SEM lost the details
                      mxComputeReportDeriv()))), silent=TRUE)
# flexmirt's LL is reported w/o prior
omxCheckCloseEnough(m2$output$fit - m2$submodels$pmodel$fitfunction$result, 33335.75, .1)
omxCheckCloseEnough(max(abs(m2$output$gradient)), 1.29, .1)
#cat(deparse(round(m2$output$confidenceIntervals,3)))
omxCheckCloseEnough(m2$output$confidenceIntervals['g1',], c(-1.687, -0.726), .01)

g1 <- mxRun(mxModel(m2, mxComputeSequence(list(
  mxComputeOnce('pmodel.fitfunction', 'gradient'),
  mxComputeReportDeriv()))), silent=TRUE)
g2 <- mxRun(mxModel(m2, mxComputeSequence(list(
  mxComputeOnce('item.fitfunction', 'gradient'),
  mxComputeReportDeriv()))), silent=TRUE)

emstat <- m2$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 18, 3)
omxCheckCloseEnough(emstat$totalMstep, 92, 5)
#omxCheckCloseEnough(emstat$semProbeCount, 108, 5)
#omxCheckCloseEnough(log(m2$output$conditionNumber), 6.12, 1)

# SEs probably wrong TODO
#cat(deparse(round(c(m2$output$standardErrors), 3)))
if (0) {
  semse <- c(0.309, 0.286, 0.232, 0.162, 0.438, 3.704, 0.248, 0.274, 0.21,  0.185, 0.226, 0.428,
             0.201, 0.227, 0.256, 0.17, 0.224, 1.547,  1.388, 0.947, 18.144, 1.216, 2.735, 10.888,
             0.52, 0.492, 9.264,  0.921, 2.031, 20.255, 0.586, 0.602, 0.348, 0.234, 0.25, 0.291 )
  omxCheckCloseEnough(c(m2$output$standardErrors), semse, .2)  # very unstable
}

# SEs probably wrong TODO
if (0) {
i1 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)
omxCheckTrue(i1$output$infoDefinite)
omxCheckCloseEnough(log(i1$output$conditionNumber), 9.7, .5)
# cat(deparse(round(c(i1$output$standardErrors), 3)))
se <- c(0.308, 0.266, 0.232, 0.216, 0.601, 4.736, 0.266, 0.301, 0.236,  0.202, 0.243,
        0.502, 0.222, 0.265, 0.295, 0.213, 0.285, 1.883,  0.278, 0.494, 5.178, 0.178,
        0.354, 1.443, 0.242, 0.267, 3.554,  0.19, 0.346, 3.115, 0.325, 0.343, 0.192,
        0.243, 0.243, 0.295)
# max(abs(c(i1$output$standardErrors) - se))
omxCheckCloseEnough(c(i1$output$standardErrors), se, .01)
}
