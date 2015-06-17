#options(error = browser)
require(OpenMx)
require(rpf)
require(ifaTools)

set.seed(2)
numItems <- 15
numPersons <- 250
slope <- 2
groupSize <- 5

items <- list()
items[1:numItems] <- rpf.drm()
correct.mat <- matrix(NA, 4, numItems)
dimnames(correct.mat) <- list(names(rpf.rparam(items[[1]])), paste("i", 1:numItems, sep=""))
correct.mat['a',] <- slope
correct.mat['b',] <- slope * seq(-1.5, 1.5, length.out = groupSize)
correct.mat['g',] <- kronecker(logit(1/2:(1+numItems/groupSize)), rep(1,groupSize))
correct.mat['u',] <- logit(1)

correct.mask <- matrix(FALSE, 4, numItems)
correct.mask[1,1] <- TRUE
correct.mask[2,] <- TRUE
correct.mask[3,] <- TRUE

mkmodel <- function() {
  maxParam <- max(vapply(items, rpf.numParam, 0))
  maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))
  
  data <- rpf.sample(numPersons, items, correct.mat)

  ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                     values=c(1, 0, NA, logit(1)), free=TRUE,
                     dimnames=list(names(rpf.rparam(items[[1]])),
                                   colnames(data)))
  ip.mat$free[1,] <- TRUE
  ip.mat$labels[1,] <- 'slope'
  ip.mat$values[3,] <- correct.mat['g',]
  ip.mat$labels[3,] <- paste0('g', 1:ncol(ip.mat))
  ip.mat$free[4,] <- FALSE
  
  m1 <- mxModel(model="model1", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items),
                mxFitFunctionML())

  m2 <- mxModel("3pl", m1,
                univariatePrior("logit-norm", paste0('g', 1:ncol(ip.mat)), correct.mat['g',]),
                mxFitFunctionMultigroup(c('model1', 'univariatePrior')),
                mxComputeEM('model1.expectation', 'scores', mxComputeNewtonRaphson()))
  m2
}

# ----------------------------------------------------------------------------

if (file.exists("inst/models/enormous/lib/stderrlib.R")) {
  source("inst/models/enormous/lib/stderrlib.R")
} else if (file.exists("lib/stderrlib.R")) {
  source("lib/stderrlib.R")
} else {
  stop("Cannot find stderrlib.R")
}

#got <- MCphase(mkmodel, reps=5, verbose=TRUE, maxCondNum=1e6)

name = "ifa-3pl-se"
recompute=FALSE
getMCdata(name, mkmodel, correct.mat[correct.mask], maxCondNum=5000, recompute=recompute)

omxCheckCloseEnough(norm(mcBias, "2"), .5525, .001)
omxCheckCloseEnough(max(abs(mcBias)), .3055, .001)
omxCheckCloseEnough(log(det(mcHessian)), 90.31, .1)

#-----------------
if (0) {
  result <- NULL
  
  set.seed(3)
  model <- mkmodel()
  model <- mxRun(mxModel(model, mxComputeSequence(list(
    mxComputeOnce('model1.expectation', 'scores'),
    mxComputeOnce('fitfunction', 'information', 'hessian'),
    mxComputeReportDeriv()
  ))))
  h1 <- model$output$hessian
  h2 <- model$output$hessian
  h3 <- model$output$hessian
  
  em <- model$compute
  fitfun <- em$mstep$fitfunction
  em$accel <- ""
  em$tolerance <- 1e-11
  for (semType in c('mr','tian')) {
    em$information <- "mr1991"
    em$infoArgs <- list(fitfunction=fitfun, semMethod=semType, semTolerance=sqrt(1e-6))
    plan <- mxComputeSequence(list(
      em,
      mxComputeHessianQuality(),
      mxComputeStandardError(),
      mxComputeReportDeriv()
    ))
    model$compute <- plan
    fit <- mxRun(model, silent=TRUE)
    result <- rbind(result, data.frame(rd=(fit$output$standardErrors - mcSE)/mcSE,
                                       algo=semType, param=names(mcSE)))
  }
  if (1) {
    em$accel <- 'ramsay1975'
    em$tolerance <- 1e-11
    em$information <- "mr1991"
    em$infoArgs <- list(fitfunction=fitfun, semMethod="agile")
    plan <- mxComputeSequence(list(
      em,
      mxComputeHessianQuality(),
      mxComputeStandardError(),
      mxComputeReportDeriv()
    ))
    if (is.null(fit)) fit <- model
    fit$compute <- plan
    # reuse the MLE, if possible
    fit <- mxRun(fit, silent=TRUE)
    result <- rbind(result, data.frame(rd=(fit$output$standardErrors - mcSE)/mcSE,
                                       algo="agile", param=names(mcSE)))
  }
  library(ggplot2)
  ggplot(result, aes(rd, param, color=algo, shape=algo)) + geom_point() + geom_vline(xintercept=0)
}

#-----------------

rda <- paste(name, "-result.rda", sep="")
load(rda)
detail <- testPhase(mkmodel, 500,
                    methods=c('re', 'estepH', 'mr', 'tian', 'agile', 'meat', 'oakes'))

if (0) {
  asem <- studyASEM(mkmodel)
  smooth <- checkSmoothness(mkmodel)
}

save(detail, asem, smooth, file=rda)

stop("done")

if (0) {
  asem <- studyASEM(mkmodel, reps = 50, targets=seq(-8.1, -3.9, .2))
  
  got <- checkSmoothness(mkmodel())
}

if (0) {
  seed <- 1
  set.seed(seed)
  m1 <- mkmodel()
  f1 <- paste("ifa-3pl-se", seed, ".csv", sep="")
  write.table(sapply(m1$data$observed, unclass)-1, f1, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  m1 <- mxRun(m1)
  covname <- "~/ifa/ifa-3pl-se/3pl-se-cov.txt"
  H <- read.csv(covname, header = FALSE, row.names=NULL)
  H <- H[,-ncol(H)]
  H <- solve(H)
  labels <- names(m1$output$estimate)
  fmOrder <- c(labels[-c(2,9,16,23,30)], labels[c(2,9,16,23,30)])
  fmPerm <- match(colnames(mcHessian), fmOrder)
  H <- H[fmPerm, fmPerm]
  H <- 2 * H
  
  m1$compute <- mxComputeSequence(list(
    mxComputeOnce('expectation','scores'),
    mxComputeOnce('fitfunction', 'information', "hessian"),
    mxComputeReportDeriv()))
  m1 <- mxRun(m1)
  

  mvn_KL_D(H, solve(H))
  mvn_KL_D(junk, solve(junk))
}

if (0) {
  require(mirt)
  seed <- 1
  set.seed(seed)
  m1 <- mkmodel()
  val <- mirt(sapply(m1$data$observed, unclass)-1, 1, "3PL", pars="values")
  val[val$name == "a1", 'value'] <- 8
  val[val$name == "a1", 'est'] <- FALSE
  gnum <- val[val$name=='g', 'parnum']
  fit <- mirt(sapply(m1$data$observed, unclass)-1, 1, "3PL", pars=val,
       constrain=list(gnum[1:6], gnum[7:12], gnum[13:18], gnum[19:24], gnum[25:30]),
       SE=TRUE, SE.type="complete")
  fit$information
}

if (0) {
  grp <- bank[[1]]$grp
  source("~/2012/sy/irtplot.R")
  
  booklet(function(item) {
    rpf.plot(grp, item, data.bins=30, basis=c(1), factor=1)
  }, colnames(correct.mat), output="4plm.pdf")
  
  rpf.plot(grp, "i1")
  rpf.plot(grp, "i6")
  rpf.plot(grp, "i11")
  rpf.plot(grp, "i16")
  grp$param <- correct.mat
}

if (0) {
  m1 <- mkmodel(1)
  if (0) {
    raw <- sapply(m1$data$observed, unclass) - 1
    write.table(raw, file="3pl1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  
}

if(1) {
  omxCheckCloseEnough(max(sapply(bank, function(t) t$condnum['meat'])), 733, 1)

  omxCheckCloseEnough(max(abs(emp$bias)), .04, .01)
  omxCheckCloseEnough(norm(check.se("meat"), "2"), 2.08, .01)
  omxCheckCloseEnough(norm(check.se("sw"), "2"), 3.298, .01)
  omxCheckCloseEnough(norm(check.se("sem"), "2"), 2.935, .01)
  
  if (0) {
    plot(emp$bias, correct.mat[correct.mask])
    hist(check.se("sem"))
    hist(check.se("meat"))
    hist(check.se("sw"))
  }
}
