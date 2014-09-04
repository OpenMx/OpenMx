#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(2)
numItems <- 16
numPersons <- 1500
slope <- 8
groupSize <- 2

items <- list()
items[1:numItems] <- rpf.drm()
correct.mat <- matrix(NA, 4, numItems)
dimnames(correct.mat) <- list(names(rpf.rparam(items[[1]])), paste("i", 1:numItems, sep=""))
correct.mat['a',] <- slope
correct.mat['b',] <- slope * seq(-1.5, 1.5, length.out = numItems)
correct.mat['g',] <- kronecker(logit(1/3:(2+numItems/groupSize)), rep(1,groupSize))
correct.mat['u',] <- logit(1)

correct.mask <- matrix(FALSE, 4, numItems)
correct.mask[1,1] <- TRUE
correct.mask[2,] <- TRUE
correct.mask[3,] <- rep(c(TRUE, rep(FALSE, groupSize-1)), numItems/groupSize)

mkmodel <- function() {
  maxParam <- max(vapply(items, rpf.numParam, 0))
  maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))
  
  data <- rpf.sample(numPersons, items, correct.mat)

  ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                     values=c(1, 0, logit(.1), logit(1)), free=TRUE,
                     dimnames=list(names(rpf.rparam(items[[1]])),
                                   colnames(data)))
  ip.mat$free[1,] <- TRUE
  ip.mat$labels[1,] <- 'slope'
  ip.mat$free[4,] <- FALSE
  ip.mat$labels[3,] <- paste('g', kronecker(1:(numItems/groupSize),rep(1,groupSize)), sep='')
#  ip.mat$labels
  
  m1 <- mxModel(model="drm", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items),
                mxFitFunctionML(),
                mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()))
  m1
}

# ----------------------------------------------------------------------------

if (file.exists("models/enormous/lib/stderrlib.R")) {
  source("models/enormous/lib/stderrlib.R")
} else if (file.exists("lib/stderrlib.R")) {
  source("lib/stderrlib.R")
} else {
  stop("Cannot find stderrlib.R")
}

#got <- MCphase(mkmodel, reps=5, verbose=TRUE, maxCondNum=1e6)

if (0) {
  # interesting model, converges but is not identified
  set.seed(18)
  m1 <- mkmodel()
  em <- m1$compute
#  em$verbose = 3L
  em$tolerance = 1e-5
  em$information = "mr1991"
  em$infoArgs = list(fitfunction="fitfunction", semMethod="agile")
  m1$compute <- mxComputeSequence(list(em,
                                       mxComputeHessianQuality()))
  m1 <- mxRun(m1)
  m1$output$conditionNumber
  stop("stopped")
}

name = "ifa-3pl-se"
getMCdata(name, mkmodel, correct.mat[correct.mask], maxCondNum=5000)

omxCheckCloseEnough(norm(mcBias, "2"), .39367, .001)
omxCheckCloseEnough(max(abs(mcBias)), .1694, .001)
omxCheckCloseEnough(log(det(mcHessian)), 90.89, .1)

detail <- testPhase(mkmodel, 500,
                    methods=c('re', 'estepH', 'mr', 'tian', 'agile', 'meat', 'oakes'))

asem <- studyASEM(mkmodel)
smooth <- checkSmoothness(mkmodel)

rda <- paste(name, "-result.rda", sep="")
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
  fit@information
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
