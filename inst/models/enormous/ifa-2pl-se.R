#options(error = browser)
require(OpenMx)
require(rpf)

numItems <- 5  # any multiple of 5
numPersons <- 1000

spec <- list()
spec[1:numItems] <- rpf.grm()
correct <- matrix(NA, 2, numItems)
dimnames(correct) <- list(names(rpf.rparam(spec[[1]])),
                          paste("i", 1:numItems, sep=""))
correct['b',] <- seq(-1.5, 1.5, length.out = 5)
correct['a',] <- seq(.5, 4, length.out=numItems)
correct['b',] <- correct['b',] * correct['a',]

mkmodel <- function() {
  maxParam <- max(vapply(spec, rpf.numParam, 0))
  maxOutcomes <- max(vapply(spec, function(i) i$outcomes, 0))
  
  data <- rpf.sample(numPersons, spec, correct)

  ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                     values=c(1,0), free=TRUE,
                     dimnames=list(names(rpf.rparam(spec[[1]])),
                                   colnames(data)))
  
  m1 <- mxModel(model="drm", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=spec),
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

#got <- MCphase(mkmodel, reps=5, verbose=TRUE)

name <- paste("ifa-2pl", numItems, "-se", sep="")
getMCdata(name, mkmodel, correct, maxCondNum=1e7, recompute=FALSE)

if (numItems == 5) {
  omxCheckCloseEnough(norm(mcBias, "2"), .8423, .001)
  omxCheckCloseEnough(max(abs(mcBias)), .665, .001)
  omxCheckCloseEnough(log(det(mcHessian)), 35.195, .1)
} else if (numItems == 20) {
  omxCheckCloseEnough(norm(mcBias, "2"), 0.3872, .001)
  omxCheckCloseEnough(max(abs(mcBias)), 0.2757, .001)
  omxCheckCloseEnough(log(det(mcHessian)), 182.08, .1)
}

detail <- testPhase(mkmodel, 500,
                    methods=c('estepH', 'meat', 're', 'mr', 'tian', 'agile', 'oakes'))
asem <- studyASEM(mkmodel)
smooth <- checkSmoothness(mkmodel)

rda <- paste(name, "-result.rda", sep="")
save(detail, asem, smooth, file=rda)

stop("done")

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
