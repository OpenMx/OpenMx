#options(error = browser)
require(OpenMx)
require(rpf)

numItems <- 20
numPersons <- 2000

spec <- list()
spec[1:numItems] <- rpf.grm(outcomes = 3)
correct <- matrix(NA, 3, numItems)
dimnames(correct) <- list(names(rpf.rparam(spec[[1]])), paste("i", 1:numItems, sep=""))
correct['a',] <- seq(.5,4, length.out = numItems)
correct['b1',] <- seq(-1.5,1.5,length.out = numItems/5L)
correct['b1',] <- correct['b1',]
correct['b2',] <- correct['b1',] - 1e-1
correct['b1',] <- correct['b1',] * correct['a',]
correct['b2',] <- correct['b2',] * correct['a',]

if (0) {
  sapply(1:numItems, function(ix) {
    rpf.prob(spec[[1]], correct[,ix], -correct['b1',ix] / correct['a',ix])
  })
  plot.icc <- function(item, param, width=3) {
    require(ggplot2)
    require(reshape2)
    pm <- t(rpf.prob(item, param, seq(-width, width, .1)))
    icc <- as.data.frame(melt(pm, varnames=c("theta",'category')))
    icc$theta <- seq(-width, width, .1)
    icc$category <- as.factor(icc$category - 1)
    ggplot(icc, aes(theta, value)) +
      geom_line(aes(color=category, linetype=category)) +
      ylim(0,1) + xlim(-width,width) + labs(y="Probability", x="Theta")
  }
  plot.icc(spec[[1]], correct[,2])
}

mkmodel <- function() {
  maxParam <- max(vapply(spec, rpf.numParam, 0))
  maxOutcomes <- max(vapply(spec, function(i) i$outcomes, 0))
  
  data <- rpf.sample(numPersons, spec, correct)

  ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                     values=c(1,.5,-.5), free=TRUE,
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

#got <- MCphase(mkmodel, reps=5, verbose=TRUE, maxCondNum = NA)

name <- paste("ifa-grm", numItems, "-se", sep="")
getMCdata(name, mkmodel, correct, maxCondNum = 1e6, recompute=FALSE)

omxCheckCloseEnough(norm(mcBias, "2"), 0.22297, .001)
omxCheckCloseEnough(max(abs(mcBias)), 0.11124, .001)
omxCheckCloseEnough(log(det(mcHessian)), 368.81, .1)

rda <- paste(name, "-result.rda", sep="")
load(rda)
detail <- testPhase(mkmodel, 500,
                    methods=c('re', 'estepH', 'mr','tian', 'agile', 'meat', 'oakes'))
if (0) {
  asem <- studyASEM(mkmodel)
  smooth <- checkSmoothness(mkmodel)
}

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
