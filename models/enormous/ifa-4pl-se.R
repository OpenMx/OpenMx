#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(2)
numItems <- 30
numPersons <- 2500
slope <- 8

items <- list()
items[1:numItems] <- rpf.drm()
correct.mat <- sapply(items, rpf.rparam)
correct.mat['a',] <- slope
correct.mat['b',] <- correct.mat['b',] * correct.mat['a',] / 2
correct.mat['g',] <- correct.mat['g',kronecker(1:5,rep(1,6))]
correct.mat['u',] <- correct.mat['u',kronecker(1:5,rep(1,6))]
colnames(correct.mat) <- paste("i", 1:numItems, sep="")

correct.mask <- matrix(FALSE, 4, numItems)
correct.mask[2,] <- TRUE
correct.mask[3,] <- rep(c(TRUE, rep(FALSE, 5)),5)
correct.mask[4,] <- rep(c(TRUE, rep(FALSE, 5)),5)

mkmodel <- function(seed) {
  set.seed(seed)
  
  maxParam <- max(vapply(items, rpf.numParam, 0))
  maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))
  
  data <- rpf.sample(numPersons, items, correct.mat)

  ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                     values=c(slope, 0, logit(.1), logit(.9)), free=TRUE)
  colnames(ip.mat) <- colnames(data)
  ip.mat$free[1,] <- FALSE
  ip.mat$labels[3,] <- paste('g', kronecker(1:5,rep(1,6)), sep='')
  ip.mat$labels[4,] <- paste('u', kronecker(1:5,rep(1,6)), sep='')
#  ip.mat$labels
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
  rownames(m.mat) <- "f1"
  cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE,
                      dimnames=list("f1","f1"))
  
  m1 <- mxModel(model="drm",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw", sort=FALSE),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, ItemParam="ItemParam"),
                mxFitFunctionML())
  m1
}

replicate <- function(seed) {
  m1 <- mxModel(mkmodel(seed),
                mxComputeSequence(steps=list(
                  mxComputeEM('expectation', 'scores',
                              mxComputeNewtonRaphson(),
                              tolerance=1e-5, information="mr1991",
                              infoArgs=list(fitfunction='fitfunction', semForcePD=TRUE)),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  
  m1 <- mxRun(m1, silent=TRUE)
  
  i1 <- mxModel(m1,
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', 'information', "meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i1 <- mxRun(i1, silent=TRUE)
  
  i2 <- mxModel(m1,
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', 'information', "sandwich"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i2 <- mxRun(i2, silent=TRUE)
  
  got <- cbind(m1$output$estimate,
               m1$output$standardErrors,
               i1$output$standardErrors,
               i2$output$standardErrors)
  colnames(got) <- c("est", "sem", "meat", "sw")
  
  sem.condnum <- m1$output$conditionNumber
  output <- list(condnum=c(sem=ifelse(is.null(sem.condnum), NA, sem.condnum),
                 meat=i1$output$conditionNumber,
                 sw=i2$output$conditionNumber),
       got=got,
       em=m1$compute$steps[[1]]$output)
  if (0) {
    output <- c(output,
                grp=list(spec=m1$expectation$ItemSpec,
                param=m1$ItemParam$values,
                mean=0,
                cov=matrix(1,1,1),
                scores=m1$expectation$output$scores,
                data=m1$data$observed))
  }
  output
}

if (0) {
  m1 <- mxModel(mkmodel(1),
                mxComputeSequence(steps=list(
                  mxComputeEM('expectation',
                              mxComputeNewtonRaphson(freeSet='ItemParam'),
                              mxComputeOnce('fitfunction', freeSet=c("mean","cov"), fit=TRUE),
                              information=TRUE, info.method="hessian",  semDebug=TRUE,
                              verbose=3L, agileMaxIter=3L, noiseTarget=.1),
                  mxComputeHessianQuality())))
  
  m1 <- mxRun(m1, silent=TRUE)
  print(m1$output$conditionNumber)
 stop("done") 
}

bank <- list()
for (repl in 1:500) {
  bank[[repl]] <- replicate(repl)
  print(repl)
}
if (0) {
  save(bank, file="ifa-bifactor-se.rda")
}

est <- sapply(bank, function (t) t$got[,'est'])
diff <- apply(est, 2, function(c) (c - correct.mat[correct.mask]))
emp <- list(bias=apply(diff, 1, mean), se=apply(diff, 1, sd))

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

check.se <- function(type) {
  se <- sapply(bank, function (t) t$got[,type])
  mask <- apply(se, 2, function(c) all(!is.na(c)))
  rd <- apply(se[,mask], 2, function (c) (c - emp$se)/emp$se)
  apply(rd, 1, mean)
}

if (1) {
  omxCheckCloseEnough(log(max(sapply(bank, function(t) t$condnum['meat']))), 7, 1)

  omxCheckCloseEnough(max(abs(emp$bias)), .04108, .001)
  omxCheckCloseEnough(norm(check.se("meat"), "2"), .2209, .01)
  omxCheckCloseEnough(norm(check.se("sem"), "2"), 2.366, .01)
  omxCheckCloseEnough(norm(check.se("sw"), "2"), 2.735, .01)
  
  if (0) {
    hist(check.se("meat"))
    hist(check.se("sem"))
    hist(check.se("sw"))
  }
}
