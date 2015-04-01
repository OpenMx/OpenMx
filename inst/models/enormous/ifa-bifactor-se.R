# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(3)
numItems <- 20
numPersons <- 1500
maxDim <- 2

items <- list()
items[1:numItems] <- rpf.grm(factors=maxDim)
correct.mat <- sapply(items, rpf.rparam)
correct.mat[3,1:10] <- seq(-2,2,length.out=10)
correct.mat[3,11:20] <- seq(-2,2,length.out=10)

mkmodel <- function(seed) {
  set.seed(seed)
  
  maxParam <- max(vapply(items, rpf.numParam, 0))
  maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))
  
  design <- matrix(c(rep(1L,numItems),
                     rep(2L,numItems/2), rep(3L, numItems/2)), byrow=TRUE, nrow=2)
  
  data <- rpf.sample(numPersons, items, correct.mat, design)

  ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                     values=c(1.414, 1, 0), free=TRUE)
  colnames(ip.mat) <- colnames(data)
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
  colnames(m.mat) <- paste("f", 1:3, sep="")
  cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE,
                      dimnames=list(colnames(m.mat), colnames(m.mat)))
  
  m1 <- mxModel(model="bifactor",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, design=design, ItemParam="ItemParam",
                                  qpoints=31, qwidth=5),
                mxFitFunctionML())
  m1
}

replicate <- function(seed) {
  m1 <- mxModel(mkmodel(seed),
                mxComputeSequence(steps=list(
                  mxComputeEM('expectation', 'scores',
                              mxComputeNewtonRaphson(), tolerance=1e-5,
                              information="mr1991", infoArgs=list(fitfunction='fitfunction')),
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
  list(condnum=c(sem=ifelse(is.null(sem.condnum), NA, sem.condnum),
                 meat=i1$output$conditionNumber,
                 sw=i2$output$conditionNumber),
       got=got,
       em=m1$compute$steps[[1]]$output)
}

if (0) {
  m1 <- mxModel(mkmodel(69),
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
diff <- apply(est, 2, function(c) (c - correct.mat))
emp <- list(bias=apply(diff, 1, mean), se=apply(diff, 1, sd))

check.se <- function(type) {
  se <- sapply(bank, function (t) t$got[,type])
  mask <- apply(se, 2, function(c) all(!is.na(c)))
  rd <- apply(se[,mask], 2, function (c) (c - emp$se)/emp$se)
  apply(rd, 1, mean)
}

if(1) {
  omxCheckCloseEnough(max(sapply(bank, function(t) t$condnum['meat'])), 130, 2)

  omxCheckCloseEnough(max(abs(emp$bias)), .0294, .001)
  omxCheckCloseEnough(norm(check.se("sem"), "2"), .284, .01)
  omxCheckCloseEnough(norm(check.se("meat"), "2"), .26, .01)
  omxCheckCloseEnough(norm(check.se("sw"), "2"), .339, .01)

  #hist(check.se("sem"))
  #hist(check.se("meat"))
  #hist(check.se("sw"))
  
#  omxCheckCloseEnough(fivenum(sapply(bank, function (t) t$em$semProbeCount) / length(c(correct.mat))),
#                      c(2.45, 2.65, 2.73, 2.90, 3.00), .01)

  omxCheckCloseEnough(sum(sapply(bank, function (t) any(is.na(t$got[,'sem'])))), 0)
  #cor(t(sapply(bank, function(t) t$condnum)), use="complete.obs")
  if (0) {
    which(sapply(bank, function (t) any(is.na(t$got[,'sem']))))
  }
}
