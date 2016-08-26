#options(error = utils::recover)
library(OpenMx)
library(rpf)

mcar <- function(data, pct) {
	size <- prod(dim(data))
	erase <- rep(TRUE, size * pct)
	mask <- c(erase, rep(FALSE, size - length(erase)))[order(runif(size))]
	shaped.mask <- array(mask, dim=dim(data))
	data[shaped.mask] <- NA
	data <- data[apply(is.na(data), 1, sum) != dim(data)[2],]  # remove when all items are missing
	data
}

set.seed(8)

numItems <- 5
i1 <- rpf.nrm(3, T.c=diag(2))
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1, version=1)
}
correct.mat <- simplify2array(correct)
correct.mat[2,] <- 1
correct.mat[3,] <- 0

good.data <- rpf.sample(500, items, correct.mat)
data <- mcar(good.data, 1/3)
for (cx in 1:length(data)) attr(data[[cx]], 'mxFactor') <- TRUE
#head(data)

if (0) {
  # for flexMIRT, write CSV
  data.fm <- sapply(data, unclass)-1
  data.fm[is.na(data.fm)] <- -9
  write.table(data.fm, file="md.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ip.mat <- mxMatrix(name="item", nrow=5, ncol=numItems,
                   values=c(1,1,0,0,0),
                   free=c(TRUE,FALSE,FALSE,TRUE,TRUE))
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c('f1', paste('a',1:2,sep=""), paste('c',1:2,sep=""))

if (1) {
#  fm <- read.flexmirt("~/irt/ifa-missingdata/ifa-md-prm.txt")
  fm.est <- structure(c(0.906661, 1, 0, -0.66474, 0.523485, 0.916341, 1,  0, -3.285, -0.882019, 0.73849, 1, 0, -1.14314, -0.0300753, 0.617796,  1, 0, -0.58211, 1.4276, 2.50623, 1, 0, 0.541075, 2.1527), .Dim = c(5L,  5L))
  fm.est.mat <- mxMatrix(name="item", nrow=5, ncol=numItems, values=fm.est)
  colnames(fm.est.mat) <- colnames(data)
  rownames(fm.est.mat) <- c('f1', paste('a',1:2,sep=""), paste('c',1:2,sep=""))
  cModel <- mxModel(model="test3", fm.est.mat,
                    mxData(observed=data, type="raw"),
                    mxExpectationBA81(ItemSpec=items),
                    mxFitFunctionML(),
                    mxComputeSequence(steps=list(
                      mxComputeOnce('fitfunction', 'fit'))))
  cModel <- mxRun(cModel)
  omxCheckCloseEnough(cModel$fitfunction$result, 2733.844, .001)
}

if (1) {
  m2 <- mxModel(model="test3", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items, debugInternal=TRUE),
                mxFitFunctionML(),
		mxComputeSequence(list(
		    mxComputeOnce('expectation', 'scores'),
		    mxComputeReportExpectation()
		)))
  m2 <- mxRun(m2)
  omxCheckCloseEnough(fivenum(m2$expectation$debug$patternLikelihood),
                      c(-6.20, -4.41, -3.62, -2.077, -1.), .01)
  omxCheckCloseEnough(match(fivenum(m2$expectation$debug$patternLikelihood), m2$expectation$debug$patternLikelihood),
                      c(82L, 24L, 9L, 5L, 7L))
  omxCheckCloseEnough(sum(m2$expectation$debug$em.expected), 1667, .1)
}

plan <- mxComputeSequence(list(
  mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()),
  mxComputeOnce('fitfunction', 'gradient'),
  mxComputeReportDeriv()))

m2 <- mxModel(model="test3", ip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=items),
              mxFitFunctionML(),
              plan)

# 	m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
# 	m2 <- mxOption(m2, "Verify level", '-1')
# m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

omxCheckCloseEnough(max(abs(m2$output$gradient)), 0, .028)

grp <- as.IFAgroup(m2)

if (0) {
  # only includes rows without missingness!
  got <- rpf.SitemFit(grp, method="pearson")
  sapply(got, function (r) r$pval)
  got <- rpf.SitemFit(grp, method="rms")
  sapply(got, function (r) r$pval)
}

emstat <- m2$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 24, 1)
omxCheckCloseEnough(emstat$totalMstep, 70, 5)
omxCheckCloseEnough(m2$output$evaluations, 148, 5)

omxCheckCloseEnough(m2$output$minimum, 2733.845, .01)
got <- cor(c(m2$matrices$item$values[c(1,4,5),]),
           c(correct.mat[c(1,4,5),]))
omxCheckCloseEnough(got, .9824, .01)

omxCheckTrue(all(abs(m2$matrices$item$values[c(1,4,5),] - fm.est[c(1,4,5),]) < .025))

omxCheckCloseEnough(summary(m2)$informationCriteria['AIC:','par'], 2763.844, .01)
omxCheckCloseEnough(summary(m2)$informationCriteria['BIC:','par'], 2826.943, .01)

if (0) {
  require(mirt)
  rdata <- sapply(data, unclass)-1
  pars <- mirt(rdata, 1, itemtype="gpcm", D=1, quadpts=49, SE=TRUE)
  # Iteration: 64, Log-Lik: -1366.922, Max-Change: 0.00009 (old)    -2 * -1366.922 = 2733.844
}

refModels <- mxRefModels(m2, run=TRUE)
omxCheckCloseEnough(refModels[['Independence']]$output$fit, 2926.20, .01)

# -------------- MAP -------------
require(mvtnorm)
mapFn <- function(grp, row, theta) {
  dims <- 1
  mean <- grp$mean
  if (is.null(mean)) mean <- rep(0, dims)
  cov <- grp$cov
  if (is.null(cov)) cov <- diag(dims)
  prob <- sapply(1:ncol(grp$param), function(ix) {
    rpf.logprob(grp$spec[[ix]], grp$param[,ix], theta)[unclass(grp$data[row,ix])]
  })
#  print(prob)
  sum(prob, na.rm=TRUE) + dmvnorm(theta, mean, cov, log = TRUE)
}
if (0) {
  mapScore <- matrix(NA, nrow(grp$data), 2)
  for (row in 1:nrow(grp$data)) {
    got <- nlm(function(x) -mapFn(grp, row, x), 0, hessian = TRUE)
    mapScore[row,] <- c(got$estimate, sqrt(solve(got$hessian)))
  }
}
#plot(Vectorize(function(x) mapFn(grp, 6, x)), -5,5)
