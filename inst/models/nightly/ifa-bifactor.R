# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(5)

numItems <- 20
numPersons <- 1000
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(factors=maxDim)
	correct[[ix]] <- rpf.rparam(items[[ix]], version=1)
	correct[[ix]][[4]] <- logit(0)   # no guessing, for now
	correct[[ix]][[5]] <- logit(1)   # upper cound
}
correct.mat <- simplify2array(correct)

# unpack factor structure

correct.mat <- rbind(correct.mat[1:2,], a3=0, correct.mat[3:5,])
correct.mat['a3',11:20] <- correct.mat['a2',11:20]
correct.mat['a2',11:20] <- 0
maxDim <- 3
items[1:numItems] <- rpf.drm(factors=maxDim)

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))

theta <- t(rmvnorm(numPersons, mean=rnorm(3, sd=.25)))
data <- rpf.sample(theta, items, correct.mat)

if (0) {
  rdata <- sapply(data, unclass)-1
  # for flexMIRT, write CSV
  write.table(rdata, file="ifa-bifactor.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 1, 0, logit(0), logit(1)),
		   free=c(rep(TRUE, 4), FALSE, FALSE))
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c(paste('f', 1:3, sep=""), 'b', 'g', 'u')
ip.mat$values['f3', 1:10] <- 0
ip.mat$free['f3', 1:10] <- FALSE
ip.mat$values['f2', 11:20] <- 0
ip.mat$free['f2', 11:20] <- FALSE

m1 <- mxModel(model="bifactor",
          ip.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(ItemSpec=items, qpoints=29, verbose=0L, debugInternal=TRUE),
	      mxFitFunctionML(),
	      mxComputeSequence(list(
		  mxComputeOnce('expectation', 'scores'),
		  mxComputeReportExpectation()
	      )))
m1 <- mxRun(m1)

omxCheckCloseEnough(fivenum(m1$expectation$debug$patternLikelihood),
                    c(-15.75, -14.96, -13.78, -12.34, -3.59), 1e-2)
omxCheckCloseEnough(sum(m1$expectation$debug$em.expected), 20000, 1)
omxCheckCloseEnough(fivenum(m1$expectation$debug$em.expected),
                    c(0, 0, 1.8e-06, 0.0034365, 43.2895967), 1e-4)

m1 <- mxModel(m1,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', 'scores'),
                mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian')),
                mxComputeReportDeriv()
              )))
m1 <- mxRun(m1)
omxCheckCloseEnough(m1$fitfunction$result, 2*11850.68, .01)
omxCheckCloseEnough(fivenum(m1$output$gradient), 2*c(-369.32879, -14.47296, 13.1165, 50.07066, 323.04627 ), .01)
omxCheckCloseEnough(fivenum(m1$output$hessian[m1$output$hessian != 0]),
                    2*c(-53.666201, -7.6857353, -6.0121325, 89.8735155, 192.6600613 ), 1e-4)

m1 <- mxModel(m1,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=items, qpoints=29),
              mxComputeEM('expectation', 'scores',
                          mxComputeNewtonRaphson(freeSet='item')))

m1 <- mxRun(m1, silent=TRUE)

emstat <- m1$compute$output
omxCheckCloseEnough(emstat$EMcycles, 42, 2)
omxCheckCloseEnough(emstat$totalMstep, 104, 10)
omxCheckCloseEnough(m1$output$evaluations, 211, 5)

#print(correct.mat)
#print(m1$matrices$item$values)
omxCheckCloseEnough(m1$output$minimum, 20859.87, .01)
mask <- is.finite(correct.mat)
got <- cor(c(m1$matrices$item$values[mask]), c(correct.mat[mask]))
omxCheckCloseEnough(got, .977, .01) 

grp <- as.IFAgroup(m1)

scores <- EAPscores(grp)
omxCheckCloseEnough(cor(c(scores[,1]), c(theta[1,])), .758, .01)
omxCheckCloseEnough(cor(c(scores[,2]), c(theta[2,])), .781, .01)
omxCheckCloseEnough(cor(c(scores[,3]), c(theta[3,])), .679, .01)

omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 2*scores[,5]), 933, 5)
omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 3*scores[,5]), 1000, 3)

i1 <- mxModel(m1,
              mxComputeSequence(steps=list(
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)

#cat(deparse(round(i1$output$standardErrors,3)))
se <- c(0.195, 0.275, 0.129, 0.12, 0.123, 0.083, 0.337, 0.197,  0.138, 0.157, 0.149,
        0.126, 0.262, 0.224, 0.194, 0.18, 0.215,  0.148, 0.223, 0.314, 0.286, 0.135,
        0.136, 0.103, 0.246, 0.362,  0.139, 0.121, 0.118, 0.079, 0.123, 0.138, 0.095,
        0.141, 0.159,  0.111, 0.161, 0.136, 0.094, 0.144, 0.154, 0.096, 0.141, 0.147,
        0.101, 0.154, 0.219, 0.169, 0.212, 0.176, 0.172, 0.152, 0.197,  0.125, 0.179,
        0.178, 0.107, 0.151, 0.134, 0.101)
omxCheckCloseEnough(log(i1$output$conditionNumber), 3.3, .2)
omxCheckCloseEnough(c(i1$output$standardErrors), se, .001)

## i2 <- mxModel(m1,
##               mxComputeSequence(steps=list(
##                 mxComputeOnce('fitfunction', 'information', "sandwich"),
##                 mxComputeStandardError(),
##                 mxComputeHessianQuality())))
## i2 <- mxRun(i2, silent=TRUE)
## omxCheckCloseEnough(log(i2$output$conditionNumber), 3.6, .2)

## swse <- c(0.255, 0.333, 0.13, 0.133, 0.145, 0.082, 0.464, 0.291,  0.17, 0.161, 0.199, 0.12,
##           0.283, 0.259, 0.171, 0.208, 0.208,  0.152, 0.276, 0.282, 0.255, 0.14, 0.145, 0.104,
##           0.309, 0.382,  0.139, 0.128, 0.138, 0.08, 0.13, 0.141, 0.097, 0.169, 0.17, 0.113,
##           0.189, 0.156, 0.096, 0.166, 0.175, 0.095, 0.17, 0.189, 0.105,  0.171, 0.209, 0.165,
##           0.236, 0.261, 0.168, 0.171, 0.22, 0.129,  0.205, 0.206, 0.11, 0.158, 0.163, 0.103)

## #cat(deparse(round(i2$output$standardErrors,3)))
## omxCheckCloseEnough(c(i2$output$standardErrors), swse, .01)
