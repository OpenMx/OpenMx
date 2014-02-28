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
	correct[[ix]] <- rpf.rparam(items[[ix]])
	correct[[ix]][[4]] <- 0   # no guessing, for now
	correct[[ix]][[5]] <- 1   # upper cound
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

design <- matrix(c(rep(1L,numItems),
		   rep(2L,numItems/2), rep(3L, numItems/2)), byrow=TRUE, nrow=2)

theta <- t(rmvnorm(numPersons, mean=rnorm(3, sd=.25)))
data <- rpf.sample(theta, items, correct, design)

if (0) {
  rdata <- sapply(data, unclass)-1
  # for flexMIRT, write CSV
  write.table(rdata, file="ifa-bifactor.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 0, 0, 1),
		   free=c(rep(TRUE, 3), FALSE, FALSE))

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE)

m1 <- mxModel(model="bifactor",
          ip.mat, m.mat, cov.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(mean="mean", cov="cov", debugInternal=TRUE,
			    ItemSpec=items, design=design, ItemParam="ItemParam", qpoints=29),
	      mxFitFunctionML(),
	      mxComputeOnce('expectation', 'scores'))
m1 <- mxRun(m1)

omxCheckCloseEnough(sum(m1@expectation@debug$patternLikelihood), -12629.4, .1)
omxCheckCloseEnough(fivenum(m1@expectation@debug$patternLikelihood),
                    c(-15.7575854, -14.9684791, -14.0992631, -12.3467773, -3.5902924 ), 1e-4)
omxCheckCloseEnough(sum(m1@expectation@debug$em.expected), 20000, 1)
omxCheckCloseEnough(fivenum(m1@expectation@debug$em.expected),
                    c(0, 0, 1.8e-06, 0.0034365, 43.2895967), 1e-4)

m1 <- mxModel(m1,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', 'scores'),
                mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian'))
              )))
m1 <- mxRun(m1)
omxCheckCloseEnough(m1@fitfunction@result, 2*11850.68, .01)
omxCheckCloseEnough(fivenum(m1@output$gradient), 2*c(-369.32879, -14.47296, 13.1165, 50.07066, 323.04627 ), .01)
omxCheckCloseEnough(fivenum(m1@output$hessian[m1@output$hessian != 0]),
                    2*c(-53.666201, -7.6857353, -6.0121325, 89.8735155, 192.6600613 ), 1e-4)

m1 <- mxModel(m1,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=items,
                                design=design,
                                ItemParam="ItemParam",
                                qpoints=29, scores="full"),
              mxComputeEM('expectation', 'scores',
			  mxComputeNewtonRaphson(free.set='ItemParam'),
			  mxComputeOnce('fitfunction', 'fit', free.set=c("mean","cov"))
				 ))

m1 <- mxRun(m1, silent=TRUE)
#print(correct.mat)
#print(m1@matrices$ItemParam@values)
omxCheckCloseEnough(m1@output$minimum, 20859.87, .01)
got <- cor(c(m1@matrices$ItemParam@values), c(correct.mat))
omxCheckCloseEnough(got, .977, .01)
scores <- m1@expectation@output$scores
omxCheckCloseEnough(cor(c(scores[,1]), c(theta[1,])), .758, .01)
omxCheckCloseEnough(cor(c(scores[,2]), c(theta[2,])), .781, .01)
omxCheckCloseEnough(cor(c(scores[,3]), c(theta[3,])), .679, .01)

omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 2*scores[,5]), 933, 5)
omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 3*scores[,5]), 1000, 2)

i1 <- mxModel(m1,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)

#cat(deparse(round(i1@output$standardErrors,3)))
se <- c(0.195, 0.275, 0.129, 0.12, 0.123, 0.083, 0.336, 0.196,  0.137, 0.157,
        0.148, 0.126, 0.262, 0.223, 0.194, 0.18, 0.215,  0.148, 0.223, 0.314,
        0.286, 0.135, 0.135, 0.103, 0.246, 0.361,  0.139, 0.121, 0.118, 0.079,
        0.123, 0.138, 0.095, 0.141, 0.159,  0.111, 0.161, 0.136, 0.094, 0.144, 0.155,
        0.096, 0.141, 0.147,  0.101, 0.154, 0.219, 0.169, 0.212, 0.177, 0.172, 0.152, 0.197,
        0.125, 0.179, 0.178, 0.107, 0.151, 0.135, 0.101)
omxCheckCloseEnough(i1@output$conditionNumber, 59, 1)
omxCheckCloseEnough(c(i1@output$standardErrors), se, .001)

i2 <- mxModel(m1,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', 'information', "sandwich"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i2 <- mxRun(i2, silent=TRUE)
omxCheckCloseEnough(i2@output$conditionNumber, 173, 1)

swse <- c(0.254, 0.333, 0.13, 0.133, 0.145, 0.082, 0.462, 0.289,  0.169, 0.161,
          0.199, 0.12, 0.282, 0.259, 0.171, 0.208, 0.208,  0.152, 0.275, 0.281,
          0.255, 0.14, 0.145, 0.104, 0.309, 0.382,  0.139, 0.128, 0.138, 0.08,
          0.131, 0.141, 0.097, 0.169, 0.17,  0.113, 0.189, 0.156, 0.096, 0.166,
          0.175, 0.095, 0.17, 0.189,  0.105, 0.171, 0.21, 0.165, 0.236, 0.262,
          0.168, 0.171, 0.22,  0.129, 0.205, 0.206, 0.11, 0.159, 0.163, 0.103)

#cat(deparse(round(i2@output$standardErrors,3)))
omxCheckCloseEnough(c(i2@output$standardErrors), swse, .001)
