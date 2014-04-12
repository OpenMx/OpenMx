#options(error = utils::recover)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 10
i1 <- rpf.drm(multidimensional=TRUE)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
  correct[[ix]][3] <- logit(0)
  correct[[ix]][4] <- logit(1)
}
correct.mat <- simplify2array(correct)
correct.mat[2,] <- correct.mat[2,] * -correct.mat[1,]

ability <- rnorm(500)
data <- rpf.sample(ability, items, correct.mat)

ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                   values=c(1,0, logit(0), logit(1)),
                   free=c(TRUE, TRUE, FALSE, FALSE))

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="drm1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(debugInternal=TRUE,
                ItemSpec=items, ItemParam="itemParam",
                mean="mean", cov="cov", qpoints=31),
              mxFitFunctionML(),
	      mxComputeOnce('expectation', 'scores'))
m2 <- mxRun(m2)
omxCheckCloseEnough(sum(m2$expectation$debug$patternLikelihood), -2032.9, .1)
omxCheckCloseEnough(fivenum(m2$expectation$debug$patternLikelihood),
                    c(-7.5454472, -7.3950031, -7.3950031, -6.9391761, -3.5411989), .001)
omxCheckCloseEnough(sum(m2$expectation$debug$em.expected), 5000, .01)
omxCheckCloseEnough(fivenum(m2$expectation$debug$em.expected),
                    c(0, 5.86e-05, 0.0687802, 7.1582354, 74.1583248), .01)

em.ex <- array(c(m2$expectation$debug$em.expected), dim=c(2,31,20))
em.tbl <- rbind(apply(em.ex[1,,], 2, sum)[1:numItems],
                apply(em.ex[2,,], 2, sum)[1:numItems])
omxCheckCloseEnough(apply(sapply(data, unclass)-1, 2, table), em.tbl, .01)

testDeriv <- mxModel(m2,
	      mxComputeIterate(list(
		  mxComputeOnce('expectation', 'scores'),
		  mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian', 'ihessian')),
      mxComputeReportDeriv()
		  )))
testDeriv <- mxRun(testDeriv)
omxCheckCloseEnough(testDeriv$fitfunction$result, 2*3221.826, .01)
omxCheckCloseEnough(fivenum(testDeriv$output$gradient), 2*c(-128.034, -8.294, 10.7, 25.814, 107.966), .01)
omxCheckCloseEnough(fivenum(testDeriv$output$hessian[testDeriv$output$hessian != 0]),
                    2*c(6.559, 6.559, 32.976, 83.554, 107.714), .01)
omxCheckCloseEnough(solve(testDeriv$output$hessian), testDeriv$output$ihessian, 1e-2)

m2 <- mxModel(m2,
              mxData(observed=data, type="raw"),  # got sorted, add it again unsorted
              mxExpectationBA81(
                ItemSpec=items, ItemParam="itemParam",
                mean="mean", cov="cov",
                qpoints=31,
                scores="full"),
	      mxComputeEM('expectation', 'scores',
	                  mxComputeNewtonRaphson()))

# 	m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
# 	m2 <- mxOption(m2, "Verify level", '-1')
# m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

emstat <- m2$compute$output
omxCheckCloseEnough(emstat$EMcycles, 12, 1)
omxCheckCloseEnough(emstat$totalMstep, 33, 5)

#print(m2$matrices$itemParam$values)
#print(correct.mat)
omxCheckCloseEnough(m2$fitfunction$result, 6216.272, .01)
got <- cor(c(m2$matrices$itemParam$values[1:2,]),
           c(correct.mat[1:2,]))
omxCheckCloseEnough(got, .988, .01)
scores <- m2$expectation$output$scores
omxCheckCloseEnough(scores[1:5,1], c(0.6783773, 0.2848123, -0.3438632, -0.1026575, -1.0820213), .001)
omxCheckCloseEnough(scores[1:5,2], c(0.6769653, 0.6667262, 0.6629124, 0.6624804, 0.6796952), 1e-4)
omxCheckCloseEnough(scores[,1], as.vector(ability), 3.5*max(scores[,2]))
omxCheckCloseEnough(cor(c(scores[,1]), ability), .737, .01)

#mxOption(NULL, 'loglikelihoodScale', -2)
i1 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', 'scores'),
                mxComputeOnce('fitfunction', 'information', "hessian"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)

#cat(deparse(round(i1$output$standardErrors,3)))
se <- c(0.11, 0.102, 0.141, 0.131, 0.109, 0.097, 0.118, 0.099,  0.095, 0.092, 0.124,
        0.112, 0.105, 0.095, 0.118, 0.108, 0.102,  0.094, 0.111, 0.11)
omxCheckCloseEnough(c(i1$output$standardErrors), se, .01)

i1 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)
se <- c(0.166, 0.111, 0.253, 0.17, 0.171, 0.104, 0.199, 0.11,  0.138, 0.095,
        0.195, 0.128, 0.158, 0.102, 0.192, 0.123, 0.149,  0.099, 0.153, 0.114)
omxCheckCloseEnough(c(i1$output$standardErrors), se, .001)

i2 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('fitfunction', 'information', "sandwich"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i2 <- mxRun(i2, silent=TRUE)

omxCheckCloseEnough(i2$output$conditionNumber, 11, .1)
#cat(deparse(round(i2$output$standardErrors,3)))
swse <- c(0.161, 0.109, 0.238, 0.161, 0.166, 0.104, 0.194,  0.109, 0.131,
          0.094, 0.213, 0.13, 0.167, 0.101, 0.188, 0.12,  0.157, 0.098, 0.157, 0.113)
omxCheckCloseEnough(c(i2$output$standardErrors), swse, .001)
