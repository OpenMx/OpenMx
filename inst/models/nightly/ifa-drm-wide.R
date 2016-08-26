#options(error = utils::recover)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 1024 + 512
i1 <- rpf.drm(multidimensional=TRUE)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1, version=1)
  correct[[ix]][3] <- logit(0)
  correct[[ix]][4] <- logit(1)
}
correct.mat <- simplify2array(correct)

numPersons <- 2500
data <- rpf.sample(numPersons, items, correct.mat)

if (0) {
  write.table(sapply(data, unclass)-1, "drm-wide.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  q()
}

ip.mat <- mxMatrix(name="item", nrow=4, ncol=numItems,
                   values=c(1,0, logit(0), logit(1)),
                   free=c(TRUE, TRUE, FALSE, FALSE))
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c('f1', 'b', 'g', 'u')

m2 <- mxModel(model="drm1", ip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=items, debugInternal=TRUE),
              mxFitFunctionML(),
	      mxComputeSequence(list(
		  mxComputeOnce('expectation', 'scores'),
		  mxComputeReportExpectation()
	      )))
m2 <- mxRun(m2)

# cat(deparse(fivenum(round(m2$expectation$debug$patternLikelihood,2))))
omxCheckCloseEnough(fivenum(m2$expectation$debug$patternLikelihood),
                    c(-1066.98, -1050.605, -996.76, -893.97, -351.4), .01)
omxCheckCloseEnough(sum(m2$expectation$debug$em.expected), numItems * numPersons, .1)

m2 <- mxModel(m2,
              mxExpectationBA81(ItemSpec=items),
              mxComputeEM('expectation', 'scores',
                          mxComputeNewtonRaphson(), verbose=0L))

m2 <- mxRun(m2)

#print(m2$matrices$item$values)
#print(correct.mat)
omxCheckCloseEnough(m2$fitfunction$result, 4045796.23, .1)

omxCheckCloseEnough(cor(c(m2$item$values[1:2,]), c(correct.mat[1:2,])), 1, .01)

emstat <- m2$compute$output
omxCheckCloseEnough(emstat$EMcycles, 45, 3)
omxCheckCloseEnough(emstat$totalMstep, 119, 5)

print(m2$output$backendTime)
