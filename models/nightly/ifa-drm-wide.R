#options(error = utils::recover)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 1024
i1 <- rpf.drm(multidimensional=TRUE)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
  correct[[ix]][3] <- 0
  correct[[ix]][4] <- 1
}
correct.mat <- simplify2array(correct)

ability <- rnorm(500)
data <- rpf.sample(ability, items, correct.mat)

if (0) {
  write.table(sapply(data, unclass)-1, "drm-wide.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  q()
}

ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                   values=c(1,0,0, 1),
                   free=c(TRUE, TRUE, FALSE, FALSE))

eip.mat <- mxAlgebra(itemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="drm1", ip.mat, m.mat, cov.mat, eip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec=items, EItemParam="EItemParam",ItemParam="itemParam",
                mean="mean", cov="cov"),
              mxFitFunctionML(),
              mxComputeOnce('expectation', context='EM')
              )
m2 <- mxRun(m2)
omxCheckCloseEnough(fivenum(m2@expectation@patternLikelihood),
                    c(-712.0873, -701.8445, -664.7972, -596.581, -225.9732), .01)
omxCheckCloseEnough(sum(m2@expectation@em.expected), 512000, .1)

m2 <- mxModel(m2,
              mxData(observed=data, type="raw"),
	      mxComputeIterate(steps=list(
				 mxComputeOnce('expectation', context='EM'),
				 mxComputeNewtonRaphson(free.set='itemParam'),
				 mxComputeOnce('expectation'),
				 mxComputeOnce('fitfunction', free.set=c("mean","cov"))
				 )))
m2 <- mxRun(m2)

#print(m2@matrices$itemParam@values)
#print(correct.mat)
# Matches flexMIRT but maybe unstable
omxCheckCloseEnough(m2@fitfunction@result, 537810.8, .01)
