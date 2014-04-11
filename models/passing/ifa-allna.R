# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)

numItems <- 4
maxDim <- 2

items <- list()
items[1:numItems] <- rpf.grm(factors=maxDim)
correct.mat <- sapply(items, rpf.rparam)

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))

design <- matrix(c(rep(1L,numItems),
		   rep(2L,numItems/2), rep(3L, numItems/2)), byrow=TRUE, nrow=2)

data <- rpf.sample(4, items, correct.mat, design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 0), free=TRUE)
m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE)

mkmodel <- function(data) {
  mxModel(model="bifactor",
          ip.mat, m.mat, cov.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(mean="mean", cov="cov", debugInternal=TRUE,
                            ItemSpec=items, design=design, ItemParam="ItemParam", qpoints=29),
          mxFitFunctionML(),
          mxComputeOnce('expectation', 'scores'))
}

tdata <- data
tdata[1,] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "Data row 1 has no information about ability 1")

tdata <- data
tdata[2,1:2] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "Data row 2 has no information about ability 2")

tdata <- data
tdata[3,3:4] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "Data row 3 has no information about ability 3")

