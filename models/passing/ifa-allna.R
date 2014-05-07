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
colnames(ip.mat) <- colnames(data)
m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
colnames(m.mat) <- paste("f", 1:3, sep="")
cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE)
dimnames(cov.mat) <- list(paste("f", 1:3, sep=""), paste("f", 1:3, sep=""))

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
omxCheckError(mxRun(mkmodel(tdata)), "1:Data row 1 has no information about ability 1
2:Data row 1 has no information about ability 2
3:Data row 1 has no information about ability 3")

tdata <- data
tdata[2,1:2] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "Data row 2 has no information about ability 2")

tdata <- data
tdata[3,3:4] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "Data row 3 has no information about ability 3")

#-------------------------------------------------

mcar <- function(ret, mcar) {
  size <- prod(dim(ret))
  mask <- rep(FALSE, size)
  mask[sample.int(size, size * mcar)] <- TRUE
  shaped.mask <- array(mask, dim=dim(ret))
  ret[shaped.mask] <- NA
  ret
}

set.seed(1)

numItems <- 12
items <- list()
items[1:numItems] <- rpf.grm()
correct.mat <- sapply(items, rpf.rparam)

maxParam <- max(vapply(items, rpf.numParam, 0))

slicen <- 50
data <- rpf.sample(slicen, items, correct.mat)
for (m in seq(.1, .9, .05)) {
  data <- rbind(data, mcar(rpf.sample(slicen, items, correct.mat), m))
}

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems, values=correct.mat)
colnames(ip.mat) <- colnames(data)

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0)
colnames(m.mat) <- 'f1'
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, dimnames=list('f1','f1'))

result <- expand.grid(ips=0:numItems, v=NA)
for (r in 1:nrow(result)) {
  m1 <- mxModel(model="perScore", ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov", ItemSpec=items, ItemParam="ItemParam",
                                  naAction="pass", scores="full", minItemsPerScore=result$ips[r]),
                mxComputeOnce('expectation'))
  fit1 <- mxRun(m1, silent=TRUE)
  v <- var(fit1$expectation$output$scores[,1], na.rm=TRUE)
  result$v[r] <- v
}

c1 <- coef(lm(v ~ ips, result))
#cat(deparse(round(c1,4)))
omxCheckCloseEnough(c1, c(0.5763, 0.0144), .001)

# ------------------------------

m1 <- mxModel(model="perScore", ip.mat, m.mat, cov.mat,
	      mxData(observed=data, type="raw"),
	      mxExpectationBA81(mean="mean", cov="cov", ItemSpec=items, ItemParam="ItemParam",
				naAction="pass", scores="full", minItemsPerScore=as.integer(numItems+1)),
	      mxComputeOnce('expectation'))
omxCheckError(mxRun(m1, silent=TRUE),
	      "perScore.expectation: minItemsPerScore (=13) cannot be larger than the number of items (=12)")
