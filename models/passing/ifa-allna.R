#options(error = browser)
require(OpenMx)
require(rpf)

numItems <- 4
maxDim <- 3

items <- list()
items[1:numItems] <- rpf.grm(factors=maxDim)
correct.mat <- sapply(items, rpf.rparam, version=1)
correct.mat['a3',1:2] <- 0
correct.mat['a2',3:4] <- 0

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i$outcomes, 0))

data <- rpf.sample(4, items, correct.mat)

ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 1, 0), free=TRUE)
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c(paste("f", 1:3, sep=""), 'b')
ip.mat$values['f3',1:2] <- 0
ip.mat$values['f2',3:4] <- 0

mkmodel <- function(data) {
  mxModel(model="bifactor",
          ip.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(items, qpoints=29, debugInternal=TRUE),
          mxFitFunctionML(),
          mxComputeOnce('expectation', 'scores'))
}

tdata <- data
tdata[1,] <- NA
omxCheckError(mxRun(mkmodel(tdata)), "You have missing data. You must set minItemsPerScore")

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
correct.mat <- sapply(items, rpf.rparam, version=1)

slicen <- 50
data <- rpf.sample(slicen, items, correct.mat)
for (m in seq(.1, .9, .05)) {
  data <- rbind(data, mcar(rpf.sample(slicen, items, correct.mat), m))
}

dimnames(correct.mat) <- list(c('f1', 'b'), colnames(data))

result <- expand.grid(ips=0:numItems, v=NA)

for (r in 1:nrow(result)) {
  grp <- list(spec=items,
              param=correct.mat,
              data=data,
              minItemsPerScore=result$ips[r])

  sc <- EAPscores(grp, naAction = "pass")
  
  v <- var(sc[,'f1'], na.rm=TRUE)
  result$v[r] <- v
}

c1 <- coef(lm(v ~ ips, result))
#cat(deparse(round(c1,4)))
omxCheckCloseEnough(c1, c(0.5763, 0.0144), .001)

# ------------------------------

grp <- list(spec=items,
            param=correct.mat,
            data=data,
            minItemsPerScore=ncol(correct.mat)+1)

omxCheckError(EAPscores(grp, naAction = "pass"),
	      "minItemsPerScore (=13) cannot be larger than the number of items (=12)")
