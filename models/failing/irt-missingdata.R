library(OpenMx)
library(rpf)

mcar <- function(data, pct) {
	erase <- rep(TRUE, length(data) * pct)
	mask <- c(erase, rep(FALSE, length(data) - length(erase)))[order(runif(length(data)))]
	shaped.mask <- array(mask, dim=dim(data))
	data[shaped.mask] <- NA
	data
}

set.seed(8)

numItems <- 5
i1 <- rpf.gpcm(3)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
}
correct.mat <- t(simplify2array(correct))

good.data <- rpf.sample(250, items, correct)
data <- mcar(good.data, 1/3)
data <- data[apply(is.na(data), 1, sum) != numItems,]  # remove all missing
colnames(data) <- paste("item", 1:dim(data)[2], sep='')

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=c(rep(mxLookupIRTItemModelID("gpcm1"), numItems),
		 rep(3, numItems),  # outcomes
		 rep(1, numItems)), # dim
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam",
                   values=c(rep(1,numItems), rep(0,numItems*2)),
                   nrow=3, ncol=numItems, byrow=TRUE, free=TRUE)

m2 <- mxModel(model="test3", ip.mat, spec,
              mxMatrix(name="A", values=c(1,0,0), nrow=1, ncol=3, free=FALSE),
              mxAlgebra(name="prior",
                        sum(dlnorm(omxSelectRows(itemParam, A), 0, .5, 1))),
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="itemParam",
                ItemPrior="prior"),
              mxFitFunctionBA81())
m2 <- mxRun(m2)

got <- cor(c(m2@matrices$itemParam@values),
           c(t(correct.mat)))
omxCheckCloseEnough(got, .9424, .01)
