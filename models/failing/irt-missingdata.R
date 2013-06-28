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
#head(data)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", nrow=3, ncol=numItems,
                   values=c(1,0,0),
		   lbound=c(1e-6, -1e6, -1e6),
                   free=TRUE)

m2 <- mxModel(model="test3", ip.mat, spec,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="itemParam"),
              mxFitFunctionBA81())
m2 <- mxOption(m2, "Analytic Gradients", 0)
if (1) {
	m2 <- mxOption(m2, "Analytic Gradients", 'yes')
	m2 <- mxOption(m2, "Verify level", '-1')
}
m2 <- mxOption(m2, "Calculate Hessian", "No")
m2 <- mxOption(m2, "Standard Errors", "No")
m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

got <- cor(c(m2@matrices$itemParam@values),
           c(t(correct.mat)))
omxCheckCloseEnough(got, .894, .01)
