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
i1 <- rpf.nrm(3, T.c=diag(2))
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
}
correct.mat <- simplify2array(correct)
correct.mat[2,] <- 1
correct.mat[3,] <- 0

good.data <- rpf.sample(250, items, correct.mat)
data <- mcar(good.data, 1/3)
#head(data)

spec <- mxMatrix(name="ItemSpec", nrow=19, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", nrow=5, ncol=numItems,
                   values=c(1,1,0,0,0),
                   free=c(TRUE,FALSE,FALSE,TRUE,TRUE))

eip.mat <- mxAlgebra(itemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="test3", ip.mat, spec, m.mat, cov.mat, eip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec="ItemSpec",
                                EItemParam="EItemParam"),
              mxFitFunctionBA81(ItemParam="itemParam"),
              mxComputeIterate(steps=list(
                mxComputeOnce("EItemParam"),
                mxComputeOnce('expectation', context='EM'),
				   mxComputeNewtonRaphson(free.set='itemParam'),
#                mxComputeGradientDescent(free.set='itemParam'),
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction'))))

	m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
	m2 <- mxOption(m2, "Verify level", '-1')
m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

got <- cor(c(m2@matrices$itemParam@values[c(1,4,5),]),
           c(correct.mat[c(1,4,5),]))
omxCheckCloseEnough(got, .936, .01)
