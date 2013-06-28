# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(mirt)
require(OpenMx)
require(rpf)
library(stringr)

data(SAT12)
data <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
df <- list()
for (cx in 1:dim(data)[2]) {
	df[[cx]] <- factor(data[,cx], levels=c(0,1))
}
df <- as.data.frame(df)
colnames(df) <- str_replace_all(colnames(data), "\\.", "")

numItems <- 32
numPersons <- dim(data)[1]
maxDim <- 2

items <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(dimensions=maxDim)
}

maxParam <- max(vapply(items, function(i) i@numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@numOutcomes, 0))

design <- matrix(c(rep(1,numItems),
		   1+c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)), byrow=TRUE, nrow=2)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=c(rep(mxLookupIRTItemModelID("drm"), numItems),
		 rep(2, numItems),
		 rep(2, numItems)),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(
		     rep(1.414,numItems),
		     rep(1,numItems),
		     rep(0,numItems),
		     rep(0,numItems)),
         free=c(rep(TRUE, numItems*3),
                rep(FALSE, numItems*1)),
         byrow=TRUE)

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

m1 <- mxModel(model="sat12",
          spec, design,
          ip.mat,
          mxMatrix(name="A", values=c(1,1,0,0), nrow=1, ncol=maxParam, free=FALSE),
          mxAlgebra(name="prior",
                    sum(dlnorm(omxSelectRows(ItemParam, A), 0, .5, 1))),
          mxData(observed=df, type="raw"),
          mxExpectationBA81(
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam",
	     ItemPrior="prior"),
          mxFitFunctionBA81()
)

m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

m1 <- mxRun(m1, silent=TRUE)
print(m1@matrices$ItemParam@values)
