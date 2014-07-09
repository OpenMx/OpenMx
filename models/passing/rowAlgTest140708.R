library(OpenMx)

set.seed(159)
xdat <- data.frame(a=rnorm(10, mean=4.2), b=1:10) # Make data set
xdat$a[3] <- NA
xdat$b[5] <- NA

xdat

cmod <- mxModel(
  name='Estimation Row Model with Missingness',
  mxData(observed=xdat, type='raw'),
  mxMatrix(values=.75, ncol=2, nrow=1, free=TRUE, name='M'),
  mxAlgebra(omxSelectCols(M, existenceVector), name='fM'),
  mxAlgebra((filteredDataRow-fM)%^%2, name='rowAlgebra'),
  mxAlgebra(sum(rowResults), name='reduceAlgebra'),
  mxFitFunctionRow(
    rowAlgebra='rowAlgebra',
    reduceAlgebra='reduceAlgebra',
    dimnames=c('a', 'b'))
)

cmodFit <- mxRun(cmod)
if (0) {
	cmodFit$M$values
	colMeans(xdat, na.rm=T)
}
omxCheckCloseEnough(as.vector(mxEval(M, cmodFit)), as.vector(colMeans(xdat, na.rm=T)), epsilon=10^(-5))
