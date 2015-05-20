library(OpenMx)

set.seed(159)
xdat <- data.frame(a=rnorm(10, mean=4.2), b=1:10) # Make data set
xdat$a[3] <- NA
xdat$b[5] <- NA

xdat

cmod <- mxModel(
  name='Estimation Row Model with Missingness',
  mxData(observed=xdat, type='raw'),
  mxMatrix(values=c(4,5), ncol=2, nrow=1, free=TRUE, name='M'),
  mxAlgebra(omxSelectCols(M, existenceVector), name='fM'),
  mxAlgebra((filteredDataRow-fM)%^%2, name='rowAlgebra'),
  mxAlgebra(sum(rowResults), name='reduceAlgebra'),
  mxAlgebra(sum(log(M)), name="sumlog"),
  mxMatrix(nrow=0, ncol=0, name="M0"),
  mxConstraint(M0 == M0, name="empty"),
  mxFitFunctionRow(
    rowAlgebra='rowAlgebra',
    reduceAlgebra='reduceAlgebra',
    dimnames=c('a', 'b'))
)

cmodFit <- omxCheckWarning(mxRun(cmod),"Constraint 'Estimation Row Model with Missingness.empty' evaluated to a 0x0 matrix and will have no effect")
if (0) {
	cmodFit$M$values
	colMeans(xdat, na.rm=T)
}

omxCheckCloseEnough(as.vector(mxEval(M, cmodFit)), as.vector(colMeans(xdat, na.rm=T)), epsilon=1e-4)
omxCheckCloseEnough(as.vector(mxEval(sumlog, cmodFit)), sum(log(colMeans(xdat, na.rm=T))), epsilon=1e-4)
