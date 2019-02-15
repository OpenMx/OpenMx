library(OpenMx)
set.seed(190127)

mu1 <- -1.5
mu2 <- 1.5

N <- 200
x <- matrix(c(rnorm(N/2,mu1,1),
              rnorm(N/2,mu2,1)),ncol=1,dimnames=list(NULL,"x"))
data4mx <- mxData(observed=x,type="raw")

class1 <- mxModel("Class1",
	mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=-.5,name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=4,name="Sigma"),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames="x"),
	mxFitFunctionML(vector=TRUE))

class2 <- mxRename(class1, "Class2")
class2$Mu$values[1,1] <- .5

mm <- mxModel(
	"Mixture", data4mx, class1, class2,
	mxAlgebra((1-Posteriors) * Class1.fitfunction, name="PL1"),
	mxAlgebra(Posteriors * Class2.fitfunction, name="PL2"),
	mxAlgebra(PL1 + PL2, name="PL"),
	mxAlgebra(PL2 / PL,  recompute='onDemand',
	          initial=matrix(runif(N,.4,.6), nrow=N, ncol = 1), name="Posteriors"),
	mxAlgebra(-2*sum(log(PL)), name="FF"),
	mxFitFunctionAlgebra(algebra="FF"),
	mxComputeEM(
	  estep=mxComputeOnce("Mixture.Posteriors"),
	  mstep=mxComputeGradientDescent(fitfunction="Mixture.fitfunction")))

mmfit <- mxRun(mm)
summary(mmfit)

confusion <- table(mmfit$Posteriors$result < .1, c(rep(TRUE,N/2),rep(FALSE,N/2)))
print(confusion)
omxCheckCloseEnough(sum(diag(confusion)), 200, 15)

omxCheckCloseEnough(coef(mmfit)[c(1,3)], c(mu1,mu2), .4)
omxCheckCloseEnough(coef(mmfit)[c(2,4)], rep(1,2), .3)

omxCheckCloseEnough(mmfit$output$fit, 539.1599, .01)
