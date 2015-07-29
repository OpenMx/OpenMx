require(OpenMx)
options(mxCondenseMatrixSlots=FALSE)
require(mvtnorm)


#Generate data:
set.seed(476)
A <- matrix(0,1000,1000)
A[lower.tri(A)] <- runif(499500, -0.025, 0.025)
A <- A + t(A)
diag(A) <- runif(1000,0.95,1.05)
y <- t(rmvnorm(1,sigma=A*0.5))
y <- y + rnorm(1000,sd=sqrt(0.5))
x <- rnorm(1000)

dat <- cbind(y,x)
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T)
gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))
plan <- mxComputeSequence(freeSet = c("Va","Ve"),steps=list(
	mxComputeNewtonRaphson(fitfunction="fitfunction"),
	mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
	mxComputeStandardError(),
	mxComputeReportDeriv()
))

testmod <- mxModel(
	"GREML_1GRM_1trait",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix("Iden",nrow=1000,name="I"),
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	ge,
	gff,
	plan
)

testrun <- mxRun(testmod)
#Test for regressions in how GREML handles its analytic derivatives:
omxCheckTrue(testrun$output$hessian[1,2]>0)
omxCheckCloseEnough(testrun$output$gradient,c(0,0),epsilon=0.15)



#Diagonalize the problem:
eigenA <- eigen(A)
yrot <- t(eigenA$vectors) %*% y
xrot <- t(eigenA$vectors) %*% cbind(1,x)
datrot <- cbind(yrot,xrot)
colnames(datrot) <- c("y","x0","x1")
testmod2 <- mxModel(
	"GREMLtest_1GRM_1trait_diagonalized",
	mxData(observed = datrot, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars=c("x0","x1"), addOnes=F),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	mxMatrix("Iden",nrow=1000,name="I"),
	mxMatrix("Diag",nrow=1000,free=F,values=eigenA$values,name="A"),
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	mxAlgebra(Va/(Va+Ve), name="h2"),
	gff,
	plan
)
testrun2 <- mxRun(testmod2)
omxCheckCloseEnough(testrun2$output$gradient,c(0,0),1e-3)
omxCheckCloseEnough(testrun$output$hessian,testrun2$output$hessian,epsilon=1)
omxCheckCloseEnough(testrun$output$estimate,testrun2$output$estimate,0.001)
omxCheckCloseEnough(testrun$output$standardErrors,testrun2$output$standardErrors,0.0001)
omxCheckCloseEnough(testrun$expectation$b,testrun2$expectation$b,1e-5)
omxCheckCloseEnough(testrun$expectation$bcov,testrun2$expectation$bcov,1e-6)
