require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  
require(mvtnorm)


#Generate data:
set.seed(476)
A1 <- matrix(0,100,100)  
A1[lower.tri(A1)] <- runif(4950, -0.025, 0.025)
A1 <- A1 + t(A1)
diag(A1) <- runif(100,0.95,1.05)
A2 <- matrix(0,100,100)  
A2[lower.tri(A2)] <- runif(4950, -0.025, 0.025)
A2 <- A2 + t(A2)
diag(A2) <- runif(100,0.95,1.05)
y <- t(rmvnorm(1,sigma=A1*0.25)+rmvnorm(1,sigma=A2*0.25))  
y <- y + rnorm(100,sd=sqrt(0.5))
#y[100] <- NA
x <- rnorm(100) 
dat <- cbind(y,x)
colnames(dat) <- c("y","x")
dat[100,1] <- NA #<--Note that the last row of the dataset is being set to NA.

testmod <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	plan <- mxComputeSequence(freeSet = c("Va1","Va2","Ve"),steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
)
testrun <- mxRun(testmod) #<--Should run without error.

dat[,1] <- y
dat[90:99,1] <- NA
testmod2 <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	plan <- mxComputeSequence(freeSet = c("Va1","Va2","Ve"),steps=list(
		mxComputeNewtonRaphson(fitfunction="fitfunction"),
		mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
		mxComputeStandardError(),
		mxComputeReportDeriv()
	)),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uni"),
	#Note that the derivatives of V are all MxAlgebras:
	mxAlgebra(vec2diag(Uni),name="I"),
	mxAlgebra(A1%x%1,name="dV_dva1"),
	mxAlgebra(A2%x%1,name="dV_dva2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="dV_dva1",va2="dV_dva2",ve="I"))
)
testrun2 <- mxRun(testmod2) #<--Should run without error.






