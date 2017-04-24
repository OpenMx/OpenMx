require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)  #<--Saves memory
require(mvtnorm)

#Generate data:
set.seed(476)
A <- matrix(0,1000,1000)  #<--Empty GRM
A[lower.tri(A)] <- runif(499500, -0.025, 0.025)
A <- A + t(A)
diag(A) <- runif(1000,0.95,1.05) #<--GRM now complete
y <- t(rmvnorm(1,sigma=A*0.5))  #<--Phenotype 'y' has a "population" variance of 1 and h2 of 0.5 
y <- y + rnorm(1000,sd=sqrt(0.5))
#Merge variables into data matrix:
dat <- t(y)
colnames(dat) <- paste("y",1:1000,sep="") #<--Column names

xpec <- mxExpectationNormal(covariance="V", means="Mu", dimnames=paste("y",1:1000,sep=""))

ff <- mxFitFunctionML()

mxdat <- mxData(observed = dat, type="raw", sort=FALSE)

#We will create some of the necessary objects inside the mxModel() statement.  We mainly want to avoid creating 
#more copies of the GRM than we need to:
testmod <- mxModel(
	"GREML_1GRM_1trait", #<--Model name
	#1x1 matrix containing residual variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "ve", lbound = 0.0001, 
					 name = "Ve"),
	#1x1 matrix containing additive-genetic variance component:
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = var(y)/2, labels = "va", name = "Va"),
	#1000x1000 identity matrix--the "relatedness matrix" for the residuals:
	mxMatrix("Iden",nrow=1000,name="I"),
	#The GRM:
	mxMatrix("Symm",nrow=1000,free=F,values=A,name="A"),
	#The model-expected covariance matrix:
	mxAlgebra((A%x%Va) + (I%x%Ve), name="V"),
	#An MxAlgebra for the heritability:
	mxAlgebra(Va/(Va+Ve), name="h2"),
	mxCI("h2"), #<--Request confidence interval for heritability
	mxdat, #<--MxData object
	xpec,
	ff,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=mean(y),labels="mu",name="meanmtx"),
	mxMatrix(type="Full",nrow=1,ncol=1000,free=F,values=1,name="Unit"),
	mxAlgebra(meanmtx %x% Unit, name="Mu")
)

testrun <- mxRun(testmod)

omxCheckEquals(testrun$output$status$code, 0)
