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
x <- rnorm(1000) #<--Covariate 'x' is actually independent of the phenotype.
#Merge variables into data matrix:
dat <- cbind(y,x)
colnames(dat) <- c("y","x") #<--Column names

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T)

gff <- mxFitFunctionGREML(dV=c(va="A",ve="I"))

plan <- mxComputeSequence(steps=list(
	mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer"),useGradient=TRUE),
	mxComputeReportExpectation()
))

mxdat <- mxData(observed = dat, type="raw", sort=FALSE)

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
	mxdat, #<--MxData object
	ge, #<--GREML expectation
	gff, #<--GREML fitfunction
	plan #<--Custom compute plan
)
testrun <- mxRun(testmod)

testmod2 <- testmod
testmod2$compute <- mxComputeSequence(steps=list(
	mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer"),useGradient=FALSE),
	mxComputeReportExpectation()
))
testrun2 <- mxRun(testmod2)

#One thing that definitely should be smaller with analytic gradient is the number of fitfunction evaluations:
if(testrun$output$evaluations >= testrun2$output$evaluations){stop("models/passing/AnalyticGradientTest.R failed")}
