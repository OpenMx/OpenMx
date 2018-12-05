require(OpenMx)
#Load and inspect data:
data(LongitudinalOverdispersedCounts)
str(longData)
str(wideData)

olsmod <- lm(y~tiem+x1+x2+x3, data=longData)
summary(olsmod)

#Let's try a generalized linear model (GLM).  We'll use the quasi-Poisson quasilikelihood function to see
#how well the y variable is approximated by a Poisson distribution (conditional on time and covariates):
glm.mod <- glm(y~tiem+x1+x2+x3, data=longData, family="quasipoisson")
#The estimate of the dispersion parameter should be about 1.0 if the data are conditionally Poisson.  
#We can see that it is actually greater than 2, indicating overdispersion:
summary(glm.mod) #<--We can get good start values for the regression coefficients from here.

#Start MxModel: ###

#MxMatrix to contain time-invariant covariates, as defintion variables:
Xdef <- mxMatrix(type="Full",nrow=1, ncol=3, free=F, labels=c("data.x1","data.x2","data.x3"), name="Xdef")

#MxMatrix to contain times (assessment waves), as definition variables:
timeVec <- mxMatrix(type="Full",nrow=4,ncol=1,free=F,labels=c("data.t0","data.t1","data.t2","data.t3"),
										name="timeVec")

#Note that the time variable is NA for participants who are missing their y score for that assessment.
#OpenMx will work around the NA's on the endogenous variable, y, but it does not tolerate NA's among the 
#definition variables.  We need to set the NA's on the time variable to a "pseudo-missing" value--something
#extreme that will throw off the results to alert us that the invalid definition-variable score was used.
wideData[,"t1"][is.na(wideData[,"y1"])] <- -999
wideData[,"t2"][is.na(wideData[,"y2"])] <- -999
wideData[,"t3"][is.na(wideData[,"y3"])] <- -999

#We're now ready to make our mydata object:
mydat <- mxData(observed = wideData, type="raw")

#a unit MxMatrix:
U <- mxMatrix(type="Unit",nrow=4,ncol=1,name="U")

#an MxMatrix of regression coefficients:
Beta <- mxMatrix(type="Full",nrow=5,ncol=1,free=T,
								 values=c(0.5,0.3,0.2,-0.1,0.8), #<--Start values from the GLM above
								 labels=c("intrcpt","betax1","betax2","betax3","betatime"),
								 name="Beta")
if(mxOption(NULL,"Default optimizer")=="SLSQP"){
	Beta$lbound[1:5,1] <- -20
	Beta$ubound[1:5,1] <- 20
}

#This MxMatrix contains the dispersion parameter.  In our case, the conditional residual variance will be
#modeled as the dispersion parameter times the conditional mean of the response variable.
Disper <- mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,lbound=0.0001,labels="dispersion",name="Disper")

#This MxMatrix will serve as the residual correlation matrix.  We will use an unstructured correlation matrix,
#with reasonable box constraints on the parameters:
R <- mxMatrix(type="Stand",nrow=4,ncol=4,free=T,values=0.1,#labels="r",
							labels=c("r12","r13","r14","r23","r24","r34"),
							lbound=0,ubound=0.9999,name="R")

#For each participant, this MxAlgebra puts the regression constant, the time-invariant covariates, and the time 
#variable together in a "long-format" matrix:
rowX <- mxAlgebra(cbind(U,(U%*%Xdef),timeVec), name="rowX")

#This MxAlgebra yields the conditional mean response, given the covariates and wave of assessment.  
#Notice that we are using a loglinear link: the conditional mean equals an exponentiated linear composite of 
#regressors, or equivalently, the log of the conditional mean equals a linear composite of predictors:
yhatAlg <- mxAlgebra( exp(t(rowX%*%Beta)), name="yhat")

#Vector of conditional standard deviations:
S <- mxAlgebra(vec2diag(sqrt(yhat%x%Disper)), name="S")

#Conditional covariance matrix:
V <- mxAlgebra(S %*% R %*% S, name="V")

expec <- mxExpectationNormal(covariance="V", means="yhat", dimnames=c("y0","y1","y2","y3"))
fitfunc <- mxFitFunctionML()

#Put everything together into our MxModel:
mygeemodel <- mxModel(
	"GEE",
	Beta,Disper,expec,fitfunc,mydat,R,rowX,S,timeVec,U,V,Xdef,yhatAlg
)

mygeemodel <- mxOption(mygeemodel, 'Optimality tolerance', 1.)

#The normal-theory standard errors are going to be wrong, but that doesn't matter for the purposes of
#debugging mxSE():
mygeerun <- mxRun(mygeemodel)
summary(mygeerun)

#All three of these algebras depend upon definition variables:
omxCheckTrue(!any(is.na(mxSE(yhat,mygeerun,defvar.row=1))))
