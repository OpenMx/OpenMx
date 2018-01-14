#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# This script (by Rob K.) demonstrates how to conduct a generalized-estimating-equations analysis in OpenMx.
# The response variables are overdispersed counts, assessed longitudinally over four waves of assessment.  The time
# variable is wave of assessment (0, 1, 2, and 3), and there are three time-invariant covariates.  Robust standard
# errors for the regression coefficients are obtained via the "sandwich estimator."

# GEE is a "semi-parametric" regression method.  A GEE model consists of three things: (1) a link function,
# which relates the conditional mean of the response variable to a linear composite of the regressors;
# (2) a variance function, which relates the residual variance to the conditional mean; and (3) a correlation
# structure that approximates the residual dependence among clustered/repeated observations.  Usually, GEE is
# intended for inference about the fixed effects (i.e., the regression coefficients) when the random effects 
# (variances and correlations) are not of interest.  GEE does not rely upon any particular distributional 
# assumptions

require(OpenMx)
#Load and inspect data:
data(LongitudinalOverdispersedCounts)
str(longData)
str(wideData)

#Let's try ordinary least-squares (OLS) regression:
olsmod <- lm(y~tiem+x1+x2+x3, data=longData)
#We will see in the diagnostic plots that the residuals are poorly approximated by normality, and are 
#heteroskedastic.  We also know that the residuals are not independent of one another, because we have 
#repeated-measures data:
plot(olsmod)
#In the summary, it looks like all of the regression coefficients are significantly different from zero, but we
#know that because the assumptions of OLS regression are violated that we should not trust its results:
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

#We will use the multivariate-normal fitfunction and expectation out of mere convenience, and not because we
#actually think multivariate normality reasonably describes the data:
expec <- mxExpectationNormal(covariance="V", means="yhat", dimnames=c("y0","y1","y2","y3"))
fitfunc <- mxFitFunctionML()

#Put everything together into our MxModel:
mygeemodel <- mxModel(
	"GEE",
	Beta,Disper,expec,fitfunc,mydat,R,rowX,S,timeVec,U,V,Xdef,yhatAlg
)
#We do not actually want OpenMx's Hessian and standard errors, because they are based on a likelihood
#(multivariate normal) that we do not actually think is realistic.  Instead, we will later obtain standard
#errors from the "sandwich estimator":
mygeemodel <- mxOption(mygeemodel,"Calculate Hessian","No")
mygeemodel <- mxOption(mygeemodel,"Standard Errors","No")

#Run our model:
mygeerun <- mxRun(mygeemodel)
#Model summary:
summary(mygeerun)

#OpenMx developer tests--compare results to those from package 'gee': ###
#Regression coefficients:
omxCheckCloseEnough(
	mygeerun$output$estimate[c("intrcpt","betax1","betax2","betax3","betatime")],
	c(0.5084808,0.3078977,0.1914246,-0.1113108,0.8389867),
	epsilon=0.05)
#Dispersion parameter:
omxCheckCloseEnough(mygeerun$output$estimate["dispersion"],2.178853,epsilon=0.1)
#Correlations:
omxCheckCloseEnough(
	mygeerun$output$estimate[c("r12","r13","r14","r23","r24","r34")],
	c(0.3226068,0.2133836,0.1502026,0.1834196,0.1224658,0.05675176),
	epsilon=0.05)

#Now, we obtain the sandwich estimator of the repeated-sampling covariance matrix of the regression coefficients:
Bread <- matrix(0,5,5)#<--The "bread" of the sandwich.
Meat <- matrix(0,5,5)#<--The "meat" of the sandwich
#We loop through the 1000 participants, and calculate the cumulative sums of the contributions of each to the 
#Bread and Meat:
for(i in 1:nrow(wideData)){
	ycurr <- mydat$observed[i,5:8] #<--response observations (y variables) for row i.
	pres <- which(!is.na(ycurr))#<--A "presence vector," indicating which y's are non-missing in row i.
	ycurr <- matrix(ycurr[pres],nrow=1) #<--Filter out any NAs among the y's.
	yhatcurr <- matrix(mxEval(yhat,mygeerun,T,defvar.row=i)[,pres],nrow=1)#<--Contional mean for row i.
	rowXcurr <- mxEval(rowX,mygeerun,T,defvar.row=i)[pres,]#<--value of MxAlgebra 'rowX' for row i.
	#D is the derivative of the conditional mean 'yhat' with respect to the regression coefficients.  It has
	#the form it has here because we are using a loglinear link function.  If we were using a different link
	#function, it would have a different form:
	D <- vec2diag(yhatcurr) %*% rowXcurr
	Vinvcurr <- solve(mxEval(V,mygeerun,T,defvar.row=i)[pres,pres]) #<--Inverse covariance matrix for row i.
	Bread <- Bread + t(D)%*%Vinvcurr%*%D
	Meat <- Meat + t(D)%*%Vinvcurr%*%t(ycurr-yhatcurr)%*%(ycurr-yhatcurr)%*%Vinvcurr%*%D
}
#Finally, the sandwich:
( sammich <- solve(Bread)%*%Meat%*%solve(Bread) )
dimnames(sammich) <- list(c("intrcpt","betax1","betax2","betax3","betatime"),
													c("intrcpt","betax1","betax2","betax3","betatime"))
#Regression coefficients with standard errors:
cbind(b=mygeerun$output$estimate[c("intrcpt","betax1","betax2","betax3","betatime")],SE=sqrt(diag(sammich)))

#Again, we do not actually think the multivariate-normal likelihood is actually reasonable for the data, so we
#do not use it for a likelihood-ratio test.  Intead, we'll do a Wald chi-square test, for the significance 
#of the non-intercept regression coefficients:
nonIntrcptLabels <- c("betax1","betax2","betax3","betatime")
Waldchi2 <- mygeerun$output$estimate[nonIntrcptLabels] %*% solve(sammich[nonIntrcptLabels,nonIntrcptLabels]) %*%
	mygeerun$output$estimate[nonIntrcptLabels]
pchisq(Waldchi2,df=4,lower.tail=F)#<--p-value



#SEs of regression coefficients:
omxCheckCloseEnough(
	sqrt(diag(sammich)),
	c(0.026086291,0.008078779,0.008583959,0.008247869,0.009411659),
	epsilon=0.005)


#Test imxRobustSE:
mygeemodel2 <- mxOption(mygeemodel,"Calculate Hessian","Yes")
mygeemodel2 <- mxOption(mygeemodel2,"Standard Errors","Yes")

mygeerun2 <- mxRun(mygeemodel2)
rse <- imxRobustSE(mygeerun2)
omxCheckCloseEnough(sqrt(diag(sammich)),rse[1:5],epsilon=0.005)
