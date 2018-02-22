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


require(OpenMx)

if (mxOption(NULL, "Default optimizer") == "SLSQP") {
	mxOption(NULL, "Optimality tolerance", "6e-10")
}

set.seed(654684)
data(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("factor")
nVar <- length(manifests)
nFac <- length(latents)

factorModel <- mxModel("One Factor ML",
											 mxData(cov(demoOneFactor), type="cov", means=colMeans(demoOneFactor), numObs=500),
											 mxMatrix("Full", 1, nVar, free=T, values=colMeans(demoOneFactor), name="M"),
											 mxMatrix("Full", nVar, nFac, free=T, values=.2, ubound=1, name="A"),
											 mxMatrix("Diag", nVar, nVar, free=T, values=1, lbound=.0001, ubound=10, name="D"),
											 mxAlgebra(A%*%t(A) + D, name="C"),
											 mxAlgebra(sqrt(diag2vec(C)),name="P"),
											 mxFitFunctionML(),mxExpectationNormal("C", "M", dimnames=manifests),
											 mxCI(c("P"))
)

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=TRUE, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","SE","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=TRUE, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","SE","HQ","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=FALSE, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","SE","HQ","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=FALSE, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","SE","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorModel <- mxOption(factorModel,"Standard Errors","No")

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=FALSE, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(!length(factorFitCI$output$standardErrors))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(!length(factorFitCI$output$standardErrors))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=F, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","SE","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","SE","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorModel <- mxOption(factorModel,"Calculate Hessian","No")

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=FALSE, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(!length(factorFitCI$output$standardErrors))
omxCheckTrue(!length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=F)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(!length(factorFitCI$output$standardErrors))
omxCheckTrue(!length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","SE","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=F, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","SE","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorModel <- mxOption(factorModel,"Standard Errors","Yes")

factorFitCI <- omxCheckWarning(
	mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=FALSE, checkHess=F),
	'the "Standard Errors" option is enabled and the "Calculate Hessian" option is disabled, which may result in poor-accuracy standard errors'
)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","SE","HQ","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
if (mxOption(NULL, "Default optimizer") == 'NPSOL') {
  omxCheckTrue(length(factorFitCI$output$standardErrors))
  omxCheckTrue(length(factorFitCI$output$calculatedHessian))
}
omxCheckTrue(!(factorFitCI$compute$.persist))

omxCheckWarning(
	factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=F),
	'the "Standard Errors" option is enabled and the "Calculate Hessian" option is disabled, which may result in poor-accuracy standard errors'
)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","SE","HQ","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
if (mxOption(NULL, "Default optimizer") == 'NPSOL') {
  omxCheckTrue(length(factorFitCI$output$standardErrors))
}
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=F, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","ND","SE","RD","RE"))
omxCheckTrue(!length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(!length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorFitCI <- mxTryHard(factorModel, fit2beat=-2863, extraTries=3, intervals=T, checkHess=T)
omxCheckEquals(names(factorFitCI$compute$steps),c("GD","CI","ND","SE","RD","RE"))
omxCheckTrue(length(factorFitCI$compute$steps$CI$output$detail))
omxCheckTrue(length(factorFitCI$output$confidenceIntervals))
omxCheckTrue(length(factorFitCI$output$confidenceIntervalCodes))
omxCheckTrue(length(factorFitCI$output$standardErrors))
omxCheckTrue(length(factorFitCI$output$calculatedHessian))
omxCheckTrue(!(factorFitCI$compute$.persist))

factorModel <- mxOption(factorModel,"Calculate Hessian","Yes")

#Test that mxTryHard() will run with other non-default arguments:

#Response variable y:
y <- matrix(c(3.86735442126894,3.21609311807948,1.6681246111281,3.54171497693329,3.02206567904312,2.40194706094571,
							4.00354871075935,3.50175679256405,3.92466558500823,4.26190374144865,2.68136931922448,3.81633303744976,
							7.63868421858379,4.42267196614715,8.92858721732324,5.23749096355528,5.79233361162228,4.68998250109233,
							4.6116111993746,5.56825200215576,3.96879794025251),dimnames=list(NULL,"y"))
#Predictor variable x:
x <- matrix(c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,2),dimnames=list(NULL,"x"))
#Data matrix:
dat <- cbind(x,y)
Laplace_rgsn_mod1 <- mxModel(
	"LaplaceReg",
	mxData(dat, type="raw"),
	mxMatrix( type="Full", nrow=1, ncol=1, 
						free=T, 
						values=5, 
						name="A", labels="a" ), #<--Intercept
	mxMatrix( type="Full", nrow=1, ncol=1, 
						free=TRUE, 
						values=0, 
						name="B", labels="b" ), #<--Slope
	mxMatrix( type="Full", nrow=1, ncol=1, 
						free=TRUE,
						values=5, 
						name="lambda", labels="lambdapar", lbound=0.0001), #<--Residual dispersion parameter
	mxMatrix(type="Full", nrow=1, ncol=1, free=F,
					 name="xdef", labels="data.x"), #<--x is definition variable.
	mxAlgebra(A + B*xdef, name="yhat", dimnames=list(NULL,"y")), #<--yhat
	mxAlgebra((log(2*lambda) + (abs(filteredDataRow - yhat)/lambda) 
	), name="rowAlgebra"), #<--Negative loglikelihood for 1 observation
	mxAlgebra(2*sum(rowResults), name="reduceAlgebra"), #<--Full-sample deviance
	mxFitFunctionRow(rowAlgebra='rowAlgebra',
									 reduceAlgebra='reduceAlgebra',
									 dimnames=c('y'))
)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,greenOK=T)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,loc=1.1)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,scale=10)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,initialGradientStepSize=0.001)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,initialGradientIterations=5)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1,initialTolerance=1)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1, scale=10, finetuneGradient=F)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1, jitterDistrib="runif")
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

Laplace_rgsn_fit1 <- mxTryHard(model=Laplace_rgsn_mod1, jitterDistrib="rcauchy", exhaustive=T, 
															 showInits=T, finetuneGradient=F)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

omxCheckError(mxTryHard(model=Laplace_rgsn_mod1, scale=-0.25, finetuneGradient=F),
							"negative value for argument 'scale'")



