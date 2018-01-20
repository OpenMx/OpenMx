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

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=FALSE)
gff <- mxFitFunctionGREML(dV=c(ve="I"))
plan <- mxComputeSequence(steps=list(
  mxComputeNewtonRaphson(fitfunction="fitfunction"),
  mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
  mxComputeStandardError(),
			      mxComputeReportDeriv(),
			      mxComputeReportExpectation()
))

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = dat, type="raw", sort=FALSE),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  ge,
  gff,
  plan
)

testrun <- mxRun(testmod)

omxCheckCloseEnough(testrun$output$estimate[1],var(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$expectation$b,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$expectation$bcov,var(dat[,1])/100,epsilon=10^-5)
omxCheckCloseEnough(testrun$output$standardErrors[1],sqrt((2*testrun$output$estimate^2)/100),epsilon=10^-3)

testrunsumm <- summary(testrun)
omxCheckEquals(testrunsumm$numObs,1)
omxCheckEquals(testrunsumm$estimatedParameters,2)
omxCheckEquals(testrunsumm$observedStatistics,100)
omxCheckEquals(testrunsumm$degreesOfFreedom,98)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrunsumm$GREMLfixeff$name,"x")
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se,sqrt(var(dat[,1])/100),epsilon=10^-5)


gff2 <- mxFitFunctionGREML()

testmod2 <- mxModel(
  "GREMLtest",
  mxData(observed = dat, type="raw", sort=FALSE),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  ge,
  gff2
)

testrun2 <- mxRun(testmod2)

omxCheckCloseEnough(testrun2$output$estimate[1],var(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$expectation$b,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$expectation$bcov,var(dat[,1])/100,epsilon=10^-5)

testrun2summ <- summary(testrun2)
omxCheckEquals(testrun2summ$numObs,1)
omxCheckEquals(testrun2summ$estimatedParameters,2)
omxCheckEquals(testrun2summ$observedStatistics,100)
omxCheckEquals(testrun2summ$degreesOfFreedom,98)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrun2summ$GREMLfixeff$name,"x")
omxCheckCloseEnough(testrun2summ$GREMLfixeff$coeff,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$se,sqrt(var(dat[,1])/100),epsilon=10^-5)


#Test GREML expectation used with FIML fitfunction:
testmod3 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge,
	mxFitFunctionML()
)

testrun3 <- mxRun(testmod3)

omxCheckCloseEnough(testrun3$output$estimate[1],var(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun3$expectation$b,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun3$expectation$bcov,(var(dat[,1]))/100,epsilon=10^-5)
#MLfit from testrun should match fit of testrun3:
omxCheckCloseEnough(testrun3$output$minimum, testrun$output$fit, epsilon=10^-2)

testrun3summ <- summary(testrun3)
#FIML fitfunction doesn't know how to tell the frontend that a GREML analysis involves only 1 observation:
omxCheckError( 
	omxCheckEquals(testrun3summ$numObs,1),
	"'100' and '1' are not equal")
omxCheckEquals(testrun3summ$estimatedParameters,2)
omxCheckEquals(testrun3summ$observedStatistics,100)
omxCheckEquals(testrun3summ$degreesOfFreedom,98)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrun3summ$GREMLfixeff$name,"x")
omxCheckCloseEnough(testrun3summ$GREMLfixeff$coeff,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun3summ$GREMLfixeff$se,sqrt(var(dat[,1])/100),epsilon=10^-5)
