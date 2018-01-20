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
dat[42,1] <- NA
dat[57,2] <- NA
start <- mxGREMLDataHandler(data=dat,Xvars="x",yvars="y",addOnes = F)
omxCheckEquals(colnames(start$yX)[1],"y")
omxCheckEquals(colnames(start$yX)[2],"x")

testmod <- mxModel(
  mxData(observed = start$yX, type = "raw"),
  mxExpectationGREML(V="V",dataset.is.yX=TRUE,casesToDropFromV=start$casesToDrop),
  mxFitFunctionGREML(),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)

testrun <- mxRun(testmod)

omxCheckCloseEnough(testrun$output$estimate[1],
                    var(dat[-testrun$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$expectation$b,
                    mean(dat[-testrun$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$expectation$bcov,
                    var(dat[-testrun$expectation$casesToDrop,1])/98,epsilon=10^-5)
omxCheckEquals(dim(testrun$V$result),c(100,100))

testrunsumm <- summary(testrun)
omxCheckEquals(testrunsumm$numObs,1)
omxCheckEquals(testrunsumm$estimatedParameters,2)
omxCheckEquals(testrunsumm$observedStatistics,98)
omxCheckEquals(testrunsumm$degreesOfFreedom,96)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrunsumm$GREMLfixeff$name,"x")
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff,
                    mean(dat[-testrun$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se,
                    sqrt(var(dat[-testrun$expectation$casesToDrop,1])/98),epsilon=10^-5)



plan <- mxComputeSequence(steps=list(
  mxComputeNewtonRaphson(fitfunction="fitfunction"),
  mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian')),
  mxComputeStandardError(),
  mxComputeReportDeriv(),
  mxComputeReportExpectation()
))

testmod2 <- mxModel(
  plan,
  mxData(observed=start$yX, type="raw"),
  mxExpectationGREML(V="V",dataset.is.yX=TRUE,casesToDropFromV=start$casesToDrop),
  mxFitFunctionGREML(dV = c(ve="I")),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)

testrun2 <- mxRun(testmod2)

omxCheckCloseEnough(testrun2$output$estimate[1],
                    var(dat[-testrun2$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$expectation$b,
                    mean(dat[-testrun2$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$expectation$bcov,
                    var(dat[-testrun2$expectation$casesToDrop,1])/98,epsilon=10^-5)
omxCheckCloseEnough(testrun2$output$standardErrors[1],sqrt((2*testrun2$output$estimate^2)/98),epsilon=10^-3)
omxCheckEquals(dim(testrun2$V$result),c(100,100))

testrun2summ <- summary(testrun2)
omxCheckEquals(testrun2summ$numObs,1)
omxCheckEquals(testrun2summ$estimatedParameters,2)
omxCheckEquals(testrun2summ$observedStatistics,98)
omxCheckEquals(testrun2summ$degreesOfFreedom,96)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrun2summ$GREMLfixeff$name,"x")
omxCheckCloseEnough(testrun2summ$GREMLfixeff$coeff,
                    mean(dat[-testrun2$expectation$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$se,
                    sqrt(var(dat[-testrun2$expectation$casesToDrop,1])/98),epsilon=10^-5)



#Test errors and warnings from mxGREMLDataHandler():
omxCheckError(mxGREMLDataHandler(data=2, Xvars="x",yvars="y"),
              "argument 'data' must be either a matrix or dataframe")
omxCheckError(mxGREMLDataHandler(data=matrix(rnorm(200),100,2), Xvars="x", yvars = "y"),
              "data must have column names")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=2, yvars = "y"),
              "elements of argument 'Xvars' must be of type 'character' (the data column names of the covariates)")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars="x", yvars=2),
              "argument 'yvars' is not of type 'character' (the data column names of the phenotypes)")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars="x",yvars=character(0)),
              "you must specify at least one phenotype in argument 'yvars'")
omxCheckWarning(mxGREMLDataHandler(data=dat, Xvars=c("x","x"), yvars=c("y")),
              "resulting 'X' matrix has some redundant column names; is it full rank?")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=list("x","x"), yvars = c("y")),
  "conflicting number of phenotypes specified by arguments 'Xvars' and 'yvars'")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=c("x","x"), yvars=c("y","y")),
              "argument 'Xvars' must be provided as a list when argument 'yvars' is of length greater than 1")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=list("x","x"), yvars=c("y","z")),
              "'z' in argument 'yvars' is not among the data column names")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=list("x","z"), yvars = c("y","y")),
              "'z' in argument 'Xvars' is not among the data column names")
omxCheckError(mxGREMLDataHandler(data=dat, Xvars=list(c("x","x"),"x"), yvars = c("y","y"), staggerZeroes = F),
              "all phenotypes must have the same number of covariates when staggerZeroes=FALSE")

start2 <- mxGREMLDataHandler(data=dat, Xvars="x", yvars="y", addOnes=F)
testmod3 <- mxModel(
  mxData(observed=start2$yX, type="raw"),
  mxExpectationGREML(V="V",dataset.is.yX=T,casesToDropFromV=start2$casesToDrop),
  mxFitFunctionGREML(),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)
omxCheckEquals(testmod3$expectation$casesToDrop,c(42,57))

testmod4 <- mxModel(
  mxData(observed=start2$yX, type="raw"),
  mxExpectationGREML(V="V",dataset.is.yX=T),
  mxFitFunctionGREML(),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=98,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)
testrun4 <- mxRun(testmod4)
omxCheckCloseEnough(testrun4$output$estimate[1],
                    var(dat[-c(42,57),1]),epsilon=10^-5)
omxCheckCloseEnough(testrun4$expectation$b,
                    mean(dat[-c(42,57),1]),epsilon=10^-5)
omxCheckCloseEnough(testrun4$expectation$bcov,
                    var(dat[-c(42,57),1])/98,epsilon=10^-5)
omxCheckEquals(dim(testrun4$V$result),c(98,98))
