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
require(mvtnorm)

set.seed(1234)
Dat <- rmvnorm(n=100,mean=rep(0,5),sigma=matrix(c(
  1,0.1,0.2,0.3,0.4,
  0.1,1,0.3,0.4,0.5,
  0.2,0.3,1,0.4,0.3,
  0.3,0.4,0.4,1,0.2,
  0.4,0.5,0.3,0.2,1),nrow=5,byrow=T))
colnames(Dat) <- c("y1","y2","x1","x2","x3")

start <- mxGREMLDataHandler(data=Dat, Xvars=list(c("x1","x2"),c("x2","x3")), yvars=c("y1","y2"))
testmod <- mxModel(
  mxData(observed=start$yX, type="raw", sort = FALSE), #<--sort=FALSE is CRITICALLY important!
  mxExpectationGREML(V = "V",dataset.is.yX = T),
  mxFitFunctionGREML(),
  #This model is actually misspecified, since it doesn't include the cross-trait covariance...
  mxMatrix(type="Diag",nrow=200,free=T,values=1.15,labels=c(rep("v1",100),rep("v2",100)),lbound=0.0001,name="V",
           condenseSlots=T)
)
testrun <- mxRun(testmod)
#testrun <- mxRun(testrun)
testrunsumm <- summary(testrun)
lm1 <- lm(y1~x1+x2, data=as.data.frame(Dat))
lm2 <- lm(y2~x2+x3, data=as.data.frame(Dat))
omxCheckCloseEnough(testrun$output$estimate, 
                    c(sum(lm1$residuals^2)/lm1$df.residual, sum(lm2$residuals^2)/lm2$df.residual),
                    epsilon=1e-5)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrun$expectation$b, testrunsumm$GREMLfixeff$coeff)
omxCheckEquals(sqrt(diag(testrun$expectation$bcov)), testrunsumm$GREMLfixeff$se)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff[1:3], lm1$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se[1:3], summary(lm1)$coeff[,2], epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff[4:6], lm2$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se[4:6], summary(lm2)$coeff[,2], epsilon=1e-5)
omxCheckEquals(testrunsumm$GREMLfixeff$name,c("y1_1","y1_x1","y1_x2","y2_1","y2_x2","y2_x3"))

#Use the non-default blockByPheno=F:
start2 <- mxGREMLDataHandler(data=Dat, Xvars=list(c("x1","x2"),c("x2","x3")), yvars=c("y1","y2"), 
                         blockByPheno = F)
testmod2 <- mxModel(
  mxData(observed = start2$yX, type="raw", sort=FALSE),
  mxExpectationGREML("V",dataset.is.yX = T),
  mxFitFunctionGREML(),
  #This model is actually misspecified, since it doesn't include the cross-trait covariance...
  mxMatrix(type="Diag",nrow=200,free=T,values=1.15,labels=rep(c("v1","v2"),100),lbound=0.0001,name="V",
           condenseSlots=T)
)
testrun2 <- mxRun(testmod2)
testrun2summ <- summary(testrun2)
omxCheckCloseEnough(testrun2$output$estimate, 
                    c(sum(lm1$residuals^2)/lm1$df.residual, sum(lm2$residuals^2)/lm2$df.residual),
                    epsilon=1e-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$coeff[1:3], lm1$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$se[1:3], summary(lm1)$coeff[,2], epsilon=1e-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$coeff[4:6], lm2$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrun2summ$GREMLfixeff$se[4:6], summary(lm2)$coeff[,2], epsilon=1e-5)



#Check output of starter function for non-default data-structuring of y and X:
start <- mxGREMLDataHandler(data=Dat, Xvars=list(c("x1","x2"),c("x2","x3")), yvars=c("y1","y2"), 
                        blockByPheno=F)
omxCheckTrue(all(start$yX[1:4,1] == matrix(c(Dat[1,1],Dat[1,2],Dat[2,1],Dat[2,2]))))
omxCheckTrue(all(start$yX[1,-1] == c(1, Dat[1,"x1"], Dat[1,"x2"], 0, 0, 0)))
omxCheckTrue(all(start$yX[2,-1] == c(0, 0, 0, 1, Dat[1,"x2"], Dat[1,"x3"])))
omxCheckTrue(all(start$yX[3,-1] == c(1, Dat[2,"x1"], Dat[2,"x2"], 0, 0, 0)))
omxCheckTrue(all(start$yX[4,-1] == c(0, 0, 0, 1, Dat[2,"x2"], Dat[2,"x3"])))

start <- mxGREMLDataHandler(data=Dat, Xvars=list(c("x1","x2"),c("x2","x3")), yvars=c("y1","y2"), 
                        staggerZeroes=F)
omxCheckTrue(all(start$yX[1,-1]==c(1,Dat[1,"x1"],Dat[1,"x2"])))
omxCheckTrue(all(start$yX[101,-1]==c(1,Dat[1,"x2"],Dat[1,"x3"])))

start <- mxGREMLDataHandler(data=Dat, Xvars=list(c("x1","x2"),c("x2","x3")), yvars=c("y1","y2"), 
                        staggerZeroes=F, blockByPheno=F)
omxCheckTrue(all(start$yX[1:4,1] == matrix(c(Dat[1,1],Dat[1,2],Dat[2,1],Dat[2,2]))))
omxCheckTrue(all(start$yX[1,-1]==c(1,Dat[1,"x1"],Dat[1,"x2"])))
omxCheckTrue(all(start$yX[2,-1]==c(1,Dat[1,"x2"],Dat[1,"x3"])))

