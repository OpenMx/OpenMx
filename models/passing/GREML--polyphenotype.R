require(OpenMx)
library(mvtnorm)

set.seed(1234)
Dat <- rmvnorm(n=100,mean=rep(0,5),sigma=matrix(c(
  1,0.1,0.2,0.3,0.4,
  0.1,1,0.3,0.4,0.5,
  0.2,0.3,1,0.4,0.3,
  0.3,0.4,0.4,1,0.2,
  0.4,0.5,0.3,0.2,1),nrow=5,byrow=T))
colnames(Dat) <- c("y1","y2","x1","x2","x3")

start <- mxGREMLStarter("foo", data=Dat, Xdata=list(c("x1","x2"),c("x2","x3")), ydata=c("y1","y2"))
testmod <- mxModel(
  start,
  mxExpectationGREML("V","X","y"),
  #This model is actually misspeficied, since it doesn't include the cross-trait covariance...
  mxMatrix(type="Diag",nrow=200,free=T,values=1.15,labels=c(rep("v1",100),rep("v2",100)),lbound=0.0001,name="V",
           condenseSlots=T)
)
testrun <- mxRun(testmod)
testrunsumm <- summary(testrun)
lm1 <- lm(y1~x1+x2, data=as.data.frame(Dat))
lm2 <- lm(y2~x2+x3, data=as.data.frame(Dat))
omxCheckCloseEnough(testrun$output$estimate, 
                    c(sum(lm1$residuals^2)/lm1$df.residual, sum(lm2$residuals^2)/lm2$df.residual),
                    epsilon=1e-5)
#Check GREML-specific part of summary() output:
omxCheckEquals(testrun$fitfunction$info$b, testrunsumm$GREMLfixeff$coeff)
omxCheckEquals(sqrt(diag(testrun$fitfunction$info$bcov)), testrunsumm$GREMLfixeff$se)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff[1:3], lm1$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se[1:3], summary(lm1)$coeff[,2], epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$coeff[4:6], lm2$coefficients, epsilon=1e-5)
omxCheckCloseEnough(testrunsumm$GREMLfixeff$se[4:6], summary(lm2)$coeff[,2], epsilon=1e-5)
omxCheckEquals(testrunsumm$GREMLfixeff$name,c("y1_1","y1_x1","y1_x2","y2_1","y2_x2","y2_x3"))
