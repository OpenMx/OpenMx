require(OpenMx)

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
dat[42,1] <- NA
dat[57,2] <- NA
start <- mxGREMLStarter("GREMLmod",data=dat,Xdata="x",ydata = "y",addOnes = F, dropNAfromV = T)

testmod <- mxModel(
  start,
  mxExpectationGREML(V="V",X="X",y="y",fixedEffects = T),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)

testrun <- mxRun(testmod)

omxCheckCloseEnough(testrun$output$estimate[1],
                    var(dat[-testrun$fitfunction$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$fitfunction$info$b,
                    mean(dat[-testrun$fitfunction$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$fitfunction$info$bcov,
                    var(dat[-testrun$fitfunction$casesToDrop,1])/98,epsilon=10^-5)


plan <- mxComputeSequence(steps=list(
  mxComputeNewtonRaphson(freeSet=c("Ve"),fitfunction="fitfunction"),
  mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian'),freeSet=c("Ve")),
  mxComputeStandardError(freeSet=c("Ve")),
  mxComputeReportDeriv(freeSet=c("Ve"))
))

testmod2 <- mxModel(
  start,
  plan,
  mxExpectationGREML(V="V",X="X",y="y",dV = c(ve="I"),fixedEffects = T),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)

testrun2 <- mxRun(testmod2)

omxCheckCloseEnough(testrun2$output$estimate[1],
                    var(dat[-testrun$fitfunction$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$fitfunction$info$b,
                    mean(dat[-testrun$fitfunction$casesToDrop,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun2$fitfunction$info$bcov,
                    var(dat[-testrun$fitfunction$casesToDrop,1])/98,epsilon=10^-5)
omxCheckCloseEnough(testrun2$output$standardErrors[1],sqrt((2*testrun2$output$estimate^2)/98),epsilon=10^-3)
