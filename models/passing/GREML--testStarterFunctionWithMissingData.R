require(OpenMx)

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
dat[42,1] <- NA
dat[57,2] <- NA
start <- mxGREMLStarter("GREMLmod",data=dat,Xdata="x",ydata = "y",addOnes = F, dropNAfromV = T)

testmod <- mxModel(
  start,
  mxExpectationGREML(V="V",X="X",y="y"),
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

testrunsumm <- summary(testrun)
omxCheckEquals(testrunsumm$numObs,98)
omxCheckEquals(testrunsumm$estimatedParameters,2)
omxCheckEquals(testrunsumm$observedStatistics,98)
omxCheckEquals(testrunsumm$degreesOfFreedom,96)


plan <- mxComputeSequence(steps=list(
  mxComputeNewtonRaphson(freeSet=c("Ve"),fitfunction="fitfunction"),
  mxComputeOnce('fitfunction', c('fit','gradient','hessian','ihessian'),freeSet=c("Ve")),
  mxComputeStandardError(freeSet=c("Ve")),
  mxComputeReportDeriv(freeSet=c("Ve"))
))

testmod2 <- mxModel(
  start,
  plan,
  mxExpectationGREML(V="V",X="X",y="y",dV = c(ve="I")),
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

testrun2summ <- summary(testrun2)
omxCheckEquals(testrun2summ$numObs,98)
omxCheckEquals(testrun2summ$estimatedParameters,2)
omxCheckEquals(testrun2summ$observedStatistics,98)
omxCheckEquals(testrun2summ$degreesOfFreedom,96)



#Test errors and warnings from mxGREMLStarter():
omxCheckError(mxGREMLStarter(2, data=dat, Xdata="x", ydata = "y"),
              "argument 'model' must be either a character string or an MxModel object")
omxCheckError(mxGREMLStarter("A", data=2, Xdata="x", ydata = "y"),
              "argument 'data' must be either a matrix or dataframe")
omxCheckError(mxGREMLStarter("A", data=matrix(rnorm(200),100,2), Xdata="x", ydata = "y"),
              "data must have column names")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata=2, ydata = "y"),
              "elements of argument 'Xdata' must be of type 'character' (the data column names of the covariates)")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata="x", ydata = 2),
              "argument 'ydata' must be of type 'character' (the data column names of the phenotypes)")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata="x", ydata = "y", Xname=2),
              "argument 'Xname' is not of type 'character' (the name for the matrix of covariates)")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata="x", ydata = "y", Xname="2"),
              paste('The name \'2\' is illegal because it can be interpreted as a number in mxMatrix(type = ',
                    '"Full", nrow = nrow(X), ncol = ncol(X), free = F, values = X, ',
                    'dimnames = list(NULL, colnames(X)), name = Xname, condenseSlots = T)',sep="")
              )
omxCheckError(mxGREMLStarter("A", data=dat, Xdata="x", ydata = "y", yname=2),
              "argument 'yname' is not of type 'character' (the name for the column vector of phenotypes)")
omxCheckError(mxGREMLStarter(testmod2, data=dat, Xdata="x", ydata = "y", Xname = "X"),
              "already an MxMatrix or MxAlgebra named 'X' in model 'GREMLmod'")
omxCheckError(mxGREMLStarter(testmod2, data=dat, Xdata="x", ydata = "y", Xname = "Xerxes", yname="y"),
              "already an MxMatrix or MxAlgebra named 'y' in model 'GREMLmod'")
omxCheckWarning(mxGREMLStarter("A", data=dat, Xdata=c("x","x"), ydata = c("y")),
              "resulting 'X' matrix has some redundant column names; is it full rank?")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata=list("x","x"), ydata = c("y")),
  "conflicting number of phenotypes specified by arguments 'Xdata' and 'ydata'")
omxCheckError(mxGREMLStarter("A", data=dat, Xdata=c("x","x"), ydata = c("y","y")),
              "argument 'Xdata' must be provided as a list when argument 'ydata' is of length greater than 1")

start2 <- mxGREMLStarter("GREMLmod",data=dat,Xdata="x",ydata = "y",Xname = "W", yname = "Z",
                         addOnes = F, dropNAfromV = T)
testmod3 <- mxModel(
  start2,
  mxExpectationGREML(V="V",X="W",y="Z"),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V")
)
testmod3$data <- NULL
omxCheckWarning(
  mxGREMLStarter(testmod3,data = dat,Xdata = "x",ydata = "y",Xname = "X",yname = "y",addOnes = F),
  "not adding MxFitFunctionGREML because model 'GREMLmod' already contains a fitfunction")
testmod3@data <-  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("dummyData","dummyData")), type="raw")
testmod3$fitfunction <- NULL
omxCheckWarning(
  mxGREMLStarter(testmod3,data = dat,Xdata = "x",ydata = "y",Xname = "X",yname = "y",addOnes = F),
  "not adding dummy MxData object because model 'GREMLmod' already contains an MxData object")
