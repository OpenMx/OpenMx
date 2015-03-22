require(OpenMx)

omxCheckError(mxExpectationGREML(V=1),
              "argument 'V' is not of type 'character' (the name of the expected covariance matrix)")
omxCheckError(mxExpectationGREML(V = "V", X = 2),
              "argument 'X' is not of type 'character' (the name of the matrix of covariates)")
omxCheckError(mxExpectationGREML(V="V",X="X",y=3),
              "argument 'y' is not of type 'character'")

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
testmod <- mxModel(
  "GREMLtest",
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",X="X",y="y"),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "the GREML expectation function does not have a dataset associated with it in model 'GREMLtest' 
consider setting up your model with mxGREMLStarter()")

testmod <- mxModel(
  testmod,
  mxData(observed = dat, type = "raw"),
  mxMatrix("Full",1,1,F,labels="data.y",name="Z")
)
omxCheckError(mxRun(testmod),
              "definition variables are incompatible (and unnecessary) with GREML expectation")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw",numObs = 0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix(type = "Full", nrow=99, ncol=1, free=F, values=dat[-100,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxMatrix("Full",nrow=100,ncol=99,free=F,values=diag(100)[,-100],name="V",condenseSlots=T),
  mxExpectationGREML(V="V",X="X",y="y"),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "X and y matrices do not have equal numbers of rows")

testmod$X <- mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T)
omxCheckError(mxRun(testmod),
              "V matrix is not square")

testmod$V <- mxMatrix("Full",nrow=99,ncol=99,free=F,values=diag(100)[-100,-100],name="V",condenseSlots=T)
omxCheckError(mxRun(testmod),
              "y and V matrices do not have equal numbers of rows")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw",numObs = 0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",X="X",y="y"),
  mxFitFunctionGREML(casesToDrop = 101L,dropNAfromV = T)
)
omxCheckWarning(
  mxRun(testmod, suppressWarnings=TRUE),
  "casesToDrop vector in GREML fitfunction contains indices greater than the number of observations")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw",numObs = 0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix("Iden",nrow=99,name="J",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",X="X",y="y",dV = c(ve="J")),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "all derivatives of V must have the same dimensions as V")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw",numObs = 0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",X="X",y="y",dV = c(ve="J")),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "The reference 'J' does not exist.  It is used by named reference 'GREMLtest.expectation' .")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw",numObs = 0),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X",condenseSlots=T),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y",condenseSlots=T),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationGREML(V="V",X="X",y="y"),
  mxFitFunctionML()
)
omxCheckError(mxRun(testmod),
              "No covariance expectation in FIML evaluation.")

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(dat[,1],1,100,dimnames=list(NULL,paste("y",1:100,sep=""))), type="raw"),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
  mxMatrix("Zero",1,100,name="Zm"),
  mxAlgebra(I %x% Ve,name="V"),
  mxExpectationNormal(covariance="V",means="Zm", dimnames=paste("y",1:100,sep="")),
  mxFitFunctionGREML()
)
omxCheckError(mxRun(testmod),
              "GREML fitfunction is currently only compatible with GREML expectation")
