require(OpenMx)

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

ge <- mxExpectationGREML(V="V",X="X", y="y")
gff <- mxFitFunctionGREML()

testmod <- mxModel(
  "GREMLtest",
  mxData(observed = matrix(as.double(NA),1,1,dimnames = list("asdf","asdf")), type="raw"),
  mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
  mxMatrix("Iden",nrow=100,name="I"),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,2], name = "X"),
  mxMatrix(type = "Full", nrow=100, ncol=1, free=F, values=dat[,1], name="y"),
  mxAlgebra(I %x% Ve,name="V"),
  ge,
  gff
)

testrun <- mxRun(testmod)

omxCheckCloseEnough(testrun$output$estimate[1],var(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$fitfunction$info$b,mean(dat[,1]),epsilon=10^-5)
omxCheckCloseEnough(testrun$fitfunction$info$bcov,var(dat[,1])/100,epsilon=10^-5)
