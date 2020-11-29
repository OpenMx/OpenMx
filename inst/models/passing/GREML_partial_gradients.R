require(OpenMx)
options(mxCondenseMatrixSlots=TRUE)
# mxOption(NULL,"Default optimizer","NPSOL")
# mxOption(NULL,"Number of threads",2)
# mxOption(NULL,"Print level",20)
# mxOption(NULL,"Print file",3)
# mxOption(NULL,"Verify level",3)
require(mvtnorm)


#Generate data:
set.seed(476)
A1 <- matrix(0,100,100)
A1[lower.tri(A1)] <- runif(4950, -0.025, 0.025)
A1 <- A1 + t(A1)
diag(A1) <- runif(100,0.95,1.05)
A2 <- matrix(0,100,100)
A2[lower.tri(A2)] <- runif(4950, -0.025, 0.025)
A2 <- A2 + t(A2)
diag(A2) <- runif(100,0.95,1.05)
y <- t(rmvnorm(1,sigma=A1*0.25)+rmvnorm(1,sigma=A2*0.25))
y <- y + rnorm(100,sd=sqrt(0.5))
y[100] <- NA
x <- rnorm(100)
dat <- cbind(y,x)
colnames(dat) <- c("y","x")

plan <- mxComputeSequence(
	steps=list(
		mxComputeGradientDescent(engine=mxOption(NULL,"Default optimizer")),
		mxComputeOnce('fitfunction', c('gradient','hessian')),
		mxComputeStandardError(),
		mxComputeHessianQuality(),
		mxComputeReportDeriv(),
		mxComputeReportExpectation()
	))


test0 <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2",ve="I"))
)
test0 <- mxRun(test0)
summary(test0)
omxCheckCloseEnough(test0$output$fit, 280.3646874, 1e-4)
omxCheckCloseEnough(coef(test0), c(0.5668669,-0.8059083,1.1844126), 1e-4)
omxCheckCloseEnough(test0$output$standardErrors[,1], c(1.2176641,0.8521144,0.8803859), 1e-4)
omxCheckCloseEnough(test0$output$gradient, c(-2.539113e-07, -3.300839e-07, -3.607828e-07), 1e-4)
omxCheckCloseEnough(test0$expectation$b[,1], c(-0.003919201,0.043639875), 1e-4)
omxCheckCloseEnough(vech(test0$expectation$bcov), c(0.0092116302,0.0003739306,0.0087043503), 1e-4)


test1sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2"))
)
test1sa <- mxRun(test1sa)
summary(test1sa)
omxCheckCloseEnough(test0$output$fit, test1sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test1sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test1sa$output$standardErrors,1e-7)
omxCheckCloseEnough(test0$output$gradient,test1sa$output$gradient,2e-6)
omxCheckCloseEnough(test0$expectation$b,test1sa$expectation$b,1e-8)
omxCheckCloseEnough(test0$expectation$bcov,test1sa$expectation$bcov,1e-9)


test1num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",va2="A2"),autoDerivType="numeric")
)
test1num <- mxRun(test1num)
summary(test1num)
omxCheckCloseEnough(test0$output$fit, test1num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test1num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test1num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test1num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test1num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test1num$expectation$bcov,1e-7)


test2sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va2="A2",ve="I"))
)
test2sa <- mxRun(test2sa)
summary(test2sa)
omxCheckCloseEnough(test0$output$fit, test2sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test2sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test2sa$output$standardErrors,1e-8)
omxCheckCloseEnough(test0$output$gradient,test2sa$output$gradient,1e-6)
omxCheckCloseEnough(test0$expectation$b,test2sa$expectation$b,1e-9)
omxCheckCloseEnough(test0$expectation$bcov,test2sa$expectation$bcov,1e-10)


test2num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va2="A2",ve="I"),autoDerivType="numeric")
)
test2num <- mxRun(test2num)
summary(test2num)
omxCheckCloseEnough(test0$output$fit, test2num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test2num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test2num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test2num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test2num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test2num$expectation$bcov,1e-7)


test3sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",ve="I"))
)
test3sa <- mxRun(test3sa)
summary(test3sa)
omxCheckCloseEnough(test0$output$fit, test3sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test3sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test3sa$output$standardErrors,1e-7)
omxCheckCloseEnough(test0$output$gradient,test3sa$output$gradient,2e-6)
omxCheckCloseEnough(test0$expectation$b,test3sa$expectation$b,1e-8)
omxCheckCloseEnough(test0$expectation$bcov,test3sa$expectation$bcov,1e-10)


test3num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1",ve="I"),autoDerivType="numeric")
)
test3num <- mxRun(test3num)
summary(test3num)
omxCheckCloseEnough(test0$output$fit, test3num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test3num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test3num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test3num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test3num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test3num$expectation$bcov,1e-7)


test4sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(ve="I"))
)
test4sa <- mxRun(test4sa)
summary(test4sa)
omxCheckCloseEnough(test0$output$fit, test4sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test4sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test4sa$output$standardErrors,1e-8)
omxCheckCloseEnough(test0$output$gradient,test4sa$output$gradient,1e-6)
omxCheckCloseEnough(test0$expectation$b,test4sa$expectation$b,1e-9)
omxCheckCloseEnough(test0$expectation$bcov,test4sa$expectation$bcov,1e-10)


test4num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(ve="I"),autoDerivType="numeric")
)
test4num <- mxRun(test4num)
summary(test4num)
omxCheckCloseEnough(test0$output$fit, test4num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test4num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test4num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test4num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test4num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test4num$expectation$bcov,1e-7)


test5sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1"))
)
test5sa <- mxRun(test5sa)
summary(test5sa)
omxCheckCloseEnough(test0$output$fit, test5sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test5sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test5sa$output$standardErrors,1e-8)
omxCheckCloseEnough(test0$output$gradient,test5sa$output$gradient,1e-6)
omxCheckCloseEnough(test0$expectation$b,test5sa$expectation$b,1e-9)
omxCheckCloseEnough(test0$expectation$bcov,test5sa$expectation$bcov,1e-10)


test5num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va1="A1"),autoDerivType="numeric")
)
test5num <- mxRun(test5num)
summary(test5num)
omxCheckCloseEnough(test0$output$fit, test5num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test5num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test5num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test5num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test5num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test5num$expectation$bcov,1e-7)


test6sa <- mxModel(
	"GREMLtest",
	plan,
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va2="A2"))
)
test6sa <- mxRun(test6sa)
summary(test6sa)
omxCheckCloseEnough(test0$output$fit, test6sa$output$fit, 1e-12)
omxCheckCloseEnough(coef(test0),coef(test6sa),1e-7)
omxCheckCloseEnough(test0$output$standardErrors,test6sa$output$standardErrors,1e-8)
omxCheckCloseEnough(test0$output$gradient,test6sa$output$gradient,1e-6)
omxCheckCloseEnough(test0$expectation$b,test6sa$expectation$b,1e-9)
omxCheckCloseEnough(test0$expectation$bcov,test6sa$expectation$bcov,1e-10)


test6num <- mxModel(
	"GREMLtest",
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values =0.5, labels = "ve", lbound = 0.0001,
					 name = "Ve"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.25, labels = "va1", name = "Va1"),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 0.20, labels = "va2", name = "Va2"),
	mxData(observed = dat, type="raw", sort=FALSE),
	mxExpectationGREML(V="V",yvars="y", Xvars="x", addOnes=T),
	mxMatrix("Iden",nrow=100,name="I"),
	mxMatrix("Symm",nrow=100,free=F,values=A1,name="A1"),
	mxMatrix("Symm",nrow=100,free=F,values=A2,name="A2"),
	mxAlgebra((A1%x%Va1) + (A2%x%Va2) + (I%x%Ve), name="V"),
	mxFitFunctionGREML(dV=c(va2="A2"),autoDerivType="numeric")
)
test6num <- mxRun(test6num)
summary(test6num)
omxCheckCloseEnough(test0$output$fit, test6num$output$fit, 5e-9)
omxCheckCloseEnough(coef(test0),coef(test6num),1e-5)
omxCheckCloseEnough(test0$output$standardErrors,test6num$output$standardErrors,0.05)
omxCheckCloseEnough(test0$output$gradient,test6num$output$gradient,5e-4)
omxCheckCloseEnough(test0$expectation$b,test6num$expectation$b,5e-6)
omxCheckCloseEnough(test0$expectation$bcov,test6num$expectation$bcov,1e-7)
