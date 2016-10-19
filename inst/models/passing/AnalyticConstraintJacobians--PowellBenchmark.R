library(OpenMx)
#mxOption(NULL,"Default optimizer","NPSOL")

powellmod1 <- mxModel(
	"PowellBenchmarkNoJacobians",
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),labels=paste("x",1:5,sep=""),name="X"),
	mxAlgebra( exp(prod(X)), name="powellfunc"),
	mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]), 
						 name="objgrad", 
						 dimnames=list(NULL,paste("x",1:5,sep="")) ),
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1"),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2"),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3"),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
)
powellrun1 <- mxRun(powellmod1)
powellrun1$fitfunction$result
powellrun1$output$iterations
powellrun1$output$evaluations

powellmod2 <- mxModel(
	"PowellBenchmarkWithJacobians",
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),labels=paste("x",1:5,sep=""),name="X"),
	mxAlgebra( exp(prod(X)), name="powellfunc"),
	mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5],
									 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]), 
						 name="objgrad", 
						 dimnames=list(NULL,paste("x",1:5,sep="")) ),
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1",jac="jac1",linear=FALSE),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2",jac="jac2",linear=FALSE),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3",jac="jac3",linear=FALSE),
	mxAlgebra(cbind(2*X[1,1],2*X[1,2],2*X[1,3],2*X[1,4],2*X[1,5]),name="jac1",
						dimnames=list(NULL,paste("x",1:5,sep=""))
	),
	mxAlgebra(cbind(0,X[1,3],X[1,2],-5*X[1,5],-5*X[1,4]),name="jac2",
						dimnames=list(NULL,paste("x",1:5,sep=""))),
	mxAlgebra(cbind(3*X[1,1]^2, 3*X[1,2]^2, 0, 0, 0),name="jac3",dimnames=list(NULL,paste("x",1:5,sep=""))),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
)
powellrun2 <- mxRun(powellmod2)
powellrun2$fitfunction$result
powellrun2$output$iterations
powellrun2$output$evaluations

#Right now, only NPSOL knows how to use analytic Jacobians:
if(mxOption(NULL,"Default optimizer")=="NPSOL"){
  #Analytic Jacobians should, if nothing else, cut down on the number of fitfunction evaluations:
	omxCheckEquals(omxGreaterThan(powellrun1$output$evaluations,powellrun2$output$evaluations),1)
}
