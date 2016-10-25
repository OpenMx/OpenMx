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
# powellmod1 <- mxRun(powellmod1,onlyFrontend=T)
# powellmod1$compute$steps$GD$verbose <- 3L
# powellmod1$compute$.persist <- TRUE
powellrun1 <- mxRun(powellmod1)
powellrun1$fitfunction$result
powellrun1$output$iterations
powellrun1$output$evaluations
summary(powellrun1)


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
# powellmod2 <- mxRun(powellmod2,onlyFrontend=T)
# powellmod2$compute$steps$GD$verbose <- 3L
# powellmod2$compute$.persist <- TRUE
powellrun2 <- mxRun(powellmod2)
powellrun2$fitfunction$result
powellrun2$output$iterations
powellrun2$output$evaluations
summary(powellrun2)


#Right now, only NPSOL knows how to use analytic Jacobians:
if(mxOption(NULL,"Default optimizer")=="NPSOL"){
  #Analytic Jacobians should, if nothing else, cut down on the number of fitfunction evaluations:
	omxCheckEquals(omxGreaterThan(powellrun1$output$evaluations,powellrun2$output$evaluations),1)
}


#This third model tests the backend's parameter mapping for analytic Jacobians:
powellmod3 <- mxModel(
	"PowellBenchmarkWithJacobians",
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-2,labels="a",name="A"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="b",name="B"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="c",name="C"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="d",name="D"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="e",name="E"),
	mxAlgebra( exp(A*B*C*D*E), name="powellfunc"),
	mxAlgebra( cbind(powellfunc*B*C*D*E,
									 powellfunc*A*C*D*E,
									 powellfunc*A*B*D*E,
									 powellfunc*A*B*C*E,
									 powellfunc*A*B*C*D),
						 name="objgrad", 
						 dimnames=list(NULL,c(letters[1:5])) ),
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1",linear=FALSE),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2",linear=FALSE),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3",linear=FALSE),
	mxAlgebra(
		cbind(2*E,2*D,2*C,2*B,2*A),
		name="jac1",
		dimnames=list(NULL,c("e","d","c","b","a"))
	),
	mxAlgebra(cbind(C,B,-5*E,-5*D,0),
						name="jac2",
						dimnames=list(NULL,c("b","c","d","e","a")) ),
	mxAlgebra(
		cbind(0,0,0,3*B^2,3*A^2),
		name="jac3",dimnames=list(NULL,c("c","e","d","b","a"))),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
)
# powellmod3 <- mxRun(powellmod3,onlyFrontend=T)
# powellmod3$compute$steps$GD$verbose <- 3L
# powellmod3$compute$.persist <- TRUE
powellrun3 <- mxRun(powellmod3)
omxCheckCloseEnough(powellrun2$fitfunction$result, powellrun3$fitfunction$result)
omxCheckEquals(powellrun2$output$iterations,powellrun3$output$iterations)
omxCheckEquals(powellrun2$output$evaluations,powellrun3$output$evaluations)
omxCheckCloseEnough(coef(powellrun2),coef(powellrun3))
summary(powellrun3)


#This model is like the last one, except that it provides no Jacobian elements for parameter b:
powellmod4 <- mxModel(
	"PowellBenchmarkWithJacobians",
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-2,labels="a",name="A"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="b",name="B"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="c",name="C"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="d",name="D"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="e",name="E"),
	mxAlgebra( exp(A*B*C*D*E), name="powellfunc"),
	mxAlgebra( cbind(powellfunc*B*C*D*E,
									 powellfunc*A*C*D*E,
									 powellfunc*A*B*D*E,
									 powellfunc*A*B*C*E,
									 powellfunc*A*B*C*D),
						 name="objgrad", 
						 dimnames=list(NULL,c(letters[1:5])) ),
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1",linear=FALSE),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2",linear=FALSE),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3",linear=FALSE),
	mxAlgebra(
		cbind(2*E,2*D,2*C,2*A),
		name="jac1",
		dimnames=list(NULL,c("e","d","c","a"))
	),
	mxAlgebra(cbind(B,-5*E,-5*D,0),
						name="jac2",
						dimnames=list(NULL,c("c","d","e","a")) ),
	mxAlgebra(
		cbind(0,0,0,3*A^2),
		name="jac3",dimnames=list(NULL,c("c","e","d","a"))),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
)
# powellmod4 <- mxRun(powellmod4,onlyFrontend=T)
# powellmod4$compute$steps$GD$verbose <- 3L
# powellmod4$compute$.persist <- TRUE
powellrun4 <- mxRun(powellmod4)
omxCheckCloseEnough(powellrun3$fitfunction$result, powellrun4$fitfunction$result,1e-7)
powellrun4$output$iterations
powellrun4$output$evaluations
omxCheckCloseEnough(coef(powellrun3),coef(powellrun4),1e-7)
summary(powellrun4)


#This model is like powellmod3, except "b" is mislabeled "z" in the Jacobian colnames,
#which should cause the backend to behave just as it did for powellmod4:
powellmod5 <- mxModel(
	"PowellBenchmarkWithJacobians",
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-2,labels="a",name="A"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="b",name="B"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=2,labels="c",name="C"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="d",name="D"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=-1,labels="e",name="E"),
	mxAlgebra( exp(A*B*C*D*E), name="powellfunc"),
	mxAlgebra( cbind(powellfunc*B*C*D*E,
									 powellfunc*A*C*D*E,
									 powellfunc*A*B*D*E,
									 powellfunc*A*B*C*E,
									 powellfunc*A*B*C*D),
						 name="objgrad", 
						 dimnames=list(NULL,c(letters[1:5])) ),
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1",linear=FALSE),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2",linear=FALSE),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3",linear=FALSE),
	mxAlgebra(
		cbind(2*E,2*D,2*C,2*B,2*A),
		name="jac1",
		dimnames=list(NULL,c("e","d","c","z","a"))
	),
	mxAlgebra(cbind(C,B,-5*E,-5*D,0),
						name="jac2",
						dimnames=list(NULL,c("z","c","d","e","a")) ),
	mxAlgebra(
		cbind(0,0,0,3*B^2,3*A^2),
		name="jac3",dimnames=list(NULL,c("c","e","d","z","a"))),
	mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
)
# powellmod5 <- mxRun(powellmod5,onlyFrontend=T)
# powellmod5$compute$steps$GD$verbose <- 3L
# powellmod5$compute$.persist <- TRUE
powellrun5 <- mxRun(powellmod5)
omxCheckCloseEnough(powellrun4$fitfunction$result, powellrun5$fitfunction$result)
omxCheckEquals(powellrun4$output$iterations,powellrun5$output$iterations)
omxCheckEquals(powellrun4$output$evaluations,powellrun5$output$evaluations)
omxCheckCloseEnough(coef(powellrun4),coef(powellrun5))
summary(powellrun5)
