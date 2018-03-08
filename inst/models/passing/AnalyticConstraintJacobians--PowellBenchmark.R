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


library(OpenMx)
#mxOption(NULL,"Default optimizer","SLSQP")

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
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1",jac="jac1" ),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2",jac="jac2" ),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3",jac="jac3" ),
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


if(mxOption(NULL,"Default optimizer")=="NPSOL"){
  #Analytic Jacobians should, if nothing else, cut down on the number of fitfunction evaluations:
	omxCheckEquals(omxGreaterThan(powellrun1$output$evaluations,powellrun2$output$evaluations),1)
	
	#At the solution, bounds should not be active, and equality constraints should be satisfied:
	omxCheckEquals(powellrun1$compute$steps$GD$output$istate,c(0,0,0,0,0,3,3,3))
	omxCheckEquals(powellrun2$compute$steps$GD$output$istate,c(0,0,0,0,0,3,3,3))
	
	#The numerical and analytic Jacobians should agree closely:
	omxCheckCloseEnough(a=powellrun1$compute$steps$GD$output$constraintJacobian,b=powellrun2$compute$steps$GD$output$constraintJacobian,
											epsilon=1e-5)
	
	#Check naming of constraint-related information:
	omxCheckEquals(
		names(powellrun1$output$constraintFunctionValues),
		c("PowellBenchmarkNoJacobians.c1[1,1]","PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(
		rownames(powellrun1$output$constraintJacobian),
		c("PowellBenchmarkNoJacobians.c1[1,1]","PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun1$output$constraintJacobian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(
		names(powellrun1$output$LagrangeMultipliers),
		c("x1.bound","x2.bound","x3.bound","x4.bound","x5.bound","PowellBenchmarkNoJacobians.c1[1,1]",
			"PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(
		names(powellrun1$output$istate),
		c("x1.bound","x2.bound","x3.bound","x4.bound","x5.bound","PowellBenchmarkNoJacobians.c1[1,1]",
			"PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(
		names(powellrun2$output$constraintFunctionValues),
		c("PowellBenchmarkWithJacobians.c1[1,1]","PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(
		rownames(powellrun2$output$constraintJacobian),
		c("PowellBenchmarkWithJacobians.c1[1,1]","PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun2$output$constraintJacobian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(
		names(powellrun2$output$LagrangeMultipliers),
		c("x1.bound","x2.bound","x3.bound","x4.bound","x5.bound","PowellBenchmarkWithJacobians.c1[1,1]",
			"PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(
		names(powellrun2$output$istate),
		c("x1.bound","x2.bound","x3.bound","x4.bound","x5.bound","PowellBenchmarkWithJacobians.c1[1,1]",
			"PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
} else if(mxOption(NULL,"Default optimizer") %in% c("CSOLNP","SLSQP")){
	#Analytic Jacobians should, if nothing else, cut down on the number of fitfunction evaluations:
	omxCheckEquals(omxGreaterThan(powellrun1$output$evaluations,powellrun2$output$evaluations),1)
	
	#At the solution, equality constraints should be satisfied within feasibility tolerance:
	omxCheckCloseEnough(powellrun1$compute$steps$GD$output$constraintFunctionValues,c(0,0,0),
											as.numeric(mxOption(NULL,"Feasibility tolerance")))
	omxCheckCloseEnough(powellrun2$compute$steps$GD$output$constraintFunctionValues,c(0,0,0),
											as.numeric(mxOption(NULL,"Feasibility tolerance")))
	
	#The numerical and analytic Jacobians should agree closely:
	omxCheckCloseEnough(a=powellrun1$compute$steps$GD$output$constraintJacobian,b=powellrun2$compute$steps$GD$output$constraintJacobian,
											epsilon=1e-5)
	
	#Check naming of constraint-related information:
	omxCheckEquals(
		names(powellrun1$output$constraintFunctionValues),
		c("PowellBenchmarkNoJacobians.c1[1,1]","PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(
		rownames(powellrun1$output$constraintJacobian),
		c("PowellBenchmarkNoJacobians.c1[1,1]","PowellBenchmarkNoJacobians.c2[1,1]","PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun1$output$constraintJacobian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(
		names(powellrun1$output$LagrangeMultipliers),
		c("PowellBenchmarkNoJacobians.c1[1,1]","PowellBenchmarkNoJacobians.c2[1,1]",
			"PowellBenchmarkNoJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun1$output$LagrHessian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(rownames(powellrun1$output$LagrHessian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(
		names(powellrun2$output$constraintFunctionValues),
		c("PowellBenchmarkWithJacobians.c1[1,1]","PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(
		rownames(powellrun2$output$constraintJacobian),
		c("PowellBenchmarkWithJacobians.c1[1,1]","PowellBenchmarkWithJacobians.c2[1,1]","PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun2$output$constraintJacobian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(
		names(powellrun2$output$LagrangeMultipliers),
		c("PowellBenchmarkWithJacobians.c1[1,1]","PowellBenchmarkWithJacobians.c2[1,1]",
			"PowellBenchmarkWithJacobians.c3[1,1]")
	)
	omxCheckEquals(colnames(powellrun2$output$LagrHessian), c("x1","x2","x3","x4","x5"))
	omxCheckEquals(rownames(powellrun2$output$LagrHessian), c("x1","x2","x3","x4","x5"))
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
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1" ),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2" ),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3" ),
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
omxCheckCloseEnough(powellrun2$fitfunction$result, powellrun3$fitfunction$result,1e-12)
if(mxOption(NULL,"Default optimizer")=="NPSOL"){
	omxCheckEquals(powellrun2$output$iterations,powellrun3$output$iterations)
	omxCheckEquals(powellrun2$output$evaluations,powellrun3$output$evaluations)
	omxCheckCloseEnough(powellrun2$compute$steps$GD$output$constraintJacobian, powellrun3$compute$steps$GD$output$constraintJacobian)
}
omxCheckCloseEnough(coef(powellrun2),coef(powellrun3),1e-7)
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
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1" ),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2" ),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3" ),
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
omxCheckCloseEnough(coef(powellrun3),coef(powellrun4),1e-6)
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
	mxConstraint( (A^2)+(B^2)+(C^2)+(D^2)+(E^2) - 10 == 0, name="c1",jac="jac1" ),
	mxConstraint(B*C - 5*D*E == 0, name="c2",jac="jac2" ),
	mxConstraint(A^3 + B^3 + 1 == 0, name="c3",jac="jac3" ),
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


#For the sake of comparison, this model uses Jacobians but no gradient:
powellmod6 <- mxModel(
	"PowellBenchmarkWithJacobiansButNoGradient",
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),labels=paste("x",1:5,sep=""),name="X"),
	mxAlgebra( exp(prod(X)), name="powellfunc"),
	# mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]), 
	# 					 name="objgrad", 
	# 					 dimnames=list(NULL,paste("x",1:5,sep="")) ),
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1",jac="jac1"),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2",jac="jac2" ),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3",jac="jac3" ),
	mxAlgebra(cbind(2*X[1,1],2*X[1,2],2*X[1,3],2*X[1,4],2*X[1,5]),name="jac1",
						dimnames=list(NULL,paste("x",1:5,sep=""))
	),
	mxAlgebra(cbind(0,X[1,3],X[1,2],-5*X[1,5],-5*X[1,4]),name="jac2",
						dimnames=list(NULL,paste("x",1:5,sep=""))),
	mxAlgebra(cbind(3*X[1,1]^2, 3*X[1,2]^2, 0, 0, 0),name="jac3",dimnames=list(NULL,paste("x",1:5,sep=""))),
	mxFitFunctionAlgebra(algebra="powellfunc")
)
# powellmod6 <- mxRun(powellmod6,onlyFrontend=T)
# powellmod6$compute$steps$GD$verbose <- 3L
# powellmod6$compute$.persist <- TRUE
powellrun6 <- mxRun(powellmod6)
powellrun6$fitfunction$result
powellrun6$output$iterations
powellrun6$output$evaluations
summary(powellrun6)


#Finally, this model uses no analytic derivatives at all:
powellmod7 <- mxModel(
	"PowellBenchmarkNoJacobians",
	mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),labels=paste("x",1:5,sep=""),name="X"),
	mxAlgebra( exp(prod(X)), name="powellfunc"),
	# mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5],
	# 								 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]), 
	# 					 name="objgrad", 
	# 					 dimnames=list(NULL,paste("x",1:5,sep="")) ),
	mxConstraint(sum(X%^%2) - 10 == 0, name="c1"),
	mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2"),
	mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3"),
	mxFitFunctionAlgebra(algebra="powellfunc")
)
# powellmod7 <- mxRun(powellmod7,onlyFrontend=T)
# powellmod7$compute$steps$GD$verbose <- 3L
# powellmod7$compute$.persist <- TRUE
powellrun7 <- mxRun(powellmod7)
powellrun7$fitfunction$result
powellrun7$output$iterations
powellrun7$output$evaluations
summary(powellrun7)

tbl <- data.frame(
	c("Yes","Yes","No","No"),c("No","Yes","Yes","No"),
	c(powellrun1$output$evaluations,powellrun2$output$evaluations,powellrun6$output$evaluations,
		powellrun7$output$evaluations), row.names=c("powellrun1","powellrun2","powellrun6","powellrun7"),
	stringsAsFactors=F
	)
colnames(tbl) <- c("Gradient?","Jacobians?","Fitfunction evaluations")
print(tbl)

if(mxOption(NULL,"Default optimizer") == "SLSQP"){
	#With GDsearch, Nelder-Mead can get a good solution even though none of its initial vertices is feasible,
	#but it has to be held to a strict feasibility tolerance:
	foo <- mxComputeNelderMead(
		iniSimplexType="smartRight", xTolProx=1e-12, fTolProx=1e-8, eqConstraintMthd="GDsearch", 
		nudgeZeroStarts=F)
	plan <- omxDefaultComputePlan()
	plan$steps <- list(foo, plan$steps$RE)
	nmpowell <- mxModel(
		"PowellBenchmarkWithJacobians",
		plan,
		mxMatrix(type="Full",nrow=1,ncol=5,free=T,values=c(-2,2,2,-1,-1),labels=paste("x",1:5,sep=""),name="X"),
		mxAlgebra( exp(prod(X)), name="powellfunc"),
		mxAlgebra( cbind(powellfunc*X[1,2]*X[1,3]*X[1,4]*X[1,5],
										 powellfunc*X[1,1]*X[1,3]*X[1,4]*X[1,5],
										 powellfunc*X[1,1]*X[1,2]*X[1,4]*X[1,5],
										 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,5],
										 powellfunc*X[1,1]*X[1,2]*X[1,3]*X[1,4]), 
							 name="objgrad", 
							 dimnames=list(NULL,paste("x",1:5,sep="")) ),
		#Nelder-Mead benefits from analytic Jacobians if using eqConstraintMthd="GDsearch":
		mxConstraint(sum(X%^%2) - 10 == 0, name="c1",jac="jac1" ),
		mxConstraint(X[1,2]*X[1,3]-5*X[1,4]*X[1,5] == 0, name="c2",jac="jac2" ),
		mxConstraint(X[1,1]^3 + X[1,2]^3 + 1 == 0, name="c3",jac="jac3" ),
		mxAlgebra(cbind(2*X[1,1],2*X[1,2],2*X[1,3],2*X[1,4],2*X[1,5]),name="jac1",
							dimnames=list(NULL,paste("x",1:5,sep=""))
		),
		mxAlgebra(cbind(0,X[1,3],X[1,2],-5*X[1,5],-5*X[1,4]),name="jac2",
							dimnames=list(NULL,paste("x",1:5,sep=""))),
		mxAlgebra(cbind(3*X[1,1]^2, 3*X[1,2]^2, 0, 0, 0),name="jac3",dimnames=list(NULL,paste("x",1:5,sep=""))),
		mxFitFunctionAlgebra(algebra="powellfunc",gradient="objgrad")
	)
	nmpowell <- mxOption(nmpowell,"Feasibility tolerance",0.001)
	nmprun <- mxRun(nmpowell)
	nmprun$compute$steps[[1]]$output$constraintFunctionValues
	omxCheckCloseEnough(coef(powellrun2),coef(nmprun),0.02)
}
