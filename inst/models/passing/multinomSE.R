library(OpenMx)
#mxOption(NULL,"Default optimizer","SLSQP")
library(MASS)

# #Custom compute plan, to calculate Hessian despite constraints:
# plan <- omxDefaultComputePlan()
# plan$steps <- list(GD=plan$steps$GD)#, ND=plan$steps$ND, SE=plan$steps$SE, RD=plan$steps$RD, RE=plan$steps$RE)

#MxModel to estimate multinomial proportions from frequencies (something that can be done by hand):
m1 <- mxModel(
	"MultinomialWithLinearConstraints",
	#plan,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
	mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
	mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
						 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m1run <- mxRun(m1)
summary(m1run)
omxCheckCloseEnough(coef(m1run),c(0.43,0.22,0.2,0.15),1e-7)
#Estimated proportions:
p <- round(coef(m1run),2)

#Basis for the nullspace of the constraint Jacobian:
U <- Null(t(m1run$output$constraintJacobian))
#Repeated-sampling covariance matrix:
( cov1 <- U%*%solve(t(U)%*%(m1run$output$hessian/2)%*%U)%*%t(U) )
omxCheckCloseEnough(cov1, m1run$output$vcov)
#Analytic value of covariance matrix, for comparison:
( cov2 <- (diag(p)-outer(p,p))/100 )
omxCheckCloseEnough(cov2, m1run$output$vcov, 1e-8)
#^^^This test will fail if this script is run single-threaded.

#Compare to the "covariance matrix" you'd get directly from the Hessian (which is wrong):
solve(m1run$output$hessian/2)
#Compare to the "covariance matrix" you get from the optimizer's BFGS approximation to the Hessian (which is even wronger):
if(mxOption(NULL,"Default optimizer")=="NPSOL"){
	print(solve( t(m1run$output$hessianCholesky) %*% m1run$output$hessianCholesky/2 ))
} else{
	print(solve(m1run$output$LagrHessian/2))
}

