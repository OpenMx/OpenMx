require(OpenMx)

gb <- function(m,verbose=FALSE,N=500){
	Sigma <- mxGetExpected(m,"covariance")
	A <- m$A$values
	C <- matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y")))
	SigmaInv <- solve(Sigma)
	edA <- -1*matrix(c(0,1,0,0),2,2)
	I_A <- solve(diag(2) - A)
	if(verbose){print("I_A:"); print(I_A)}
	firstPart <- -1*I_A%*%edA%*%Sigma
	if(verbose){print("firstPart:"); print(firstPart)}
	secondPart <- matrix(0,2,2)
	thirdPart <- -1*Sigma%*%t(edA)%*%t(I_A)
	if(verbose){print("thirdPart:"); print(thirdPart)}
	Der <- firstPart+secondPart+thirdPart
	CinvDer_trace <- sum(diag(SigmaInv%*%Der))
	if(verbose){message(paste0("CinvDer_trace: ",CinvDer_trace))}
	secondTerm <- sum(diag((N-1)/N*SigmaInv%*%C%*%SigmaInv%*%Der))
	if(verbose){message(paste0("secondTerm: ",secondTerm))}
	#return((N)*(CinvDer_trace - secondTerm))
	covPart <- (N)*(CinvDer_trace - secondTerm)
	if(verbose){message(paste0("covPart: ",covPart))}

	Nu <- c(1.0,1.5) - mxGetExpected(m,"means")
	if(verbose){print("Nu:");print(Nu)}
	dMu_dtheta <- c(0,mxEval(m1,m,T)[1])
	if(verbose){print("dMu_dtheta:");print(dMu_dtheta)}
	CinvNu <- SigmaInv %*% t(Nu)
	if(verbose){print("CinvNu:");print(CinvNu)}
	dMu_dtheta_CinvNu <- dMu_dtheta %*% CinvNu
	if(verbose){print("dMu_dtheta_CinvNu:");print(dMu_dtheta_CinvNu)}
	meanscorrection <- -N*dMu_dtheta_CinvNu
	if(verbose){message(paste0("meanscorrection: ",meanscorrection))}
	return(covPart + 2*meanscorrection)
}

# Only S paths free ####

plan <- mxComputeSequence(list(mxComputeOnce("fitfunction",c("fit","gradient")),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","Yes")
m1a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.1,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m1a <- mxRun(m1a)

plan3 <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=F,hessian=F),mxComputeReportDeriv(),mxComputeReportExpectation()))
mxOption(NULL,"Analytic gradients","No")
m1n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.1,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m1n <- mxRun(m1n)

omxCheckCloseEnough(m1a$output$gradient-m1n$output$gradient, c(0,0), 5e-7)


# Only A paths free ####

mxOption(NULL,"Analytic gradients","Yes")
m2a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m2a <- mxRun(m2a)
gb(m2a,T)
m2a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m2n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=1.0,free=F),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m2n <- mxRun(m2n)
m2n$output$gradient

omxCheckCloseEnough(m2a$output$gradient-m2n$output$gradient, 0, 5e-7)


# Only M paths free ####

mxOption(NULL,"Analytic gradients","Yes")
m3a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.1,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m3a <- mxRun(m3a)
m3a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m3n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.5,labels="b",free=F),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.1,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m3n <- mxRun(m3n)
m3n$output$gradient

omxCheckCloseEnough(m3a$output$gradient-m3n$output$gradient, c(0,0), 5e-7)


# A & M paths free ####

mxOption(NULL,"Analytic gradients","Yes")
m4a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m4a <- mxRun(m4a)
#gb(m4a,T)
m4a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m4n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m4n <- mxRun(m4n)
m4n$output$gradient

# omxCheckCloseEnough(m4a$output$gradient-m4n$output$gradient, c(0,0,0), 5e-7)


## A paths & m2 free ####

mxOption(NULL,"Analytic gradients","Yes")
m5a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(1.0,0.1),free=c(F,T)),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m5a <- mxRun(m5a)
#gb(m4a)
m5a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m5n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(1.0,0.1),free=c(F,T)),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m5n <- mxRun(m5n)
m5n$output$gradient

omxCheckCloseEnough(m5a$output$gradient-m5n$output$gradient, c(0,0), 5e-8)


## A paths & m1 free ####

mxOption(NULL,"Analytic gradients","Yes")
m6a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(0.2,1.0),free=c(T,F)),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m6a <- mxRun(m6a)
#gb(m6a)
m6a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m6n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=c(0.998,0.7485),labels=c("v1","v2"),free=F),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=c(0.2,1.0),free=c(T,F)),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m6n <- mxRun(m6n)
m6n$output$gradient


# All paths free ####

mxOption(NULL,"Analytic gradients","Yes")
m7a <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.5,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m7a <- mxRun(m7a)
m7a$output$gradient

mxOption(NULL,"Analytic gradients","No")
m7n <- mxModel(
	"TwoByTwo",
	type="RAM",
	plan3,
	manifestVars = c("x","y"),
	mxPath(from=c("x","y"),arrows=2, values=0.5,labels=c("v1","v2"),free=T),
	mxPath(from="x",to="y",arrows=1,values=0.1,labels="b",free=T),
	mxPath(from="one",to=c("x","y"),arrows=1,labels=c("m1","m2"),values=0.2,free=T),
	mxData(matrix(c(1,.5,.5,1),2,2,dimnames=list(c("x","y"), c("x","y"))), type="cov", means=c(x=1.0,y=1.5), numObs=500)
)
m7n <- mxRun(m7n)
m7n$output$gradient
