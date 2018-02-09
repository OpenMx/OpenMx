library(mvtnorm)
library(Matrix)
library(OpenMx)
set.seed(180122)

#Number of simulees (participants):
N <- 1000

#True parameter values, for data generation:
truevals <- c(
	exposed_tau=qnorm(0.5), #<--Threshold #1 for exposure; 0
	va1=0.5,
	ve1=0.5,
	s_lat=0.5,
	va2=0.5,
	ve2=0.25
)

exposed_lat_mz <- rmvnorm(n=N/4,mean=c(0,0),sigma=(truevals["va1"]*matrix(c(1,1,1,1),nrow=2))+(truevals["ve1"]*diag(1,2)))
exposed_lat_dz <- rmvnorm(n=N/4,mean=c(0,0),sigma=(truevals["va1"]*matrix(c(1,0.5,0.5,1),nrow=2))+(truevals["ve1"]*diag(1,2)))
exposed_mz <- (exposed_lat_mz > truevals["exposed_tau"]) + 0.0
exposed_dz <- (exposed_lat_dz > truevals["exposed_tau"]) + 0.0

dep_mz <- truevals["s_lat"]*exposed_lat_mz + 
	rmvnorm(n=N/4,mean=c(0,0),sigma=(truevals["va2"]*matrix(c(1,1,1,1),nrow=2))+(truevals["ve2"]*diag(1,2)))
dep_dz <- truevals["s_lat"]*exposed_lat_dz + 
	rmvnorm(n=N/4,mean=c(0,0),sigma=(truevals["va2"]*matrix(c(1,0.5,0.5,1),nrow=2))+(truevals["ve2"]*diag(1,2)))

twindata2 <- data.frame(
	exposed1=mxFactor(c(exposed_mz[,1],exposed_dz[,1]),levels=c(0,1)),
	dep1=c(dep_mz[,1],dep_dz[,1]),
	exposed2=mxFactor(c(exposed_mz[,2],exposed_dz[,2]),levels=c(0,1)),
	dep2=c(dep_mz[,2],dep_dz[,2]),
	zyg=c(rep(1,N/4),rep(0.5,N/4))
)

#Make non-missing dependence conditional on exposure:
twindata2$dep1[twindata2$exposed1==0] <- NA
twindata2$dep2[twindata2$exposed2==0] <- NA

twinmod2 <- mxModel(
	"Twin2",
	mxData(observed=twindata2, type="raw", sort=FALSE),
	mxMatrix(
		type="Full",nrow=1,ncol=4,free=c(F,T,F,T),
		values=c(0,mean(c(twindata2$dep1,twindata2$dep2),na.rm=T),0,mean(c(twindata2$dep1,twindata2$dep2),na.rm=T)),
		labels=c("mu_exposed","mu_dep","mu_exposed","mu_dep"),name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=F,labels="data.zyg",name="Zyg"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.5,labels="va1",name="Va1"),
	mxAlgebra(1-Va1,name="Ve1"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.1,labels="s",name="S"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,
					 values=var(c(twindata2$dep1, twindata2$dep2),na.rm=T)/2,
					 labels="va2",name="Va2"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,
					 values=var(c(twindata2$dep1, twindata2$dep2),na.rm=T)/2,
					 labels="ve2",name="Ve2"),
	mxAlgebra(rbind(
		cbind(Va1,S*Va1),
		cbind(S*Va1,Va1*(S^2)+Va2)
	),name="A"),
	mxAlgebra(rbind(
		cbind(Ve1,S*Ve1),
		cbind(S*Ve1,Ve1*(S^2)+Ve2)
	),name="E"),
	mxAlgebra(A+E, name="V"),
	mxAlgebra(
		rbind(
			cbind(V, Zyg%x%A),
			cbind(Zyg%x%A, V)
		), name="Cov"
	),
	mxMatrix(
		type="Full",nrow=1,ncol=2,free=T,
		values=qnorm(mean(rbind(exposed_mz,exposed_dz),na.rm=T)),
		labels="tau",name="Tau"),
	mxExpectationNormal(covariance="Cov",means="Mu",dimnames=c("exposed1","dep1","exposed2","dep2"),thresholds="Tau",
											threshnames=c("exposed1","exposed2")),
	mxFitFunctionML(jointConditionOn="continuous")
)
#Works:
twinmod2 <- mxModel(twinmod2,
	mxComputeSequence(list(
		mxComputeOnce('fitfunction','fit'),
		mxComputeNumericDeriv(checkGradient=FALSE),
		mxComputeReportDeriv())))
twinmod3 <- twinmod2

twinmod2 <- mxRun(twinmod2)

twinmod3$fitfunction <- mxFitFunctionML(jointConditionOn="ordinal")
twinmod3 <- mxRun(twinmod3)

print(twinmod2$output$fit - twinmod3$output$fit)
omxCheckCloseEnough(twinmod2$output$fit, twinmod3$output$fit, 1e-3)
omxCheckCloseEnough(cor(twinmod2$output$gradient, twinmod3$output$gradient), 1, 1e-6)
omxCheckCloseEnough(cor(vech(twinmod2$output$hessian), vech(twinmod3$output$hessian)), 1, 1e-5)
