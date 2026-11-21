#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
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
if(mxOption(NULL,"Default optimizer")=="CSOLNP"){stop("SKIP")}
#mxOption(NULL,"Default optimizer","SLSQP")
#mxOption(NULL,"Number of threads",2)

library(mvtnorm)

set.seed(170209)
Rtrue <- matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3,byrow=T)
latdat <- rmvnorm(n=5000,mean=rep(0,3),sigma=Rtrue)
latdat[,1] <- ifelse(latdat[,1]>1.64,1,0)
latdat[,2] <- ifelse(latdat[,2]>1.64,1,0)
latdat[,3] <- ifelse(latdat[,3]>1.64,1,0)
table(apply(latdat,1,function(x){paste(x[1],x[2],x[3],sep="")}))

obsdat <- data.frame(
	y1=mxFactor(latdat[,1],levels=c(0,1)),
	y2=mxFactor(latdat[,2],levels=c(0,1)),
	y3=mxFactor(latdat[,3],levels=c(0,1))
)

mxdat <- mxData(observed=obsdat,type="raw")

safeT <- mxConstraint(-1%x%Sigma<Zilch,name="safety")

m1 <- mxModel(
	"mod1",
	mxdat,
	mxMatrix(type="Lower",nrow=3,free=T,values=c(1,1e-7,1e-7,1,1e-7,1),labels=c("s11","s21","s31","s22","s32","s33"),name="S", lbound=-10, ubound=10),
	mxMatrix(type="Zero",nrow=1,ncol=3,name="Mu"),
	mxMatrix(type="Unit",nrow=3,ncol=1,name="ONE"),
	mxMatrix(type="Zero",nrow=3,ncol=3,name="Zilch"),
	mxMatrix(type="Full",nrow=1,ncol=3,free=T,values=1.64,labels=c("tau1","tau2","tau3"),name="Tau",ubound=5),
	mxAlgebra(S%*%t(S),name="Sigma"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="Sigma",means="Mu",dimnames=c("y1","y2","y3"),thresholds="Tau",threshnames=c("y1","y2","y3")),
	safeT,
	mxConstraint(diag2vec(Sigma)==ONE,name="identifying"),
	mxConstraint(2%x%diag2vec(Sigma)==2%x%ONE,name="extra")
)

#All three of the following warning tests should pass if all three expected warnings are being raised:
omxCheckWarning(mxRun(m1),"Ignoring constraint 'mod1.identifying[1]' because it is redundant (|0.000000| < 0.000000)")
omxCheckWarning(mxRun(m1),"Ignoring constraint 'mod1.identifying[2]' because it is redundant (|0.000000| < 0.000000)")
m2 <- omxCheckWarning(mxRun(m1),"Ignoring constraint 'mod1.identifying[3]' because it is redundant (|0.000000| < 0.000000)")

omxCheckEquals(m2$output$effectiveNumEqualityMxConstraints,3)
omxCheckEquals(nrow(m2$output$constraintJacobian),12)
s2 <- summary(m2)
print(s2)
omxCheckEquals(s2$constraints[[1]],3)
omxCheckEquals(names(s2$constraints),"the set of linearly independent equality MxConstraints")
omxCheckEquals(s2$observedStatistics,15003)
