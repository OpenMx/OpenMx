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


require(OpenMx)
require(mvtnorm)

set.seed(21204)
Sigma <- matrix(0.5,4,4)
diag(Sigma) <- 1
Dat <- rmvnorm(n=26,mean=rep(0,4),sigma=Sigma)
Dat <- matrix(t(Dat))
colnames(Dat) <- "y"
famid <- sort(rep(LETTERS,4))
CRM <- outer(famid,famid,"==") + 0
ge <- mxExpectationGREML(V="V", yvars="y")

Cmod <- mxModel(
	"GREMLtest",
	mxData(observed = Dat, type="raw", sort=FALSE),
	ge,
	mxFitFunctionML(),
	mxMatrix(type = "Diag", nrow = 104, ncol=104, free=T, values = 0.5, labels = "ve", 
					 lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.5,labels="vc",name="Vc"),
	mxMatrix(type="Full",nrow=104,ncol=104,free=F,values=CRM,name="C",condenseSlots = T),
	mxAlgebra( (Vc%x%C) + Ve, name="V")
)
Crun <- mxRun(Cmod)
Csumm <- summary(Crun)

#Compare to lme4 REML results:
omxCheckCloseEnough(Csumm$GREMLfixeff$coeff[1],-0.01313328,1e-5)
omxCheckCloseEnough(Csumm$GREMLfixeff$se[1],0.1533185,1e-4)
omxCheckCloseEnough(Crun$output$estimate,c(0.3797,0.5162),1e-3)

#Intercept should equal grand mean:
omxCheckCloseEnough(Csumm$GREMLfixeff$coeff[1],mean(Dat),1e-12)
