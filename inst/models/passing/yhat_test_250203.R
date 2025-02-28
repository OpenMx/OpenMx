#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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

set.seed(1234)
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")

# `REML=FALSE`, 'yhat' provided: ####

ge <- mxExpectationGREML(V="V",yvars="y",REML=FALSE,yhat="foo")
gff <- mxFitFunctionGREML(autoDerivType="numeric")

testmod <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge,
	gff
)

testrun <- mxRun(testmod)
summary(testrun)
omxCheckCloseEnough(testrun$output$estimate[1],var(dat[,"y"])*99/100,1e-6)
omxCheckCloseEnough(testrun$output$estimate[2],mean(dat[,"y"]),1e-6)
m <- lm(dat[,"y"]~1)
omxCheckCloseEnough(testrun$output$fit,-2*logLik(m),1e-4)
omxCheckCloseEnough(testrun$output$standardErrors[1],sqrt((2*(var(dat[,"y"])*99/100)^2)/100),1e-5)
omxCheckCloseEnough(testrun$output$standardErrors[2],sqrt(var(dat[,"y"])*99/10000),1e-6)

# `REML=FALSE`, 'yhat' empty: ####

ge2 <- mxExpectationGREML(V="V",yvars="y",Xvars="x",addOnes=F,REML=F)

testmod2 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge2,
	gff
)
testrun2 <- mxRun(testmod2)
( sm <- summary(testrun2) )
omxCheckCloseEnough(testrun2$output$estimate[1],var(dat[,"y"])*99/100,1e-7)
omxCheckCloseEnough(sm$GREMLfixeff$coeff[1],mean(dat[,"y"]),1e-7)
omxCheckCloseEnough(testrun2$output$fit,-2*logLik(m),1e-4)
omxCheckCloseEnough(testrun2$output$standardErrors[1],sqrt((2*(var(dat[,"y"])*99/100)^2)/100),1e-5)
omxCheckCloseEnough(sm$GREMLfixeff$se[1],sqrt(var(dat[,"y"])*99/10000),1e-7)

# `REML=TRUE`: ####

ge3 <- mxExpectationGREML(V="V",yvars="y",Xvars="x",addOnes=F)

testmod3 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge3,
	gff
)
testrun3 <- mxRun(testmod3)
( sm <- summary(testrun3) )
omxCheckCloseEnough(testrun3$output$estimate[1],var(dat[,"y"]),1e-7)
omxCheckCloseEnough(sm$GREMLfixeff$coeff[1],mean(dat[,"y"]),1e-7)
omxCheckCloseEnough(testrun3$output$standardErrors[1],sqrt((2*(var(dat[,"y"]))^2)/100),1e-3)
omxCheckCloseEnough(sm$GREMLfixeff$se[1],sqrt(var(dat[,"y"])/100),1e-7)

####
#### Re-run the 3 models, now with semi-analytic derivatives: ####
####
gff4 <- mxFitFunctionGREML(autoDerivType="semiAnalyt")

# `REML=FALSE`, 'yhat' provided: semi-analytic derivatives Not Yet Implemented ####

testmod4 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge,
	gff4
)
testrun4 <- mxRun(testmod4)
# There is something wrong with the second derivatives...:
( sm <- summary(testrun4) )
omxCheckCloseEnough(testrun4$output$estimate[1],var(dat[,"y"])*99/100,1e-7)
omxCheckCloseEnough(testrun4$output$estimate[2],mean(dat[,"y"]),1e-7)
omxCheckCloseEnough(testrun4$output$fit,-2*logLik(m),1e-4)
# omxCheckWarning(
# 	mxRun(testmod4),
# 	"use of semi-analytic derivatives with 'yhat' is Not Yet Implemented; numeric derivatives will be used instead"
# )

# `REML=FALSE`, 'yhat' empty: ####

testmod5 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge2,
	gff4
)
testrun5 <- mxRun(testmod5)
( sm <- summary(testrun5) )
omxCheckCloseEnough(testrun5$output$estimate[1],var(dat[,"y"])*99/100,1e-7)
omxCheckCloseEnough(sm$GREMLfixeff$coeff[1],mean(dat[,"y"]),1e-7)
omxCheckCloseEnough(testrun5$output$fit,-2*logLik(m),1e-4)
omxCheckCloseEnough(testrun5$output$standardErrors[1],sqrt((2*(var(dat[,"y"])*99/100)^2)/100),1e-5)
omxCheckCloseEnough(sm$GREMLfixeff$se[1],sqrt(var(dat[,"y"])*99/10000),1e-7)
omxCheckCloseEnough(testrun2$output$fit,testrun5$output$fit,1e-5)

# `REML=TRUE`: ####

testmod6 <- mxModel(
	"GREMLtest",
	mxData(observed = dat, type="raw", sort=FALSE),
	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
	mxAlgebra(I %x% Ve,name="V"),
	ge3,
	gff4
)
testrun6 <- mxRun(testmod6)
( sm <- summary(testrun6) )
omxCheckCloseEnough(testrun6$output$estimate[1],var(dat[,"y"]),1e-7)
omxCheckCloseEnough(sm$GREMLfixeff$coeff[1],mean(dat[,"y"]),1e-7)
omxCheckCloseEnough(testrun6$output$standardErrors[1],sqrt((2*(var(dat[,"y"]))^2)/100),1e-3)
omxCheckCloseEnough(sm$GREMLfixeff$se[1],sqrt(var(dat[,"y"])/100),1e-7)
omxCheckCloseEnough(testrun3$output$fit,testrun6$output$fit,1e-5)

# Analytic derivatives with explicit means model: ####

# testmod7 <- mxModel(
# 	"GREMLtest",
# 	mxData(observed = dat, type="raw", sort=FALSE),
# 	mxMatrix(type = "Full", nrow = 1, ncol=1, free=T, values = 2, labels = "ve", lbound = 0.0001, name = "Ve"),
# 	mxMatrix(type="Full",nrow=100,ncol=1,name="foo",free=T,values=0.12345,labels="bar"),
# 	mxMatrix(type="Unit",nrow=100,ncol=1,name="Uno",condenseSlots=T),
# 	mxMatrix("Iden",nrow=100,name="I",condenseSlots=T),
# 	mxAlgebra(I %x% Ve,name="V"),
# 	ge,
# 	mxFitFunctionGREML(dV=c(ve="I"),dyhat=c(bar="Uno"))
# )
# testrun7 <- mxRun(testmod7)
