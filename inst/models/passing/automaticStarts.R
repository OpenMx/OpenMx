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

#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-01-28 19:49:13
# Filename: automaticStarts.R
# Purpose: Test automatic starting values with mxAutoStart()
#------------------------------------------------------------------------------

library(OpenMx)
data(demoOneFactor)


latents  = c("G") # the latent factor
manifests = names(demoOneFactor) # manifest variables to be modeled

m1 <- mxModel("One Factor", type = "RAM", 
	manifestVars = manifests, latentVars = latents, 
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2, values=-.2),
	mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
	mxPath(from = "one", to = manifests),
	mxData(demoOneFactor, type = "raw")
)

start <- Sys.time()
m1s <- mxAutoStart(m1)
stop <- Sys.time()
m1sr <- mxRun(m1s)

badStart <- coef(m1)

mxOption(NULL, "Number of Threads", 1) # otherwise error message changes

errmsg <- "The job for model 'One Factor' exited abnormally with the error message: fit is not finite (The continuous part of the model implied covariance (loc2) is not positive definite in data 'One Factor.data' row 111. Detail:
covariance =  matrix(c(    # 5x5
 -0.19, 0.01, 0.01, 0.01, 0.01
, 0.01, -0.19, 0.01, 0.01, 0.01
, 0.01, 0.01, -0.19, 0.01, 0.01
, 0.01, 0.01, 0.01, -0.19, 0.01
, 0.01, 0.01, 0.01, 0.01, -0.19), byrow=TRUE, nrow=5, ncol=5)
)"
warnmsg <- "In model 'One Factor' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()"

omxCheckWarning( omxCheckError( mxRun(m1), message=errmsg ), warnmsg)


# change start values by hand
m1 <- omxSetParameters(m1, values=0, labels=names(coef(m1)))
m1r <- mxRun(m1)


cmp <- cbind(BadStart=badStart, UserStart=coef(m1), AutoStart=coef(m1s), UserEst=coef(m1r), AutoEst=coef(m1sr))


round(cmp, 3)

omxCheckCloseEnough(cor(cmp[,-c(1:2)]), matrix(1, 3, 3), 1e-6)


cbind(AutoTime=stop-start, RunFromAutoTime=summary(m1sr)$wallTime, RunFromUserTime=summary(m1r)$wallTime)


#------------------------------------------------------------------------------
# Test for covariance data

m2 <- m1
m2$data <- mxData(cov(demoOneFactor), type='cov',
	numObs=nrow(demoOneFactor), means=colMeans(demoOneFactor))
m2s <- mxAutoStart(m2)

omxCheckCloseEnough(coef(m1s), coef(m2s), 1e-6)


#------------------------------------------------------------------------------


