#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# -----------------------------------------------------------------------
# Script by Robert M. Kirkpatrick.
# AGES Workshop, October 2017, Richmond, Virginia, USA.

library(OpenMx) #library(OpenMx,lib.loc=.libPaths()[2])
# Only CSOLNP gets status code 0 with every bootstrap replication:
if(mxOption(NULL,"Default optimizer") != 'CSOLNP') stop("SKIP")
mxVersion()
#mxOption(NULL,"Default optimizer","CSOLNP")

set.seed(12345)

trail <- c(rep(1,4000), rep(2,6000))
#The sample is drawn from two different normal distributions.  40% of the sample is drawn from a normal distribution 
#with a mean of 1, and 60%, from one with a mean of 2 (both distributions have a variance of 1):
trailN <- sapply(trail, function(v) rnorm(1, mean=v))

#This is the null model.  It has only one free parameter, a mean.  It's a "mixture" of one class:
mix1 <- mxModel(
	"mix1",
	type="RAM", #<--This model is being set up via RAM-path specification.
	mxData(data.frame(ob=trailN), "raw"),
	manifestVars=c("ob"),
	mxPath("one", "ob", value=mean(trailN), free=TRUE), #<--Path for the mean.
	mxPath("ob", arrows=2, value=1, free=FALSE)) #<--In this scenario, we know that the true variance is 1.
mix1Fit <- mxRun(mix1)
summary(mix1Fit)

classes <- list()

#This is apparently how you set up the classes for a mixture model in OpenMx:
for (cl in 1:2) {
	classes[[cl]] <- mxModel(
		paste0("class", cl), type="RAM",
		manifestVars=c("ob"),
		mxPath("one", "ob", value=cl, free=TRUE),
		mxPath("ob", arrows=2, value=1, free=FALSE),
		mxFitFunctionML(vector=TRUE))
}

#This is the alternative model.  It's a mixture of two classes.  It has three free parameters--a mean for one class,
#a mean for the other, and a mixture proportion:
mix2 <- mxModel(
	"mix2", classes,
	mxData(data.frame(ob=trailN), "raw"),
	mxMatrix(values=1, nrow=1, ncol=2, free=c(FALSE,TRUE), name="weights"),
	mxExpectationMixture(paste0("class",1:2), scale="sum"),
	mxFitFunctionML())

mix2Fit <- mxTryHardOrig(mix2)
#Note that the free parameter(s) for the mixture proportion(s) are unnormalized weights:
summary(mix2Fit)
#Normalized proportions:
mix2Fit$expectation$output$weights

#LRT chi-square test gives a p-value on the order of 1e-54, saying we should reject the null hypothesis
#of only one class:
mxCompare(base=mix2Fit, comparison=mix1Fit)
#Bootstrap LRT test (slow!):
blrt <- mxCompare(base=mix2Fit, comparison=mix1Fit, boot=T, replications=100)
#^^^Not sure about warning--all status codes look good (see below).
#Bootstrap LRT also says we should reject the null of one mixture class, albeit with a much less extreme p-value:
blrt
omxCheckEquals(blrt$df,c(9997,9999))
omxCheckCloseEnough(blrt$diffLL[2],246.9217,0.0002)

( tbl1 <- table(attr(blrt@results,"bootData")[[1]]$statusCode) )
( tbl2 <- table(attr(blrt@results,"bootData")[[2]]$statusCode) )
( p <- mean( (attr(blrt@results,"bootData")[[1]]$fit - attr(blrt@results,"bootData")[[2]]$fit) > 246.9217) )
omxCheckEquals(tbl1[1],100)
omxCheckEquals(tbl2[1],100)
omxCheckCloseEnough(p,0.01,0.002)
