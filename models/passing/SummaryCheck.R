#
#   Copyright 2007-2014 The OpenMx Project
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
# Date: 2014.06.10
# Filename: SummaryCheck.R
# Purpose: Actually test that pieces of MxSummary() do what they're supposed 
#  to do.
# Procedure: If you're going to muck around with some part of summary,
#  first make sure this file tests that part.  If this file does not test the
#  part of summary that you want to change, then you must add to this test file
#  before modifying summary.  Once the test has been added, muck around with
#  summary to your heart's content.
# TODO add multigroup case to make sure this works
# TODO add model with constraint to make sure it adds back DoF
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load package and data

require(OpenMx)

data(demoTwoFactor)

#------------------------------------------------------------------------------
# Specify raw data model

vpf <- 5
varnames <- names(demoTwoFactor)
latnames <- c("F1", "F2")

fl <- mxPath(from=rep(latnames, each=vpf), to=varnames, values=.8, labels=paste("load_", varnames, sep=""), free=TRUE)
fc <- mxPath(latnames, connect="unique.pairs", arrows=2, free=c(FALSE, TRUE, FALSE), values=c(1, .5, 1), labels=c("varF1", "covF1F2", "varF2"))
mm <- mxPath(from="one", to=varnames, free=TRUE, values=colMeans(demoTwoFactor), labels=paste("mean_", varnames, sep=""))
rv <- mxPath(from=varnames, arrows=2, free=TRUE, values=.2, labels=paste("resid_", varnames, sep=""))

raw <- mxModel("Raw Test Model to Check MxSummary",
	type="RAM", manifestVars=varnames, latentVars=latnames,
	fl, fc, mm, rv,
	mxData(demoTwoFactor, "raw")
)

raw.fit <- mxRun(raw)

#------------------------------------------------------------------------------
# Specify Saturated raw data model

sat.fit <- omxSaturatedModel(raw, run=TRUE)


#------------------------------------------------------------------------------
# Specify cov data models
cov <- mxModel("Covariance Test Model to Check MxSummary",
	type="RAM", manifestVars=varnames, latentVars=latnames,
	fl, fc, rv,
	mxData(cov(demoTwoFactor), "cov", numObs=nrow(demoTwoFactor))
)

cov.fit <- mxRun(cov)

covm <- mxModel("Covariance and Means Test Model to Check MxSummary",
	type="RAM", manifestVars=varnames, latentVars=latnames,
	fl, fc, rv, mm,
	mxData(cov(demoTwoFactor), "cov", numObs=nrow(demoTwoFactor), means=colMeans(demoTwoFactor))
)

covm.fit <- mxRun(covm)


#------------------------------------------------------------------------------
# Get summaries of models

raw.sum <- summary(raw.fit)
sat.sum <- summary(sat.fit[[1]])
raws.sum <- summary(raw.fit, SaturatedLikelihood=sat.fit[[1]])
cov.sum <- summary(cov.fit)
covs.sum <- summary(covm.fit)


#------------------------------------------------------------------------------
# Compare Fit statistics to what they should be.

# Things to check for raw w/o sat, raw w/ sat, cov, cov & means, multigroup
#	number of observed statistics
#	TODO add multigroup case to make sure this works
#	degrees of freedom
#	-2LL, sat -2LL, numObs
#	chi-square, chi dof, chi p
#	AIC df AIV, param, BIC df, BIC param, BIC sample size
#	CFI, TLI, RMSEA



#	number of observed statistics
omxCheckEquals(raw.sum$observedStatistics, 5000)
omxCheckEquals(raws.sum$observedStatistics, 5000)
omxCheckEquals(cov.sum$observedStatistics, 55)
omxCheckEquals(covs.sum$observedStatistics, 65)

#	TODO add multigroup case to make sure this works


#	degrees of freedom
# TODO add model with constraint to make sure it adds back DoF
ep <- 31
omxCheckEquals(raw.sum$degreesOfFreedom, 5000-ep)
omxCheckEquals(raws.sum$degreesOfFreedom, 5000-ep)
omxCheckEquals(sat.sum$degreesOfFreedom, 5000-sat.sum$estimatedParameters)
omxCheckEquals(raws.sum$satDoF, 5000-sat.sum$estimatedParameters)
omxCheckEquals(cov.sum$degreesOfFreedom, 55-(ep-10)) #no means estimated
omxCheckEquals(covs.sum$degreesOfFreedom, 65-ep)


#	-2LL, sat -2LL, numObs
omxCheckWithinPercentError(raw.sum$Minus2LogLikelihood, 9236.675, percent=1e-4)
omxCheckTrue(is.na(raw.sum$SaturatedLikelihood))
omxCheckWithinPercentError(raws.sum$SaturatedLikelihood, 9186.911, percent=1e-4)
omxCheckEquals(raw.sum$numObs, 500)
omxCheckEquals(raws.sum$numObs, 500)
omxCheckEquals(cov.sum$numObs, 500)
omxCheckEquals(covs.sum$numObs, 500)


#	chi-square, chi dof, chi p
omxCheckTrue(is.na(raw.sum$Chi))
omxCheckTrue(!is.na(raws.sum$Chi))
#cov.sum$Chi
#covs.sum$Chi

# ChiDoF = df - sat.df, i.e. df = obsStat-ep, sat.df = obsStat-sat.ep, so ChiDoF = sat.ep - ep
chi.df <- 34
#omxCheckTrue(is.na(raw.sum$ChiDoF)) # is still populated, but it reported as NA when satLik is NA
omxCheckEquals(raws.sum$ChiDoF, chi.df)
omxCheckEquals(cov.sum$ChiDoF, chi.df)
omxCheckEquals(covs.sum$ChiDoF, chi.df)

omxCheckTrue(is.na(raw.sum$p))
omxCheckCloseEnough(raws.sum$p, pchisq(raws.sum$Chi, raws.sum$ChiDoF, lower.tail=FALSE), 1e-8)
omxCheckCloseEnough(cov.sum$p, pchisq(cov.sum$Chi, cov.sum$ChiDoF, lower.tail=FALSE), 1e-8)
omxCheckCloseEnough(covs.sum$p, pchisq(covs.sum$Chi, covs.sum$ChiDoF, lower.tail=FALSE), 1e-8)


#	AIC df, AIC param, BIC df, BIC param, BIC sample size
#	CFI, TLI, RMSEA


# RMSEA
omxCheckTrue(is.na(raw.sum$RMSEA))
omxCheckTrue(!is.na(raws.sum$RMSEA))
omxCheckCloseEnough(raws.sum$RMSEA, 0.03045151, .001)
omxCheckCloseEnough(raws.sum$RMSEACI, c(0, 0.05064688), .001)
omxCheckCloseEnough(cov.sum$RMSEA, 0.03035523, .001)
omxCheckCloseEnough(cov.sum$RMSEACI, c(0, 0.05056999), .001)

# CFI
omxCheckTrue(is.na(raw.sum$CFI))
omxCheckCloseEnough(cov.sum$CFI, 0.9978234, .001)

# TLI
omxCheckTrue(is.na(raw.sum$TLI))
omxCheckCloseEnough(cov.sum$TLI, 0.9971192, .001)

# Information Criteria
omxCheckCloseEnough(raw.sum$informationCriteria['AIC:','par'], 9298.67, .01)
omxCheckCloseEnough(raw.sum$informationCriteria['BIC:','par'], 9429.32, .01)

omxCheckCloseEnough(sat.sum$informationCriteria['AIC:','par'], 9316.91, .01)
omxCheckCloseEnough(sat.sum$informationCriteria['BIC:','par'], 9590.86, .01)

omxCheckCloseEnough(raws.sum$informationCriteria['AIC:','par'], 9298.67, .01)
omxCheckCloseEnough(raws.sum$informationCriteria['BIC:','par'], 9429.32, .01)

omxCheckCloseEnough(cov.sum$informationCriteria['AIC:','par'], 91.66, .01)
omxCheckCloseEnough(cov.sum$informationCriteria['BIC:','par'], 180.17, .01)

omxCheckCloseEnough(covs.sum$informationCriteria['AIC:','par'], 111.66, .01)
omxCheckCloseEnough(covs.sum$informationCriteria['BIC:','par'], 242.31, .01)


