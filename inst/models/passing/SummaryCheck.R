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

sat.fit <- mxRefModels(raw, run=TRUE)

cmp <- mxCompare(sat.fit[['Saturated']], raw.fit)
omxCheckCloseEnough(cmp$diffLL[2], 49.764, 1e-2)
omxCheckCloseEnough(cmp$diffdf[2], 34)
omxCheckCloseEnough(cmp$p[2], 0.03961044, 1e-3)

cmp <- mxCompare(raw.fit, raw.fit)
omxCheckTrue(is.na(cmp$p[2]))

cmp <- omxCheckWarning(mxCompare(raw.fit, sat.fit[['Saturated']]),
		       "Model 'Raw Test Model to Check MxSummary' has more degrees of freedom than Saturated Raw Test Model to Check MxSummary which means that the models need to be compared in the opposite order")
omxCheckTrue(is.na(cmp$p[2]))

#------------------------------------------------------------------------------
# Specify a multiple group model


data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests),
      mxData(demoOneFactor, type="raw"))
manifests <- manifests[-1]
factorModelLess <- mxModel("OneFactorLess",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from = 'one', to = manifests),
      mxData(demoOneFactor, type="raw"))
mg <- mxModel("bob", factorModel, factorModelLess, mxFitFunctionMultigroup(c("OneFactor", "OneFactorLess")))


mg.fit <- mxRun(mg)

mg.sat <- mxRefModels(mg.fit, run=TRUE)

omxCheckEquals(dim(mxEval(satCov, mg.sat[[1]]$`Saturated OneFactor`))[1], 5)
omxCheckEquals(dim(mxEval(satCov, mg.sat[[1]]$`Saturated OneFactorLess`))[1], 4)
# the saturated model for OneFactorLess should be only on the variables used x2:x5, not all the variables x1:x5.

refStats <- lapply(mg.sat, function(model) {
	list(model$output$Minus2LogLikelihood, summary(model)$degreesOfFreedom)
})
refStats$Independence[[1]] <- mg.fit$output$fit - 100
ign <- omxCheckWarning(summary(mg.fit, refModels=refStats),
		"Your model may be mis-specified (and fit worse than an independence model), or you may be using the wrong independence model, see ?mxRefModels")

#------------------------------------------------------------------------------
# Specify a multiple group cov model

data(demoOneFactor)
latents  = c("G")
manifests = names(demoOneFactor)

m1 <- mxModel("model1", type = "RAM",
    manifestVars = manifests, latentVars = latents,
    mxPath(from = latents, to = manifests),
    mxPath(from = manifests, arrows = 2),
    mxPath(from = latents, arrows = 2, free = F, values = 1.0),
    mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)
m2 <- mxModel("model2", type = "RAM",
    manifestVars = manifests, latentVars = latents,
    mxPath(from = latents, to = manifests),
    mxPath(from = manifests, arrows = 2),
    mxPath(from = latents, arrows = 2, free = F, values = 1.0),
    mxData(cov(demoOneFactor), type = "cov", numObs = 500)
)

m3 = mxModel("bob", m1, m2,
    mxFitFunctionMultigroup(c("model1.fitfunction","model2.fitfunction"))
)


m3r = mxRun(m3)
ref <- mxRefModels(m3r, run=TRUE)

omxCheckError(summary(m3r, SaturatedLikelihood=ref),
              "List of length two (illegal argument) passed to 'SaturatedLikelihood' argument of summary function. You probably meant to use the refModels argument instead.")

#TODO add checking of fit stats for this model


#------------------------------------------------------------------------------
# Specify cov data models
cov <- mxModel("Covariance Test Model to Check MxSummary",
	type="RAM", manifestVars=varnames, latentVars=latnames,
	fl, fc, rv,
	mxData(cov(demoTwoFactor), "cov", numObs=nrow(demoTwoFactor))
)

cov.fit <- mxRun(cov)
cov.sat <- mxRefModels(cov.fit, run=TRUE)

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
rawr.sum <- summary(raw.fit, refModels=sat.fit)
cov.sum <- summary(cov.fit)
covs.sum <- summary(covm.fit)
mg.sum <- summary(mg.fit, refModels=mg.sat)


#------------------------------------------------------------------------------
# Compare Fit statistics to what they should be.

# Things to check for raw w/o sat, raw w/ sat, cov, cov & means, multigroup
#	number of observed statistics
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
omxCheckEquals(mg.sum$observedStatistics, 5*500 + 4*500)



#	degrees of freedom
# TODO add model with constraint to make sure it adds back DoF
ep <- 31
omxCheckEquals(raw.sum$degreesOfFreedom, 5000-ep)
omxCheckEquals(raws.sum$degreesOfFreedom, 5000-ep)
omxCheckEquals(sat.sum$degreesOfFreedom, 5000-sat.sum$estimatedParameters)
omxCheckEquals(raws.sum$satDoF, 5000-sat.sum$estimatedParameters)
omxCheckEquals(cov.sum$degreesOfFreedom, 55-(ep-10)) #no means estimated
omxCheckEquals(covs.sum$degreesOfFreedom, 65-ep)
omxCheckEquals(mg.sum$degreesOfFreedom, 5*500+4*500-27)


#	-2LL, sat -2LL, numObs
omxCheckWithinPercentError(raw.sum$Minus2LogLikelihood, 9236.675, percent=1e-4)
omxCheckTrue(is.na(raw.sum$SaturatedLikelihood))
omxCheckWithinPercentError(raws.sum$SaturatedLikelihood, 9186.911, percent=1e-4)
omxCheckWithinPercentError(cov.sum$SaturatedLikelihood, 7.520532, percent=1e-4)
omxCheckWithinPercentError(covs.sum$SaturatedLikelihood, 7.520532, percent=1e-4)
omxCheckEquals(raw.sum$numObs, 500)
omxCheckEquals(raws.sum$numObs, 500)
omxCheckEquals(cov.sum$numObs, 500)
omxCheckEquals(covs.sum$numObs, 500)
omxCheckEquals(mg.sum$numObs, 1000)


#	chi-square, chi dof, chi p
omxCheckTrue(is.na(raw.sum$Chi))
omxCheckTrue(!is.na(raws.sum$Chi))
omxCheckCloseEnough(cov.sum$Chi, 49.75399, epsilon=1e-3)
omxCheckCloseEnough(covs.sum$Chi, 49.75399, epsilon=1e-3)


# ChiDoF = df - sat.df, i.e. df = obsStat-ep, sat.df = obsStat-sat.ep, so ChiDoF = sat.ep - ep
chi.df <- 34
#omxCheckTrue(is.na(raw.sum$ChiDoF)) # is still populated, but it reported as NA when satLik is NA
omxCheckEquals(raws.sum$ChiDoF, chi.df)
omxCheckEquals(cov.sum$ChiDoF, chi.df)
omxCheckEquals(covs.sum$ChiDoF, chi.df)
omxCheckEquals(mg.sum$ChiDoF, 7)

omxCheckCloseEnough(mg.sum$Chi, 13.1811, 0.001)

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
omxCheckCloseEnough(mg.sum$RMSEA, 0.02971556, .001)
omxCheckCloseEnough(mg.sum$RMSEACI, c(0, 0.0582353), .001)

# CFI
omxCheckTrue(is.na(raw.sum$CFI))
omxCheckCloseEnough(cov.sum$CFI, 0.9978234, .001)
omxCheckCloseEnough(rawr.sum$CFI, 0.997814, .001)
omxCheckCloseEnough(mg.sum$CFI, 0.999077, .001)

# TLI
omxCheckTrue(is.na(raw.sum$TLI))
omxCheckCloseEnough(cov.sum$TLI, 0.9971192, .001)
omxCheckCloseEnough(rawr.sum$TLI, 0.9971067, .001)
omxCheckCloseEnough(mg.sum$TLI, 0.9978904, .001)

# Information Criteria
omxCheckCloseEnough(raw.sum$informationCriteria['AIC:','par'], 9298.67, .01)
omxCheckCloseEnough(raw.sum$informationCriteria['BIC:','par'], 9429.32, .01)

omxCheckCloseEnough(sat.sum$informationCriteria['AIC:','par'], 9316.91, .01)
omxCheckCloseEnough(sat.sum$informationCriteria['BIC:','par'], 9590.86, .01)

omxCheckCloseEnough(raws.sum$informationCriteria['AIC:','par'], 9298.67, .01)
omxCheckCloseEnough(raws.sum$informationCriteria['BIC:','par'], 9429.32, .01)

omxCheckCloseEnough(cov.sum$informationCriteria['AIC:','par'], 91.75, .01)
omxCheckCloseEnough(cov.sum$informationCriteria['BIC:','par'], 180.26, .01)

omxCheckCloseEnough(covs.sum$informationCriteria['AIC:','par'], 111.75, .01)
omxCheckCloseEnough(covs.sum$informationCriteria['BIC:','par'], 242.40, .01)

omxCheckCloseEnough(mg.sum$informationCriteria[c(1:4,6)], c(-6936.958, -28889.347, 2063.042, 2195.552, 2109.798), .01)

# -----------------------------------

library(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G1")
fit1 <- mxRun(mxModel(model="One Factor", type="RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxPath(from = latents, to=manifests),
                      mxPath(from = manifests, arrows = 2),
                      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
                      mxPath(from = 'one', manifests, free=FALSE),
                      mxData(cov(demoOneFactor), colMeans(demoOneFactor),
                             type = "cov", numObs = 500)
))

fit2 <- mxRun(mxModel(model="One Factor", type="RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxPath(from = latents, to=manifests, values=0, free=c(F,T,T,T,T)),
                      mxPath(from = manifests, arrows = 2),
                      mxPath(from = latents, arrows = 2, free = FALSE, values = 1.0),
                      mxPath(from = 'one', manifests, free=FALSE),
                      mxData(cov(demoOneFactor), colMeans(demoOneFactor),
                             type = "cov", numObs = 500)
))

got <- mxCompare(fit1, fit2, boot=T, replications = 10)
omxCheckCloseEnough(got[2,'p'], 0, .01)
