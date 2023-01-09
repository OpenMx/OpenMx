#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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
# Date: 2022-09-07
# Filename: WLSCompare.R
# Purpose: Test WLS behavior with mxCompare on single and multiple groups
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load packages

library(OpenMx)


#------------------------------------------------------------------------------
# Generate Data for 1 Factor Model

p <- 4
N <- 3200 #3200
man <- paste0('x', 1:4)


#--------------------------------------
# Stupidly simple model, but Mplus can't handle it
m1t <- mxModel(model="OneFOrd", type='RAM',
    latentVars=c('F'),
    manifestVars=man,
    mxPath('F', arrow=2, values=1.5, free=TRUE, labels='FVar'),
    mxPath(man, arrows=2, values=0.5, free=TRUE,
        labels=paste0('resid', 1:p)),
    mxPath('F', man, labels=paste0('load', 1:p),
        free=c(FALSE, rep(TRUE, p-1)), values=c(1, 1 + (2:p - p/2 - 1)/5)),
    mxPath('one', 'F', values=0.5, free=TRUE, labels='FMean'),
    mxPath('one', man, values=0, free=FALSE),
    mxThreshold(man, nThresh=3, free=c(TRUE, FALSE, TRUE),
        labels=outer(man, paste0('thr', 1:3), FUN=paste0))
    )

# set.seed(37)
# ds1 <- mxGenerateData(m1t, nrow=N)


#--------------------------------------
# Altered model so Mplus can deal
m1t <- mxModel(model="OneFOrd", type='RAM',
    latentVars=c('F'),
    manifestVars=man,
    # Esimate factor variacne
    mxPath('F', arrow=2, values=1.5, free=TRUE, labels='FVar'),
    # Fix residual variances to 1
    mxPath(man, arrows=2, values=1, free=FALSE,
        labels=paste0('resid', 1:p)),
    # Fix first factor loading
    mxPath('F', man, labels=paste0('load', 1:p),
        free=c(FALSE, rep(TRUE, p-1)), values=c(1, 1 + (2:p - p/2 - 1)/5)),
    # Fix factor mean
    mxPath('one', 'F', values=0, free=FALSE, labels='FMean'),
    # Fix residual means (intercepts) to 0
    mxPath('one', man, values=0, free=FALSE),
    # Estimate all thresholds
    mxThreshold(man, nThresh=3, free=c(TRUE, TRUE, TRUE),
        labels=t(outer(man, paste0('thr', 1:3), FUN=paste0)))
    )

set.seed(37)
ds1 <- mxGenerateData(m1t, nrow=N)

# write.table(ds1, file='dataWls1.dat', row.names=FALSE, col.names=FALSE)


#------------------------------------------------------------------------------
# Generate Data for tau-equivalent 1 Factor Model

m2t <- omxSetParameters(m1t, labels=paste0('load', 1:p), newlabels='load',
    values=1, free=FALSE, name="OneFOrdTau")

set.seed(41)
ds2 <- mxGenerateData(m2t, nrow=N)


#------------------------------------------------------------------------------
# Generate Data for 2 Factor Model

# TODO

#m3t <- mxModel(model="TwoFOrd")


#------------------------------------------------------------------------------
# Fit null and alternative models when null is false

# Data from "congeneric" factor model
# (i.e., general 1 factor model with various loadings)
# Fit congeneric and tau-equivalent models
# Null is false

m1 <- mxModel(m1t, mxData(ds1, 'raw'), mxFitFunctionWLS('DWLS'))
m2 <- omxSetParameters(m1, labels=paste0('load', 1:p), newlabels='load',
    values=1, free=FALSE, name="OneFOrdTau")

m1r <- mxRun(m1)
m2r <- mxRun(m2)

(cmp <- mxCompare(m1r, m2r))


#------------------------------------------------------------------------------
# Fit null and alternative models when null is true

# Data from tau-equivalent factor model
# Fit congeneric and tau-equivalent models
# Null is true

n1 <- mxModel(m1, mxData(ds2, 'raw'))
n2 <- mxModel(m2, mxData(ds2, 'raw'))

n1r <- mxRun(n1)
n2r <- mxRun(n2)

(ncmp <- mxCompare(n1r, n2r))


#------------------------------------------------------------------------------
# What is asymCov in relation to useWeight?

#plot(
#    diag(solve(diag(diag(m1r$data$observedStats$asymCov))))/N/N,
#    diag(m1r$data$observedStats$useWeight))
#abline(a=0, b=1)

# Check relation with omxCheckCloseEnough()
ascov <- diag(solve(diag(diag(m1r$data$observedStats$asymCov))))/N/N
aswei <- diag(m1r$data$observedStats$useWeight)
omxCheckCloseEnough(ascov, aswei, 1e-10)


#------------------------------------------------------------------------------
# Specify model in lavaan

runLavaan <- FALSE
if(runLavaan){
require(lavaan)

lavFull <- '
    # Fix first factor loadings
    f1 =~ 1*x1 + x2 + x3 + x4
    # Estimate factor variance
    f1 ~~ varF1*f1
    # Fix factor mean to 0
    f1 ~ 0*1
    # Fix intercepts to zero
    x1 ~ 0*1
    x2 ~ 0*1
    x3 ~ 0*1
    x4 ~ 0*1
    # Fix residual variances to 1
    x1 ~~ 1*x1
    x2 ~~ 1*x2
    x3 ~~ 1*x3
    x4 ~~ 1*x4
    # Estimate thresholds
    x1 | t1
    x1 | t2
    x1 | t3
    x2 | t1
    x2 | t2
    x2 | t3
    x3 | t1
    x3 | t2
    x3 | t3
    x4 | t1
    x4 | t2
    x4 | t3'

lavTau <- '
    # Fix first factor loadings
    f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
    # Estimate factor variance
    f1 ~~ varF1*f1
    # Fix factor mean to 0
    f1 ~ 0*1
    # Fix intercepts to zero
    x1 ~ 0*1
    x2 ~ 0*1
    x3 ~ 0*1
    x4 ~ 0*1
    # Fix residual variances to 1
    x1 ~~ 1*x1
    x2 ~~ 1*x2
    x3 ~~ 1*x3
    x4 ~~ 1*x4
    # Estimate thresholds
    x1 | t1
    x1 | t2
    x1 | t3
    x2 | t1
    x2 | t2
    x2 | t3
    x3 | t1
    x3 | t2
    x3 | t3
    x4 | t1
    x4 | t2
    x4 | t3'


lmo <- lavaan(lavFull, data=ds1, parameterization='theta', estimator='WLS')
lmt <- lavaan(lavTau, data=ds1, parameterization='theta', estimator='WLS')
lno <- lavaan(lavFull, data=ds2, parameterization='theta', estimator='WLS')
lnt <- lavaan(lavTau, data=ds2, parameterization='theta', estimator='WLS')


lmo <- lavaan(lavFull, data=ds1, parameterization='theta', estimator='WLS')
summary(lmo)$test$standard$stat # chi
lmo <- lavaan(lavFull, data=ds1, parameterization='theta', estimator='WLSMV', test='mean.var.adjusted')
summary(lmo)$test$standard$stat # fit
summary(lmo)$test$mean.var.adjusted$scaling.factor # mvadjust
summary(lmo)$test$mean.var.adjusted$df # dfstar
summary(lmo)$test$mean.var.adjusted$stat # chimv
lmo <- lavaan(lavFull, data=ds1, parameterization='theta', estimator='WLSMV', test='Satorra.Bentler')
summary(lmo)$test$satorra.bentler$scaling.factor # madjust
summary(lmo)$test$satorra.bentler$stat # chim

}

# Hard coded results from lavaan 2023-01-09, packageVersion('lavaan') 0.6.12
lavo <- c(chi=2.027596, fit=0.9178031, mvadjust=0.4546067,
    dfstar=1.998749, chimv=2.018895,
    madjust=0.4543222, chim=2.020159)

opmo <- c(chi=m1r$output$chi, fit=m1r$output$fit, mvadjust=m1r$output$chiMVadjust,
    dfstar=m1r$output$chiDoFstar, chimv=m1r$output$chiMV,
    madjust=m1r$output$chiMadjust, chim=m1r$output$chiM)

omxCheckCloseEnough(opmo, lavo, .0005)


# lmo <- lavaan(lavFull, data=ds1, parameterization='theta', estimator='WLSM')
# lmt <- lavaan(lavTau, data=ds1, parameterization='theta', estimator='WLSM')
# anova(lmo, lmt, method='satorra.bentler.2001')
# cmp
omxCheckCloseEnough(cmp$SBchisq[2], 52.236, .02)


# lno <- lavaan(lavFull, data=ds2, parameterization='theta', estimator='WLSM')
# lnt <- lavaan(lavTau, data=ds2, parameterization='theta', estimator='WLSM')
# anova(lno, lnt, method='satorra.bentler.2001')
# ncmp
omxCheckCloseEnough(ncmp$SBchisq[2], 1.7838, .02)



#------------------------------------------------------------------------------
# Check on analytic gradients for WLS

mgrad <- mxOption(m1, 'Major iterations', 1)
mrgrad <- suppressWarnings(mxRun(mgrad))

hi <- mxRun(mxModel(mrgrad, mxComputeJacobian()))
J <- hi$compute$output$jacobian

W <- mrgrad$data$observedStats$useWeight
obs <- c(
    c(mrgrad$data$observedStats$thresholds),
    vechs(mrgrad$data$observedStats$cov))
mi <- mxGetExpected(mrgrad, 'standVector')

# Hand check fit function computation
omxCheckCloseEnough(
    as.numeric(mrgrad$fitfunction$result),
    t(obs - mi) %*% W %*% (obs - mi),
    1e-10)

# Analytic and numeric gradients
anagrad <- -2 * t(J) %*% W %*% (obs - mi)
numgrad <- mrgrad$output$gradient

# plot(cbind(anagrad, numgrad))
# text(anagrad, numgrad, labels=names(numgrad))

omxCheckCloseEnough(anagrad, numgrad, 1e-5)
omxCheckCloseEnough(sum(abs(anagrad - numgrad)), 0, 1e-5)


#------------------------------------------------------------------------------
# Another form of  RAM gradient

# I need to account for the quadratic function from the
#  parameters to the summary statistics in the model.
#  That is, there's another step in the chain rule:
# RAM -> MVN -> summary stats -> fit function

# RAM gradiet/Jacobian
amat <- mrgrad$A$values
smat <- mrgrad$S$values
fmat <- mrgrad$F$values
imat <- diag(1, nrow=nrow(amat))
uinv <- solve(imat - amat)
dA <- matrix(0, nrow=nrow(amat), ncol=ncol(amat))
dA[2, 5] <- 1
dS <- matrix(0, nrow=nrow(amat), ncol=ncol(amat))
D <- -uinv %*% (imat - dA) %*% uinv
B <- D %*% smat %*% t(uinv)
dC <- fmat %*% (B + uinv %*% dS %*% t(uinv) + t(B)) %*% t(fmat)


#------------------------------------------------------------------------------
# End