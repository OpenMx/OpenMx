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
# Date: 2024.03.39
# Filename: test-ModelIdentification2.R
# Purpose: Create further and more intense checks of the model identification
#  checking function.
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# Load OpenMx

library(OpenMx)
library(testthat)
context("Model Identification 2")


#------------------------------------------------------------------------------
# Outline of Tests

# X 1 RAM Factor Model
#       Works
# X 2 Latent growth curve model
#       Works
# X 3 Latent growth curve model with definition variable loadings
#       Works but
# TODO Check dimensions of Jacobian to verify that mxCheckID only used single values?
# X 4 Single group ACE model using definition variable for relatedness
#       Works
# X 5 Factor model with factor loading set to b0 + b1*data.sex
#       Works in theory but fails because defvar sorting
# What if first and last rows have the same defvar values?
# X 6 Multiple group (normal) model with definition variable on means
#       Works
# X 7 Factor model that uses mxConstraint for fixing factor variance
#       Works but
#       Requires data to do a data-independent identification check
# X 8 Multiple group (normal) model that uses mxConstraint to set
#       variance to b0 + b1*sex
#       Fails.  As designed?


#------------------------------------------------------------------------------
# TODO Add checks for models with definition variables but are not identified

# Factor model with b0+b1*data.sex loading but data has only one value of sex.
#  This model is not identified.  The b1 parameter should show as non-identified.



#------------------------------------------------------------------------------
# TODO Add empirical identification functionality
# TODO Add empirical identification checks


#------------------------------------------------------------------------------
# Check 1 RAM Factor Model

mv <- letters[23:26]
c1 <- mxModel('RAM Factor Model', type='RAM', manifestVars=mv, latentVars='H',
              mxPath('H', mv, values=c(.2, .3, .4, .5),
                     labels=paste0('load', 1:4)),
              mxPath('H', arrows=2, values=1, free=FALSE, labels='hvar'),
              mxPath(mv, arrows=2, values=c(.8, .9, 1.1, 1.2)),
              mxPath('one', mv, values=2:5))
c1id <- mxCheckIdentification(c1)
omxCheckTrue(c1id$status)


#------------------------------------------------------------------------------
# Check 2 Latent growth curve model

data(myLongitudinalData)
mld <- myLongitudinalData
mv <- names(mld)

lv <- c('Intercept', 'Slope')
mld[, mv] <- matrix(1:5, nrow=nrow(mld), ncol=5, byrow=TRUE)
c2 <- mxModel('Latent growth curve model', type='LISREL',
              manifestVars=list(ex=mv), latentVars=list(ex=lv),
              mxData(mld, 'raw'),
              mxPath(lv, arrows=2, connect='unique.pairs', values=c(1, .5, 1)),
              mxPath(lv[1], mv, values=1, free=FALSE),
              mxPath(lv[2], mv, values=-2:2, free=FALSE),
              mxPath(mv, arrows=2, values=.1, labels='evar'),
              mxPath('one', lv))
# c2r <- mxRun(c2)
# summary(c2r)
c2id <- mxCheckIdentification(c2)
omxCheckTrue(c2id$status)


#------------------------------------------------------------------------------
# Check 3 Latent growth curve model with definition variable loadings

data(myLongitudinalData)
mld <- myLongitudinalData
mv <- names(mld)
dv <- paste0('t', 1:5)
lv <- c('L0', 'L1')
mld[, dv] <- matrix(1:5, nrow=nrow(mld), ncol=5, byrow=TRUE)
c3 <- mxModel('defvar LGM', type='RAM', manifestVars=mv, latentVars=lv,
              mxData(mld, 'raw'),
              mxPath(lv, arrows=2, connect='unique.pairs', values=c(1, .5, 1)),
              mxPath(lv[1], mv, values=1, free=FALSE),
              mxPath(lv[2], mv, labels=paste0('data.', dv), free=FALSE),
              mxPath(mv, arrows=2, values=.1, labels='evar'),
              mxPath('one', lv))
# c3r <- mxRun(c3)
# summary(c3r)
c3id <- mxCheckIdentification(c3)
omxCheckTrue(c3id$status)
# Check dimensions of Jacobian to verify that mxCheckID only used single values
omxCheckTrue(nrow(c3id$jacobian) == 20)


#------------------------------------------------------------------------------
# Check 4 Single group ACE model using definition variable for relatedness

data(twinData)
twinVar <- names(twinData)
selVars <- c('ht1', 'ht2')
mzdzData <- subset(twinData, zyg %in% c(1, 3), c(selVars, 'zyg'))
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]
lv <- c('A1', 'A2', 'C1', 'C2', 'E1', 'E2')

c4 <- mxModel('defvar ACE', type='LISREL', manifestVars=list(ex=selVars),
              latentVars=list(ex=lv),
              mxPath(lv, arrows=2, values=1, free=FALSE),
              mxPath('one', selVars, labels='m'),
              mxPath(lv, selVars, values=.2,
                     labels=rep(c('a', 'c', 'e'), each=2)),
              mxPath(lv[1], lv[2], arrows=2, labels='data.RCoef', free=FALSE),
              mxPath(lv[3], lv[4], arrows=2, values=1, free=FALSE),
              mxData(mzdzData, 'raw'))
c4r <- mxRun(c4)
summary(c4r)

c4id <- mxCheckIdentification(c4r)
omxCheckTrue(c4id$status)


#------------------------------------------------------------------------------
# Check 5 Factor model with factor loading set to b0 + b1*data.sex

data(HS.ability.data)
mv <- c('visual', 'cubes', 'flags', 'straight')
selv <- c('id', 'Gender', 'agey', mv)
c5d <- HS.ability.data[,selv]
c5d$sex <- as.numeric(c5d$Gender) - 1

loadLab <- paste0('loads[', 1:4, ',1]')
c5 <- mxModel(type='RAM', manifestVars=mv, latentVars='V',
              mxPath('V', mv, values=c(.2, .3, .4, .5), labels=loadLab, free=FALSE),
              mxPath('V', arrows=2, values=1, free=FALSE, labels='hvar'),
              mxPath(mv, arrows=2, values=c(.8, .9, 1.1, 1.2)),
              mxPath('one', mv, values=2:5),
              mxMatrix('Full', 4, 2, values=c(.2, -.9), free=TRUE, name='B',
                       byrow=TRUE),
              mxMatrix('Full', 2, 1, values=c(1, NA), labels=c(NA, 'data.sex'),
                       name='x'),
              mxAlgebra(B %*% x, 'loads'),
              mxData(c5d, 'raw')
              )
c5r <- mxRun(c5)

summary(c5r)

# What if first and last rows have the same defvar values?
c5id <- mxCheckIdentification(c5r)
omxCheckTrue(c5id$status)


#------------------------------------------------------------------------------
# Check 6 Multiple group (normal) model with definition variable on means

data(example2)

mzd <- example2[example2$Zygosity %in% 'MZ', ]
dzd <- example2[example2$Zygosity %in% 'DZ', ]

#---------------------------------------
# Each group with its own data

c6mz <- mxModel('MZ', mxData(mzd, 'raw'),
                mxMatrix('Full', 1, 2, values=c(1, -1), free=TRUE,
                         labels=c('b0m', 'b1m'), name='BMZ'),
                mxMatrix('Full', 1, 1, values=.5, free=TRUE, labels='var',
                         name='V'),
                mxAlgebra(b0m + b1m*(data.TwinNum - 1), name='M'),
                mxExpectationNormal(covariance='V', means='M', dimnames='X'),
                mxFitFunctionML())
c6dz <- mxModel('DZ', mxData(dzd, 'raw'),
                mxMatrix('Full', 1, 2, values=c(1, -1), free=TRUE,
                         labels=c('b0m', 'b1m'), name='BDZ'),
                mxMatrix('Full', 1, 1, values=.5, free=TRUE, labels='var',
                         name='V'),
                mxAlgebra(b0m + b1m*(data.TwinNum - 1), name='M'),
                mxExpectationNormal(covariance='V', means='M', dimnames='X'),
                mxFitFunctionML())
c6 <- mxModel('mg', c6mz, c6dz, mxFitFunctionMultigroup(c('MZ', 'DZ')))
c6r <- mxRun(c6)
summary(c6r)

omxCheckTrue(imxHasDefinitionVariable(c6))
omxManifestModelByParameterJacobian(c6, defvar.row=c(MZ=14, DZ=201))
omxManifestModelByParameterJacobian(c6, defvar.row=c(MZ=14, DZ=199))

c6id <- mxCheckIdentification(c6)
omxCheckTrue(c6id$status)
# Evaluating the definition variable at multiple values
#  when there are multiple groups too

#---------------------------------------
# Inherited data from super-model


#------------------------------------------------------------------------------
# Check 7 Factor model that uses mxConstraint for fixing factor variance

mv <- letters[23:26]
c7 <- mxModel('Constraint RAM Factor Model', type='RAM',
              manifestVars=mv, latentVars='H',
              mxPath('H', mv, values=c(.2, .3, .4, .5),
                     labels=paste0('load', 1:4)),
              mxPath('H', arrows=2, values=1, labels='hvar'),
              mxPath(mv, arrows=2, values=c(.8, .9, 1.1, 1.2)),
              mxPath('one', mv, values=2:5), mxData(data.frame(w=NA, x=NA, y=NA, z=NA), 'raw'),
              mxConstraint(hvar == 1))
c7id <- mxCheckIdentification(c7)
# Error, says it needs a data set
# Checking ID for models with constraints now requires a data set
# Desired behavior?
omxCheckTrue(c7id$status)


#------------------------------------------------------------------------------
# Check 8 Multiple group (normal) model that uses mxConstraint to set
#     variance to b0 + b1*sex

data(example2)

mzd <- example2[example2$Zygosity %in% 'MZ', ]
dzd <- example2[example2$Zygosity %in% 'DZ', ]

c8mz <- mxModel('MZ', mxData(mzd, 'raw'),
                mxMatrix('Full', 1, 2, values=c(1, -1), free=TRUE,
                         labels=c('b0m', 'b1m'), name='BMZ'),
                mxMatrix('Full', 1, 1, values=.5, free=TRUE, labels='var',
                         name='V'),
                mxAlgebra(b0m + b1m*(data.TwinNum - 1), name='M'),
                mxMatrix('Full', 1, 1, free=TRUE, labels='mf', name='MF'),
                mxConstraint(M == MF),
                mxExpectationNormal(covariance='V', means='MF', dimnames='X'),
                mxFitFunctionML())
c8dz <- mxModel('DZ', mxData(dzd, 'raw'),
                mxMatrix('Full', 1, 2, values=c(1, -1), free=TRUE,
                         labels=c('b0m', 'b1m'), name='BDZ'),
                mxMatrix('Full', 1, 1, values=.5, free=TRUE, labels='var',
                         name='V'),
                mxAlgebra(b0m + b1m*(data.TwinNum - 1), name='M'),
                mxMatrix('Full', 1, 1, free=TRUE, labels='df', name='DF'),
                mxConstraint(M == DF),
                mxExpectationNormal(covariance='V', means='M', dimnames='X'),
                mxFitFunctionML())
c8 <- mxModel('mg', c8mz, c8dz, mxFitFunctionMultigroup(c('MZ', 'DZ')))
c8r <- mxRun(c8)
summary(c8r)
# Hmm ... it looks like constraints with definition variables might not be honored

#c8id <- mxCheckIdentification(c8)
# Error: foreign function call
#omxCheckTrue(c8id$status)


#------------------------------------------------------------------------------
# End
