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
# Date: 2011.09.22
# Filename: LISRELFactorRegression_Matrix.R
# Purpose: Create a test for the mxExpectationLISREL function.  The LISREL
#  objective is still under active development and probably has some memory
#  leaks, but this models estimates correctly and is much fast/easier than the
#  RObjective implementation: demo/RObjectiveLISRELFactorRegression.R
#------------------------------------------------------------------------------

# Revision History:
#   Thu Sep 22 20:56:35 Central Daylight Time 2011 -- Michael Hunter created file
#   Mon Mar 05 08:41:32 Central Standard Time 2012 -- Michael Hunter added comparison to published estimated parameters
#   Tue Mar 06 12:52:01 Central Standard Time 2012 -- Michael Hunter added OpenMx copyright and moved from models/failing to models/passing and renamed from LISRELObjectiveTest_Matrix.R to LISRELFactorRegression_Matrix.R
#   Fri Mar 30 03:10:06 Central Daylight Time 2012 -- Michael Hunter added dimnames to observed data
#   Sat Apr 07 20:06:13 Central Daylight Time 2012 -- Michael Hunter added dimnames to factor loading matrices
#


#------------------------------------------------------------------------------
# Define the data to be fit
obsCor <- matrix(c(
    1.00, 0.43, 0.52, 0.54, 0.83,
    0.43, 1.00, 0.67, 0.45, 0.45,
    0.52, 0.67, 1.00, 0.63, 0.56,
    0.54, 0.45, 0.63, 1.00, 0.52,
    0.83, 0.45, 0.56, 0.52, 1.00),
    nrow=5, ncol=5, byrow=TRUE
)

# Define the exogenous and endogenous variables
endInd <- c(1, 5) #Indexes of the endogenous/dependent variables
exoInd <- c(2:4) #Indexes of the exogenous/independent variables

# Correlation matrix rearranged into separate blocks of
# endogenous and exogenous variables
gObsCor <- cbind(rbind(obsCor[endInd, endInd], obsCor[exoInd, endInd]), rbind(obsCor[endInd, exoInd], obsCor[exoInd, exoInd]) )
rownames(gObsCor) <- c(paste('y', 1:length(endInd), sep=''), paste('x', 1:length(exoInd), sep='') )
colnames(gObsCor) <- rownames(gObsCor)



#------------------------------------------------------------------------------
# Specify the model

require(OpenMx)

numManExo <- length(exoInd)
numManEnd <- length(endInd)
numLatExo <- 1
numLatEnd <- 2

manExoNam <- paste('x', 1:numManExo, sep='')
manEndNam <- paste('y', 1:numManEnd, sep='')
latExoNam <- paste('xi', 1:numLatExo, sep='')
latEndNam <- paste('eta', 1:numLatEnd, sep='')


# Specify the nine LISREL matrices, no model for means
lx <-     mxMatrix(
        name='LX',
        nrow=numManExo,
        ncol=numLatExo,
        free=T,
        values=c(.8, .8, .8),
        labels=c('lam1', 'lam2', 'lam3'),
        dimnames=list(manExoNam, latExoNam)
        )
ly <-     mxMatrix(
        name='LY',
        type='Iden',
        free=F,
        nrow=numManEnd,
        ncol=numLatEnd,
        dimnames=list(manEndNam, latEndNam)
        )
be <-     mxMatrix(
        name='BE',
        nrow=numLatEnd,
        ncol=numLatEnd,
        free=c(F, T, F, F),
        values=c(0, .8, 0, 0),
        labels=c(NA, 'bet', NA, NA)
        )
ga <-     mxMatrix(
        name='GA',
        nrow=numLatEnd,
        ncol=numLatExo,
        free=T,
        values=c(.8, 1.8),
        labels=c('gam1', 'gam2')
        )
ps <-     mxMatrix(
        name='PS',
        nrow=numLatEnd,
        ncol=numLatEnd,
        free=c(T, F, F, T),
        values=c(.8, 0, 0, 1.8),
        labels=c('zet1sq', NA, NA, 'zet2sq')
        )
ph <-     mxMatrix(
        name='PH',
        nrow=numLatExo,
        ncol=numLatExo,
        free=F,
        values=1
        )
te <-    mxMatrix(
        name='TE',
        nrow=numManEnd,
        ncol=numManEnd,
        free=F,
        values=0
        )
td <-     mxMatrix(
        name='TD',
        type='Symm',
        nrow=numManExo,
        ncol=numManExo,
        free=c(T, F, T, F, T, F, T, F, T),
        values=c(.8, 0, .1, 0, .8, 0, .1, 0, .8),
        labels=c('del1', NA, 'del13', NA, 'del2', NA, 'del13', NA, 'del3')
        )
th <-     mxMatrix(
        name='TH',
        nrow=numManExo,
        ncol=numManEnd,
        free=F,
        values=0
        )


Jor82Ex1 <- mxModel(
    name='LISREL Test from 1982 Paper with mxExpectationLISREL',
    mxData(observed=gObsCor, type='cor', numObs=1000),
    lx, ly, be, ga, ph, ps, td, te, th,
    mxExpectationLISREL(LX=lx$name, LY=ly$name, BE=be$name, GA=ga$name, PH=ph$name, PS=ps$name, TD=td$name, TE=te$name, TH=th$name),
    mxFitFunctionML()
)



#------------------------------------------------------------------------------
# Fit the model

# Uncomment the following lines when debugbing
#ex1Run <- mxRun(Jor82Ex1, onlyFrontend=TRUE) # This runs fine.
#Jor82Ex1 <- mxOption(Jor82Ex1, "Calculate Hessian", "No")
#Jor82Ex1 <- mxOption(Jor82Ex1, "Standard Errors", "No")
#Jor82Ex1 <- mxOption(Jor82Ex1, "Major iterations", 0)

ex1Run <- omxCheckWarning(mxRun(Jor82Ex1),
                          paste("OpenMx does not yet correctly handle mxData(type='cor')",
                                'standard errors and fit statistics.',
                                'See Steiger (1980), "Tests for comparing elements of a correlation matrix".'))

summary(ex1Run)


#------------------------------------------------------------------------------
# Published Results

TDp <- matrix(c(0.64, 0, -0.31, 0, 0.52, 0, -0.31, 0, 0.65), nrow=3, ncol=3)
LXp <- matrix(c(0.76, 0.85, 0.76), nrow=3, ncol=1)
GAp <- matrix(c(0.63, 0.21), nrow=2, ncol=1)
BEp <- matrix(c(0, 0.70, 0, 0), nrow=2, ncol=2)
PSp <- diag(c(0.78, 0.53))


#------------------------------------------------------------------------------
# Fitted Results

# Covariance matrices were scaled as correlation matrices
#  with the standard deviations down the diagonal
PSf <- cov2cor(ex1Run$PS$values)
diag(PSf) <- sqrt(diag(ex1Run$PS$values))

TDf <- cov2cor(ex1Run$TD$values)
diag(TDf) <- sqrt(diag(ex1Run$TD$values))


# The other matrices are reported raw.
LXf <- ex1Run$LX$values
GAf <- ex1Run$GA$values
BEf <- ex1Run$BE$values

#------------------------------------------------------------------------------
# Compare Fitted and Published Results

# Note: Because the published paper only reports
#  two decimals of precision, the comparison can only
#  be made down to two decimal places (i.e. epsilon=0.01)

omxCheckCloseEnough(TDf, TDp, epsilon=0.01)
omxCheckCloseEnough(LXf, LXp, epsilon=0.01)
omxCheckCloseEnough(GAf, GAp, epsilon=0.01)
omxCheckCloseEnough(BEf, BEp, epsilon=0.01)
omxCheckCloseEnough(PSf, PSp, epsilon=0.01)

#------------------------------------------------------------------------------
# End
