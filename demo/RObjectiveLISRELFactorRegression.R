#
#   Copyright 2007-2012 The OpenMx Project
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
# Author: Michael D Hunter
# Filename: RObjectiveLISRELFactorRegression.R
# Date: 2011.06.30

# Purpose: To demonstrate an example of LISREL model
#  specification in OpenMx as an mxRObjective.
#  Later there will be a model of type LISREL.

# Revision History:
#   Michael Hunter -- 2011.06.30 Created File
#   Michael Hunter -- 2011.07.21 Added extensive comments, contributed to svn
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# The source of the model and the data presented here is
# Joreskog, K. G., & Sorbom, D.  (1982).  Recent developments in
#  structural equation modeling.  Journal of Marketing Research,
#  19(4), p. 404-416.
# The BibTeX is
#@article{ ,
#    Author = {Karl G {J\"{o}reskog} and Dag S\"{o}rbom},
#    Journal = {Journal of Marketing Research},
#    Month = {November},
#    Number = {4},
#    Pages = {404-416},
#    Title = {Recent Developments in Structural Equation Modeling},
#    Volume = {19},
#    Year = {1982}
#    }

# The data occurs on page 409 in Table 2
# The model is example 1, on pages 409-411.


#------------------------------------------------------------------------------
# Load the OpenMx package
library(OpenMx)


#------------------------------------------------------------------------------
# On LISREL
# 
# LISREL requires specifying 9 matrices.  These are
# LambdaX, LambdaY, ThetaDelta, ThetaEpsilon, ThetaDeltaEpsilon,
#  Beta, Gamma, Psi, and Phi
# 
# The extended LISREL model (a model including means) requires 
#  specifying 4 additional matrices, for a total of 13.  These are
#  Alpha, Kappa, TauX, and TauY.
# 
# Here is a brief description of the matrices
# LambdaX (abbreviated LX): matrix of exogenous factor loadings
# LambdaY (abbreviated LY): matrix of endogenous factor loadings
# ThetaDelta (abbreviated TD): residual covariance matrix of exogenous manifest variables (often diagonal)
# ThetaEpsilon (abbreviated TE): residual covariance matrix of endogenous manifest variables (often diagonal)
# ThetaDeltaEpsilon (abbreviated TH): residual covariance matrix of exogenous manifest variables with endogenous manifest variables (often zero)
# Beta (abbreviated BE): matrix of regression weights of the endogenous latent variables predicting themselves
# Gamma (abbreviated GA): matrix of regression weights of the exogenous latent variables predicting the endogenous latent variables
# Psi (abbreviated PS): residual covariance matrix of the endogenous latent variables
# Phi (abbreviated PH): covariance matrix of the exogenous latent variables (factor covariance matrix)
# Alpha (abbreviated AL): relates to means of endogenous latent variables
# Kappa (abbreviated KA): means of the exogenous latent variables
# TauX (abbreviated TX): residual means of the exogenous manifest variables
# TauY (abbreviated TY): residual means of the endogenous manifest variables


#------------------------------------------------------------------------------
# Define the function to perform the LISREL optimization


#--------------------------------------
# Author: Michael D Hunter
# Filename: funcLISREL.R
# Date: 2011.06.30
# Purpose: To define a number of helper functions to allow
#  OpenMx to more easily estimate models under the LISREL
#  specification.  This implementation uses an mxRObjective
#  function.  Later these functions will be translated and
#  moved to the C backend for an mxModel of type LISREL
#  as opposed to RAM.
#--------------------------------------


#--------------------------------------
# Take in the Lambda matrices of factor loadings
#  and put them together in blocks to form one
#  larger Lambda matrix

lisrelBlockLambda <- function(LambdaY, LambdaX){
    
    upper <- matrix(0, nrow=nrow(LambdaY), ncol=ncol(LambdaX))
    lower <- matrix(0, nrow=nrow(LambdaX), ncol=ncol(LambdaY))
    
    ret <- rbind(cbind(LambdaY, upper),
                 cbind(lower, LambdaX)
                )
    return(ret)
}


#--------------------------------------
# Take in the Theta matrices of measurement error covariances
#  and put them together in blocks to form one larger
#  Theta matrix

lisrelBlockTheta <- function(ThetaEpsilon, ThetaDelta, ThetaDeltaEpsilon){
    
    ret <- rbind(cbind(ThetaEpsilon, t(ThetaDeltaEpsilon)),
                 cbind(ThetaDeltaEpsilon, ThetaDelta)
            )
    return(ret)
}


#--------------------------------------
# Take in the regression matrices (Beta and Gamma),
#  the covariance matrix of the exogenous variables (Phi),
#  and the latent structural error covariance matrix (Psi)
#  then block them together to form a matrix I'm calling W

lisrelBlockW <- function(Beta, Gamma, Phi, Psi){
    
    I <- diag(1, dim(Beta)[1])
    A <- solve(I - Beta)
    
    AG <- A %*% Gamma
    
    block1 <- AG %*% Phi %*% t(AG) + A %*% Psi %*% t(A)
    block2 <- AG %*% Phi
    
    ret <- rbind(cbind(block1, block2),
                 cbind(t(block2), Phi)
                )
    return(ret)
}


#--------------------------------------
# Take all of the LISREL matrices and block them and
#  produce the model implied (aka expected) covariance
#  matrix (here called miCov)

lisrelCov <- function(LambdaX, LambdaY, ThetaEpsilon, ThetaDelta, ThetaDeltaEpsilon, Beta, Gamma, Phi, Psi){
    
    Lambda <- lisrelBlockLambda(LambdaY, LambdaX)
    Theta <- lisrelBlockTheta(ThetaEpsilon, ThetaDelta, ThetaDeltaEpsilon)
    W <- lisrelBlockW(Beta, Gamma, Phi, Psi)
    
    miCov <- Lambda %*% W %*% t(Lambda) + Theta
    
    return(miCov)
}


#--------------------------------------
# Take the model implied covariance matrix and the observed covariance
#  matrix, and produce a value that can be minimized to produce the
#  maximum likelihood estimates of the parameters
# Objective Function

lisrelML <- function(miCov, obsCov){
    ret <- suppressWarnings(log(det(miCov))) + tr(obsCov %*% solve(miCov))
    return(ret)
}


#--------------------------------------
# Rescale the output of lisrelML to official log likelihood units
# x is the output of lisrelML

lisrelMLRescale <- function(x, obsCov){
    ret <- x - suppressWarnings(log(det(obsCov))) - ncol(obsCov)
    return(ret)
}


#--------------------------------------
# Produce the Unweighted Least Squares (ULS) value that can be minimized
#  to make ULS estimates of the parameters
# Objective Function

lisrelULS <- function(miCov, obsCov){
    ret <- tr((obsCov-miCov)^2)/2
    return(ret)
}


#--------------------------------------
# Take in an mxModel object and produce a function value to be minimized
#  to produce maximum likelihood estimates of the free parameters.
# mxRObjective function

lisrelCovMx <- function(model, state){
    expectedCov <- lisrelCov(
        model$LambdaX@values,
        model$LambdaY@values,
        model$ThetaEpsilon@values,
        model$ThetaDelta@values,
        model$ThetaDeltaEpsilon@values,
        model$Beta@values,
        model$Gamma@values,
        model$Phi@values,
        model$Psi@values)
    ret <- lisrelML(expectedCov, model@data@observed)
    return(ret)
}


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


#------------------------------------------------------------------------------
# Specify the model

numManExo <- length(exoInd)
numManEnd <- length(endInd)
numLatExo <- 1
numLatEnd <- 2

Jor82Ex1 <- mxModel(
    name='LISREL 5 Model Example 1 from 1982 Paper',
    mxData(observed=gObsCor, type='cor', numObs=1000), # ADJUST HERE
    mxMatrix(
        name='LambdaY',
        type='Iden',
        free=F,
        nrow=numManEnd,
        ncol=numLatEnd
        ),
    mxMatrix(
        name='LambdaX',
        nrow=numManExo,
        ncol=numLatExo,
        free=T,
        values=c(.8, .8, .8),
        labels=c('lam1', 'lam2', 'lam3')
        ),
    mxMatrix(
        name='Beta',
        nrow=numLatEnd,
        ncol=numLatEnd,
        free=c(F, T, F, F),
        values=c(0, .8, 0, 0),
        labels=c(NA, 'bet', NA, NA)
        ),
    mxMatrix(
        name='Gamma',
        nrow=numLatEnd,
        ncol=numLatExo,
        free=T,
        values=c(.8, 1.8),
        labels=c('gam1', 'gam2')
        ),
    mxMatrix(
        name='Psi',
        nrow=numLatEnd,
        ncol=numLatEnd,
        free=c(T, F, F, T),
        values=c(.8, 0, 0, 1.8),
        labels=c('zet1sq', NA, NA, 'zet2sq')
        ),
    mxMatrix(
        name='Phi',
        nrow=numLatExo,
        ncol=numLatExo,
        free=F,
        values=1
        ),
    mxMatrix(
        name='ThetaEpsilon',
        nrow=numManEnd,
        ncol=numManEnd,
        free=F,
        values=0
        ),
    mxMatrix(
        name='ThetaDelta',
        type='Symm',
        nrow=numManExo,
        ncol=numManExo,
        free=c(T, F, T, F, T, F, T, F, T),
        values=c(.8, 0, .1, 0, .8, 0, .1, 0, .8),
        labels=c('del1', NA, 'del13', NA, 'del2', NA, 'del13', NA, 'del3')
        ),
    mxMatrix(
        name='ThetaDeltaEpsilon',
        nrow=numManExo,
        ncol=numManEnd,
        free=F,
        values=0
        ),
    mxRObjective(lisrelCovMx)
)

#------------------------------------------------------------------------------
# Fit the model
ex1Run <- mxRun(Jor82Ex1)

summary(ex1Run)



#------------------------------------------------------------------------------
# Fitted Results

# The paper (Joreskog & Sorbom, 1982) reports correlation sacross
#  variables and standard deviations within variables
#  for the residual covariance matrices.
# This takes a covariance matrix (the matrix estimated,
#  turns it into a correlation matrix with the square root of
#  the variances (i.e. the standard deviations) down the diagonal.
TDf <- cov2cor(ex1Run$ThetaDelta@values)
diag(TDf) <- sqrt(diag(ex1Run$ThetaDelta@values))

PSf <- cov2cor(ex1Run$Psi@values)
diag(PSf) <- sqrt(diag(ex1Run$Psi@values))

# The other matrices are reported raw.
LXf <- ex1Run$LambdaX@values
GAf <- ex1Run$Gamma@values
BEf <- ex1Run$Beta@values


#------------------------------------------------------------------------------
# Published Results

TDp <- matrix(c(0.64, 0, -0.31, 0, 0.52, 0, -0.31, 0, 0.65), nrow=3, ncol=3)
LXp <- matrix(c(0.76, 0.85, 0.76), nrow=3, ncol=1)
GAp <- matrix(c(0.63, 0.21), nrow=2, ncol=1)
BEp <- matrix(c(0, 0.70, 0, 0), nrow=2, ncol=2)
PSp <- diag(c(0.78, 0.53))

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
