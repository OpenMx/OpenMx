#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2011.09.22
# Filename: LISRELObjectiveTest_Matrix.R
# Purpose: Create a test for the mxLISRELObjective function.  Currently,
#  (Thu Sep 22 20:24:04 Central Daylight Time 2011) this only tests the
#  front end because there is no back end yet.
#------------------------------------------------------------------------------

# Revision History:
#   Thu Sep 22 20:56:35 Central Daylight Time 2011 -- Michael Hunter created file
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


#------------------------------------------------------------------------------
# Specify the model

require(OpenMx)

numManExo <- length(exoInd)
numManEnd <- length(endInd)
numLatExo <- 1
numLatEnd <- 2


lx <-     mxMatrix(
        name='LX',
        nrow=numManExo,
        ncol=numLatExo,
        free=T,
        values=c(.8, .8, .8),
        labels=c('lam1', 'lam2', 'lam3')
        )
ly <-     mxMatrix(
        name='LY',
        type='Iden',
        free=F,
        nrow=numManEnd,
        ncol=numLatEnd
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
    name='LISREL Test from 1982 Paper with mxLISRELObjective',
    mxData(observed=gObsCor, type='cor', numObs=1000),
    lx, ly, be, ga, ph, ps, td, te, th,
    mxLISRELObjective(LX=lx@name, LY=ly@name, BE=be@name, GA=ga@name, PH=ph@name, PS=ps@name, TD=td@name, TE=te@name, TH=th@name)
)



#------------------------------------------------------------------------------
# Fit the model
#ex1Run <- mxRun(Jor82Ex1, onlyFrontend=TRUE) # This runs fine.
Jor82Ex1 <- mxOption(Jor82Ex1, "Calculate Hessian", "No")
Jor82Ex1 <- mxOption(Jor82Ex1, "Standard Errors", "No")
Jor82Ex1 <- mxOption(Jor82Ex1, "Major iterations", 0)
ex1Run <- mxRun(Jor82Ex1)

summary(ex1Run)



