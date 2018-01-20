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
# Date: 2011.04.09
# Filename: LISRELExoEndoOnly.R
# Purpose: Create a test for the mxExpectationLISREL function using only
#  exogenous or only endogenous variables.  This test was created based on
#  models/passing/LISRELFactorRegressionWithMeans_Matrix*.R.
#------------------------------------------------------------------------------

# Revision History:
# Mon Apr 09 19:16:11 Central Daylight Time 2012 -- Michael Hunter created test from LISRELFactorRegressionWithMeans_Matrix*.R
# 


#--------------------------------------------------------------------
# Load OpenMx

require(OpenMx)


#--------------------------------------------------------------------
# Read in and set up the data

IndManExo <- 1:8
IndManEnd <- 9:12

# The data
data(latentMultipleRegExample1)

# Rearange Columns to separate exogenous and endogenous variables
rawlisdat <- latentMultipleRegExample1[, c(IndManEnd, IndManExo)]
rawlisy <- latentMultipleRegExample1[, IndManEnd]
rawlisx <- latentMultipleRegExample1[, IndManExo]


# Take covariance and means
covlisdat <- cov(rawlisdat)
mealisdat <- colMeans(rawlisdat)


# Number of manifest and latent exogenous and endogenous variables
numLatExo <- 2
numLatEnd <- 1
numManExo <- 8
numManEnd <- 4


# Dimnames
LatExo <- paste('xi', 1:numLatExo, sep='')
LatEnd <- paste('eta', 1:numLatEnd, sep='')
ManExo <- names(rawlisdat)[(numManEnd+1):(numManEnd+numManExo)]
ManEnd <- names(rawlisdat)[1:numManEnd]

#--------------------------------------------------------------------
# Specify the 13 extended LISREL matrices


lx <- mxMatrix("Full", numManExo, numLatExo,
	free=c(F,T,T,T,F,F,F,F,F,F,F,F,F,T,T,T),
	values=c(1, .2, .2, .2, 0, 0, 0, 0, 0, 0, 0, 0, 1, .2, .2, .2),
	labels=c( paste('l', 1, 1:4, sep=''), rep(NA, 8),  paste('l', 2, 5:8, sep='')),
	name='LX',
	dimnames=list(ManExo, LatExo)
) #DONE

ly <- mxMatrix("Full", numManEnd, numLatEnd,
	free=c(F,T,T,T),
	values=c(1, .2, .2, .2),
	labels= paste('l', 3, 9:12, sep=''),
	name='LY',
	dimnames=list(ManEnd, LatEnd)
) #DONE

be <- mxMatrix("Zero", numLatEnd, numLatEnd, name='BE', dimnames=list(LatEnd, LatEnd)) #DONE

ga <- mxMatrix("Full", numLatEnd, numLatExo,
	free=T,
	values=.2,
	labels=c('b13', 'b23'),
	name='GA',
	dimnames=list(LatEnd, LatExo)
) #DONE

ph <- mxMatrix("Symm", numLatExo, numLatExo,
	free=c(T,T,T),
	values=c(.8, .3, .8),
	labels=c('varF1', 'covF1F2', 'varF2'),
	name='PH',
	dimnames=list(LatExo, LatExo)
) #DONE

ps <- mxMatrix("Symm", numLatEnd, numLatEnd,
	free=T,
	values=.8,
	labels='varF3',
	name='PS',
	dimnames=list(LatEnd, LatEnd)
) #DONE

td <- mxMatrix("Diag", numManExo, numManExo,
	free=T,
	values=.8,
	labels=paste('d', 1:8, sep=''),
	name='TD',
	dimnames=list(ManExo, ManExo)
) #DONE

te <- mxMatrix("Diag", numManEnd, numManEnd,
	free=T,
	values=.8,
	labels=paste('e', 9:12, sep=''),
	name='TE',
	dimnames=list(ManEnd, ManEnd)
) #DONE

th <- mxMatrix("Zero", numManExo, numManEnd, name='TH', dimnames=list(ManExo, ManEnd)) #DONE

tx <- mxMatrix("Full", numManExo, 1,
	free=T,
	values=.1,
	labels=paste('m', 1:8, sep=''),
	name='TX',
	dimnames=list(ManExo, "TXMeans")
) #DONE

ty <- mxMatrix("Full", numManEnd, 1,
	free=T,
	values=.1,
	labels=paste('m', 9:12, sep=''),
	name='TY',
	dimnames=list(ManEnd, "TYMeans")
) #DONE

ka <- mxMatrix("Zero", numLatExo, 1, name='KA', dimnames=list(LatExo, "KAMeans")) #DONE

al <- mxMatrix("Zero", numLatEnd, 1, name='AL', dimnames=list(LatEnd, "ALMeans")) #DONE



#--------------------------------------------------------------------
# Define Endogenous-only the model

ymod <- mxModel(
	name='LISREL Endogenous Model with Means',
	mxData(observed=rawlisy, type='raw'),
	ly, be, ps, te, ty, al,
	mxExpectationLISREL(
		LY=ly$name,
		BE=be$name,
		PS=ps$name,
		TE=te$name,
		TY=ty$name,
		AL=al$name
	),
	mxFitFunctionML()
)


#--------------------------------------------------------------------
# Run the Endogenous-only model


# Uncomment the following lines when debugbing
#ymodRun <- mxRun(ymod, onlyFrontend=TRUE) # This runs fine.
#ymod <- mxOption(ymod, "Calculate Hessian", "No")
#ymod <- mxOption(ymod, "Standard Errors", "No")
#ymod <- mxOption(ymod, "Major iterations", 1)


ymodRun <- mxRun(ymod)
summary(ymodRun)



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Define the model

xmod <- mxModel(
	name='LISREL Exogenous Model with Means',
	mxData(observed=rawlisx, type='raw'),
	lx, ph, td, tx, ka,
	mxExpectationLISREL(
		LX=lx$name,
		PH=ph$name,
		TD=td$name,
		TX=tx$name,
		KA=ka$name
	),
	mxFitFunctionML()
)


#--------------------------------------------------------------------
# Run the Exogenous-only model


# Uncomment the following lines when debugbing
#xmodRun <- mxRun(xmod, onlyFrontend=TRUE) # This runs fine.
#xmod <- mxOption(xmod, "Calculate Hessian", "No")
#xmod <- mxOption(xmod, "Standard Errors", "No")
#xmod <- mxOption(xmod, "Major iterations", 1)


xmodRun <- mxRun(xmod)
summary(xmodRun)





#--------------------------------------------------------------------
# Create RAM models that mirror the LISREL ones

xman <- names(rawlisx)

xrmod <- mxModel(
	mxData(observed=rawlisx, type='raw'),
	name="RAM Exogenous",
	type="RAM",
	manifestVars=xman,
	latentVars=c("f1", "f2"),
	mxPath("f1", xman[1:4], free=c(F, T, T, T), values=c(1, .2, .2, .2)),
	mxPath("f2", xman[5:8], free=c(F, T, T, T), values=c(1, .2, .2, .2)),
	mxPath(c("f1", "f2"), connect="unique.pairs", free=T, arrows=2, values=c(.8, .3, .8)),
	mxPath(xman, arrows=2, free=T, values=.8),
	mxPath("one", c(xman, "f1", "f2"), free=c(rep(T, 8), F, F), values=c(rep(.1, 8), 0, 0))
)

xrmodRun <- mxRun(xrmod)





yman <- names(rawlisy)

yrmod <- mxModel(
	mxData(observed=rawlisy, type='raw'),
	name="RAM Endogenous",
	type="RAM",
	manifestVars=yman,
	latentVars="f1",
	mxPath("f1", yman[1:4], free=c(F, T, T, T), values=c(1, .2, .2, .2)),
	mxPath("f1", free=T, arrows=2, values=.8),
	mxPath(yman, arrows=2, free=T, values=.8),
	mxPath("one", c(yman, "f1"), free=c(rep(T, 4), F), values=c(rep(.1, 4), 0))
)

yrmodRun <- mxRun(yrmod)
summary(yrmodRun)




#--------------------------------------------------------------------
# Compare the estimate parameters in LISREL and RAM

# Check the exogenous only model
omxCheckCloseEnough(xmodRun$output$estimate, xrmodRun$output$estimate[c(1:6, 15:17, 7:14, 18:25)], epsilon=0.001)


# Check the endogenoug only model
omxCheckCloseEnough(ymodRun$output$estimate, yrmodRun$output$estimate[c(1:3, 8, 4:7, 9:12)], epsilon=0.001)



#--------------------------------------------------------------------
#require(rbenchmark)


#benchmark(mxRun(xmod), mxRun(xrmod), replications=10)
# LISREL and RAM models here take about the same amount of time
# LISREL is acutally .9% faster

#benchmark(mxRun(ymod), mxRun(yrmod), replications=10)
# LISREL takes about 9% longer than RAM models here

# For the combined model, LISREL is about 10% slower.



#--------------------------------------------------------------------
# End

