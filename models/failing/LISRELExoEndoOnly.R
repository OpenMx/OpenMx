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
# Author: Michael D. Hunter
# Date: 2011.04.09
# Filename: LISRELExoEndoOnly.R
# Purpose: Create a test for the mxLISRELObjective function using only
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
	mxLISRELObjective(
		LY=ly@name,
		BE=be@name,
		PS=ps@name,
		TE=te@name,
		TY=ty@name,
		AL=al@name
	)
)


#--------------------------------------------------------------------
# Run the Endogenous-only model


# Uncomment the following lines when debugbing
#ymodRun <- mxRun(ymod, onlyFrontend=TRUE) # This runs fine.
ymod <- mxOption(ymod, "Calculate Hessian", "No")
ymod <- mxOption(ymod, "Standard Errors", "No")
ymod <- mxOption(ymod, "Major iterations", 1)


ymodRun <- mxRun(ymod)
summary(ymodRun)



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Define the model

xmod <- mxModel(
	name='LISREL Exogenous Model with Means',
	mxData(observed=rawlisx, type='raw'),
	lx, ph, td, tx, ka,
	mxLISRELObjective(
		LX=lx@name,
		PH=ph@name,
		TD=td@name,
		TX=tx@name,
		KA=ka@name
	)
)


#--------------------------------------------------------------------
# Run the Exogenous-only model


# Uncomment the following lines when debugbing
#xmodRun <- mxRun(xmod, onlyFrontend=TRUE) # This runs fine.
xmod <- mxOption(xmod, "Calculate Hessian", "No")
xmod <- mxOption(xmod, "Standard Errors", "No")
xmod <- mxOption(xmod, "Major iterations", 1)


xmodRun <- mxRun(xmod)
summary(xmodRun)



#--------------------------------------------------------------------
# Compare the estimate parameters to known values

# Recal that The known parameters are taken from model 3 (AKA 
#  threeLatentMultipleReg1 and threeLatentMultipleReg1Out) in
#  models/passing/IntroSEM-ThreeLatentMultipleRegTest1.R.



expectedParam <- c(
 0.809021, 1.148334, 1.305356, 0.805289, 1.232038, 1.189698, 
 0.87486, 0.786604, 0.701819, 0.990778, 0.456683, 1.92904, 
 0.153399, 1.316878, 0.770605, 1.133502, 1.067425, 1.063866,
 1.053663, 0.847609, 0.761789, 1.186354, 1.009936, 0.942481, 
 0.971017, 0.880311, 0.956065,  0.060813, 0.037374, -0.048679,
-0.013194, 0.193348, 0.220002, 0.256785, 0.171892, 0.166191, 
 0.23483, 0.173025, 0.157615
)



#omxCheckCloseEnough(lmodRun@output$estimate, expectedParam, epsilon=0.001)



#--------------------------------------------------------------------
# End

