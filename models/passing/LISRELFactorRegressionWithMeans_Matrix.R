#
#   Copyright 2007-2013 The OpenMx Project
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
# Date: 2011.03.22
# Filename: LISRELFactorRegressionWithMeans_Matrix.R
# Purpose: Create a test for the mxExpectationLISREL function.  This test includes
#  a model for the means.  The known parameters are taken from model 3 (AKA 
#  threeLatentMultipleReg1 and threeLatentMultipleReg1Out) in
#  models/passing/IntroSEM-ThreeLatentMultipleRegTest1.R.  I re-estimated the
#  parameters of that model with ML instead of FIML to get the comparison
#  estimates.
#------------------------------------------------------------------------------

# Revision History:
# Thu Mar 22 21:21:39 Central Daylight Time 2012 -- Michael Hunter created complete file
# Sat Mar 24 01:24:23 Central Daylight Time 2012 -- Michael Hunter added dimnames to all matrices
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
# Define the model

lmod <- mxModel(
	name='LISREL Factor Regression Model with Means',
	mxData(observed=covlisdat, type='cov', means=mealisdat, numObs=nrow(rawlisdat)),
	#mxData(observed=rawlisdat, type='raw'),
	lx, ly, be, ga, ph, ps, td, te, th, tx, ty, ka, al,
	imxExpectationLISREL(LX=lx@name, LY=ly@name, BE=be@name,
		GA=ga@name, PH=ph@name, PS=ps@name, TD=td@name,
		TE=te@name, TH=th@name, TX=tx@name, TY=ty@name,
		KA=ka@name, AL=al@name),
	mxFitFunctionML()
)


#--------------------------------------------------------------------
# Run the model


# Uncomment the following lines when debugbing
#lmodRun <- mxRun(lmod, onlyFrontend=TRUE) # This runs fine.
#lmod <- mxOption(lmod, "Calculate Hessian", "No")
#lmod <- mxOption(lmod, "Standard Errors", "No")
#lmod <- mxOption(lmod, "Major iterations", 0)


lmodRun <- mxRun(lmod)
summary(lmodRun)


#--------------------------------------------------------------------
# Compare the estimate parameters to known values

# Recal that The known parameters are taken from model 3 (AKA 
#  threeLatentMultipleReg1 and threeLatentMultipleReg1Out) in
#  models/passing/IntroSEM-ThreeLatentMultipleRegTest1.R.  I re-estimated the
#  parameters of that model with ML instead of FIML to get the comparison
#  estimates.


expectedParam <- c(
 0.80902164,  1.14833409,  1.30535700,  0.80529009,  1.23203934,  1.18969961, 
 0.87485965,  0.78660379,  0.70181917,  0.99077779,  0.45668438,  1.93873470, 
 0.15417010,  1.32349163,  0.77447746,  1.13919808,  1.07278935,  1.06921199, 
 1.05895855,  0.85186877,  0.76561617,  1.19231559,  1.01501106,  0.94721747, 
 0.97589600,  0.88473398,  0.96086846,  0.06081483,  0.03737652, -0.04867645, 
-0.01319064,  0.19334869,  0.22000314,  0.25678636,  0.17189360,  0.16619428, 
 0.23483328,  0.17302787,  0.15761761 
)


omxCheckCloseEnough(lmodRun@output$estimate, expectedParam, epsilon=0.001)


#--------------------------------------------------------------------
# End

