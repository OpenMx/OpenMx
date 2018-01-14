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
# Date: 2014.12.17
# Filename: ModelIdentification.R
# Purpose: Check the model identification checking function.
#------------------------------------------------------------------------------




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
	free=c(T,T,T,T,F,F,F,F,F,F,F,F,F,T,T,T),
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
#--------------------------------------------------------------------
# Model identification

id <- mxCheckIdentification(xmod)

omxCheckEquals(id$status, FALSE)
omxCheckEquals(id$non_identified_parameters, c("l11", "l12", "l13", "l14", "varF1", "covF1F2"))

#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Add constraint, now I don't know what to do
xmod2 <- mxModel(xmod, mxConstraint(TX[1,1] == TX[2,1], name='conman'))

omxCheckError(mxCheckIdentification(xmod2), "Whoa Nelly.  I found an MxConstraint in your model.  I just cannot work under these conditions. I will be in my trailer until you reparameterize your model without using mxConstraint().")

