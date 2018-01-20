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
# Date: 2014-12-05
# Filename: LisrelTypeCheck.R
# Purpose: Check functionality of the LISREL type while implementing pieces of
#  it.  Currently, the only working function is to create an MxModel with
#  the correct zero matrices, having accurately processed the manifestVars and
#  latentVars arguments.
#------------------------------------------------------------------------------

# The next TODO item is to start implementing removeEntriesLISREL in the
#  MxLISRELModel.R file.

#------------------------------------------------------------------------------
# load needed package(s)
require(OpenMx)


#------------------------------------------------------------------------------


mod1 <- mxModel(name="A type LISREL model",
	manifestVars=list(exogenous='x1', endogenous='y1'),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL"
)

# make sure partial matching works
mod2 <- mxModel(name="partial matching LISREL model",
	manifestVars=list(exo='x1', end='y1'),
	latentVars=list(exo='ksi1', end='eta1'),
	type="LISREL"
)

# exogenous only
mod3 <- mxModel(name="exogenous only LISREL model",
	manifestVars=list(exo='x1'),
	latentVars=list(exo='ksi1'),
	type="LISREL"
)

# endogenous only
mod4 <- mxModel(name="endogenous LISREL model",
	manifestVars=list(end='y1'),
	latentVars=list(end='eta1'),
	type="LISREL"
)


#------------------------------------------------------------------------------
# Now try adding entries


data(demoTwoFactor)

# data
mod1a <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw')
)

# means
mod1b <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw'),
	mxPath(from="one", to='x1', values=1, free=TRUE),
	mxPath(from="one", to='y3', values=3, free=TRUE),
	mxPath(from="one", to=c('x2', 'y2'), values=2, free=TRUE, labels='TwoValue'),
	mxPath(from="one", to=c("ksi1", "eta1"), values=7, free=FALSE)
)

# latent exo paths
mod1c <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw'),
	mxPath(from='ksi1', to='x1', values=1, free=TRUE, labels='load1'), #LX
	mxPath(from='ksi1', to='eta1', values=2, free=TRUE, labels='reg1', ubound=20), #GA
	mxPath(from='ksi1', arrow=2, values=3, free=FALSE, labels='lvar1') #PH
)

# latent endo paths
mod1d <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw'),
	mxPath(from='eta1', to='y1', values=1, free=TRUE, labels='load1'), #LY
	mxPath(from='eta1', arrow=2, values=3, free=FALSE, labels='lvar1'), #PS
	mxPath(from='eta1', arrow=1, values=3, free=FALSE, labels='lvar1') #BE, N.B. bc 1 dimensional, BE path overwrites PS path.
)


# manifest exo paths
mod1e <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw'),
	mxPath(from='x1', to='x1', values=1, arrows=2, free=TRUE, labels='resid1'), #TD
	mxPath(from='x1', to='y3', values=2, arrows=2, free=TRUE, labels='residX1') #TH
)

# manifest enod paths
mod1f <- mxModel(name="A type LISREL model",
	manifestVars=list(exo=names(demoTwoFactor)[1:5], endo=names(demoTwoFactor)[6:10]),
	latentVars=list(exogenous='ksi1', endogenous='eta1'),
	type="LISREL",
	mxData(demoTwoFactor, 'raw'),
	mxPath(from='y4', to='y2', values=1, arrows=2, free=TRUE, labels='residC3'), #TE
	mxPath(from='y2', to='x5', values=2, arrows=2, free=TRUE, labels='residX2') #TH
)

#------------------------------------------------------------------------------

