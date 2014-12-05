#
#   Copyright 2007-2014 The OpenMx Project
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

# The next TODO item is to start implementing addEntriesLISREL in the
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

