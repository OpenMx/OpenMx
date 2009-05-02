#
#   Copyright 2007-2009 The OpenMx Project
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

# ----------------------------------
# Read data from disk.

# myData <- read.csv("myFile")
# myManifest <- names(myData)

# ----------------------------------
# Define 100 indicators and five first level latent variables.

myManifest <- sprintf("%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")

ex2Model <- mxModel(type = "RAM", latentVars = myLatent, manifestVars = myManifest)

# ----------------------------------
# Attach five latent variables to the 100 indicators in simple factor structure.
paths <- list()
for (i in 1:5) {
    j <- i*20
    paths <- c(paths, mxPath(from = myLatent[i], 
                            to = myManifest[(j - 19) : j], 
                            arrows = 1,
                            free = c(FALSE, rep(TRUE, 19)), 
                            start = c(1, rep(0.75, 19))))
}
ex2Model <- mxModel(ex2Model, paths)

# ----------------------------------
# Allow the factors to freely covary.

ex2Model <- mxModel(ex2Model, 
                 mxPath(from = myLatent,
                 		all = TRUE,
                        arrows = 2,
                        free = TRUE, 
                        start = 1))
                        
# ----------------------------------

print(ex2Model)

# ----------------------------------
# Create a second model from the first model.
# Remove the free factor covariances and create a second order common factor.

ex3Model <- mxModel(ex2Model, latentVars = "G6",
			mxPath(from=myLatent, all=TRUE, free=FALSE, start=0), # Remove free covs
			mxPath(from=myLatent, free=TRUE, start=1), # Add back in free vars
			mxPath(from="G6", to=myLatent, 
				all=TRUE, free=TRUE, start=.75, arrows=1), # Add 2nd order factor
			mxPath(from="G6", free=FALSE, start=1, arrows=2) # Add free var
           )

# ----------------------------------

print(ex3Model)


