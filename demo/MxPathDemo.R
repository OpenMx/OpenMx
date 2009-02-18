# ----------------------------------
# Read data from disk.

# myData <- read.csv("myFile")
# myManifest <- names(myData)

# ----------------------------------
# Define 100 indicators and five first level latent variables.

myManifest <- sprintf("%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")

ex2Model <- mxModel()

# ----------------------------------
# Attach five latent variables to the 100 indicators in simple factor structure.

for (i in 1:5) {
    j <- i*20
    ex2Model <- mxModel(ex2Model, latentVars = myLatent, manifestVars = myManifest, 
                     mxPath(from=myLatent[i], 
                            to=myManifest[(j-19):j], 
                            arrows=1,
                            free=c(FALSE,rep(TRUE, 19)), 
                            startVal=c(1,rep(0.75,19))))
}

# ----------------------------------
# Allow the factors to freely covary.

ex2Model <- mxModel(ex2Model, 
                 mxPath(from=myLatent,
                 		all=TRUE,
                        arrows=2,
                        free=TRUE, 
                        startVal=1))
                        
# ----------------------------------
# Convert the model to Open Mx (this step is needed?) and print it.

ex2Model <- omxConvertPathModel(ex2Model)

print(ex2Model)

# ----------------------------------
# Create a second model from the first model.
# Remove the free factor covariances and create a second order common factor.

ex3Model <- mxModel(ex2Model, latentVars = "G6",
           mxPath(from=myLatent, all=TRUE, free=FALSE, startVal=0), # Remove free covs
           mxPath(from=myLatent, free=TRUE, startVal=1), # Add back in free vars
           mxPath(from="G6", to=myLatent, 
                  all=TRUE, free=TRUE, startVal=.75, arrows=1), # Add 2nd order factor
           mxPath(from="G6", free=FALSE, startVal=1, arrows=2) # Add free var
           )

# ----------------------------------
# Convert the model to Open Mx (this step is needed?) and print it.

ex3Model <- omxConvertPathModel(ex3Model)

print(ex3Model)


