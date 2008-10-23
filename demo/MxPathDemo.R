
#myData <- read.csv("myFile")
#myManifest <- names(myData)

myManifest <- sprintf("%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")

ex2Model <- mxModel()

for (i in 1:5) {
    j <- i*20
    ex2Model <- mxModel(ex2Model, latentVars = myLatent, manifestVars = myManifest, 
                     mxPath(from=myLatent[i], 
                            to=myManifest[(j-19):j], 
                            arrows=1,
                            free=c(FALSE,rep(TRUE, 19)), 
                            startVal=c(1,rep(0.75,19))))
}

ex2Model <- mxModel(ex2Model, 
                 mxPath(from=myLatent,
                 		all=TRUE,
                        arrows=2,
                        free=TRUE, 
                        startVal=1))
                             
ex3Model <- mxModel(ex2Model,
           mxPath(from=myLatent, all=TRUE, free=FALSE, startVal=0),
           mxPath(from=myLatent, free=TRUE, startVal=1),
           mxPath(from="G6", to=myLatent, 
                  all=TRUE, free=TRUE, startVal=.75),
           mxPath(from="G6", free=FALSE, startVal=1),
           latentVars = "G6")

print(ex3Model)


