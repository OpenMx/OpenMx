
#myData <- read.csv("myFile")
#myManifest <- names(myData)

myManifest <- as.character(c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")

ex2Model <- new("MxPathModel")

ex2Model@latentVars <- myLatent
ex2Model@manifestVars <- myManifest


for (i in 1:5) {
    j <- i*20
    ex2Model <- mxAddPath(ex2Model, 
                     mxCreatePath(from=myLatent[i], 
                                 to=myManifest[(j-19):j], 
                                 arrows=1,
                                 free=c(FALSE,rep(TRUE, 19)), 
                                 startVal=c(1,rep(0.75,19))))
}

ex2Model <- mxAddPath(ex2Model, 
                 mxCreatePath(from=myLatent,
                 			 all=TRUE,
                             arrows=2,
                             free=TRUE, 
                             startVal=1))
                             
ex3Model <- mxAddPath(ex2Model,
           mappend(mxCreatePath(from=myLatent, all=TRUE, free=FALSE, startVal=0),
                   mxCreatePath(from=myLatent, free=TRUE, startVal=1),
                   mxCreatePath(from="G6", to=myLatent, 
                      			all=TRUE, free=TRUE, startVal=.75),
                   mxCreatePath(from="G6", free=FALSE, startVal=1)))

ex3Model@latentVars <- append(ex3Model@latentVars, "G6")

print(ex3Model)


