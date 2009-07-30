require(OpenMx)

#Prepare Data
twinData <- read.table("ozbmi.data", header=T, na.strings = ".")
twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
#dimnames(twinData) <- list(NULL, twinVars)
summary(twinData)
selVars <- c('bmi1','bmi2')
mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
cov(mzfData)
cov(dzfData)

#Fit ACE Model with RawData and Matrices Input
twinACEModel <- mxModel("twinACE", 
    mxMatrix("Full", 1, 2, T, c(0,0), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanMZ"), 
    mxMatrix("Full", 1, 2, T, c(0,0), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanDZ"), 
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, name="X"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, name="Y"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, name="Z"),
    mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
    mxAlgebra(X * t(X), name="A"),
    mxAlgebra(Y * t(Y), name="C"),
    mxAlgebra(Z * t(Z), name="E"), 
    mxAlgebra(rbind (cbind(A + C + E, A + C),
                     cbind(A + C    , A + C + E)), 
                     dimnames = list(selVars, selVars), 
                     name="aceCovMZ"),
    mxAlgebra(rbind (cbind(A + C + E  , h %x% A + C),
                     cbind(h %x% A + C, A + C + E)), 
                     dimnames = list(selVars, selVars),
                     name="aceCovDZ"),
    mxModel("MZ",
        mxData(mzfData, type="raw"), 
        mxFIMLObjective("twinACE.aceCovMZ", "twinACE.expMeanMZ")),
    mxModel("DZ", 
        mxData(dzfData, type="raw"), 
        mxFIMLObjective("twinACE.aceCovDZ", "twinACE.expMeanDZ")),
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)
