# SCRIPT: UnivariateTwinAnalysisModerated_Matrix.R
# Author: 
# History:  Sat 26 Sep 2009 14:07:23 BST
#    changed to use data instead of read.table (tb)
# OpenMx: http://www.openmx.virginia.com
##########################################
require(OpenMx)
#Prepare Data
data("myTwinData", package="OpenMx") # data are in the file "myTwinData.txt" in the package's data folder
myTwinData[myTwinData=="."]=NA # na.strings="."
twinData <- subset(myTwinData, !is.na(age))
summary(twinData)
selVars <- c('bmi1','bmi2')
mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2,age)))
dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2,age)))
cov(mzfData, use="complete.obs")
cov(dzfData, use="complete.obs")
summary(mzfData)

#Fit ACE Model with RawData and Matrices Input
twinACEModel <- mxModel("ACE", 
    mxMatrix("Full", nrow=1, ncol=1, free=T, values=0.3352, label="a", name="X"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=-.6, label="c", name="Y"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=-.6, label="e", name="Z"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.01, label="amod", name="Xmod"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.01, label="cmod", name="Ymod"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.01, label="emod", name="Zmod"),
    mxBounds(c("c", "e"), NA, 0),
    mxBounds(c("a"), 0, NA),

    mxModel("MZ",
        mxData(mzfData, type="raw"), 
        mxMatrix("Full", nrow=1, ncol=1, free=F, label="data.age", name="Age"),
        mxAlgebra( ACE.X + Age * ACE.Xmod, name="XplusXmod"),
        mxAlgebra( ACE.Y + Age * ACE.Ymod, name="YplusYmod"),
        mxAlgebra( ACE.Z + Age * ACE.Zmod, name="ZplusZmod"),
        mxAlgebra(XplusXmod %*% t(XplusXmod), name="A"),
        mxAlgebra(YplusYmod %*% t(YplusYmod), name="C"),
        mxAlgebra(ZplusZmod %*% t(ZplusZmod), name="E"), 

        mxAlgebra(rbind (cbind(A + C + E, A + C),
                         cbind(A + C    , A + C + E)), 
                         dimnames = list(selVars, selVars), name="expCovMZ"),

        mxMatrix("Full", 1, 2, T, c(20,20), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanMZ"), 
        mxFIMLObjective("expCovMZ", "expMeanMZ")
    ),

    mxModel("DZ", 
        mxData(dzfData, type="raw"), 
        mxMatrix("Full", nrow=1, ncol=1, free=F, label="data.age", name="Age"),
        mxAlgebra( ACE.X + Age * ACE.Xmod, name="XplusXmod"),
        mxAlgebra( ACE.Y + Age * ACE.Ymod, name="YplusYmod"),
        mxAlgebra( ACE.Z + Age * ACE.Zmod, name="ZplusZmod"),
        mxAlgebra(XplusXmod %*% t(XplusXmod), name="A"),
        mxAlgebra(YplusYmod %*% t(YplusYmod), name="C"),
        mxAlgebra(ZplusZmod %*% t(ZplusZmod), name="E"), 
        mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
        
        mxAlgebra(rbind (cbind(A + C + E  , h %x% A + C),
                         cbind(h %x% A + C, A + C + E)), dimnames = list(selVars, selVars), name="expCovDZ"),
        
        mxMatrix("Full", 1, 2, T, c(20,20), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanDZ"), 
        mxFIMLObjective("expCovDZ", "expMeanDZ")),
 
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin")
)

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

#Check results against hard-coded Mx results
MZc    <- mxEval(expCovMZ,  twinACEFit$MZ)
DZc    <- mxEval(expCovDZ,  twinACEFit$DZ)
M      <- mxEval(expMeanMZ, twinACEFit$MZ)
X      <- mxEval(ACE.X, twinACEFit)
Y      <- mxEval(ACE.Y, twinACEFit)
Z      <- mxEval(ACE.Z, twinACEFit)
Xmod   <- mxEval(ACE.Xmod,  twinACEFit)
Ymod   <- mxEval(ACE.Ymod,  twinACEFit)
Zmod   <- mxEval(ACE.Zmod,  twinACEFit)
LL_ACE <- mxEval(objective, twinACEFit)


#Run Mx scripts with omxOriginalMx (commented out for now)
#mymatrices1 <- omxOriginalMx("mx-scripts/UnivariateTwinAnalysisModerated.mx", "temp-files")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
#Mx.A <- mymatrices1$A1.1
#Mx.C <- mymatrices1$C1.1
#Mx.E <- mymatrices1$E1.1
#Mx.M <- mymatrices1$M2.1
#Mx.LL_ACE <- mymatrices1$F4.1;

#Mx answers hard-coded
#1: Heterogeneity Model
Mx.X    <-  0.3352
Mx.Y    <- -0.9489
Mx.Z    <- -0.1013
Mx.Xmod <-  0.0174
Mx.Ymod <-  0.0332
Mx.Zmod <- -0.0131
Mx.M <- matrix(c(21.3810, 21.3810),1,2)
Mx.LL_ACE <- 4013.616


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(X,Mx.X,.001)
omxCheckCloseEnough(Y,Mx.Y,.001)
omxCheckCloseEnough(Z,Mx.Z,.001)
omxCheckCloseEnough(Xmod,Mx.Xmod,.001)
omxCheckCloseEnough(Ymod,Mx.Ymod,.001)
omxCheckCloseEnough(Zmod,Mx.Zmod,.001)
omxCheckCloseEnough(M,Mx.M,.001)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)

#
#
##Run AE model
#twinAEModel <- mxModel(twinACEModel, 
#    mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y")
#    )
#twinAEFit <- mxRun(twinAEModel)
#
#MZc <- mxEval(expCovMZ, twinAEFit)
#DZc <- mxEval(expCovDZ, twinAEFit)
#A <- mxEval(A, twinAEFit)
#C <- mxEval(C, twinAEFit)
#E <- mxEval(E, twinAEFit)
#V <- (A + C + E)
#a2 <- A / V
#c2 <- C / V
#e2 <- E / V
#AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
#LL_AE <- mxEval(objective, twinAEFit)
#
#LRT_ACE_AE <- LL_AE - LL_ACE
#
##Print relevant output
#ACEest
#AEest
#LRT_ACE_AE
#