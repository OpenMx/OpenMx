setwd("~/Applications/bin/OpenMx/trunk/demo/ExamplesH")
require(OpenMx)

#Prepare Data
twinData <- read.table("ozbmi.data", header=T, na.strings=".")
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
    mxMatrix("Full", 1, 2, T, c(20,20), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanMZ"), 
    mxMatrix("Full", 1, 2, T, c(20,20), labels= c("mean","mean"), dimnames=list(NULL, selVars), name="expMeanDZ"), 
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
    mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
    mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
    mxAlgebra(X * t(X), name="A"),
    mxAlgebra(Y * t(Y), name="C"),
    mxAlgebra(Z * t(Z), name="E"), 
    mxAlgebra(rbind (cbind(A + C + E, A + C),
                     cbind(A + C    , A + C + E)), dimnames = list(selVars, selVars), name="expCovMZ"),
    mxAlgebra(rbind (cbind(A + C + E  , h %x% A + C),
                     cbind(h %x% A + C, A + C + E)), dimnames = list(selVars, selVars), name="expCovDZ"),
    mxModel("MZ",
        mxData(mzfData, type="raw"), 
        mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMeanMZ")),
    mxModel("DZ", 
        mxData(dzfData, type="raw"), 
        mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMeanDZ")),
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEvaluate(expCovMZ, twinACEFit)
DZc <- mxEvaluate(expCovDZ, twinACEFit)
a <- mxEvaluate(A, twinACEFit)
c <- mxEvaluate(C, twinACEFit)
e <- mxEvaluate(E, twinACEFit)
v <- (a+c+e)
a2 <- a/v
c2 <- c/v
e2 <- e/v
ACEest <- rbind(cbind(a,c,e),cbind(a2,c2,e2))
LL_ACE <- mxEvaluate(objective, twinACEFit)

#Run Mx scripts
mymatrices1 <- omxOriginalMx("mx-scripts/univACE.mx", "temp-files")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
Mx.a <- mymatrices1$A1.1
Mx.c <- mymatrices1$C1.1
Mx.e <- mymatrices1$E1.1
Mx.M <- mymatrices1$M2.1
Mx.LL_ACE <- mymatrices1$F4.1;

#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(a,Mx.a,.001)
omxCheckCloseEnough(c,Mx.c,.001)
omxCheckCloseEnough(e,Mx.e,.001)
omxCheckCloseEnough(c(mean,mean),Mx.M,.001)


#Run AE model
twinAEModel <- mxModel(twinACEModel, 
    mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y")
    )
twinAEFit <- mxRun(twinAEModel)

MZc <- mxEvaluate(expCovMZ, twinAEFit)
DZc <- mxEvaluate(expCovDZ, twinAEFit)
a <- mxEvaluate(A, twinAEFit)
c <- mxEvaluate(C, twinAEFit)
e <- mxEvaluate(E, twinAEFit)
v <- (a+c+e)
a2 <- a/v
c2 <- c/v
e2 <- e/v
AEest <- rbind(cbind(a,c,e),cbind(a2,c2,e2))
LL_AE <- mxEvaluate(objective, twinAEFit)

LRT_ACE_AE <- LL_AE-LL_ACE

ACEest
AEest
LRT_ACE_AE
