require(OpenMx)

#Prepare Data
twinData <- read.table("myTwinData.txt", header=T, na.strings=".")
twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
#dimnames(twinData) <- list(NULL, twinVars)
summary(twinData)
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))
colMeans(mzfData,na.rm=TRUE)
colMeans(dzfData,na.rm=TRUE)
cov(mzfData,use="complete")
cov(dzfData,use="complete")

#Fit ACE Model with RawData and Path-Style Input
share <- mxModel("share", 
    type="RAM",
    manifestVars=selVars,
    latentVars=aceVars,
    mxPath(from=aceVars, arrows=2, free=FALSE, values=1),
    mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0),
    mxPath(from="one", to=selVars, arrows=1, free=TRUE, values=20, labels= c("mean","mean")),
    mxPath(from=c("A1","C1","E1"), to="bmi1", arrows=1, free=TRUE, values=.6, label=c("a","c","e")),
    mxPath(from=c("A2","C2","E2"), to="bmi2", arrows=1, free=TRUE, values=.6, label=c("a","c","e")),
    mxPath(from="C1", to="C2", arrows=2, free=FALSE, values=1)
    )
twinACEModel <- mxModel("twinACE", 
    mxModel(share,
        mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1),
        mxData(mzfData, type="raw"), 
        type="RAM", name="MZ"),
    mxModel(share, 
        mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=.5),
        mxData(dzfData, type="raw"), 
        type="RAM", name="DZ"),
    mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
    mxAlgebraObjective("twin"))

#Run ACE model
twinACEFit <- mxRun(twinACEModel)

MZc <- mxEval(MZ.covariance, twinACEFit)
DZc <- mxEval(DZ.covariance, twinACEFit)
M <- mxEval(MZ.means, twinACEFit)
A <- mxEval(a*a, twinACEFit)
C <- mxEval(c*c, twinACEFit)
E <- mxEval(e*e, twinACEFit)
V <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
LL_ACE <- mxEval(objective, twinACEFit)


#Mx answers hard-coded
#1: Heterogeneity Model
Mx.A <- 0.6173023
Mx.C <- 5.595822e-14
Mx.E <- 0.1730462
Mx.M <- matrix(c(21.39293, 21.39293),1,2)
Mx.LL_ACE <- 4067.663


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LL_ACE,Mx.LL_ACE,.001)
omxCheckCloseEnough(A,Mx.A,.001)
omxCheckCloseEnough(C,Mx.C,.001)
omxCheckCloseEnough(E,Mx.E,.001)
omxCheckCloseEnough(M,t(Mx.M),.001)


#Run AE model
twinAEModel <- mxModel(twinACEModel, type="RAM",
    manifestVars=selVars,
    latentVars=aceVars,
    mxPath(from=c("A1","C1","E1"), to="bmi1", arrows=1, free=c(T,F,T), values=c(.6,0,.6), label=c("a","cfixed","e")),
    mxPath(from=c("A2","C2","E2"), to="bmi2", arrows=1, free=c(T,F,T), values=c(.6,0,.6), label=c("a","cfixed","e"))
    )
twinAEFit <- mxRun(twinAEModel)

MZc <- mxEval(MZ.covariance, twinAEFit)
DZc <- mxEval(DZ.covariance, twinAEFit)
M <- mxEval(MZ.means, twinAEFit)
A <- mxEval(a*a, twinAEFit)
C <- mxEval(c*c, twinAEFit)
E <- mxEval(e*e, twinAEFit)
V <- (A + C + E)
a2 <- A / V
c2 <- C / V
e2 <- E / V
AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(objective, twinAEFit)

LRT_ACE_AE <- LL_AE - LL_ACE

#Print relevant output
ACEest
AEest
LRT_ACE_AE
