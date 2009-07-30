require(OpenMx)

#Simulate Data
require(MASS)
#group 1
set.seed(200)
rs=.5
xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
set.seed(200)
#group 2
rs=.4
xy2 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))

#Print Descriptive Statistics
selVars <- c('X','Y')
summary(xy1)
cov(xy1)
dimnames(xy1) <- list(NULL, selVars)
summary(xy2)
cov(xy2)
dimnames(xy2) <- list(NULL, selVars)

#Fit Heterogeneity Model
bivHetModel <- mxModel("bivHet",
    mxModel("group1",
        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=2, 
            free=c(T,T,F,T), 
            values=c(1,.2,0,1),
            labels=c("vX1", "cXY1", "zero", "vY1"),
            dimnames=list(selVars, selVars), 
            name="Chol1"
        ), 
        mxAlgebra(
            Chol1 %*% t(Chol1), 
            name="EC1", 
            dimnames=list(selVars, selVars)
        ), 
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=2, 
            free=T, 
            values=c(0,0), 
            labels=c("mX1", "mY1"), 
            dimnames=list(NULL, selVars), 
            name="EM1"
        ), 
        mxData(
            xy1, 
            type="raw"
        ), 
        mxFIMLObjective(
            "EC1", 
            "EM1")
        ),
    mxModel("group2",
        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=2, 
            free=c(T,T,F,T), 
            values=c(1,.2,0,1),
            labels=c("vX2", "cXY2", "zero", "vY2"),
            dimnames=list(selVars, selVars), 
            name="Chol2"
        ), 
        mxAlgebra(
            Chol2 %*% t(Chol2), 
            name="EC2", 
            dimnames=list(selVars, selVars)
        ), 
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=2, 
            free=T, 
            values=c(0,0), 
            labels=c("mX2", "mY2"), 
            dimnames=list(NULL, selVars), 
            name="EM2"
        ), 
        mxData(
            xy2, 
            type="raw"
        ), 
        mxFIMLObjective(
            "EC2", 
            "EM2")
        ),
    mxAlgebra(
            group1.objective + group2.objective, 
            name="h12"
        ),
    mxAlgebraObjective("h12")
    )
bivHetFit <- mxRun(bivHetModel)
    EM1Het <- mxEvaluate(group1.EM1, bivHetFit)
    EM2Het <- mxEvaluate(group2.EM2, bivHetFit)
    EC1Het <- mxEvaluate(group1.EC1, bivHetFit)
    EC2Het <- mxEvaluate(group2.EC2, bivHetFit)
    LLHet <- mxEvaluate(objective, bivHetFit)

#Fit Homnogeneity Model
bivHomModel <- bivHetModel
    bivHomModel[['group2.Chol2']]@labels <- bivHomModel[['group1.Chol1']]@labels
    bivHomModel[['group2.EM2']]@labels <- bivHomModel[['group1.EM1']]@labels
bivHomFit <- mxRun(bivHomModel)
    EM1Hom <- mxEvaluate(group1.EM1, bivHomFit)
    EM2Hom <- mxEvaluate(group2.EM2, bivHomFit)
    EC1Hom <- mxEvaluate(group1.EC1, bivHomFit)
    EC2Hom <- mxEvaluate(group2.EC2, bivHomFit)
    LLHom <- mxEvaluate(objective, bivHomFit)

    Chi= LLHom-LLHet
    LRT= rbind(LLHet,LLHom,Chi)
    LRT


#Mx answers hard-coded
#1: Heterogeneity Model
Mx.EM1Het <- matrix(c(0.03211284, -0.004889846),1,2)
Mx.EC1Het <- matrix(c(1.0092856, 0.4813512, 0.4813512, 0.9935414),2,2)
Mx.EM2Het <- matrix(c(0.03341992, -0.007112054),1,2)
Mx.EC2Het <- matrix(c(1.012324, 0.3799160, 0.379916, 0.9956605),2,2)
Mx.LLHet <- 10944.873

#2: Homogeneity Model
Mx.EM1Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.EC1Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.EM2Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.EC2Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.LLHom <- 10954.368


#OpenMx summary
cov <- rbind(cbind(EC1Het,EC2Het),cbind(EC1Hom,EC2Hom))
mean <- rbind(cbind(EM1Het, EM2Het),cbind(EM1Hom,EM2Hom))
like <- rbind(cbind(LLHet),cbind(LLHom))
cov; mean; like

#old Mx summary
Mx.cov <- rbind(cbind(Mx.EC1Het,Mx.EC2Het),cbind(Mx.EC1Hom,Mx.EC2Hom))
Mx.mean <- rbind(cbind(Mx.EM1Het, Mx.EM2Het),cbind(Mx.EM1Hom,Mx.EM2Hom))
Mx.like <- rbind(cbind(Mx.LLHet),cbind(Mx.LLHom))
Mx.cov; Mx.mean; Mx.like


#Compare OpenMx results to Mx results (LL: likelihood; EC: expected covariance, EM: expected means)
omxCheckCloseEnough(LLHet,Mx.LLHet,.001)
omxCheckCloseEnough(EC1Het,Mx.EC1Het,.001)
omxCheckCloseEnough(EM1Het,Mx.EM1Het,.001)
omxCheckCloseEnough(EC2Het,Mx.EC2Het,.001)
omxCheckCloseEnough(EM2Het,Mx.EM2Het,.001)

omxCheckCloseEnough(LLHom,Mx.LLHom,.001)
omxCheckCloseEnough(EC1Hom,Mx.EC1Hom,.001)
omxCheckCloseEnough(EM1Hom,Mx.EM1Hom,.001)
omxCheckCloseEnough(EC2Hom,Mx.EC2Hom,.001)
omxCheckCloseEnough(EM2Hom,Mx.EM2Hom,.001)


