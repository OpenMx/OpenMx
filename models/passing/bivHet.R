setwd("~/Applications/bin/OpenMx/trunk/demo/ExamplesH")
require(OpenMx)

#Simulate Data
require(MASS)
set.seed(200); rs=.5; xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2)); group=1;
set.seed(200); rs=.4; xy2 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2)); group=2;
selVars <- c('X','Y')
summary(xy1); cov(xy1); dimnames(xy1) <- list(NULL, selVars)
summary(xy2); cov(xy2); dimnames(xy2) <- list(NULL, selVars)

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
    		name="Chol1"), 
    	mxAlgebra(
    		Chol1 %*% t(Chol1), 
    		name="EC1", 
    		dimnames=list(selVars, selVars)), 
    	mxMatrix(
    		type="Full", 
    		nrow=1, 
    		ncol=2, 
    		free=T, 
    		values=c(0,0), 
            labels=c("mX1", "mY1"), 
    		dimnames=list(NULL, selVars), 
    		name="EM1"), 
    	mxData(
    		xy1, 
    		type="raw"), 
    	mxFIMLObjective(
    		"EC1", 
    		"EM1")),
	mxModel("group2",
    	mxMatrix(
    		type="Full", 
    		nrow=2, 
    		ncol=2, 
    		free=c(T,T,F,T), 
    		values=c(1,.2,0,1),
    		labels=c("vX2", "cXY2", "zero", "vY2"),
    		dimnames=list(selVars, selVars), 
    		name="Chol2"), 
    	mxAlgebra(
    		Chol2 %*% t(Chol2), 
    		name="EC2", 
    		dimnames=list(selVars, selVars)), 
    	mxMatrix(
    		type="Full", 
    		nrow=1, 
    		ncol=2, 
    		free=T, 
    		values=c(0,0), 
            labels=c("mX2", "mY2"), 
    		dimnames=list(NULL, selVars), 
    		name="EM2"), 
    	mxData(
    		xy2, 
    		type="raw"), 
    	mxFIMLObjective(
    		"EC2", 
    		"EM2")),
	mxAlgebra(group1.objective + group2.objective, name="h12"),
	mxAlgebraObjective("h12"))
	
bivHetFit <- mxRun(bivHetModel)
EM1 <- bivHetFit[['EM1']]@values
EM2 <- bivHetFit[['EM2']]@values
EC1 <- mxEvaluate(Chol1 %*% t(Chol1),bivHetFit)
EC2 <- mxEvaluate(Chol2 %*% t(Chol2),bivHetFit)
LLHet <- mxEvaluate(objective,bivHetFit);

bivHomModel <- bivHetModel
bivHomModel[['group2.EC2']]@labels <- bivHomModel[['group1.EC1']]@labels
bivHomModel[['group2.EM2']]@labels <- bivHomModel[['group1.EM1']]@labels

bivHomFit <- mxRun(bivHomModel)
EM1 <- bivHomFit[['EM1']]@values
EM2 <- bivHomFit[['EM2']]@values
EC1 <- mxEvaluate(Chol1 %*% t(Chol1),bivHomFit)
EC2 <- mxEvaluate(Chol2 %*% t(Chol2),bivHomFit)
LLHom <- mxEvaluate(objective,bivHomFit);

Chi= LLhet-LLhom;
LRT= rbind(LLhet,LLhom,Chi); LRT


