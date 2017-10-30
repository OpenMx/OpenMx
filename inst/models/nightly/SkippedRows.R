library(OpenMx)
library(mvtnorm)
library(MASS)

set.seed(1)
nVar <- 3
Data1<-rmvnorm(2000, mean=rep(0,nVar), sigma=50*diag(nVar))
selVars<-paste0("Var",1:nVar)
colnames(Data1)<-selVars
Data1 <- as.data.frame(Data1)
Data1$group <- sample.int(nrow(Data1), replace=TRUE)

for (condOn in c('ordinal', 'continuous')) {
	m1 <- mxModel("m1",
		      mxMatrix("Full", 1,1, labels="data.group", name="ForceSomewhatRowwise"),
		      mxMatrix("Symm", nVar, nVar, values=diag(.1, nVar, nVar), name="expCov"),
		      mxMatrix("Full", 1, nVar, unlist(Data1[1,1:nVar]), free=TRUE, name="expMean"),
		      mxData(Data1, type="raw"),
		      mxExpectationNormal("expCov", "expMean", dimnames=selVars),
		      mxFitFunctionML(jointConditionOn=condOn))
	m1$expCov$free[,] <- m1$expCov$values != 0
#	m1$fitfunction$verbose <- 2L
	# Need to retry lots of times because it's hard to take the
	# gradient with the same number of skipped rows.
	m1fit <- mxTryHard(m1, extraTries = 39)
	print(summary(m1fit))
	omxCheckCloseEnough(m1fit$output$fit, 40719.05, .01)
}

# -----------------------------------------------------------------------
# Mixture distribution model: 50:50 mixture of two normal distributions
#  that have different means but same variance.
#
# (Mike Neale's twisted idea)
# -----------------------------------------------------------------------

nClass <-2
meanDiff <- 50
Data1<-mvrnorm(100, mu=rep(0,nVar), Sigma=diag(nVar))
Data2<-mvrnorm(100, mu=rep(meanDiff,nVar), Sigma=diag(nVar))
Data <- as.data.frame(rbind(Data1,Data2))
selVars<-paste0("Var",1:nVar)
names(Data)<-selVars
labFun <- function(name="matrix",nrow=1,ncol=1,lower=FALSE,symmetric=FALSE)
            {
                matlab <- matrix(paste(rep(name, each=nrow*ncol), rep(rep(1:nrow),ncol), rep(1:ncol,each=nrow),sep="_"),nrow,ncol)
                if (lower) {for (i in 1:(nrow-1)){for (j in (i+1):ncol){matlab[i,j]<-NA}}}
				if(symmetric){for (i in 1:(nrow-1)){for (j in (i+1):ncol){matlab[i,j]<-matlab[j,i]}}}
                return(matlab)
            }

# Build simple single-variable model
g1Model <-  mxModel("group1",
                mxMatrix("Symm", nVar, nVar, diag(2*pi,nVar,nVar), labels=labFun("expCov",nrow=nVar,ncol=nVar,symmetric=T), name="expCov",free=T),
                mxMatrix("Full", 1, nVar, values=rep(meanDiff/2+1,3), free=T, name="expMean"),
                mxData(Data, type="raw"), 
                mxExpectationNormal("expCov", "expMean",dimnames=selVars), mxFitFunctionML(vector=T)
                )

# Can repeat for more classes in a loop as needed
g2Model <-  mxModel(g1Model, name="group2", 
                mxMatrix("Full", 1, nVar, values=rep(meanDiff/2-1,3), free=T, name="expMean")
                )

mixtureModel <- mxModel("mixture", g1Model, g2Model,
                    mxMatrix(type="Full", nrow=nClass, ncol=1, values=1:nClass, free=c(F,T), lbound=1e-6, name="pRaw"),
					mxAlgebra(pRaw/sum(pRaw),name="p"),
                    mxAlgebra(-2*sum(log(cbind(group1.objective, group2.objective) %*% p )), name="min2LL"), 
                    mxFitFunctionAlgebra("min2LL")
                )

mixtureModelFit <- mxRun(mixtureModel)

print(summary(mixtureModelFit))

omxCheckCloseEnough(mixtureModelFit$group1$expCov$values, diag(3), .15)
omxCheckCloseEnough(mixtureModelFit$group1$expMean$values[1,], rep(50,3), .22)
omxCheckCloseEnough(mixtureModelFit$group2$expMean$values[1,], rep(0,3), .22)
omxCheckCloseEnough(mixtureModelFit$pRaw$values[2,1], 1, 1e-4)
