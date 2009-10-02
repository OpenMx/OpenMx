require(OpenMx)

# -----------------------------------------------------------------------
# Program: BivariateCorrelation20090925.R  
#  Author: Hermine Maes
#    Date: Wed Sep 25 11:45:52 EDT 2009
#
# Revision History
#   Hermine Maes -- Wed Sep 25 11:45:52 EDT 2009 BivariateCorrelation20090925.R
# -----------------------------------------------------------------------

# Simulate Data: two standardized variables X & Y with correlation of .5
# -----------------------------------------------------------------------
require(MASS)
set.seed(200)
rs=.5
xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
testData <- xy
selVars <- c('X','Y')
dimnames(testData) <- list(NULL, selVars)
summary(testData)
cov(testData)


# Fit Saturated Model with Raw Data and Matrix-style Input
# -----------------------------------------------------------------------
bivCorModel <- mxModel("bivCor",
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=TRUE, 
        values=c(0,0), 
        name="expMean"
    ), 
    mxMatrix(
        type="Lower", 
        nrow=2, 
        ncol=2, 
        free=TRUE,
        values=.5, 
        dimnames=list(selVars, selVars), 
        name="Chol"
    ), 
    mxAlgebra(
        expression=Chol %*% t(Chol), 
        name="expCov", 
    ), 
    mxData(
        observed=testData, 
        type="raw"
    ), 
    mxFIMLObjective(
        covariance="expCov", 
        means="expMean",
        dimnames=selVars)
    )

# Run Model and Generate Output
# -----------------------------------------------------------------------
bivCorFit <- mxRun(bivCorModel)
EM <- mxEval(expMean, bivCorFit)
EC <- mxEval(expCov, bivCorFit)
LL <- mxEval(objective, bivCorFit)


# Specify SubModel testing Covariance=Zero
# -----------------------------------------------------------------------
bivCorModelSub <-mxModel(bivCorModel,
    mxMatrix(
        type="Diag", 
        nrow=2, 
        ncol=2, 
        free=TRUE, 
        name="Chol", 
        dimnames=list(selVars, selVars)
    )
)

# Run Model and Generate Output
# -----------------------------------------------------------------------
bivCorFitSub <- mxRun(bivCorModelSub)
EMs <- mxEval(expMean, bivCorFitSub)
ECs <- mxEval(expCov, bivCorFitSub)
LLs <- mxEval(objective, bivCorFitSub)
Chi= LLs-LL;
LRT= rbind(LL,LLs,Chi); LRT


# Mx Answers of Saturated Model Hard-coded
# -----------------------------------------------------------------------
Mx.EM <- matrix(c(0.03211656, -0.004883885), 1, 2)
Mx.EC <- matrix(c(1.0092853, 0.4813504, 0.4813504, 0.9935390), 2, 2)
Mx.LL <- 5415.772

# Compare OpenMx Results to Mx Results 
# -----------------------------------------------------------------------
#LL: likelihood; EC: expected covariance, EM: expected means
omxCheckCloseEnough(LL,Mx.LL,.001)
omxCheckCloseEnough(EC,Mx.EC,.001)
omxCheckCloseEnough(EM,Mx.EM,.001)
