require(OpenMx)
library(MASS)
# Definition Variable Test 3
# Author : Mike Neale
# History: tbates: September 5: fixed typo in values, rearranged for clarity, added expected results
# Date   : July 29 2009

#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has means of 1 and 2 for x and y
#The group with a definition value of 0 has means af zero for x and y 
#The definition variable is then used to define a mean deviation of the group with definition value 1

# make some data!
set.seed(200)
n = 500
Sigma <- matrix(c(1,.5,.5,1),2,2)
group1<-mvrnorm(n=n, mu=c(1,2), Sigma) # use mvrnorm from MASS package
group2<-mvrnorm(n=n, mu=c(0,0), Sigma)

y<-rbind(group1,group2)           # Bind both groups together by rows
dimnames(y)[2]<-list(c("x","y")); # Add names
def    <-rep(c(1,0),each=n);      # Add a definition variable 2n in length for group status
selvars<-c("x","y")               # Make a selection variables object
# write data to a file for the mx script to read (not necessary for running in R)
# write.table(cbind(y,def),file="temp-files/xydefmeans.rec",col.names=FALSE,row.names=FALSE)

# Three sets of covariance paths: 
#  "cov" for the zero relationship group
#  "def" for the definition variable, 
#   and "beta" for estimating difference between groups' covariances
# One common mean vector, "M"

#define the model, including a FIML objective function, which will optimize the matrix S
defMeansModel<-mxModel("Definition Means via Paths", type="RAM",
		manifestVars=c("x","y"),
		latentVars  ="DefDummy",
    mxPath(from=c("x","y"),                arrows=2, free= TRUE,  values=1,  labels=c("Varx","Vary")),     # variances
    mxPath(from="x",        to="y",        arrows=2, free= TRUE,  values=.1, labels=c("Covxy")),           # covariances
    mxPath(from="one",      to=c("x","y"), arrows=1, free= TRUE,  values=1,  labels=c("meanx","meany")),   # means
    mxPath(from="one",      to="DefDummy", arrows=1, free= FALSE, values=1,  labels="data.def"),           # defn value
    mxPath(from="DefDummy", to=c("x","y"), arrows=1, free= TRUE,  values=1,  labels=c("beta_1","beta_2")), # beta weights
    mxData(observed=data.frame(y,def), type="raw")
)

# Run the model
defMeansFit<-mxRun(defMeansModel)
defMeansFit@matrices
defMeansFit@algebras
# mxAlgebra 'Z' 
#      [,1] [,2]      [,3]
# [1,]    1    0 0.9653445
# [2,]    0    1 1.9713524
# [3,]    0    0 1.0000000
# 
# mxAlgebra 'covariance' 
#           x         y
# x 1.0537625 0.4952057
# y 0.4952057 0.9626567
# 
# mxAlgebra 'means' 
#               x          y
# [1,] 0.05725722 0.03859851

# Compare OpenMx estimates to summary statistics from raw data, remembering to knock off 1 and 2 from group 1's
# data, so as to estimate variance of combined sample without the mean correction.
# First we compute some summary statistics from the data
ObsCovs        <- cov(rbind(group1 - rep(c(1,2), each=n), group2))
ObsMeansGroup1 <- c(mean(group1[,1]), mean(group1[,2]))
ObsMeansGroup2 <- c(mean(group2[,1]), mean(group2[,2]))

# Second we extract the parameter estimates and matrix algebra results from the model
Sigma <- mxEval(S[1:2,1:2], defMeansFit)
Mu    <- mxEval(M[1:2], defMeansFit)
beta  <- mxEval(A[1:2,3], defMeansFit)

# Third, we check to see if things are more or less equal
omxCheckCloseEnough(ObsCovs,Sigma,.01)
omxCheckCloseEnough(ObsMeansGroup1,as.vector(Mu+beta),.001)
omxCheckCloseEnough(ObsMeansGroup2,as.vector(Mu),.001)