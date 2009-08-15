require(OpenMx)
library(MASS)
#Definition Variable Test 3
#Author: Mike Neale
#Date: July 29 2009

#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has means of 1 and 2 for x and y
#The group with a definition value of 0 has means af zero for x and y 
#The definition variable is then used to define a mean deviation of the group with definition value 1

#make some data!
set.seed(200)
n = 500

Sigma <- matrix(c(1,.5,.5,1),2,2)
group1<-mvrnorm(n=n, c(1,2), Sigma)
group2<-mvrnorm(n=n, c(0,0), Sigma)

#put them both together, add a definition variable, and make an selection variables object
y<-rbind(group1,group2)
dimnames(y)[2]<-list(c("x","y"))
def<-rep(c(1,0),each=n)
selvars<-c("x","y")
# write data to a file for the mx script to read (not necessary for running in R)
# write.table(cbind(y,def),file="temp-files/xydefmeans.rec",col.names=F,row.names=F)

# Three covariance model matrices: 
#  "cov" for the zero relationship group
#  "def" for the definition variable, 
#   and "beta" for estimating difference between groups' covariances
# One common mean vector, "M"

#define the model, including a FIML objective function, which will optimize the matrix S
defMeansModel<-mxModel("Definition Means via Paths", 
    type="RAM",
    mxData(
        observed=data.frame(y,def), 
        type="raw"),
    manifestVars=c("x","y"),
    latentVars="DefDummy",
    mxPath(from=c("x","y"), 
        arrows=2,
        free=TRUE,
        values=c(1,.1,1),
        labels=c("Varx","Vary")
    ), # variances
    mxPath(from="x", to="y",
        arrows=2,
        free=TRUE,
        values=c(.1),
        labels=c("Covxy")
    ), # covariances
    mxPath(from="one",
        to=c("x","y","DefDummy"),
        arrows=1,
        free=c(TRUE,TRUE,FALSE),
        values=c(1,1,1),
        labels =c("meanx","meany","data.def")
    ),
    mxPath(from="DefDummy",
        to=c("x","y"),
        arrows=1,
        free=c(TRUE,TRUE),
        values=c(1,1),
        labels =c("beta_1","beta_2")
    )
)

# Run the model
defMeansFit<-mxRun(defMeansModel)
defMeansFit@matrices
defMeansFit@algebras

# Compare OpenMx estimates to summary statistics from raw data, remembering to knock off 1 and 2 from group 1's
# data, so as to estimate variance of combined sample without the mean correction.

# First we compute some summary statistics from the data
ObsCovs <- cov(rbind(group1 - rep(c(1,2), each=n), group2))
ObsMeansGroup1 <- c(mean(group1[,1]), mean(group1[,2]))
ObsMeansGroup2 <- c(mean(group2[,1]), mean(group2[,2]))

# Second we extract the parameter estimates and matrix algebra results from the model
Sigma <- mxEval(S[1:2,1:2], defMeansFit)
Mu <- mxEval(M[1:2], defMeansFit)
beta <- mxEval(A[1:2,3], defMeansFit)

# Third, we check to see if things are more or less equal
omxCheckCloseEnough(ObsCovs,Sigma,.01)
omxCheckCloseEnough(ObsMeansGroup1,as.vector(Mu+beta),.001)
omxCheckCloseEnough(ObsMeansGroup2,as.vector(Mu),.001)

