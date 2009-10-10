# -----------------------------------------------------------------------
# Program: DefinitionMeans_PathRaw.R  
#  Author: Mike Neale
#    Date: 08 01 2009 
#
# Definition Means model to estimate moderation effect of measured variable 
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has means of 1 and 2 for x and y
#The group with a definition value of 0 has means af zero for x and y 
#The definition variable is used to define a mean deviation of the group with definition value 1

# Simulate data
# -----------------------------------------------------------------------
library(MASS)
set.seed(200)
N=500
Sigma <- matrix(c(1,.5,.5,1),2,2)
group1<-mvrnorm(N, c(1,2), Sigma) # use mvrnorm from MASS package
group2<-mvrnorm(N, c(0,0), Sigma)

y<-rbind(group1,group2)           # Bind both groups together by rows
dimnames(y)[2]<-list(c("x","y")); # Add names
def    <-rep(c(1,0),each=N);      # Add a definition variable 2n in length for group status
selVars<-c("x","y")               # Make a selection variables object

#Define model
# -----------------------------------------------------------------------
defMeansModel <- mxModel("Definition Means -- Path Specification", 
	type="RAM",
	manifestVars=selVars,
	latentVars  ="DefDummy",
	# variances
    mxPath(
    	from=c("x","y"), 
    	arrows=2, 
    	free= TRUE, 
    	values=1,  
    	labels=c("Varx","Vary")
    ),
    # covariances  
    mxPath(
    	from="x", 
    	to="y", 
    	arrows=2, 
    	free= TRUE, 
    	values=.1, 
    	labels=c("Covxy")
    ), 
    # means      
    mxPath(
    	from="one", 
    	to=c("x","y"), 
    	arrows=1, 
    	free= TRUE, 
    	values=1, 
    	labels=c("meanx","meany")
    ), 
    # definition value 
    mxPath(
    	from="one", 
    	to="DefDummy", 
    	arrows=1, 
    	free= FALSE, 
    	values=1, 
    	labels="data.def"
    ),    
    # beta weights
    mxPath(
    	from="DefDummy", 
    	to=c("x","y"), 
    	arrows=1, 
    	free= TRUE, 
    	values=1, 
    	labels=c("beta_1","beta_2")
    ), 
    mxData(
    	observed=data.frame(y,def), 
    	type="raw"
    )
)

# Run the model
# -----------------------------------------------------------------------
defMeansFit<-mxRun(defMeansModel)

defMeansFit@matrices
defMeansFit@algebras


# Compare OpenMx estimates to summary statistics from raw data, 
# -----------------------------------------------------------------------
# Remember to knock off 1 and 2 from group 1's data
# so as to estimate variance of combined sample without the mean correction.
# First we compute some summary statistics from the data
ObsCovs        <- cov(rbind(group1 - rep(c(1,2), each=N), group2))
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
