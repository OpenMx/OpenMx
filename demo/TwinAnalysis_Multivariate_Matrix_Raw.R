# Multivariate matrix-style twin script...
# History: tbates 24.09.2009 first draft

require(OpenMx)
#Import Data
data(twinData)
twinVar = names(twinData); twinVar
# [1] "fam"   "age"   "zyg"   "part"  "wt1"   "wt2"   "ht1"   "ht2"   "htwt1" "htwt2" "bmi1"  "bmi2" 
selVars <- c('bmi1','ht1', 'bmi2', 'ht2');  # pick out variables to be modeled (in this case two), for twin 1 and then twin 2
mzfData <- as.matrix(subset(twinData, zyg==1, selVars)) # assumes MZ F pairs are coded 1
dzfData <- as.matrix(subset(twinData, zyg==3, selVars)) # assumes MZ M pairs are coded 3
head(mzfData);

nVar = length(selVars)/2; # number of dependent variables ** per INDIVIDUAL ( so times-2 for a family)**
# Examine the raw data before going on. You might also plot them
obsMZmeans = colMeans(mzfData,na.rm=TRUE); obsMZmeans;
colMeans(dzfData,na.rm=TRUE);
cov2cor(cov(mzfData,use="complete"));
cov2cor(cov(dzfData,use="complete"));

##### Fit ACE Model  ##### 
# Define 1*1 constant of 0.5 for DZ cov A, or .25 for D
hMatrix = mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h")
qMatrix = mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.25, name="q")

# use MZT1's means for the two traits as start values for the means
meanStarts = rep(obsMZmeans[1:2],2);
	#      bmi1       ht1      bmi1       ht1 
	# 21.353328  1.629724 21.353328  1.629724 

meanLabels = paste("Trait", rep(1:nVar,2),  "mean", sep=""); # make labels for the means matrix
# [1] "Trait1mean" "Trait2mean" "Trait1mean" "Trait2mean"

# grand mean phenotypes (grand mean because labels are the same across zyg and T1 and T2 for each trait
expMZMeans = mxMatrix("Full", nrow=(nVar*2), ncol=1, free=TRUE, values=meanStarts, label=meanLabels, dimnames=list(selVars, "means"), name="expMeanMZ");
# Clone the MZ means matrix for  DZs
expDZMeans = expMZMeans; 
expDZMeans@name="expMeanDZ"; 

# Matrices for path coefficients, Labels for cells, names for the matrices
aLabels <- c(
	"a1v1",
	"a1v2", "a2v2");
cLabels <- c(
	"c1v1",
	"c1v2", "c2v2");
eLabels <- c(
	"e1v1",
	"e1v2", "e2v2");

aMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, labels=aLabels, name="a") # Additive genetic path coefficient
cMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, labels=cLabels, name="c") # Common environmental path coefficient
eMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, labels=eLabels, name="e") # Unique environmental path coefficient

twinACEModel <- mxModel("twinACE",
	hMatrix, expMZMeans, expDZMeans,
	aMatrix, cMatrix, eMatrix,
	# Multiply by each path coefficient by its inverse to get variance component
    mxAlgebra(a %*% t(a), name="A"), # additive genetic variance
    mxAlgebra(c %*% t(c), name="C"), # common environmental variance
    mxAlgebra(e %*% t(e), name="E"), # unique environmental variance

	# MZ ACE algebra
	mxAlgebra(rbind (cbind(A+C+E  , A+C),
                     cbind(A+C    , A+C+E)), dimnames = list(selVars, selVars), name="expCovMZ"),
	# DZ ACE algebra
    mxAlgebra(rbind (cbind(A+C+E  , h%x%A+C),
                     cbind(h%x%A+C, A+C+E)), dimnames = list(selVars, selVars), name="expCovDZ"),
	# Import raw data for each group, setting its objective to the respective group algebras
    mxModel("MZ", mxData(mzfData, type="raw"), mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMeanMZ")),
    mxModel("DZ", mxData(dzfData, type="raw"), mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMeanDZ")),
	# Make algebra for sum of the group objectives,    
	mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
	# Set this to be the objective for the supermodel
	mxAlgebraObjective("twin")
)

#Run ACE model
twinACEFit = mxRun(twinACEModel)

#Extract various fitted results
MZc = mxEval(expCovMZ,  twinACEFit) # Same effect as expCovMZ@matrices$twinACEFit
DZc = mxEval(expCovDZ,  twinACEFit)
M   = mxEval(expMeanMZ, twinACEFit)
A   = mxEval(A, twinACEFit)
C   = mxEval(C, twinACEFit)
E   = mxEval(E, twinACEFit)

# Calculate the genetic (and C & E) covariances
rA = cov2cor(A) # calculates rg (off diagonal)
rC = cov2cor(C) # calculates rc (off diagonal)
rE = cov2cor(E) # calculates re (off diagonal)
ACErs = data.frame(round(cbind(rA,rC,rE),2), row.names=selVars[1:nVar]); 
names(ACErs)<- paste(rep(c("A", "C", "E"),each=2), rep(1:2), sep=""); ACErs

# Calc standardised variance components
Vtot= (A+C+E)
a2  = A/Vtot
c2  = C/Vtot
e2  = E/Vtot
ACEest = data.frame(round(cbind(a2,c2,e2),2), row.names=selVars[1:nVar])
# make a nice-ish table
names(ACEest)<- paste(rep(c("A", "C", "E"),each=2), rep(1:2), sep="")

# Get the fit, and print along with standardised ACE
LL_ACE = mxEval(objective, twinACEFit)
ACEest; LL_ACE
