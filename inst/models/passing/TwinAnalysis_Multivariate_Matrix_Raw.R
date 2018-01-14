#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


# Multivariate matrix-style twin script...
# History: tbates 24.08.2009 first draft
# History: tbates 12.09.2009 
# 	Convert to use three models: model, mz, and dz
# 	Reference MZ-group A,C,E matrices in DZ algebra
# 	Compute standardized path coefficients:
# 		var = matrix(0,nVar,nVar);  diag(var)=diag(Vtot);  # variances on the diagonal, 0's elsewhere.
# 		SD =  solve(sqrt(var))   # Inverse (solve) of diagonal matrix of standard deviations: \sqrt(var)~ in oldMx-speak
# 	Improved comments
# TODO (tb) Write equivalent script in old-mx
# TODO (tb) Add output from old-mx
# TODO (tb) Add near-enough calls to verify output

require(OpenMx)

#Import Data

data(twinData)
twinVar = names(twinData); twinVar
# [1] "fam"   "age"   "zyg"   "part"  "wt1"   "wt2"   "ht1"   "ht2"   "htwt1" "htwt2" "bmi1"  "bmi2" 
selVars <- c('ht1', 'bmi1','ht2','bmi2');  # pick out variables to be modeled (in this case two), for twin 1 and then twin 2
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

# use MZT1's means for the two traits as start values for the means
meanStarts = rep(obsMZmeans[1:2],2);
	#      bmi1       ht1      bmi1       ht1 
	# 21.353328  1.629724 21.353328  1.629724 
meanLabels = paste("Trait", rep(1:nVar,2),  "mean", sep=""); # make labels for the means matrix
# [1] "Trait1mean" "Trait2mean" "Trait1mean" "Trait2mean"

# grand mean phenotypes (grand mean because labels are the same across zyg and T1 and T2 for each trait
expMZMeans = mxMatrix("Full", nrow=1, ncol=(nVar*2), free=TRUE, values=meanStarts, label=meanLabels, dimnames=list("means", selVars), name="expMeanMZ");
# Clone the MZ means matrix for  DZs
expDZMeans = expMZMeans; 
expDZMeans$name="expMeanDZ"; 

# Matrices for path coefficients
lb <- matrix(NA, nVar, nVar); diag(lb) <- 0
aMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, lbound=lb, name="a") # Additive genetic path coefficient
cMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, lbound=lb, name="c") # Common environmental path coefficient
eMatrix = mxMatrix("Lower", nrow=nVar, ncol=nVar, free=TRUE, values=.5, lbound=lb, name="e") # Unique environmental path coefficient

# Make the mz group: define variance as square of path coefficients, define algebra of twin covariance, 
# import mz data, and set objective to model the covariance and means observed in these data

mzGroup <- mxModel("mz",
	hMatrix, aMatrix, cMatrix, eMatrix,
	expMZMeans,
	# Multiply by each path coefficient by its inverse to get variance component
  mxAlgebra(a %*% t(a), name="A"), # additive genetic variance
  mxAlgebra(c %*% t(c), name="C"), # common environmental variance
  mxAlgebra(e %*% t(e), name="E"), # unique environmental variance
	# MZ ACE algebra
	mxAlgebra(rbind (cbind(A+C+E  , A+C),
                   cbind(A+C    , A+C+E)), dimnames = list(selVars, selVars), name="expCov"),
	# Import raw data setting and set objective
	mxData(mzfData, type="raw"), 
	mxExpectationNormal("expCov", "expMeanMZ"),
	mxFitFunctionML()
)

# make the dz group: much simpler - just has its data and means, and refers to the ACE matrices 
# in mz group so that both groups are fitted of the same parameters
dzGroup <- mxModel("dz",
	hMatrix, 
	expDZMeans,
	 mxAlgebra(rbind (cbind(mz.A+ mz.C+ mz.E  , h %x% mz.A + mz.C),
	                  cbind(h %x% mz.A + mz.C, mz.A + mz.C + mz.E)), dimnames = list(selVars, selVars), name="expCov"),
	mxData(dzfData, type="raw"), 
	mxExpectationNormal("expCov", "expMeanDZ"),
	mxFitFunctionML()
)

# Combine the mz and dz groups in a supermodel which can have as its objective maximising the likelihood ofboth groups simultaneously.
model = mxModel("ACE", mzGroup, dzGroup,
	mxFitFunctionMultigroup(c('mz', 'dz'))
)

#Run ACE model
fit = mxRun(model)

#Extract various fitted results
MZc = mxEval(mz.expCov,  fit); # Same effect as expCovMZ$matrices$fit
DZc = mxEval(dz.expCov,  fit);
GEc = mxGetExpected(fit, 'covariance')
omxCheckCloseEnough(list(MZc, DZc), GEc, 1e-10)
M   = mxEval(mz.expMeanMZ, fit);
GEM = mxGetExpected(fit, 'means')
omxCheckCloseEnough(M, GEM[[1]], 1e-10)
A   = mxEval(mz.A, fit);
C   = mxEval(mz.C, fit);
E   = mxEval(mz.E, fit);
Vtot= A+C+E;

a   = mxEval(mz.a, fit);
c   = mxEval(mz.c, fit);
e   = mxEval(mz.e, fit);

# Calc standardised variance components
var = diag(diag(Vtot), nrow=nrow(Vtot));  # variances on the diagonal, 0's elsewhere.
SD    <- solve(sqrt(var))   # Inverse (solve) of diagonal matrix of standard deviations: \sqrt(var)~ in oldMx-speak

# standardized _path_ coefficients ready to be stacked together
a2 <- SD %*% a ; # Standardized path coefficients
c2 <- SD %*% c ;
e2 <- SD %*% e ;
ACEest = data.frame(round(cbind(a2,c2,e2),3), row.names=selVars[1:nVar])
# make a table
names(ACEest)<- paste(rep(c("A", "C", "E"),each=2), rep(1:2), sep="");
print(ACEest);
#          A1    A2 C1 C2     E1    E2
# ht1   0.938 0.000  0  0  0.345 0.000
# bmi1 -0.030 0.883  0  0 -0.145 0.445

# Get the fit, and print along with standardised ACE
LL_ACE = mxEval(fitfunction, fit); print(LL_ACE)

#           [,1]
# [1,] -1500.460
