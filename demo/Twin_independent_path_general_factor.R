# SCRIPT: Twin_independent_path_general_factor.R
# Author: Timothy  Bates tim.bates@ed.ac.uk
# History:  version 1: Fri Sep 18 12:56:23 BST 2009
# Rhelp: http://www.statmethods.net
# OpenMx: http://www.openmx.virginia.com
##########################################
# The diagram of this model is as follows: (not including variance and means) (open in omnigraffle or other .dot capable program)
# digraph G {
# 	Ag->height;
# 	Ag->weight;
# 	Ag->BMI;
# 
# 	Cg->height;
# 	Cg->weight;
# 	Cg->BMI;

# 	Eg->height;
# 	Eg->weight;
# 	Eg->BMI;

# 	Asp1->height;
# 	Asp2->weight;
# 	Asp3->BMI;
# 
# 	Csp1->height;
# 	Csp2->weight;
# 	Csp3->BMI;
# 
# 	Esp1->height;
# 	Esp2->weight;
# 	Esp3->BMI;
# }
require(OpenMx); 

# Get Data
data(twinData)
mzData  = subset(twinData, zyg==1, selVars) # assumes MZ pairs are coded 1
dzData  = subset(twinData, zyg==3, selVars) # assumes DZ pairs are coded 3

# SELVARS must be all values FOR twin1, THEN all values for twin 2....
selVars = c('ht1', 'wt1', 'bmi1',
            'ht2', 'wt2', 'bmi2');
nFac = 1; # 1 general factor
nSib = 2; # two siblings (twins) per family
nVar = length(selVars)/nSib; # number of dependent variables per INDIVIDUAL

obsMZmeans = colMeans(mzData,na.rm=TRUE); # calculate means from the data to give accurate start values
meanStarts = rep(obsMZmeans[1:nVar],2);   # ?rep to learn about this very helpful but slightly opaque command
meanLabels = paste("Trait", rep(1:nVar,2),  "mean", sep=""); # make labels for the means matrix: same label for each trait in each twin
# [1] "Trait1mean" "Trait2mean" "Trait3mean" "Trait1mean" "Trait2mean" "Trait3mean"

##### Fit ACE Model  #####
share = mxModel("all",
	# Variance and means model shared across all groups, and in a group whose name we won't wish to change
	mxMatrix("Full", nrow=nVar, ncol=nFac, free=T, values=.6, name="a_c"), # column of general Additive genetic path coefficient
	mxMatrix("Full", nrow=nVar, ncol=nFac, free=T, values=.6, name="c_c"), # column of general Common environmental path coefficient
	mxMatrix("Full", nrow=nVar, ncol=nFac, free=T, values=.6, name="e_c"), # column of general Unique environmental path coefficient

	mxMatrix("Diag", nrow=nVar, ncol=nVar, free=T, values=.6, name="a"), # Additive genetic path coefficient
	mxMatrix("Diag", nrow=nVar, ncol=nVar, free=T, values=.6, name="c"), # Common environmental path coefficient
	mxMatrix("Diag", nrow=nVar, ncol=nVar, free=T, values=.6, name="e"), # Unique environmental path coefficient

	# Multiply by each path coefficient by its inverse to get variance component
	mxAlgebra(a_c %*% t(a_c) + a %*% t(a), name="A"), # Total genetic variance is (general variance + specific)
	mxAlgebra(c_c %*% t(c_c) + c %*% t(c), name="C"), # common environmental variance
	mxAlgebra(e_c %*% t(e_c) + e %*% t(e), name="E"), # unique environmental variance
	mxMatrix("Full", nrow=1, ncol=(nVar*2), free=TRUE, values=meanStarts, label=meanLabels, dimnames=list( "means",selVars), name="expMean")
)

mzModel = mxModel(name="MZ", 
	# MZ ACE algebra: note how this refers to the A C and E matrices in the shared "all" group
	mxAlgebra(rbind (cbind(all.A+all.C+all.E, all.A+all.C),
                   cbind(all.A+all.C,   all.A+all.C+all.E) ), dimnames = list(selVars, selVars), name="expCov"),
	mxData(mzData, type="raw"),
	mxFIMLObjective("expCov", "all.expMean")
)

dzModel = mxModel(name="DZ", hMatrix,
	mxAlgebra(rbind (cbind( all.A+all.C+all.E , .5 %x% all.A+all.C),
                   cbind(.5 %x% all.A+all.C,  all.A+all.C+all.E)  ), dimnames = list(selVars, selVars), name="expCov"),
	mxData(dzData, type="raw"),
	mxFIMLObjective("expCov", "all.expMean")
)

# Build and Run ACE model
model = mxModel("twinACE",share, mzModel, dzModel, mxAlgebra(MZ.objective + DZ.objective, name="twin"), mxAlgebraObjective("twin") )
fit   = mxRun(model);

a   = mxEval(all.a, fit); # specific path coefficients
c   = mxEval(all.c, fit);
e   = mxEval(all.e, fit);

a_c   = mxEval(all.a_c, fit); # general path coefficients
c_c   = mxEval(all.c_c, fit);
e_c   = mxEval(all.e_c, fit);

A   = mxEval(all.A, fit);# Variances
C   = mxEval(all.C, fit);
E   = mxEval(all.E, fit);
Vtot=A+C+E;			# total variance
I     <- matrix(0,nVar,nVar);  diag(I)=diag(matrix(1,nVar,nVar));
SD    <- solve(sqrt(I*Vtot))   # inverse of diagonal matrix of standard deviations  (same as "(\sqrt(I.Vtot))~"

logLikelihood = mxEval(objective, fit);
print(logLikelihood[1,1]);

# standardized _path_ coefficients ready to be stacked together
a2 <- SD %*% a ; # Standardized path coefficients
c2 <- SD %*% c ;
e2 <- SD %*% e ;

StandardizedEstimates  = data.frame(cbind(a2,c2,e2), row.names=selVars[1:nVar]);
names(StandardizedEstimates)= paste(rep(c("A", "C", "E"), each=nVar), rep(1:nVar), sep="");
print(round(StandardizedEstimates,accuracy)); 

a_c2 <- SD %*% a_c ; # Standardized estimate of common A factor
c_c2 <- SD %*% c_c ;
e_c2 <- SD %*% e_c ;

StandardizedEstimates  = data.frame(cbind(a_c2,c_c2,e_c2), row.names=selVars[1:nVar]);
names(StandardizedEstimates)= c("gA", "gC", "gE");
print(round(StandardizedEstimates,accuracy)); 