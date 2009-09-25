# SCRIPT: NTF_design.R - NTF design in OpenMx 
# Author: Matt Keller & Sarah Medland
# History:  Thu Sep 24 17:23:33 BST 2009
# 2009-09-24: (tb): Use three groups (change MZ and DZ family algebra so that they refer to a common set of matrices in a new "NTF" group); 
# 2009-09-24: (tb): Simplify/speed algebra using Quadratic operator and pre-calculating variance components i.e., e %*% t(e) = E etc; 
# Rhelp: http://www.statmethods.net
# OpenMx: http://www.openmx.virginia.com
##########################################
require(OpenMx); require(foreign); require(MASS);

# TODO: Add code to simulate data
# TODO: Add output from old MZ to verfiy correctness
# TODO: Verify algebra/constraint specification
#        Note that there might be a minor error somewhere in the algebra as the 
#        Var(D) is slightly underestimated and Var(S) is slightly overestimated:

# > res.mat
#                  Var.E   Var.F  Var.A  Var.S   Var.D  Cor.Sps Var.Phen
# OpenMx-Estimated 0.310   0      0.421  0.259   0.049   0.189    1.065
# Simulated        0.305   0      0.432  0.200   0.099   0.202    1.036

require(OpenMx)
#Run Gene Evolve
setwd("/temp/OMX")
source("GE-73.R")
# Mat or Sarah - any chance of putting a simulated-data file online so the script is a self-contained test?
# NuclearTwinFamSim <- as.matrix(read.table("http://http://openmx.psyc.virginia.edu/sites/default/files/NuclearTwinFamSim",header=FALSE,col.names=selVars))

selVars <- c("Famid","Tw1","Tw2","Fa","Mo","bro1","bro2","sis1","sis2","sp.tw1","sp.tw2","son1.tw1","son2.tw1","dau1.tw1","dau2.tw1", "son1.tw2","son2.tw2","dau1.tw2","dau2.tw2","Tw1.age","Tw2.age","Fa.age","Mo.age","bro1.age","bro2.age","sis1.age","sis2.age","sp.tw1.age","sp.tw2.age", "son1.tw1.age","son2.tw1.age","dau1.tw1.age","dau2.tw1.age", "son1.tw2.age","son2.tw2.age","dau1.tw2.age","dau2.tw2.age")
mzf <- as.matrix(read.table("MZF",header=FALSE,col.names=selVars))
mzf <- mzf[,2:5]
mzf[mzf==-999] <- NA
dzf <- as.matrix(read.table("DZF",header=FALSE,col.names=selVars))
dzf <- dzf[,2:5]
dzf[dzf==-999] <- NA

#Fit NTF Model with RawData and Matrices Input
ntf <- mxModel(model="NTF",
	# Matrices
	mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, values=1,  label="FamilialVar",  name="f"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.1,  label="FamilialPath", name="m"),

	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.7,  label="Env",          name="e"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.6,  label="AddGen",       name="a"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.3,  label="Sib",          name="s"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="Dominance",    name="d"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.1,  label="AMCopath",     name="mu"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1.2, label="VarPhen",      name="Vp1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="VarF",         name="x1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.25, label="CovPhenGen",   name="delta1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1,   label="VarAddGen",    name="q1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.15, label="CovFA",        name="w1"),

	#mxAlgebra section - nonlinear constraints
	mxAlgebra(e %*% t(e), name="E"),
	mxAlgebra(d %*% t(d), name="D"),
	mxAlgebra(s %*% t(s), name="S"),

	mxAlgebra((a %&% q1) + (f %&% x1) + 2 %x% a %*% w1 %*% t(f) + E + D + S, name="Vp2"),
	mxAlgebra(2 %x% m %&% Vp1 + 2 %x% (m %&% (mu %&% Vp1)), name="x2"),

	mxAlgebra(q1 %*% a + w1 %*% f, name="delta2"),
	mxAlgebra(1 + delta1 %*% mu %*% t(delta1), name="q2"),
	mxAlgebra(delta1 %*% m + delta1 %*% mu %*% Vp1 %*% t(m), name="w2"),

	#constraints - equating nonlinear constraints and parameters
	mxConstraint('Vp1',   "=", 'Vp2',   name='VpCon'),
	mxConstraint('x1',    "=", 'x2',    name='xCon'),
	mxConstraint('delta1',"=", 'delta2',name='deltaCon'),
	mxConstraint('q1',    "=", 'q2',    name='qCon'),
	mxConstraint('w1',    "=", 'w2',    name='wCon'),
	#mxAlgebra section - relative covariances
	mxAlgebra(a   %&% q1 + f %&% x1 + 2 %x% a %*% w1 %*% t(f) + D + S, name="CvMz"),
	mxAlgebra(a   %&% (q1-.5)  + .25 %x% d %*% t(d) + f %*% x1 %*% t(f) + 2 %x% a %*% w1 %*% t(f) + S, name="CvDz"),
	mxAlgebra(.5  %x% a %*% (q1 %*% a + w1 %*% f) + .5 %x% a %*% (q1 %*% a + w1 %*% f) %*% mu %*% Vp1 + m %*% Vp1 + m %*% Vp1 %&% mu, name="ParChild"),
	mxAlgebra(Vp1 %&% mu, name="CvSps")
)

mzModel <- mxModel(name = "MZNTF",
	mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values= .25, label="mean", dimnames=list(NULL, selVars), name="expMeanMz"),
	# Algebra for expected variance/covariance matrix in MZF
	mxAlgebra(expression=rbind(
		cbind(NTF.Vp1,      NTF.CvMz,     NTF.ParChild, NTF.ParChild),
		cbind(NTF.CvMz,     NTF.Vp1,      NTF.ParChild, NTF.ParChild),
		cbind(NTF.ParChild, NTF.ParChild, NTF.Vp1,      NTF.CvSps),
		cbind(NTF.ParChild, NTF.ParChild, NTF.CvSps,    NTF.Vp1)
	), 
	dimnames=list(selVars,selVars),name="expCovMz"),
	mxData(observed=mzf, type="raw"),
	mxFIMLObjective(covariance="expCovMz",means="expMeanMz")
)

dzModel <- mxModel(name = "DZNTF",
	mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values= .25, label="mean", dimnames=list(NULL, colnames(dzf)), name="expMeanDz"),
	# Algebra for expected variance/covariance matrix in MZF
	mxAlgebra(expression=rbind(
		cbind(NTF.Vp1,      NTF.CvDz,     NTF.ParChild, NTF.ParChild),
		cbind(NTF.CvDz,     NTF.Vp1,      NTF.ParChild, NTF.ParChild),
		cbind(NTF.ParChild, NTF.ParChild, NTF.Vp1,      NTF.CvSps),
		cbind(NTF.ParChild, NTF.ParChild, NTF.CvSps,    NTF.Vp1)
	), 
	dimnames=list(selVars,selVars),name="expCovDz"),
	mxData(observed=dzf, type="raw"),
	mxFIMLObjective(covariance="expCovDz",means="expMeanDz")
)

NTFModel <- mxModel(model="NucTwFam", mzModel, dzModel, NucTwFam,
	mxAlgebra(expression=MZNTF.objective + DZNTF.objective, name="ntffit"), #MZ.objective is the automatic name for the -2LL of mzModel
	mxAlgebraObjective("ntffit")
)

#Run MX
NTFModelFit <- mxRun(NTFModel)
#Look at results
summary(NTFModelFit)
res <- NTFModelFit@output$estimate
round(res,3)
#compare to simulation
res.mat <- rbind(round(c(res[1:5]^2,res[6:7]),3),round(ALL$track.changes[c('var.U','var.F','var.A','var.S','var.D','cor.spouses','var.cur.phenotype'),'data.t1'],3))
dimnames(res.mat) <- list(c('OpenMx-Estimated','Simulated'),c('Var.E','Var.F','Var.A','Var.S','Var.D','Cor.Sps','Var.Phen'))

#look at implied Covariances
round(NTFModelFit@output$algebras$MZNTF.expCovMz,3)
round(NTFModelFit@output$algebras$DZNTF.expCovDz,3)