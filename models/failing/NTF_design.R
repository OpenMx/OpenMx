# SCRIPT: NTF_design.R - NTF design in OpenMx 
# Author: Matt Keller & Sarah Medland; edited by Tim Bates
# History:  Thu Sep 24 17:23:33 BST 2009
# For description of model, see Keller, Medland, Duncan, Hatemi, Neale, Maes, & Eaves (2009) TRHG, 21, p.8 - 18.
# 2009-09-24: (tb): Use three groups (change MZ and DZ family algebra so that they refer to a common set of matrices in a new "ntf" group); 
# 2009-09-26: (tb) read data from DEMO folder	
# 2009-09-30: (tb): Simplify/speed algebra using Quadratic operator and pre-calculating variance components i.e., e %*% t(e) = E etc; 
# 2009-09-30: (tb): 2* speedup using Quadratic operator and pre-calculating variance components i.e., e %*% t(e) = E etc; 
# OpenMx: http://www.openmx.virginia.com
# TODO (tb): Add closeenough calls
# TODO (tb): Resolve question about fixing m...
##########################################

require(OpenMx);

start = proc.time()[1]; # record start time
ntf <- mxModel(model="ntf",
	# Matrices
	# NOTE: NTF design does not allow Vs & Vf to be estimated simultaneously; for identifiability, choose either m or s to be free and the other to be fixed at 0             
	# TODO This DOESN"T run if m is fixed@0
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.7,  label="Env",          name="e"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, values=1,  label="FamilialVar",  name="f"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=0.1, label="FamilialPath", name="m"),  # fix m=0 if you want Vf=0  (nope... non pos def)
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.6,  label="AddGen",       name="a"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.3,  label="Sib",          name="s"),  # fix s=0 if you want Vs=0
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="Dominance",    name="d"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.1,  label="AMCopath",     name="mu"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1.2, label="VarPhen",      name="Vp1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="VarF",         name="x1"),  # keep this parameter free, even if Vf is fixed to 0
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.25, label="CovPhenGen",   name="delta1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1,   label="VarAddGen",    name="q1"),
	mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.15, label="CovFA",        name="w1"),
  # mxAlgebra section - nonlinear constraints
  mxAlgebra(e %*% t(e), name="E"),
  mxAlgebra(d %*% t(d), name="D"),
  mxAlgebra(s %*% t(s), name="S"),

  mxAlgebra(expression= a %&% q1 + f %&% x1 + 2 %x% (a %*% w1 %*% t(f)) + E + D + S, name="Vp2"),
  mxAlgebra(expression= 2 %x% (m %&% Vp1) + 2 %x% m %*% Vp1 %*% mu %*% t(Vp1) %*% t(m), name="x2"),
  mxAlgebra(expression= (q1 %*% a) + (w1 %*% f), name="delta2"),
  mxAlgebra(expression= 1 + delta1 %&% mu, name="q2"),
  mxAlgebra(expression= delta1 %*% m + delta1 %*% mu %*% Vp1 %*% t(m), name="w2"),
  #constraints - equating nonlinear constraints and parameters
  mxConstraint("Vp1",   "=", "Vp2",   name='VpCon'),
  mxConstraint("x1",    "=", "x2",    name='xCon'),
  mxConstraint("delta1","=", "delta2",name='deltaCon'),
  mxConstraint("q1",    "=", "q2",    name='qCon'),
  mxConstraint("w1",    "=", "w2",    name='wCon'),
  # mxAlgebra section - relative covariances
  mxAlgebra(expression= (a %&% q1) + (f %&% x1) + (2 %x% a) %*% w1 %*% t(f) + D+S, name="CvMz"),
  mxAlgebra(expression= a %*% (q1-.5) %*% t(a) + (.25 %x% d) %*% t(d) + f %*% x1 %*% t(f) + (2 %x% a) %*% w1 %*% t(f) + S, name="CvDz"),
  mxAlgebra(expression= (.5 %x% a) %*% (q1 %*% a + w1 %*% f) + (.5 %x% a) %*% (q1 %*% a + w1 %*% f) %*% (mu %*% Vp1) + (m %*% Vp1) + m %*% Vp1 %*% mu %*% t(Vp1), name="ParChild"),
  mxAlgebra(expression= Vp1 %&% mu, name="CvSps")
)

mzModel <- mxModel(name = "MZNTF",
  mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values= .25, label="mean", dimnames=list(NULL, colnames(mzData)), name="expMeanMz"),
  # Algebra for expected variance/covariance matrix in MZF
  mxAlgebra(expression=rbind(
 	  cbind(ntf.Vp1,       ntf.CvMz,       ntf.ParChild, ntf.ParChild),
	  cbind(ntf.CvMz,      ntf.Vp1,        ntf.ParChild, ntf.ParChild),
	  cbind(ntf.ParChild,  ntf.ParChild,   ntf.Vp1,      ntf.CvSps),
	  cbind(ntf.ParChild,  ntf.ParChild,   ntf.CvSps,    ntf.Vp1)
  ), dimnames=list(colnames(mzData),colnames(mzData)),name="expCovMz"),
  mxData(observed=mzData, type="raw"),
  mxFIMLObjective(covariance="expCovMz",means="expMeanMz")
)

dzModel <- mxModel(ntf, name = "DZNTF",
  mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values= .25, label="mean", dimnames=list(NULL, colnames(dzData)), name="expMeanDz"),
  # Algebra for expected variance/covariance matrix in MZF
  mxAlgebra(expression=rbind(
 	  cbind(ntf.Vp1,       ntf.CvDz,       ntf.ParChild, ntf.ParChild),
	  cbind(ntf.CvDz,      ntf.Vp1,        ntf.ParChild, ntf.ParChild),
	  cbind(ntf.ParChild,  ntf.ParChild,   ntf.Vp1,      ntf.CvSps),
	  cbind(ntf.ParChild,  ntf.ParChild,   ntf.CvSps,    ntf.Vp1)
  ), dimnames=list(colnames(dzData),colnames(dzData)),name="expCovDz"),
  mxData(observed=dzData, type="raw"),
  mxFIMLObjective(covariance="expCovDz",means="expMeanDz")
)

model <- mxModel(model="NuclearTwinFamily", ntf, mzModel, dzModel,
  mxAlgebra(expression=MZNTF.objective + DZNTF.objective, name="ntffit"), #groupName.objective is the name for the -2LL of submodels
  mxAlgebraObjective("ntffit")
)

#Run MX
fit <- mxRun(model)
#Look at results
estimate <- fit@output$estimate
# estimate == estimateTarget
endTime = proc.time()[1]; (endTime -start) # print elapsed time (5 seconds on an old Intel mac laptop)

round(estimate,3)
# #compare to simulation
# estimate.mat <- rbind(round(c(estimate[1:5]^2,res[6:7]),3),round(ALL$track.changes[c('var.U','var.F','var.A','var.S','var.D','cor.spouses','var.cur.phenotype'),'data.t1'],3))
# dimnames(estimate.mat) <- list(c('OpenMx-Estimated','Simulated'),c('Var.E','Var.F','Var.A','Var.S','Var.D','Cor.Sps','Var.Phen'))
# 
#NOTE: the minor difference between simulated & estimated parameters are to be expected & are due to sampling error
# > res.mat
#                  Var.E   Var.F  Var.A  Var.S   Var.D  Cor.Sps Var.Phen
#OpenMx-Estimated 0.443   0.100   0.294   0.149     0   0.209     0.986
#Old.mx           0.442   0.099   0.294   0.149     0   0.206     0.986
#Simulated        0.433   0.100   0.304   0.197     0   0.200     1.032