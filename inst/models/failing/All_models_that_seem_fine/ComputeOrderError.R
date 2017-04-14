# SCRIPT: NTF_design.R - NTF design in OpenMx 
# Author: Matt Keller & Sarah Medland; edited by Tim Bates
# History:  Thu Sep 24 17:23:33 BST 2009
# 	2017-04-14 04:15PM Updated by tbates to use multifgroup. but this code seem to have a bug (missing ALL object),
# 	and did not check the run model, but rather the un-run model...

# For description of model, see Keller, Medland, Duncan, Hatemi, Neale, Maes, & Eaves (2009) TRHG, 21, p.8 - 18.
# 2009-09-24: (tb): Use three groups (change MZ and DZ family algebra so that they refer to a common set of matrices in a new "NTF" group); 
# 2009-09-24: (tb): Simplify/speed algebra using Quadratic operator and pre-calculating variance components i.e., e %*% t(e) = E etc; 
# Mon Sep 28 22:12:31 BST 2009: (tb) read data from DEMO folder	
# Rhelp: http://www.statmethods.net
# OpenMx: http://www.openmx.virginia.com
##########################################
require(OpenMx);

# TODO: Add output from old MZ to verify correctness
# TODO: Verify algebra/constraint specification

#NOTE: the minor difference between simulated & estimated parameters are to be expected & are due to sampling error
# > res.mat
#                  Var.E   Var.F  Var.A  Var.S   Var.D  Cor.Sps Var.Phen
#OpenMx-Estimated 0.443   0.100   0.294   0.149     0   0.209     0.986
#Old.mx           0.442   0.099   0.294   0.149     0   0.206     0.986
#Simulated        0.433   0.100   0.304   0.197     0   0.200     1.032

#Get Data
data(nuclear_twin_design_data)
selVars <-names(nuclear_twin_design_data)[1:4]
mzData <- nuclear_twin_design_data[nuclear_twin_design_data$zyg=='mz',selVars]
dzData <- nuclear_twin_design_data[nuclear_twin_design_data$zyg=='dz',selVars]

#Fit NTF Model with RawData and Matrices Input
ntf <- mxModel(model="NTF",
   # Matrices
   # NOTE: NTF design does not allow Vs & Vf to be estimated simultaneously; for identifiability, choose either m or s to be free and the other to be fixed at 0             
   mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, values=0,  label="FamilialPath", name="m"),  #fix m=0 if you want Vf=0
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.3,  label="Sib",          name="s"),  #fix s=0 if you want Vs=0
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.7,  label="Env",          name="e"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.6,  label="AddGen",       name="a"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="Dominance",    name="d"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.1,  label="AMCopath",     name="mu"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1.2, label="VarPhen",      name="Vp1"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.2,  label="VarF",         name="x1"),  #keep this parameter free, even if Vf is fixed to 0
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.25, label="CovPhenGen",   name="delta1"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=1,   label="VarAddGen",    name="q1"),
   mxMatrix(type="Full", nrow=1, ncol=1, free=TRUE, values=.15, label="CovFA",        name="w1"),
   
   #mxAlgebra section - nonlinear constraints
   mxAlgebra(e %*% t(e), name="E"),
   mxAlgebra(d %*% t(d), name="D"),
   mxAlgebra(s %*% t(s), name="S"),
   
   mxAlgebra((a %&% q1) + x1 + 2 %x% a %*% w1 + E + D + S, name="Vp2"),
   mxAlgebra(2 %x% m %&% Vp1 + 2 %x% (m %&% (mu %&% Vp1)), name="x2"),
   
   mxAlgebra(q1 %*% a + w1, name="delta2"),
   mxAlgebra(1 + delta1 %*% mu %*% t(delta1), name="q2"),
   mxAlgebra(delta1 %*% m + delta1 %*% mu %*% Vp1 %*% t(m), name="w2"),
   
   #constraints - equating nonlinear constraints and parameters
   mxConstraint(Vp1 == Vp2,   name='VpCon'),
   mxConstraint(x1 == x2,    name='xCon'),
   mxConstraint(delta1 == delta2,name='deltaCon'),
   mxConstraint(q1   == q2,    name='qCon'),
   mxConstraint(w1   == w2,    name='wCon'),
   #mxAlgebra section - relative covariances
   mxAlgebra(a   %&% q1 + x1   + 2 %x% a %*% w1 + D + S, name="CvMz"),
   mxAlgebra(a   %&% (q1-.5) + .25 %x% d %*% t(d) + x1 + 2 %x% a %*% w1 + S, name="CvDz"),
   mxAlgebra(.5  %x% a %*% (q1 %*% a + w1) + .5 %x% a %*% (q1 %*% a + w1) %*% mu %*% Vp1 + m %*% Vp1 + m %*% Vp1 %&% mu, name="ParChild"),
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
	mxData(mzData, type="raw"),
	mxFIMLObjective(covariance="expCovMz",means="expMeanMz")
)

dzModel <- mxModel(name = "DZNTF",
	mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values= .25, label="mean", dimnames=list(NULL, selVars), name="expMeanDz"),
	# Algebra for expected variance/covariance matrix in DZF
	mxAlgebra(expression=rbind(
		cbind(NTF.Vp1,      NTF.CvDz,     NTF.ParChild, NTF.ParChild),
		cbind(NTF.CvDz,     NTF.Vp1,      NTF.ParChild, NTF.ParChild),
		cbind(NTF.ParChild, NTF.ParChild, NTF.Vp1,      NTF.CvSps),
		cbind(NTF.ParChild, NTF.ParChild, NTF.CvSps,    NTF.Vp1)
	), 
	dimnames=list(selVars,selVars),name="expCovDz"),
	mxData(dzData, type="raw"),
	mxFIMLObjective(covariance="expCovDz",means="expMeanDz")
)

model <- mxModel(model="NucTwFam", mzModel, dzModel, ntf,
	mxFitFunctionMultigroup(c("MZNTF", "DZNTF"))
	# mxAlgebra(expression=MZNTF.objective + DZNTF.objective, name="ntffit"), # modelName.objective is the automatic name for the -2LL of modelName
	# mxFitFunctionAlgebra("ntffit")
)

#Run MX
fit <- mxRun(model)
#Look at results
summary(fit)
res <- model$output$estimate
round(res,3)
#compare to simulation

# ============================================================
# = The ALL object doesn't appear to exist in this script... =
# ============================================================
res.mat <- rbind(round(c(res[1:5]^2,res[6:7]),3),
		round(ALL$track.changes[c('var.U','var.F','var.A','var.S','var.D','cor.spouses','var.cur.phenotype'),'data.t1'],3))
		
dimnames(res.mat) <- list(c('OpenMx-Estimated','Simulated'),c('Var.E','Var.F','Var.A','Var.S','Var.D','Cor.Sps','Var.Phen'))

# look at implied Covariances
round(fit$output$algebras$MZNTF.expCovMz, 3)
round(fit$output$algebras$DZNTF.expCovDz, 3)
