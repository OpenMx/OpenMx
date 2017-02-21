# Program: ExtremeMultiformity_Mdd_Can.R  
# Models of Comorbidity for Multifactorial Disorders; Michael C. Neale & Kenneth S. Kendler
# A. J. Hum. Genet. 57:935-953, 1995
#
# This is a useful regression test because at starting values
# omxAllInt returns probabilities that are exactly zero.
#
# ****************
# Extreme Multiformity
# ****************

library(OpenMx)

nv                <- 1
MZtot		<- 565
DZtot		<- 641

# --------------------------------------------------------------------------------
# Fit Extreme Multiformity Model with Cell Frequencies Input, ONE overall Threshold
# --------------------------------------------------------------------------------

# Specify objects to store a, c, e path coefficients
pathAf1 <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.6, label="a11", name="af1")
pathCf1	<- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.4, label="c11", name="cf1")
pathEf1	<- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.6, label="e11", name="ef1")

pathAf2 <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.6, label="a12", name="af2")
pathCf2	<- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.4, label="c12", name="cf2")
pathEf2	<- mxMatrix( type="Full", nrow=nv, ncol=nv, free=T, values=.6, label="e12", name="ef2")

# Matrices A, C, and E compute variance components
covAf1	<- mxAlgebra( expression=a11 %*% t(a11), name="Af1")
covCf1	<- mxAlgebra( expression=c11 %*% t(c11), name="Cf1")
covEf1	<- mxAlgebra( expression=e11 %*% t(e11), name="Ef1")
covVf1	<- mxAlgebra( expression=Af1+Cf1+Ef1, name="Vf1")

covAf2  	<- mxAlgebra( expression=a12 %*% t(a12), name="Af2")
covCf2	<- mxAlgebra( expression=c12 %*% t(c12), name="Cf2")
covEf2	<- mxAlgebra( expression=e12 %*% t(e12), name="Ef2")
covVf2	<- mxAlgebra( expression=Af2+Cf2+Ef2, name="Vf2")

# Constraint on variance of ordinal variables
matI	<- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
conf1	<- mxConstraint( expression=diag2vec(Vf1)==I, name="Varf1")
conf2	<- mxConstraint( expression=diag2vec(Vf2)==I, name="Varf2")


# Matrix & Algebra for expected means vector
Meanf1	<- mxMatrix( type="Zero", nrow=1, ncol=2, name="expMeanf1" )
Meanf2  <- mxMatrix( type="Zero", nrow=1, ncol=2, name="expMeanf2" )

# Algebra for expected variance/covariance matrix in MZ & DZ pairs
covMZf1	<- mxAlgebra( expression= rbind  ( cbind(Af1+Cf1+Ef1 , Af1+Cf1)		,cbind(Af1+Cf1   , Af1+Cf1+Ef1)), name="expCovMZf1" )
covDZf1	<- mxAlgebra( expression= rbind  ( cbind(Af1+Cf1+Ef1 , 0.5%x%Af1+Cf1)	,cbind(0.5%x%Af1+Cf1 , Af1+Cf1+Ef1)),  name="expCovDZf1" ) 

covMZf2 <- mxAlgebra( expression= rbind  ( cbind(Af2+Cf2+Ef2 , Af2+Cf2)		,cbind(Af2+Cf2   , Af2+Cf2+Ef2)), name="expCovMZf2" )
covDZf2	<- mxAlgebra( expression= rbind  ( cbind(Af2+Cf2+Ef2 , 0.5%x%Af2+Cf2)	,cbind(0.5%x%Af2+Cf2 , Af2+Cf2+Ef2)),  name="expCovDZf2" ) 

# Matrix & Algebra for expected thresholds

# Matrices for thresholds two each for Factors 1 and 2, and a matrix for means
inf   	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, values=Inf, name="infinity")
Tsf1 		<- mxMatrix( type="Full", nrow=2, ncol=1, free=TRUE, values=c(.8, .2), label=c("thr1f1","deltathr2f1"), lbound=c(-Inf,0.0001),name="P" )
Tsf2 		<- mxMatrix( type="Full", nrow=2, ncol=1, free=TRUE, values=c(.6,.3), label=c("thr1f2","deltathr2f2"), lbound=c(-Inf,0.0001),name="R" )
inc     	<- mxMatrix( type="Lower", nrow=2, ncol=2, free=FALSE, values=1, name="Incrementer")
Thresf1 	<- mxAlgebra( expression=rbind(-infinity,Incrementer%*%P,infinity),name="thresholdsFactor1")
Thresf2 	<- mxAlgebra( expression=rbind(-infinity,Incrementer%*%R,infinity),name="thresholdsFactor2")

# ??????
# ??????
incrTf1 <-mxAlgebra( cbind(thresholdsFactor1, thresholdsFactor1), name="increasingThresholdsf1")
incrTf2 <-mxAlgebra( cbind(thresholdsFactor2, thresholdsFactor2), name="increasingThresholdsf2")
AllintMZf1 <- mxAlgebra(omxAllInt(expCovMZf1, expMeanf1, increasingThresholdsf1), name="AllintMZf1")
AllintMZf2 <- mxAlgebra(omxAllInt(expCovMZf2, expMeanf2, increasingThresholdsf2), name="AllintMZf2")
AllintDZf1 <- mxAlgebra(omxAllInt(expCovDZf1, expMeanf1, increasingThresholdsf1), name="AllintDZf1")
AllintDZf2 <- mxAlgebra(omxAllInt(expCovDZf2, expMeanf2, increasingThresholdsf2), name="AllintDZf2")

# Using algebras to extract each of the 6 integrals on the diagonal and below
#MZ
CMZ <- mxAlgebra(expression=AllintMZf1[1,1], name="CMZ")
DMZ <- mxAlgebra(expression=AllintMZf1[2,1], name="DMZ")
EMZ <- mxAlgebra(expression=AllintMZf1[3,1], name="EMZ")
FMZ <- mxAlgebra(expression=AllintMZf1[5,1], name="FMZ")
GMZ <- mxAlgebra(expression=AllintMZf1[6,1], name="GMZ")
HMZ <- mxAlgebra(expression=AllintMZf1[9,1], name="HMZ")

IMZ <- mxAlgebra(expression=AllintMZf2[1,1], name="IMZ")
JMZ <- mxAlgebra(expression=AllintMZf2[2,1], name="JMZ")
KMZ <- mxAlgebra(expression=AllintMZf2[3,1], name="KMZ")
LMZ <- mxAlgebra(expression=AllintMZf2[5,1], name="LMZ")
MMZ <- mxAlgebra(expression=AllintMZf2[6,1], name="MMZ")
NMZ <- mxAlgebra(expression=AllintMZf2[9,1], name="NMZ")

#DZ
CDZ <- mxAlgebra(expression=AllintDZf1[1,1], name="CDZ")
DDZ <- mxAlgebra(expression=AllintDZf1[2,1], name="DDZ")
EDZ <- mxAlgebra(expression=AllintDZf1[3,1], name="EDZ")
FDZ <- mxAlgebra(expression=AllintDZf1[5,1], name="FDZ")
GDZ <- mxAlgebra(expression=AllintDZf1[6,1], name="GDZ")
HDZ <- mxAlgebra(expression=AllintDZf1[9,1], name="HDZ")

IDZ <- mxAlgebra(expression=AllintDZf2[1,1], name="IDZ")
JDZ <- mxAlgebra(expression=AllintDZf2[2,1], name="JDZ")
KDZ <- mxAlgebra(expression=AllintDZf2[3,1], name="KDZ")
LDZ <- mxAlgebra(expression=AllintDZf2[5,1], name="LDZ")
MDZ <- mxAlgebra(expression=AllintDZf2[6,1], name="MDZ")
NDZ <- mxAlgebra(expression=AllintDZf2[9,1], name="NDZ")

# Working out the Expected Frequencies
MZexpFr	        <- mxAlgebra(rbind(
        CMZ*IMZ ,
        2.0*(CMZ*JMZ) ,
        2.0*(DMZ*IMZ) ,
        2.0*(EMZ*(IMZ+JMZ+KMZ) + DMZ*(JMZ+KMZ) + CMZ*KMZ) ,
        CMZ*LMZ ,
        2.0*(DMZ*JMZ) , 
        2.0*(EMZ*(JMZ+LMZ+MMZ) + DMZ*(LMZ+MMZ) + CMZ*MMZ) ,
        FMZ*IMZ ,
        2.0*(GMZ*(IMZ+JMZ+KMZ) +FMZ*(JMZ+KMZ) + DMZ*KMZ) ,
        HMZ + NMZ - HMZ*NMZ + GMZ*(JMZ+KMZ+LMZ+MMZ+MMZ) + EMZ*(KMZ+MMZ) + GMZ*(JMZ+LMZ+MMZ+KMZ+MMZ) + FMZ*(LMZ+MMZ+MMZ) + 
                DMZ*MMZ + EMZ*(KMZ+MMZ) + DMZ*MMZ), name="MZExpectedProb" )

# The DZ model
# Working out the Expected Frequencies
DZexpFr         <- mxAlgebra(rbind(
        CDZ*IDZ ,
        2.0*(CDZ*JDZ) ,
        2.0*(DDZ*IDZ) ,
        2.0*(EDZ*(IDZ+JDZ+KDZ) + DDZ*(JDZ+KDZ) + CDZ*KDZ) ,
        CDZ*LDZ ,
        2.0*(DDZ*JDZ) , 
        2.0*(EDZ*(JDZ+LDZ+MDZ) + DDZ*(LDZ+MDZ) + CDZ*MDZ) ,
        FDZ*IDZ ,
        2.0*(GDZ*(IDZ+JDZ+KDZ) +FDZ*(JDZ+KDZ) + DDZ*KDZ) ,
        HDZ + NDZ - HDZ*NDZ + GDZ*(JDZ+KDZ+LDZ+MDZ+MDZ) + EDZ*(KDZ+MDZ) + GDZ*(JDZ+LDZ+MDZ+KDZ+MDZ) + FDZ*(LDZ+MDZ+MDZ) + 
                DDZ*MDZ + EDZ*(KDZ+MDZ) + DDZ*MDZ), name="DZExpectedProb" )

Munit   	<- mxMatrix( type="Unit", nrow=10, ncol=1, name="mU10")
MZOfreq	<- mxAlgebra( MZExpectedProb * (mU10 %x% MZtot), name="MZExpectedFrequencies" )
DZOfreq	<- mxAlgebra( DZExpectedProb * (mU10 %x% DZtot), name="DZExpectedFrequencies" )

# DATA: cell Frequency inputs 
dataMZ	<- mxMatrix(type="Full", nrow=10, ncol=1, free=FALSE, values=c(273,96,46,23,39,9,12,23,25,19), name="MZObservedFrequencies")
dataDZ	<- mxMatrix(type="Full", nrow=10, ncol=1, free=FALSE, values=c(251,109,74,52,26,29,18,35,35,12), name="DZObservedFrequencies")

# Specify the User Defined Fitfunctions for each group: likelihood
#objMZ	<- mxAlgebra( -2 * sum(MZObservedFrequencies * log(MZExpectedFrequencies)),name="MZobj")		
#objDZ	<- mxAlgebra( -2 * sum(DZObservedFrequencies * log(DZExpectedFrequencies)),name="DZobj")		

# Specify the User Defined Fitfunctions for each group: Chi-square
objMZ	<- mxAlgebra(sum ( (  (MZObservedFrequencies - MZExpectedFrequencies) * (MZObservedFrequencies - MZExpectedFrequencies))/MZExpectedFrequencies ), name="MZobj")
objDZ	<- mxAlgebra(sum ( (  (DZObservedFrequencies - DZExpectedFrequencies) * (DZObservedFrequencies - DZExpectedFrequencies))/DZExpectedFrequencies ), name="DZobj")


# Combine Groups
pars		<- list(pathAf1, pathCf1, pathEf1, pathAf2, pathCf2, pathEf2, covAf1, covCf1, covEf1, covVf1, covAf2, covCf2, covEf2, covVf2, matI,  Meanf1, Meanf2, Tsf1,Tsf2,inc, Thresf1, Thresf2,incrTf1,incrTf2,inf, Munit)
modelMZ	<- mxModel(pars, covMZf1, covMZf2, AllintMZf1,AllintMZf2,CMZ,DMZ,EMZ,FMZ,GMZ,HMZ,IMZ,JMZ,KMZ,LMZ,MMZ,NMZ,MZexpFr, dataMZ, objMZ, MZOfreq, conf1, conf2,name="MZ")
modelDZ	<- mxModel(pars, covDZf1, covDZf2, AllintDZf1,AllintDZf2,CDZ,DDZ,EDZ,FDZ,GDZ,HDZ,IDZ,JDZ,KDZ,LDZ,MDZ,NDZ,DZexpFr, dataDZ, objDZ, DZOfreq, name="DZ")
minus2ll	<- mxAlgebra(expression=MZ.MZobj + DZ.DZobj, name="ChiSq")
obj		<- mxFitFunctionAlgebra("ChiSq")
ConfACE         <- mxCI (c ('MZ.Af1[1,1]', 'MZ.Cf1[1,1]', 'MZ.Ef1[1,1]',
                            'MZ.Af2[1,1]', 'MZ.Cf2[1,1]', 'MZ.Ef2[1,1]'))
ExtremeMultiformityModel	<-mxModel("ExtremeMultiformity", pars, modelMZ, modelDZ, minus2ll, obj,ConfACE) 

#ExtremeMultiformityRun<-mxTryHardOrdinal(ExtremeMultiformityModel, greenOK = TRUE, checkHess = FALSE, finetuneGradient=FALSE, exhaustive=TRUE, OKstatuscodes=c(0,1,5,6), wtgcsv=c("prev"),intervals=F )

ExtremeMultiformityFit		<- mxRun(ExtremeMultiformityModel, intervals=T)

(ExtremeMultiformitySumm		<- summary(ExtremeMultiformityFit))

omxCheckCloseEnough(mxEval(MZ.Af1, ExtremeMultiformityFit), 0.5434, .05)
omxCheckCloseEnough(mxEval(MZ.Af2, ExtremeMultiformityFit), 0.4425, .05)
omxCheckCloseEnough(mxEval(MZ.Cf1, ExtremeMultiformityFit), 0.2537, .05)
omxCheckCloseEnough(mxEval(MZ.Cf2, ExtremeMultiformityFit), 0, 1e-2)
omxCheckCloseEnough(mxEval(MZ.Ef1, ExtremeMultiformityFit), 0.2028, .05)
omxCheckCloseEnough(mxEval(MZ.Ef2, ExtremeMultiformityFit), 0.5574, .05)

# ******************************************************
# Get Chi Square
(chisquare	<-mxEval(objective,ExtremeMultiformityFit))

# Work out degrees of freedom 
pTable	<-data.frame(ExtremeMultiformitySumm$parameters)
(np		<-nrow(pTable))
(df		<-20-np)

# Work out p-value
#  Note: pchisq(q,df) is a distribution function in R
(pvalue	<-pchisq(chisquare,df,lower.tail=FALSE))


# -------------------------------------------------------------------------------------------------------------------------------------
# Sub1 model: delta put to extremely high value, Extreme Multiformity of first disorder only
# --------------------------------------------------------------------------------------------------------------------------------------


sub1Model	<- mxModel(ExtremeMultiformityFit, name="sub1")
sub1Model	<- omxSetParameters(sub1Model, labels="deltathr2f2", free=F, values=6)
sub1Fit	<- mxRun(sub1Model, intervals=T)
(sub1Summ	<- summary(sub1Fit))

# -------------------------------------------------------------------------------------------------------------------------------------

# ******************************************************
# Get Chi Square
(chisquare	<-mxEval(objective,sub1Fit))

# Work out degrees of freedom 
pTable	<-data.frame(sub1Summ$parameters)
(np		<-nrow(pTable))
(df		<-20-np)

# Work out p-value
#  Note: pchisq(q,df) is a distribution function in R
(pvalue	<-pchisq(chisquare,df,lower.tail=FALSE))



# -------------------------------------------------------------------------------------------------------------------------------------
# Sub2 model: delta put to extremely high value, Extreme Multiformity of second disorder only
# --------------------------------------------------------------------------------------------------------------------------------------


sub2Model	<- mxModel(ExtremeMultiformityFit, name="sub2")
sub2Model	<- omxSetParameters(sub2Model, labels="deltathr2f1", free=F, values=6)
sub2Fit	<- mxRun(sub2Model,intervals=T)
(sub2Summ	<- summary(sub2Fit))

# -------------------------------------------------------------------------------------------------------------------------------------

# ******************************************************
# Get Chi Square
(chisquare	<-mxEval(objective,sub2Fit))
omxCheckCloseEnough(chisquare, 37.3105, 1e-2)

# Work out degrees of freedom 
pTable	<-data.frame(sub2Summ$parameters)
(np		<-nrow(pTable))
(df		<-20-np)

# Work out p-value
#  Note: pchisq(q,df) is a distribution function in R
(pvalue	<-pchisq(chisquare,df,lower.tail=FALSE))
omxCheckCloseEnough(log(pvalue), -9.188, 1e-2)


# -------------------------------------------------------------------------------------------------------------------------------------
# Test that math is correct by checking whether probabilties add up to 1
# --------------------------------------------------------------------------------------------------------------------------------------


oneTableDZ <- data.frame(ExtremeMultiformityFit$output$algebras$DZ.DZExpectedProb)
oneTableMZ <- data.frame(ExtremeMultiformityFit$output$algebras$MZ.MZExpectedProb)

(sumDZprob <- colSums (oneTableDZ, na.rm = FALSE, dims = 1))
(sumMZprob <- colSums (oneTableMZ, na.rm = FALSE, dims = 1))

omxCheckCloseEnough(sumDZprob, 1, 1e-5)
omxCheckCloseEnough(sumMZprob, 1, 1e-5)
