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


##############################################################################
# Stable Trait, Autoregressive Trait, State Model with Multiple Indicators  ##
#                                                                           ##  
# There are three indicators of the construct at each of 10 waves.  However ##
# two of the waves (waves 5 and 7) are missing, so phantom variables are    ##
# used in these years. A latent occasion factor is modeled in each year     ##
# from the three indicators.  This occasion factor is decomposed into a     ##
# stable trait factor, an autoregressive trait factor, and an occasion-     ##
# specific state factor.  Stationarity is assumed, so the stabilities and   ##
# variances are constrained to be equal across waves.                       ##
#                                                                           ##
# The script was designed to handle different numbers of waves and different##
# missing waves, so it has a lot of extra stuff to make it flexible.  But   ##
# there are probably more efficient ways of doing this.                     ##
##############################################################################

data <- suppressWarnings(try(read.csv("data/fakeSTARTSMissing.csv",header=TRUE))) #With Missing Data
if (is(data, "try-error")) data <- read.csv("models/nightly/data/fakeSTARTSMissing.csv",header=TRUE)
waves <- 10 #Total number including phantom waves
indicators <- 3
phantom <- c(5,7) #List of waves that are missing
stabilitySV <- .85 #Starting value for stability
require(OpenMx)
  
#SETUP FOR PHANTOM VARIABLES
indicatorPrefix <- c(paste("ind",1:indicators,sep=""))
allOccasions <- c(paste("occ",1:waves,sep=""))
if(phantom[1]=="0") obsOccasions <- allOccasions else
  obsOccasions <- allOccasions[-phantom]
nActualWaves <- length(obsOccasions)
x <- matrix(obsOccasions,nrow=indicators,ncol=nActualWaves,byrow=T) 
indicatorNames <- c(paste(x,indicatorPrefix,sep=""))
occNamesForPaths <- as.vector(x)
allStates <- c(paste("state",1:waves,sep=""))
if(phantom[1]=="0") obsStates <- allStates else
  obsStates <- allStates[-phantom]
residualNames <- paste(indicatorNames,"res",sep="")
residualLabels <- paste(indicatorNames,"residVar",sep="")

#DEFINE MANIFEST AND LATENT VARIABLES
#Rename variables
names(data) <- indicatorNames
manifests <- names(data)
stableTraits <- c("ST")
arTraits <- c(paste("AR",1:waves,sep=""))
if(phantom[1]=="0") obsArTraits <- arTraits else
  obsArTraits <- arTraits[-phantom]

#DATA STATEMENT
startsData <- mxData(observed=data, type="raw")

#CALCULATE START VALUES
##Start value for stability is set by user (.85 as default)
##Start values for the ST and first AR variances are set to be 1/3 of average variance of the indicators
##Start values for the S and residual variances are set to be 1/9 and 2/9 of the average variance
##Start value for the AR residual is set to AR1-AR1*stability^2.

#Calculate average variance for indicators
sumVar=0
for (i in 1:(nActualWaves*indicators))  {
  sumVar <- sumVar+var(data, use="pairwise.complete.obs")[i,i]
}
avgVar <- sumVar/(nActualWaves*indicators)

##Set starting values
residSV <- 2*avgVar/9
stateSV <- avgVar/9
stSV <- avgVar/3
ar1SV <- avgVar/3
arSV <- ar1SV-ar1SV*stabilitySV*stabilitySV

##Set starting values for means to the average item mean
meanSV <- mean(colMeans(data, na.rm=TRUE))

#PATH STATEMENTS
##Variances
residVar <- mxPath(from=residualNames,
                   arrows=2,
                   #labels=residualLabels, #use if equality constraints aren't wanted
                   labels=c("resid1Var","resid2Var","resid3Var"), #constrain to be equal
                   free=TRUE,
                   values=residSV,
                   lbound=0)
stVar <- mxPath(from=stableTraits,
                arrows=2,
                labels="stableTrait",
                free=TRUE,
                values=stSV,
                lbound=0)
ar1Var <- mxPath(from=arTraits[1],
                 arrows=2,
                 labels="ar1Trait",
                 free=TRUE,
                 values=ar1SV) 
arVar <- mxPath(from=arTraits[2:waves],
                arrows=2,
                labels="arTrait",
                free=TRUE,
                values=arSV)
occasionVar <- mxPath(from=obsOccasions,
                      arrows=2,
                      labels="occTrait",
                      free=FALSE,
                      values=0)  #All occasion variance is ST, AR, or S
stateVar <- mxPath(from=obsStates,
                   arrows=2,
                   labels="States",
                   free=TRUE,
                   values=stateSV,
                   lbound=0) 

##Factor Loadings
stLoadings <- mxPath(from=stableTraits,
                     to=obsOccasions,
                     arrows=1,
                     free=FALSE,
                     values=1)
arLoadings <- mxPath(from=obsArTraits,
                     to=obsOccasions,
                     arrows=1,
                     free=FALSE,
                     values=1)
stateLoadings <- mxPath(from=obsStates,
                        to=obsOccasions,
                        arrows=1,
                        free=FALSE,
                        values=1)
residualLoadings <- mxPath(from=residualNames,
                           to=indicatorNames,
                           arrows=1,
                           free=FALSE,
                           values=1)
first <- (3*(c(1:nActualWaves))-2) #Indicator loadings to set to "1"
indicatorLoadings1 <- mxPath(from=occNamesForPaths[first],
                            to=indicatorNames[first],
                            arrows=1,
                            free=FALSE,
                            values=1)
indicatorLoadings <- mxPath(from=occNamesForPaths[-first],
                            to=indicatorNames[-first],
                            arrows=1,
                            labels=c("indLoading1","indLoading2"),
                            free=TRUE,
                            values=1)

 
##Stabilities
stabilityPath <- mxPath(from=arTraits[1:(waves-1)],
                   to=arTraits[2:waves],
                   arrows=1,
                   free=T,
                   values=stabilitySV,
                   labels="stability")

#Means
manifestMeans <- mxPath(from="one",
                        to=manifests,
                        arrows=1,
                        free=TRUE,
                        values=meanSV,
                        #labels=paste("mean",1:(nActualWaves*indicators),sep="")
                        labels=c("indicator1Mean","indicator2Mean","indicator3Mean")
                        )
stMeans <- mxPath(from="one", to=stableTraits, arrows=1, free=FALSE, values=0)
arMeans <- mxPath(from="one", to=arTraits, arrows=1, free=FALSE, values=0)
occMeans <- mxPath(from="one", to=obsOccasions, arrows=1, free=FALSE, values=0)
stateMeans <- mxPath(from="one", to=obsStates, arrows=1, free=FALSE, values=0)
residualMeans <- mxPath(from="one",to=residualNames,arrows=1,free=FALSE,values=0)
  
#Constraints:  Stationarity Constraint for Autoregressive Component
matrixAr <- mxMatrix("Full",
                     nrow=1,
                     ncol=1,
                     free=T,
                     values=arSV,
                     labels="arTrait",
                     name="cMatAr")
matrixAr1 <- mxMatrix("Full",
                      nrow=1,
                      ncol=1,
                      free=T,
                      values=ar1SV,
                      labels="ar1Trait",
                      name="cMatAr1")
matrixStability <- mxMatrix("Full",
                            nrow=1,
                            ncol=1,
                            free=T,
                            values=stabilitySV,
                            labels="stability",
                            name="cMatStab")
algebraSTARTS <- mxAlgebra(cMatAr1-cMatAr1*cMatStab*cMatStab, name="AR2Variance")
stationarityConstraint <- mxConstraint(cMatAr == AR2Variance ,name="Stationarity")

#CREATE PATTERN OF CORRELATED RESIDUALS FOR SAME INDICATOR AT DIFFERENT WAVE
#FIRST COLUMN IS FIRST RESIDUAL; SECOND COLUMN IS SECOND RESIDUAL
residCorrMatrix <- matrix(NA,1,2)
tempMatrix <- matrix(NA,1,2)
for (i in 1:indicators) {
  for (j in 1:(nActualWaves-1)) {
    for (k in (j+1):nActualWaves) {
      if(is.na(residCorrMatrix[1,1])) {
        residCorrMatrix[1,1] <- (i+(j-1)*indicators)
        residCorrMatrix[1,2] <- (i+(j-1)*indicators)+(k-j)*indicators
      } else {
        tempMatrix[1,1] <- (i+(j-1)*indicators)
        tempMatrix[1,2] <- (i+(j-1)*indicators)+(k-j)*indicators
        residCorrMatrix <- rbind(residCorrMatrix,tempMatrix)
      }
    }
  }
}

#PATH STATEMENTS FOR CORRELATED RESIDUALS
correlatedResiduals <- mxPath(from=residualNames[residCorrMatrix[,1]],
                              to=residualNames[residCorrMatrix[,2]],
                              arrows=2,
                              free=TRUE,
                              values=.05)
                                

#MODEL
STARTSM <- mxModel("STARTS",
                  type="RAM",
                  startsData,
                  manifestVars=manifests,
                  latentVars=c(stableTraits,arTraits,obsOccasions,obsStates,residualNames),
                  residVar,
                  stVar,
                  ar1Var,
                  arVar,
                  occasionVar,
                  stateVar,
                  stLoadings,
                  arLoadings,
                  stateLoadings,
                  indicatorLoadings1,
                  indicatorLoadings,
                  residualLoadings,
                  stabilityPath,
                  manifestMeans,
                  stMeans,
                  arMeans,
                  occMeans,
                  stateMeans,
                  residualMeans,
                  matrixAr,
                  matrixAr1,
                  matrixStability,
                  algebraSTARTS,
                  stationarityConstraint,
                  correlatedResiduals,
		  mxFitFunctionML(rowwiseParallel=FALSE))
startsModel <- mxRun(STARTSM)
omxCheckCloseEnough(startsModel$output$fit, 2718.410, .1)

if (.Platform$OS.type != 'windows' && detectCores() > 1) {
	omxCheckTrue(startsModel$compute$steps[['GD']]$output$maxThreads > 1)
}
