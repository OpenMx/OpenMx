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

require(OpenMx)

#Ordinal Data test, based on poly3dz.mx

# Data
nthresh1 <- 1
nthresh2 <- 12	
data <- suppressWarnings(try(read.table("models/passing/data/mddndzf.dat", na.string=".", 
	col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))))
if (is(data, "try-error")) data <- read.table("data/mddndzf.dat", na.string=".", 
	col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))
data[,c(1,3)] <- mxFactor(data[,c(1,3)], c(0 : nthresh2))
data[,c(2,4)] <- mxFactor(data[,c(2,4)], c(0 : nthresh1))

diff <- nthresh2 - nthresh1
nvar <- 4

Mx1Threshold <- rbind(
c(-1.9209, 0.3935, -1.9209, 0.3935),
c(-0.5880, 0    , -0.5880, 0    ),
c(-0.0612, 0    , -0.0612, 0    ),
c( 0.3239, 0    ,  0.3239, 0    ),
c( 0.6936, 0    ,  0.6936, 0    ),
c( 0.8856, 0    ,  0.8856, 0    ),
c( 1.0995, 0    ,  1.0995, 0    ),
c( 1.3637, 0    ,  1.3637, 0    ),
c( 1.5031, 0    ,  1.5031, 0    ),
c( 1.7498, 0    ,  1.7498, 0    ),
c( 2.0733, 0    ,  2.0733, 0    ),
c( 2.3768, 0    ,  2.3768, 0    ))

Mx1R <- rbind(
    c(1.0000,  0.2955,  0.1268, 0.0760),
    c(0.2955,  1.0000, -0.0011, 0.1869),
    c(0.1268, -0.0011,  1.0000, 0.4377),
    c(0.0760,  0.1869,  0.4377, 1.0000))

nameList <- names(data)
# Define the model
model <- mxModel("ThresholdTest", type="RAM", manifestVars=nameList)
# Variances
model <- mxModel(model, mxPath(from=nameList, connect="unique.bivariate", arrows=2, free=TRUE))
# Covariances
model <- mxModel(model, mxPath(from=nameList, connect="single", arrows=2, values=1, free=FALSE))
# Zero Means
model <- mxModel(model, mxPath(from="one", to=nameList, values=0, free=FALSE))
# Thresholds 
model <- mxModel(model, mxThreshold(nameList[c(1,3)], nThresh=nthresh2, 
                                    values=omxNormalQuantiles(nthresh2), 
                                    free=TRUE, labels=paste("neur", 1:nthresh2, sep="")))
model <- mxModel(model, mxThreshold(nameList[c(2,4)], nThresh=nthresh1, 
                                    values=1, free=TRUE, 
                                    labels=paste("mddd4l", 1:nthresh1, sep="")))

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the data to the model
model <- mxModel(model, dataMatrix)

# Run the job
modelOut <- mxRun(model)

estimates <- modelOut$output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(Thresholds, modelOut)[,1], Mx1Threshold[,1], 0.05)
omxCheckCloseEnough(mxEval(Thresholds, modelOut)[1,3], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(S, modelOut), Mx1R, 0.01)  # Note that when F == I and A == 0, S == R
omxCheckCloseEnough(modelOut$output$Minus2LogLikelihood, 4081.48, 0.2)
