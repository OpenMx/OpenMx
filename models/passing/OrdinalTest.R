#
#   Copyright 2007-2010 The OpenMx Project
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
data <- read.table("data/mddndzf.dat", na.string=".", 
	col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))
data[,1] <- mxFactor(data[,1], levels = c(0 : nthresh2))
data[,2] <- mxFactor(data[,2], levels = c(0 : nthresh1))
data[,3] <- mxFactor(data[,3], levels = c(0 : nthresh2))
data[,4] <- mxFactor(data[,4], levels = c(0 : nthresh1))

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
model <- mxModel()
model <- mxModel(model, mxMatrix("Stand", name = "R", # values=c(.2955, .1268, -.0011, .0760, .1869, .4377), 
            nrow = nvar, ncol = nvar, free=TRUE))
model <- mxModel(model, mxMatrix("Zero", name = "M", nrow = 1, ncol = nvar, free=FALSE))
model <- mxModel(model, mxMatrix("Full", 
            name="thresh", 
            # values = Mx1Threshold,
            values=cbind(
                    seq(-1.9, 1.9, length.out=nthresh2),          # t1Neur1: 12 thresholds evenly spaced from -1.9 to 1.9
                   c(rep(1, nthresh1), rep(0, diff)),               # t1mddd4l: 1 threshold at 1
                    seq(-1.9, 1.9, length.out=nthresh2),          # t2Neur1: 12 thresholds same as t1Neur1
                   c(rep(1, nthresh1), rep(0, diff))                # t2mddd4l: 1 threshold same as t1mddd4l
                    ),
            free = c(rep(c( rep(TRUE, nthresh2), 
                            rep(TRUE, nthresh1), rep(FALSE, diff)
                            ), 2)), 
            dimnames = list(c(), nameList), 
            labels = rep(c(paste("neur", 1:nthresh2, sep=""),
                        paste("mddd4l", 1:nthresh1, sep=""), rep(NA, diff))
                        )))

# Define the objective function
objective <- mxFIMLObjective(covariance="R", means="M", dimnames=nameList, thresholds="thresh")

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix)

# Run the job
model <- mxRun(model)

estimates <- model@output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(thresh, model)[,1], Mx1Threshold[,1], 0.01)
omxCheckCloseEnough(mxEval(thresh, model)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, model), Mx1R, 0.01)
omxCheckCloseEnough(model@output$Minus2LogLikelihood, 4081.48, 0.02)
