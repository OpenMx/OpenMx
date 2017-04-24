# ===========
# = HISTORY =
# ===========
# 2017-04-14 05:16PM TBATES
# No longer runs, needs data shifted to mxFactor() format.

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

# Ordinal Data test, based on poly3dz.mx

# Data

# data <- read.table("~/bin/OpenMx/inst/models/passing/data/mddndzf.dat", na.string=".", col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))
data <- read.table("../passing/data/mddndzf.dat", na.string=".", col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))
mxFactor()
nthresh1 <- 1
nthresh2 <- 12
diff     <- nthresh2-nthresh1
nvar     <- 4

Mx1Threshold <- rbind(
c(-1.9209, 0.3935, -1.9209, 0.3935),
c(-0.5880, NA    , -0.5880, NA    ),
c(-0.0612, NA    , -0.0612, NA    ),
c( 0.3239, NA    ,  0.3239, NA    ),
c( 0.6936, NA    ,  0.6936, NA    ),
c( 0.8856, NA    ,  0.8856, NA    ),
c( 1.0995, NA    ,  1.0995, NA    ),
c( 1.3637, NA    ,  1.3637, NA    ),
c( 1.5031, NA    ,  1.5031, NA    ),
c( 1.7498, NA    ,  1.7498, NA    ),
c( 2.0733, NA    ,  2.0733, NA    ),
c( 2.3768, NA    ,  2.3768, NA    ))

Mx1R <- rbind(
    c(1.0000,  0.2955,  0.1268, 0.0760),
    c(0.2955,  1.0000, -0.0011, 0.1869),
    c(0.1268, -0.0011,  1.0000, 0.4377),
    c(0.0760,  0.1869,  0.4377, 1.0000))

nameList <- names(data)
data$breaking <- rep(1, length(data$t1neur1))  ### <-- This line causes the model to break.


# Define the model
model <- mxModel('model')
model <- mxModel(model, mxMatrix("Stand", name = "R", nrow = nvar, ncol = nvar, free=TRUE))
model <- mxModel(model, mxMatrix("Zero", name = "M", nrow = 1, ncol = nvar, free=FALSE))
model <- mxModel(model, mxMatrix("Full", 
            name="thresh", 
            values=cbind(
                    seq(-1.9, 1.9, by=(3.8/(nthresh2-1))),  # t1Neur1: 12 thresholds evenly spaced from -1.9 to 1.9
                   c(rep(1, nthresh1),rep(NA, diff)),       # t1mddd4l: 1 threshold at 1
                    seq(-1.9, 1.9, by=(3.8/(nthresh2-1))),  # t2Neur1: 12 thresholds same as t1Neur1
                   c(rep(1, nthresh1),rep(NA, diff))        # t2mddd4l: 1 threshold same as t1mddd4l
                    ),
            free = c(rep(c( rep(TRUE, nthresh2), 
                            rep(TRUE, nthresh1), rep(FALSE, diff)
                            ), 2)), 
            dimnames = list(c(), nameList), 
            labels = rep(c(paste("neur", 1:nthresh2, sep=""),
                        paste("mddd4l", 1:nthresh1, sep=""), rep(NA, diff))
                     )
					)
)

# Add the objective function and the data to the model
model <- mxModel(model, 
	mxExpectationNormal(covariance="R", means="M", dimnames=nameList, thresholds="thresh"),
	mxFitFunctionML(),
	mxData(data, type='raw')
)

# Run the job
model <- mxRun(model)

estimates <- model$output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(thresh, model)[,1], Mx1Threshold[,1], 0.01)
omxCheckCloseEnough(mxEval(thresh, model)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, model), Mx1R, 0.01)
omxCheckCloseEnough(model$output$Minus2LogLikelihood, 4081.48, 0.02)
