#
#   Copyright 2007-2009 The OpenMx Project
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
data <- read.table("data/mddndzf.dat", na.string=".", col.names=c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l"))

nthresh1 <- 1
nthresh2 <- 12
diff <- nthresh2-nthresh1
nvar <- 4

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

# Define the model
model <- mxModel()
model <- mxModel(model, mxMatrix("Stand", name = "R", # values=c(.2955, .1268, -.0011, .0760, .1869, .4377), 
            nrow = nvar, ncol = nvar, free=TRUE, 
            dimnames=list(nameList, nameList)))
model <- mxModel(model, mxMatrix("Zero", name = "M", nrow = 1, ncol = nvar, free=FALSE, dimnames = list(NULL, nameList)))

# Threshold differences:
model <- mxModel(model, mxMatrix("Full", name="T", ncol = 1, nrow = nthresh1, free=T, values=c(.2, rep(.3, nthresh1-1)), labels = paste("mddd4lThreshold", 1:nthresh1, sep="")))
model <- mxModel(model, mxMatrix("Full", name="U", ncol = 1, nrow = nthresh2, free=T, values=c(-2, rep(.3, nthresh2-1)), labels=paste("Neur1Threshold", 1:nthresh2, sep="")))

# For Multiplication
model <- mxModel(model, mxMatrix("Lower", name="I1", nrow = nthresh1, ncol = nthresh1, free=F, values=1))
model <- mxModel(model, mxMatrix("Lower", name="I2", nrow = nthresh2, ncol = nthresh2, free=F, values=1))

# Algebras
model <- mxModel(model, mxAlgebra(I1%*%T, name="OneMddd4lThreshold"))
model <- mxModel(model, mxAlgebra(cbind(OneMddd4lThreshold, OneMddd4lThreshold), name ="thresh1", dimnames=list(NULL, c("t1mddd4l", "t2mddd4l"))))
model <- mxModel(model, mxAlgebra(I2%*%U, name="OneNeur1Threshold"))
model <- mxModel(model, mxAlgebra(cbind(OneNeur1Threshold, OneNeur1Threshold), name ="thresh2", dimnames=list(NULL, c("t1neur1", "t2neur1"))))

model <- mxModel(model, mxBounds(parameters=c(paste("mddd4lThreshold", 2:nthresh1, sep=""), paste("Neur1Threshold", 2:nthresh2, sep="")), min = 0, ))

# Define the objective function
objective <- mxFIMLObjective(covariance="R", means="M", thresholds=c("thresh1", "thresh2"))

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix)

# Run the job
model <- mxRun(model)

estimates <- model@output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(thresh2, model)[,1], Mx1Threshold[,1], 0.01)
omxCheckCloseEnough(mxEval(thresh1, model)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, model), Mx1R, 0.01)
omxCheckCloseEnough(model@output$Minus2LogLikelihood, 4081.48, 0.02)
