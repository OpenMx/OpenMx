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
colNames <- c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l")
data <- suppressWarnings(try(read.table("models/passing/data/mddndzf.dat", na.string=".", col.names=colNames)))
if (is(data, "try-error")) data <- read.table("data/mddndzf.dat", na.string=".", col.names=colNames)
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
model <- mxModel('model')
model <- mxModel(model, mxMatrix("Stand", name = "R", # values=c(.2955, .1268, -.0011, .0760, .1869, .4377), 
            nrow = nvar, ncol = nvar, free=TRUE, 
            dimnames=list(nameList, nameList)))
model <- mxModel(model, mxMatrix("Zero", name = "M",
			nrow = 1, ncol = nvar, free=FALSE, dimnames = list(NULL, nameList)))

# Threshold differences:
model <- mxModel(model, mxMatrix("Full", name="T", ncol = 1, nrow = nthresh1, 
		free=T, values=c(.2, rep(.3, nthresh1-1)), 
		labels = paste("mddd4lThreshold", 1:nthresh1, sep="")))
model <- mxModel(model, mxMatrix("Full", name="U", ncol = 1, nrow = nthresh2, 
		free=T, values=c(-2, rep(.3, nthresh2-1)), 
		labels=paste("Neur1Threshold", 1:nthresh2, sep="")))

# For Multiplication
model <- mxModel(model, mxMatrix("Lower", name="I1", 
			nrow = nthresh1, ncol = nthresh1, free=F, values=1))
model <- mxModel(model, mxMatrix("Lower", name="I2", 
			nrow = nthresh2, ncol = nthresh2, free=F, values=1))

# Algebras
model$OneMddd4lThreshold <- mxAlgebra(I1 %*% T)
model$thresh1 <- mxAlgebra(cbind(OneMddd4lThreshold, OneMddd4lThreshold), 
	dimnames=list(NULL, c("t1mddd4l", "t2mddd4l")))
model$OneNeur1Threshold <- mxAlgebra(I2 %*% U)
model$thresh2 <- mxAlgebra(cbind(OneNeur1Threshold, OneNeur1Threshold), 
	dimnames=list(NULL, c("t1neur1", "t2neur1")))
model$zeros <- mxMatrix('Zero', nrow = 11, ncol = 2)
model$thresholds <- mxAlgebra(cbind(rbind(thresh1, zeros), thresh2))

if (nthresh1 > 1) {
	model <- mxModel(model, mxBounds(parameters = paste("mddd4lThreshold", 2:nthresh1, sep=""), min = 0))
}
if (nthresh2 > 1) {
	model <- mxModel(model, mxBounds(parameters = paste("Neur1Threshold", 2:nthresh2, sep=""), min = 0))
}

# Define the objective function
objective <- mxExpectationNormal(covariance="R", means="M", thresholds="thresholds")

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix, mxFitFunctionML())

# Run the job
model <- mxRun(model)

estimates <- model$output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(thresh2, model)[,1], Mx1Threshold[,1], 0.045)
omxCheckCloseEnough(mxEval(thresh1, model)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, model), Mx1R, 0.01)
omxCheckCloseEnough(model$output$Minus2LogLikelihood, 4081.48, 0.2)
