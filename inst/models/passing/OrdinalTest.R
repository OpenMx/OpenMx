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

v <- 1:3
omxCheckError(mxFactor(v, levels=1:3, exclude=3), "Factor levels and exclude vector are not disjoint; both contain '3'")

v <- 1:4
omxCheckError(mxFactor(v, levels=1:3), "The following values are not mapped to factor levels and not excluded: '4'")

cf <- omxCheckError(mxFactor(sample(1:2, 10, replace=TRUE), levels=1:2,
                             labels=c("incorrect", "incorrect")),
                    "Duplicate labels and collapse=TRUE not specified: 'incorrect'")
cf <- mxFactor(sample(1:2, 10, replace=TRUE), levels=1:2,
               labels=c("incorrect", "incorrect"), collapse=TRUE)
omxCheckEquals(length(levels(cf)), 1)
omxCheckEquals(levels(cf), 'incorrect')
omxCheckTrue(all(cf == "incorrect"))

foo <- data.frame(x=c(1:3),y=c(4:6),z=c(7:9))
foo <- mxFactor(foo, c(1:9), labels=c(1,1,1,2,2,2,3,3,3), collapse=TRUE)
omxCheckTrue(all(foo == matrix(kronecker(1:3, rep(1,3)),3,3)))

v <- sample.int(50, 200, replace=TRUE)
vl <- v %% 11
mask <- !duplicated(v)
v2 <- mxFactor(v, levels=v[mask], labels=vl[mask], collapse = TRUE)
omxCheckTrue(all(v2 == vl))

#Ordinal Data test, based on poly3dz.mx

# Data
nthresh1 <- 1
nthresh2 <- 12	
cnames <- c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l")
data <- suppressWarnings(try(read.table("data/mddndzf.dat", na.string=".", col.names=cnames)))
if (is(data, "try-error")) data <- read.table("models/passing/data/mddndzf.dat", na.string=".", col.names=cnames)
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
            labels = rep(c(paste("neur", 1:nthresh2, sep=""),
                        paste("mddd4l", 1:nthresh1, sep=""), rep(NA, diff))
                        )))

# Define the objective function
objective <- mxExpectationNormal(covariance="R", means="M", dimnames=nameList, thresholds="thresh")

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix, mxFitFunctionML())

# Run the job
modelOut <- mxRun(model)

estimates <- modelOut$output$estimate

# Results from old Mx:
omxCheckCloseEnough(mxEval(thresh, modelOut)[,1], Mx1Threshold[,1], 0.03)
omxCheckCloseEnough(mxEval(thresh, modelOut)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, modelOut), Mx1R, 0.01)
omxCheckCloseEnough(modelOut$output$Minus2LogLikelihood, 4081.48, 0.08)
