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


# ---------------------------------------------------------------------
# Program: ThreeLatentMultipleRegTest1.R
#  Author: Steven M. Boker
#    Date: Sun Mar 14 15:21:00 EDT 2010
#
# This program tests variations on a latent variable multiple regression
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Mar 14 14:42:39 EDT 2010
#      Created ThreeLatentMultipleRegTest1.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(latentMultipleRegExample1)

numberFactors <- 3
indicators <- names(latentMultipleRegExample1)
numberIndicators <- length(indicators)
totalVars <- numberIndicators + numberFactors

# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- paste("F", 1:numberFactors, sep="")
loadingLabels <- paste("b_F", rep(1:numberFactors, each=numberIndicators), 
                              rep(indicators, numberFactors), sep="") 
loadingLabels

uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

latents1 <- latents[1]
indicators1 <- indicators[1:4]
loadingLabels1 <- paste("b_F1", indicators[1:4], sep="") 
latents2 <- latents[2]
indicators2 <- indicators[5:8]
loadingLabels2 <- paste("b_F2", indicators[5:8], sep="") 
latents3 <- latents[3]
indicators3 <- indicators[9:12]
loadingLabels3 <- paste("b_F3", indicators[9:12], sep="") 

threeLatentOrthoRaw1 <- mxModel("threeLatentOrthogonal",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents1, to=indicators1, 
#           arrows=1, all=TRUE,
           arrows=1, connect="all.pairs",
           free=TRUE, values=.2, 
           labels=loadingLabels1),
    mxPath(from=latents2, to=indicators2, 
#           arrows=1, all=TRUE, 
           arrows=1, connect="all.pairs", 
           free=TRUE, values=.2, 
           labels=loadingLabels2),
    mxPath(from=latents3, to=indicators3, 
#           arrows=1, all=TRUE,
           arrows=1, connect="all.pairs", 
           free=TRUE, values=.2, 
           labels=loadingLabels3),
    mxPath(from=latents1, to=indicators1[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents2, to=indicators2[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=latents3, to=indicators3[1], 
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=latentMultipleRegExample1, type="raw")
    )

threeLatentOrthoRaw1Out <- mxRun(threeLatentOrthoRaw1, suppressWarnings=TRUE)

summary(threeLatentOrthoRaw1Out)

# model 2
threeLatentObliqueRaw1 <- mxModel(threeLatentOrthoRaw1,
#    mxPath(from=latents,to=latents,all=TRUE,
    mxPath(from=latents,to=latents,connect="unique.pairs",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=latents,
           arrows=2, 
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    name="threeLatentOblique"
    )

threeLatentObliqueRaw1Out <- mxRun(threeLatentObliqueRaw1, suppressWarnings=TRUE)

summary(threeLatentObliqueRaw1Out)

# model 3
threeLatentMultipleReg1 <- mxModel(threeLatentOrthoRaw1,
    mxPath(from="F1",to="F2",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=c("F1","F2"), to="F3",
           arrows=1, 
           free=TRUE, values=.2, 
           labels=c("b1", "b2")),
    name="threeLatentMultipleReg"
    )

threeLatentMultipleReg1Out <- mxRun(threeLatentMultipleReg1, suppressWarnings=TRUE)

summary(threeLatentMultipleReg1Out)

# model 4
threeLatentMultipleReg2 <- mxModel(threeLatentMultipleReg1,
    mxPath(from="F1",to="F3",
           arrows=1, 
           free=FALSE, values=0),
    name="threeLatentMultipleReg2"
    )

threeLatentMultipleReg2Out <- mxRun(threeLatentMultipleReg2, suppressWarnings=TRUE)

summary(threeLatentMultipleReg2Out)

# model 5
threeLatentMultipleReg3 <- mxModel(threeLatentMultipleReg1,
    mxPath(from="F2",to="F3",
           arrows=1, 
           free=FALSE, values=0),
    name="threeLatentMultipleReg3"
    )

threeLatentMultipleReg3Out <- mxRun(threeLatentMultipleReg3, suppressWarnings=TRUE)

summary(threeLatentMultipleReg3Out)

#---------------------
# check values: threeLatentObliqueRaw1Out

expectVal1 <- c(0.809023, 1.148336, 1.305359, 0.80529, 1.232041, 
1.189699, 0.87486, 0.786604, 0.70182, 1.133502, 1.067424, 1.06387, 
1.053661, 0.847609, 0.761788, 1.186356, 1.009935, 0.942482, 0.971016, 
0.880312, 0.956066, 1.92902, 0.153394, 1.316872, 1.981284, 0.753374, 
3.077671, 0.060812, 0.037374, -0.04868, -0.013194, 0.193347, 
0.220002, 0.256785, 0.171893, 0.166191, 0.23483, 0.173025, 0.157615
)

expectSE1 <- c(0.073166, 0.088502, 0.096505, 0.0788, 0.113654, 0.104821, 0.056678, 
0.052296, 0.051769, 0.137988, 0.121393, 0.140895, 0.155914, 0.114258, 
0.093581, 0.164655, 0.144287, 0.134106, 0.124408, 0.109581, 0.111582, 
0.29534, 0.1298, 0.213446, 0.269646, 0.176748, 0.402252, 0.123592, 
0.107811, 0.134146, 0.147162, 0.104016, 0.089866, 0.126178, 0.119852, 
0.14156, 0.128784, 0.117827, 0.111053)

expectMin <- 7892.3

omxCheckCloseEnough(expectVal1, threeLatentObliqueRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE1, 
    as.vector(threeLatentObliqueRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentObliqueRaw1Out$output$minimum, 0.001)

#---------------------
# check values: threeLatentMultipleReg1Out

expectVal2 <- c(0.809021, 1.148334, 1.305356, 0.990778, 0.805289, 
1.232038, 1.189698, 0.456683, 0.87486, 0.786604, 0.701819, 1.133502, 
1.067425, 1.063866, 1.053663, 0.847609, 0.761789, 1.186354, 1.009936, 
0.942481, 0.971017, 0.880311, 0.956065, 1.92904, 0.153399, 1.316878, 
0.770605, 0.060813, 0.037374, -0.048679, -0.013194, 0.193348, 
0.220002, 0.256785, 0.171892, 0.166191, 0.23483, 0.173025, 0.157615
)

expectSE2 <- c(0.073152, 0.088483, 0.096479, 0.088473, 0.078783, 0.113624, 
0.104798, 0.084247, 0.056673, 0.052292, 0.051766, 0.138045, 0.121331, 
0.140846, 0.155822, 0.114265, 0.093585, 0.164694, 0.144279, 0.134095, 
0.124405, 0.109585, 0.111583, 0.295218, 0.129883, 0.213411, 0.151414, 
0.123699, 0.107898, 0.134262, 0.147251, 0.103957, 0.089827, 0.126105, 
0.119782, 0.141617, 0.128844, 0.117876, 0.111064)

expectMin <- 7892.3

omxCheckCloseEnough(expectVal2, threeLatentMultipleReg1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE2, 
    as.vector(threeLatentMultipleReg1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg1Out$output$minimum, 0.001)

#---------------------
# check values: threeLatentMultipleReg2Out

expectVal3 <- c(0.823312, 1.178891, 1.348899, 0.791461, 1.220494, 
1.188543, 0.594946, 0.869191, 0.794743, 0.699582, 1.212983, 1.076301, 
1.037146, 0.975338, 0.834923, 0.782921, 1.204745, 0.995637, 0.943946, 
1.002546, 0.841619, 0.966428, 1.849562, 0.227649, 1.32956, 2.605626, 
0.060816, 0.037378, -0.048675, -0.013189, 0.193348, 0.220002, 
0.256785, 0.171892, 0.166194, 0.234833, 0.173028, 0.157617)

expectSE3 <- c(0.077248, 0.094604, 0.104575, 0.078149, 0.113009, 0.104635, 
0.122078, 0.058502, 0.053706, 0.053269, 0.149475, 0.125404, 0.150462, 
0.169893, 0.114201, 0.095033, 0.166428, 0.144272, 0.146543, 0.13493, 
0.113343, 0.116675, 0.294892, 0.130169, 0.214737, 0.359998, 0.12363, 
0.10784, 0.134161, 0.147166, 0.103842, 0.089742, 0.125972, 0.119647, 
0.141509, 0.128744, 0.11779, 0.111004)

expectMin <- 8053.265

omxCheckCloseEnough(expectVal3, threeLatentMultipleReg2Out$output$estimate, 0.001)

omxCheckWithinPercentError(expectSE3, 
    as.vector(threeLatentMultipleReg2Out$output[['standardErrors']]), 1)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg2Out$output$minimum, 0.001)

#---------------------
# check values: threeLatentMultipleReg3Out

expectVal4 <- c(0.814491, 1.150163, 1.298783, 1.045046, 0.791358, 
1.216789, 1.198978, 0.865098, 0.774729, 0.693622, 1.138792, 1.053805, 
1.062754, 1.095613, 0.835546, 0.783526, 1.217672, 0.963405, 0.888211, 
0.982675, 0.904797, 0.965161, 1.92375, 0.246242, 1.328938, 1.030995, 
0.060817, 0.037378, -0.048675, -0.013189, 0.193347, 0.220002, 
0.256786, 0.171892, 0.166197, 0.234836, 0.17303, 0.157619)

expectSE4 <- c(0.073471, 0.089084, 0.097044, 0.094182, 0.078301, 0.113157, 
0.105603, 0.056038, 0.051684, 0.051288, 0.138987, 0.120673, 0.142134, 
0.160282, 0.115057, 0.095646, 0.16789, 0.144607, 0.133183, 0.12695, 
0.11242, 0.113455, 0.295381, 0.131575, 0.214744, 0.1797, 0.123552, 
0.107784, 0.134083, 0.147057, 0.103971, 0.089832, 0.126128, 0.1198, 
0.141564, 0.128796, 0.117815, 0.111048)

expectMin <-  7923.5

omxCheckCloseEnough(expectVal4, threeLatentMultipleReg3Out$output$estimate, 0.001)

omxCheckCloseEnough(c((c(threeLatentMultipleReg3Out$output[['standardErrors']]) - expectSE4)/expectSE4),
                    rep(0, length(expectSE4)), 0.01)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg3Out$output$minimum, 0.001)

#---------------------
# check values: threeLatentOrthoRaw1Out

expectVal5 <- c(0.822342, 1.178877, 1.353015, 0.794595, 1.218425, 
1.201487, 0.862384, 0.785135, 0.693446, 1.216073, 1.081345, 1.041502, 
0.960428, 0.840164, 0.779621, 1.219233, 0.962068, 0.902475, 1.007963, 
0.862748, 0.97278, 1.846472, 1.32432, 3.117711, 0.060816, 0.037378, 
-0.048674, -0.013187, 0.19335, 0.220003, 0.256787, 0.171894, 
0.166197, 0.234836, 0.17303, 0.157619)

expectSE5 <- c(0.077294, 0.094581, 0.104766, 0.078581, 0.113482, 0.106016, 
0.057994, 0.053189, 0.0529, 0.14973, 0.125644, 0.150833, 0.169725, 
0.115323, 0.095546, 0.167946, 0.144947, 0.146358, 0.135671, 0.114815, 
0.117761, 0.294174, 0.214511, 0.408378, 0.123742, 0.10793, 0.134322, 
0.147337, 0.104042, 0.089895, 0.126216, 0.119883, 0.141757, 0.128956, 
0.117976, 0.111161)

expectMin <- 8079.447

omxCheckCloseEnough(expectVal5, threeLatentOrthoRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE5, 
    as.vector(threeLatentOrthoRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentOrthoRaw1Out$output$minimum, 0.001)


