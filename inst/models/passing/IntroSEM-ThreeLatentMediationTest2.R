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
# Program: ThreeLatentMediationTest2.R
#  Author: Steven M. Boker
#    Date: Sun Mar 14 16:28:34 EDT 2010
#
# This program tests variations on a latent mediation model
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Mar 14 16:28:38 EDT 2010
#      Created ThreeLatentMediationTest2.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(latentMultipleRegExample2)

numberFactors <- 3
indicators <- names(latentMultipleRegExample2)
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
    mxData(observed=latentMultipleRegExample2, type="raw")
    )

threeLatentOrthoRaw1Out <- mxRun(threeLatentOrthoRaw1, suppressWarnings=TRUE)

summary(threeLatentOrthoRaw1Out)

threeLatentMediation1 <- mxModel(threeLatentOrthoRaw1,
    mxPath(from="F1",to="F3",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b13"),
    name="threeLatentMediation1"
    )

threeLatentMediation1Out <- mxRun(threeLatentMediation1, suppressWarnings=TRUE)

summary(threeLatentMediation1Out)

threeLatentMediation2 <- mxModel(threeLatentMediation1,
    mxPath(from="F2",to="F3",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b23"),
    name="threeLatentMediation2"
    )

threeLatentMediation2Out <- mxRun(threeLatentMediation2)

summary(threeLatentMediation2Out)

threeLatentMediation3 <- mxModel(threeLatentMediation2,
    mxPath(from="F1",to="F2",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b12"),
    name="threeLatentMediation3"
    )

threeLatentMediation3Out <- mxRun(threeLatentMediation3)

summary(threeLatentMediation3Out)

threeLatentMediation4 <- mxModel(threeLatentMediation3,
    mxPath(from="F1",to="F3",
           arrows=1, 
           free=FALSE, values=0),
    name="threeLatentMediation4"
    )

threeLatentMediation4Out <- mxRun(threeLatentMediation4)

summary(threeLatentMediation4Out)

#---------------------------------------
# check values: threeLatentMediation4Out

expectVal <- c(0.900257, 1.208829, 1.450517, 0.639961, 0.714244, 
1.294949, 1.144996, 0.53999, 0.925796, 0.907652, 0.857963, 1.03848, 
1.037751, 0.833618, 0.930552, 0.865279, 0.968228, 1.033988, 1.072282, 
1.360518, 1.022325, 0.759086, 0.966465, 1.742763, 0.346848, 0.80116, 
0.067507, 0.12139, 0.088572, -0.034884, -0.094191, 0.012754, 
-0.067873, -0.059442, -0.070524, -0.049504, -0.049853, -0.098781
)

expectSE <- c(0.076992, 0.088067, 0.102978, 0.068953, 0.089564, 0.121074, 
0.113108, 0.10292, 0.118615, 0.116758, 0.113076, 0.123447, 0.118805, 
0.116293, 0.146865, 0.106722, 0.105903, 0.142954, 0.134282, 0.174763, 
0.133614, 0.110377, 0.122868, 0.266351, 0.078456, 0.171579, 0.118087, 
0.110839, 0.130227, 0.151884, 0.098249, 0.086923, 0.118773, 0.111146, 
0.111238, 0.099416, 0.091558, 0.094503)

expectMin <- 7726.295

omxCheckCloseEnough(expectVal, threeLatentMediation4Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation4Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentMediation4Out$output$minimum, 0.001)

#---------------------------------------
# check values: threeLatentMediation3Out

expectVal <- c(0.899883, 1.21197, 1.447472, 0.62833, 0.481911, 0.700915, 
1.29579, 1.138492, -0.010669, 0.906008, 0.891849, 0.847204, 1.038427, 
1.038888, 0.820293, 0.945828, 0.835867, 0.973787, 0.9823, 1.049919, 
1.331469, 1.038727, 0.767558, 0.965986, 1.742814, 0.401946, 0.745863, 
0.067508, 0.121391, 0.088574, -0.034882, -0.094188, 0.012755, 
-0.067872, -0.059441, -0.070523, -0.049504, -0.049853, -0.098781
)

expectSE <- c(0.076954, 0.087963, 0.102668, 0.069253, 0.127926, 0.088107, 
0.119887, 0.11135, 0.152795, 0.115277, 0.113172, 0.110223, 0.122928, 
0.1186, 0.114618, 0.145991, 0.105993, 0.106802, 0.141538, 0.133933, 
0.172114, 0.133217, 0.109453, 0.121739, 0.266425, 0.084583, 0.157381, 
0.117906, 0.110697, 0.130002, 0.151597, 0.098095, 0.086834, 0.118547, 
0.110953, 0.11118, 0.099373, 0.091479, 0.094449)

expectMin <- 7710.615

omxCheckCloseEnough(expectVal, threeLatentMediation3Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation3Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentMediation3Out$output$minimum, 0.001)


#---------------------------------------
# check values: threeLatentMediation2Out

expectVal <- c(0.900993, 1.207067, 1.425325, 0.43279, 0.694426, 
1.355306, 1.152128, 0.078525, 0.906157, 0.891475, 0.846991, 1.016895, 
1.017922, 0.80959, 1.01297, 0.871724, 1.000945, 0.876173, 1.063462, 
1.331154, 1.038162, 0.76807, 0.966173, 1.764348, 1.054153, 0.747071, 
0.06751, 0.121392, 0.088576, -0.034879, -0.094187, 0.012757, 
-0.067869, -0.059439, -0.070522, -0.049502, -0.049852, -0.09878)

expectSE <- c(0.076472, 0.088267, 0.102919, 0.092556, 0.092707, 0.136435, 
0.119071, 0.108581, 0.115255, 0.113106, 0.11018, 0.124685, 0.119489, 
0.120868, 0.159665, 0.117298, 0.111136, 0.160984, 0.146537, 0.172013, 
0.133155, 0.109448, 0.12173, 0.268676, 0.188841, 0.15717, 0.118021, 
0.110782, 0.13014, 0.151765, 0.098149, 0.086882, 0.118609, 0.110993, 
0.109976, 0.098261, 0.09034, 0.093452)

expectMin <- 7840.327

omxCheckCloseEnough(expectVal, threeLatentMediation2Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation2Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentMediation2Out$output$minimum, 0.001)


#---------------------------------------
# check values: threeLatentMediation1Out

expectVal <- c(0.900912, 1.207724, 1.426027, 0.472404, 0.692089, 
1.351976, 1.149559, 0.907366, 0.893115, 0.848558, 1.017988, 1.019067, 
0.80839, 1.011663, 0.86774, 1.002452, 0.878396, 1.064432, 1.333989, 
1.037999, 0.766996, 0.965187, 1.763251, 1.058137, 0.743454, 0.067509, 
0.121392, 0.088575, -0.034881, -0.094187, 0.012757, -0.067868, 
-0.059439, -0.070523, -0.049503, -0.049852, -0.098781)

expectSE <- c(0.076476, 0.088205, 0.102851, 0.077261, 0.092403, 0.136048, 
0.11869, 0.115463, 0.113329, 0.110403, 0.124548, 0.119397, 0.120363, 
0.159031, 0.117181, 0.111273, 0.161261, 0.146516, 0.172106, 0.133218, 
0.109401, 0.121713, 0.268393, 0.188989, 0.156945, 0.118051, 0.110785, 
0.130146, 0.151787, 0.098155, 0.08689, 0.118623, 0.110996, 0.111222, 
0.099412, 0.091541, 0.094504)

expectMin <- 7840.856

omxCheckCloseEnough(expectVal, threeLatentMediation1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentMediation1Out$output$minimum, 0.001)


#---------------------------------------
# check values: threeLatentOrthoRaw1Out

expectVal <- c(0.901891, 1.20158, 1.427049, 0.692088, 1.351975, 
1.149558, 0.987615, 0.966986, 0.902343, 1.013963, 1.012679, 0.828678, 
0.998327, 0.86774, 1.002452, 0.878394, 1.064433, 1.459186, 0.987214, 
0.727833, 0.960052, 1.767278, 1.058139, 1.01176, 0.067512, 0.121394, 
0.088578, -0.034877, -0.094186, 0.012757, -0.067868, -0.059437, 
-0.070521, -0.049501, -0.04985, -0.098779)

expectSE <- c(0.076528, 0.088527, 0.103665, 0.092459, 0.136218, 0.118764, 
0.130232, 0.129587, 0.123389, 0.125806, 0.119984, 0.12382, 0.163071, 
0.117211, 0.111291, 0.161268, 0.14648, 0.181275, 0.136797, 0.113505, 
0.125792, 0.269576, 0.189346, 0.226953, 0.117883, 0.110644, 0.129952, 
0.151555, 0.098128, 0.086865, 0.118584, 0.110962, 0.111143, 0.099336, 
0.091474, 0.094429)

expectMin <- 7897.082

omxCheckCloseEnough(expectVal, threeLatentOrthoRaw1Out$output$estimate, 0.01)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentOrthoRaw1Out$output[['standardErrors']]), 0.01)

omxCheckCloseEnough(expectMin, threeLatentOrthoRaw1Out$output$minimum, 0.001)

