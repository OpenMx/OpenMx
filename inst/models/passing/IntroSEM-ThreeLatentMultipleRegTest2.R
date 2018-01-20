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
# Program: ThreeLatentMultipleRegTest2.R
#  Author: Steven M. Boker
#    Date: Sun Mar 14 16:28:34 EDT 2010
#
# This program tests variations on a latent variable multiple regression
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Mar 14 16:28:38 EDT 2010
#      Created ThreeLatentMultipleRegTest2.R.
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

expectVal <- c(0.899885, 1.211974, 1.447476, 0.700916, 1.295793, 
1.138494, 0.90601, 0.89185, 0.847205, 1.038429, 1.038887, 0.820293, 
0.945828, 0.835869, 0.973786, 0.9823, 1.049919, 1.331468, 1.038726, 
0.767559, 0.965987, 1.742792, 1.095046, 1.089994, 0.828187, 0.516085, 
1.139466, 0.067508, 0.121391, 0.088573, -0.034883, -0.094189, 
0.012756, -0.067871, -0.059441, -0.070524, -0.049503, -0.049853, 
-0.098781)

expectSE <- c(0.07695, 0.087954, 0.102656, 0.088136, 0.119946, 0.111449, 
0.11528, 0.113175, 0.110252, 0.122903, 0.11859, 0.114598, 0.145985, 
0.106003, 0.106815, 0.141517, 0.133942, 0.17212, 0.133215, 0.10943, 
0.121737, 0.26635, 0.162563, 0.18551, 0.158593, 0.118376, 0.235046, 
0.117896, 0.110676, 0.129977, 0.151584, 0.098098, 0.086849, 0.118565, 
0.110935, 0.11116, 0.099326, 0.091469, 0.094438)

expectMin <- 7710.615

omxCheckCloseEnough(expectVal, threeLatentObliqueRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentObliqueRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeLatentObliqueRaw1Out$output$minimum, 0.001)

#---------------------
# check values: threeLatentMultipleReg1Out

expectVal <- c(0.899885, 1.211973, 1.447476, 0.481912, 0.700916, 
1.295792, 1.138493, -0.010669, 0.906009, 0.891849, 0.847204, 
1.038428, 1.038887, 0.820293, 0.945828, 0.835868, 0.973786, 0.982301, 
1.049919, 1.331467, 1.038726, 0.767558, 0.965987, 1.742797, 1.09505, 
1.089998, 0.745861, 0.067508, 0.12139, 0.088573, -0.034883, -0.09419, 
0.012755, -0.067872, -0.059441, -0.070524, -0.049503, -0.049853, 
-0.098781)

expectSE <- c(0.0769, 0.087865, 0.102543, 0.127828, 0.088084, 0.119791, 0.111236, 
0.152625, 0.115206, 0.113084, 0.110164, 0.122883, 0.118535, 0.114589, 
0.145932, 0.10598, 0.106799, 0.141507, 0.133841, 0.171978, 0.13313, 
0.109413, 0.121718, 0.26595, 0.16229, 0.185135, 0.157228, 0.11835, 
0.111086, 0.13054, 0.152265, 0.098407, 0.087037, 0.119005, 0.111335, 
0.11135, 0.09954, 0.091682, 0.094603)

expectMin <- 7710.615

omxCheckCloseEnough(expectVal, threeLatentMultipleReg1Out$output$estimate, 0.01)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMultipleReg1Out$output[['standardErrors']]), 0.01)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg1Out$output$minimum, 0.01)

#---------------------
# check values: threeLatentMultipleReg2Out

expectVal <- c(0.90026, 1.208833, 1.450523, 0.714247, 1.294955, 
1.145001, 0.539993, 0.925794, 0.907651, 0.857963, 1.038481, 1.037751, 
0.833617, 0.930552, 0.865281, 0.968228, 1.033987, 1.072283, 1.360517, 
1.022325, 0.759087, 0.966464, 1.742743, 1.115285, 1.060584, 0.801163, 
0.067507, 0.12139, 0.088573, -0.034883, -0.09419, 0.012755, -0.067872, 
-0.059442, -0.070524, -0.049503, -0.049853, -0.098781)

expectSE <- c(0.07703, 0.088125, 0.103052, 0.089582, 0.121094, 0.113147, 
0.102949, 0.118661, 0.116812, 0.113105, 0.123476, 0.118851, 0.116296, 
0.146876, 0.106735, 0.105913, 0.142948, 0.13435, 0.17478, 0.133635, 
0.110398, 0.122877, 0.266637, 0.163601, 0.183091, 0.171665, 0.11752, 
0.110322, 0.129462, 0.150946, 0.097806, 0.086686, 0.118153, 0.110619, 
0.111, 0.099205, 0.091341, 0.09432)

expectMin <- 7726.295

omxCheckCloseEnough(expectVal, threeLatentMultipleReg2Out$output$estimate, 0.01)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMultipleReg2Out$output[['standardErrors']]), 0.01)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg2Out$output$minimum, 0.01)

#---------------------
# check values: threeLatentMultipleReg3Out

expectVal <- c(0.899897, 1.211937, 1.447488, 0.474834, 0.701118, 
1.296061, 1.138723, 0.905957, 0.891751, 0.8471, 1.038272, 1.03872, 
0.820213, 0.945435, 0.836219, 0.973649, 0.982129, 1.049802, 1.331303, 
1.038698, 0.767626, 0.96607, 1.742954, 1.094587, 1.089649, 0.746656, 
0.067506, 0.121389, 0.088571, -0.034885, -0.09419, 0.012755, 
-0.067873, -0.059442, -0.070525, -0.049504, -0.049854, -0.098782)

expectSE <- c(0.076951, 0.087964, 0.102685, 0.077557, 0.088124, 0.119953, 
0.111436, 0.115278, 0.11316, 0.110212, 0.122951, 0.118561, 0.114641, 
0.145913, 0.105921, 0.106784, 0.141498, 0.133944, 0.172141, 0.133254, 
0.109461, 0.121748, 0.266409, 0.162381, 0.185465, 0.157074, 0.117716, 
0.110489, 0.129714, 0.151265, 0.097968, 0.086774, 0.118361, 0.110787, 
0.111027, 0.09925, 0.091375, 0.094362)

expectMin <-  7710.62

omxCheckCloseEnough(expectVal, threeLatentMultipleReg3Out$output$estimate, 0.01)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMultipleReg3Out$output[['standardErrors']]), 0.01)

omxCheckCloseEnough(expectMin, threeLatentMultipleReg3Out$output$minimum, 0.01)

#---------------------
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

omxCheckCloseEnough(expectMin, threeLatentOrthoRaw1Out$output$minimum, 0.01)
