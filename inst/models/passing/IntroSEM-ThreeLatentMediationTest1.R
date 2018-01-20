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

# model 2
threeLatentMediation1 <- mxModel(threeLatentOrthoRaw1,
    mxPath(from="F2",to="F3",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b23"),
    name="threeLatentMediation1"
    )

threeLatentMediation1Out <- mxRun(threeLatentMediation1, suppressWarnings=TRUE)

summary(threeLatentMediation1Out)

# model 3
threeLatentMediation2 <- mxModel(threeLatentMediation1,
    mxPath(from="F1",to="F3",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b13"),
    name="threeLatentMediation2"
    )

threeLatentMediation2Out <- mxRun(threeLatentMediation2, suppressWarnings=TRUE)

summary(threeLatentMediation2Out)

# model 4
threeLatentMediation3 <- mxModel(threeLatentMediation2,
    mxPath(from="F2",to="F1",
           arrows=1, 
           free=TRUE, values=.2,
           labels="b12"),
    name="threeLatentMediation3"
    )

threeLatentMediation3Out <- mxRun(threeLatentMediation3, suppressWarnings=TRUE)

summary(threeLatentMediation3Out)

# model 5
threeLatentMediation4 <- mxModel(threeLatentMediation3,
    mxPath(from="F2",to="F3",
           arrows=1, 
           free=FALSE, values=0),
    name="threeLatentMediation4"
    )

threeLatentMediation4Out <- mxRun(threeLatentMediation4, suppressWarnings=TRUE)

summary(threeLatentMediation4Out)

#---------------------------------------
# check values: threeLatentMediation4Out

expectVal <- c(0.899897, 1.211936, 1.447487, 0.474833, 0.701115, 
1.296055, 1.138717, 1.004529, 0.905956, 0.891751, 0.847099, 1.038272, 
1.03872, 0.820213, 0.945435, 0.836217, 0.973649, 0.982129, 1.049803, 
1.331303, 1.038699, 0.767626, 0.966071, 0.643407, 1.089664, 0.746656, 
0.067508, 0.12139, 0.088573, -0.034884, -0.09419, 0.012755, -0.067872, 
-0.059441, -0.070524, -0.049503, -0.049853, -0.098781)

expectSE <- c(0.076951, 0.087951, 0.102676, 0.077557, 0.0881, 0.119881, 0.111364, 
0.111189, 0.115273, 0.113165, 0.110221, 0.122962, 0.118599, 0.114646, 
0.145942, 0.105904, 0.106794, 0.141494, 0.133961, 0.172149, 0.133256, 
0.109461, 0.121744, 0.125662, 0.185381, 0.157077, 0.117532, 0.110355, 
0.129547, 0.151052, 0.097896, 0.086741, 0.118288, 0.11071, 0.110982, 
0.099194, 0.091335, 0.094306)

expectMin <- 7710.62

omxCheckCloseEnough(expectVal, threeLatentMediation4Out$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation4Out$output[['standardErrors']]), 0.002)

omxCheckCloseEnough(expectMin, threeLatentMediation4Out$output$minimum, 0.002)

#---------------------------------------
# check values: threeLatentMediation3Out

expectVal <- c(0.899884, 1.211972, 1.447474, 0.481912, 0.700912, 
1.295785, 1.138487, 1.004631, -0.010669, 0.906008, 0.891848, 
0.847203, 1.038428, 1.038887, 0.820293, 0.94583, 0.835866, 0.973786, 
0.982301, 1.049919, 1.331468, 1.038727, 0.767558, 0.965987, 0.642672, 
1.090014, 0.745861, 0.067507, 0.12139, 0.088573, -0.034884, -0.09419, 
0.012755, -0.067872, -0.059441, -0.070524, -0.049504, -0.049854, 
-0.098781)

expectSE <- c(0.076937, 0.087925, 0.102634, 0.127972, 0.088123, 0.119893, 
0.111359, 0.111145, 0.152816, 0.115262, 0.113151, 0.110214, 0.122903, 
0.11856, 0.11459, 0.145961, 0.105997, 0.106805, 0.141516, 0.133879, 
0.17206, 0.133238, 0.109431, 0.121728, 0.125978, 0.185416, 0.15734, 
0.117811, 0.110551, 0.129828, 0.151418, 0.098074, 0.086838, 0.118514, 
0.110886, 0.11105, 0.099275, 0.091429, 0.094377)

expectMin <- 7710.615

omxCheckCloseEnough(expectVal, threeLatentMediation3Out$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation3Out$output[['standardErrors']]), 0.002)

omxCheckCloseEnough(expectMin, threeLatentMediation3Out$output$minimum, 0.002)


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

omxCheckCloseEnough(expectVal, threeLatentMediation2Out$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation2Out$output[['standardErrors']]), 0.002)

omxCheckCloseEnough(expectMin, threeLatentMediation2Out$output$minimum, 0.002)


#---------------------------------------
# check values: threeLatentMediation1Out

expectVal <- c(0.901892, 1.20158, 1.427048, 0.704502, 1.355613, 
1.1605, 0.465948, 0.94444, 0.924281, 0.870142, 1.013963, 1.01268, 
0.828678, 0.998326, 0.879945, 0.990166, 0.890404, 1.054125, 1.389673, 
1.009607, 0.750164, 0.965169, 1.767278, 1.045934, 0.854189, 0.067511, 
0.121393, 0.088577, -0.034878, -0.094187, 0.012757, -0.067869, 
-0.059439, -0.070522, -0.049502, -0.049852, -0.09878)

expectSE <- c(0.076491, 0.08846, 0.103591, 0.092754, 0.133509, 0.118961, 
0.099896, 0.121764, 0.120148, 0.115831, 0.125683, 0.119873, 0.123774, 
0.162982, 0.115216, 0.110034, 0.154831, 0.143633, 0.176605, 0.134085, 
0.111037, 0.12368, 0.269229, 0.186339, 0.185038, 0.117895, 0.110663, 
0.129977, 0.151571, 0.098181, 0.086905, 0.118655, 0.111025, 0.111196, 
0.099395, 0.091536, 0.094484)

expectMin <- 7867.323

omxCheckCloseEnough(expectVal, threeLatentMediation1Out$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentMediation1Out$output[['standardErrors']]), 0.002)

omxCheckCloseEnough(expectMin, threeLatentMediation1Out$output$minimum, 0.002)


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

omxCheckCloseEnough(expectVal, threeLatentOrthoRaw1Out$output$estimate, 0.002)

omxCheckCloseEnough(expectSE, 
    as.vector(threeLatentOrthoRaw1Out$output[['standardErrors']]), 0.002)

omxCheckCloseEnough(expectMin, threeLatentOrthoRaw1Out$output$minimum, 0.002)

