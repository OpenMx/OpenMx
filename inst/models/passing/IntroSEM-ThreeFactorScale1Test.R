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
# Program: ThreeFactorScale1Test.R
#  Author: Steven M. Boker
#    Date: Sun Mar 14 14:42:36 EDT 2010
#
# This program tests the number of factors in simulated data.
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Mar 14 14:42:39 EDT 2010
#      Created ThreeFactorScale1Test.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorScaleExample1)

numberFactors <- 3
indicators <- names(factorScaleExample1)
numberIndicators <- length(indicators)
totalVars <- numberIndicators + numberFactors

# ----------------------------------
# Build an Old-style RAM OpenMx single factor FIML model with fixed variance

latents <- paste("F", 1:numberFactors, sep="")
loadingLabels <- paste("b_F", rep(1:numberFactors, each=numberIndicators), rep(indicators, numberFactors), sep="") 
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

threeFactorOrthoRaw1 <- mxModel("threeFactorOrthogonal",
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
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=factorScaleExample1, type="raw")
    )

threeFactorOrthoRaw1Out <- mxRun(threeFactorOrthoRaw1)

summary(threeFactorOrthoRaw1Out)


threeFactorObliqueRaw1 <- mxModel(threeFactorOrthoRaw1,
    #mxPath(from=latents,to=latents,all=TRUE,
	mxPath(from=latents,to=latents,connect="unique.pairs",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorOblique"
    )

threeFactorObliqueRaw1Out <- mxRun(threeFactorObliqueRaw1)

summary(threeFactorObliqueRaw1Out)

threeFactorObliqueRaw2 <- mxModel(threeFactorOrthoRaw1,
    mxPath(from="F1",to="F2",
           arrows=2, 
           free=FALSE, values=1),
    mxPath(from="F1",to="F3",
           arrows=2, 
           free=TRUE, values=.3,labels="C1"),
    mxPath(from="F2",to="F3",
           arrows=2, 
           free=TRUE, values=.3,labels="C1"),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorOblique2"
    )

threeFactorObliqueRaw2Out <- mxRun(threeFactorObliqueRaw2)

summary(threeFactorObliqueRaw2Out)

threeFactorObliqueRaw3 <- mxModel(threeFactorOrthoRaw1,
    mxPath(from=latents, to=latents,
#           arrows=2, all=TRUE,
           arrows=2, connect="unique.pairs",
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorAllOne"
    )

threeFactorObliqueRaw3Out <- mxRun(threeFactorObliqueRaw3)

summary(threeFactorObliqueRaw3Out)

#------------------
# Check values: threeFactorOrthoRaw1Out

expectVal <- c(1.190861, 1.027769, 1.254965, 1.363502, 0.841032, 
0.611059, 1.176141, 0.962213, 1.001306, 0.872364, 0.732371, 0.781316, 
0.879178, 0.948706, 0.802852, 0.957057, 0.829901, 0.966251, 1.16657, 
1.046025, 1.0449, 0.930938, 0.886489, 0.857332, -0.015984, 0.044187, 
-0.045358, -0.045237, 0.042802, -0.015518, 0.025968, -0.029904, 
-0.09063, 0.066277, -0.03439, -0.102583)

expectSE <- c(0.095071, 0.091627, 0.095687, 0.104192, 0.089938, 0.087914, 
0.115636, 0.103374, 0.105377, 0.096345, 0.088913, 0.089604, 0.118024, 
0.115144, 0.118351, 0.1404, 0.114143, 0.110775, 0.190218, 0.150236, 
0.156371, 0.130233, 0.111504, 0.112966, 0.107235, 0.100174, 0.109089, 
0.118727, 0.087664, 0.081834, 0.112898, 0.099287, 0.101177, 0.09197, 
0.084346, 0.085667)

expectMin <- 7637.69

omxCheckCloseEnough(expectVal, threeFactorOrthoRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeFactorOrthoRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeFactorOrthoRaw1Out$output$minimum, 0.001)


#------------------
# Check values: threeFactorObliqueRaw1Out

expectVal <- c(1.199074, 1.024535, 1.237762, 1.375742, 0.842897, 
0.603621, 1.188015, 0.954562, 0.979856, 0.872316, 0.735397, 0.797505, 
0.85955, 0.955344, 0.845734, 0.92353, 0.82676, 0.975287, 1.138494, 
1.060691, 1.087397, 0.931022, 0.882049, 0.831771, 0.536646, 0.370669, 
0.260529, -0.015986, 0.044185, -0.045361, -0.04524, 0.042801, 
-0.015519, 0.025965, -0.029905, -0.09063, 0.066277, -0.03439, 
-0.102583)

expectSE <- c(0.094315, 0.091308, 0.095556, 0.103031, 0.08778, 0.086482, 
0.111739, 0.100458, 0.104029, 0.094893, 0.088243, 0.088802, 0.11479, 
0.114309, 0.117201, 0.134742, 0.109258, 0.109441, 0.17717, 0.142814, 
0.152609, 0.12695, 0.110274, 0.111154, 0.066984, 0.078801, 0.088475, 
0.107079, 0.100044, 0.10893, 0.118549, 0.087645, 0.081832, 0.112886, 
0.099263, 0.101146, 0.091946, 0.084319, 0.085644)

expectMin <- 7575.204

omxCheckCloseEnough(expectVal, threeFactorObliqueRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeFactorObliqueRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw1Out$output$minimum, 0.001)




#------------------
# Check values: threeFactorObliqueRaw2Out

expectVal <- c(1.195198, 1.004327, 1.214682, 1.348694, 0.54545, 
0.37662, 0.778512, 0.603039, 0.979279, 0.87233, 0.735412, 0.797955, 
0.868832, 0.996346, 0.90234, 0.997225, 1.239722, 1.197804, 1.9438, 
1.608228, 1.088527, 0.930999, 0.882027, 0.831054, 0.379883, -0.015985, 
0.044185, -0.04536, -0.045239, 0.042801, -0.015519, 0.025966, 
-0.029904, -0.090631, 0.066277, -0.034391, -0.102583)

expectSE <- c(0.094148, 0.091872, 0.095946, 0.103538, 0.089709, 0.085627, 
0.114201, 0.102037, 0.104034, 0.094897, 0.088254, 0.088808, 0.113821, 
0.117047, 0.11859, 0.136206, 0.129157, 0.12209, 0.205295, 0.167271, 
0.152618, 0.126959, 0.110283, 0.111156, 0.077844, 0.107162, 0.100119, 
0.109026, 0.118651, 0.087666, 0.081832, 0.112911, 0.099296, 0.101198, 
0.091996, 0.084361, 0.085684)

expectMin <- 7685.782

omxCheckCloseEnough(expectVal, threeFactorObliqueRaw2Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeFactorObliqueRaw2Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw2Out$output$minimum, 0.001)

#------------------
# Check values: threeFactorObliqueRaw3Out

expectVal <- c(1.181873, 1.00225, 1.186484, 1.354415, 0.545864, 
0.380693, 0.773758, 0.608262, 0.380386, 0.393442, 0.361358, 0.430867, 
0.900508, 1.000514, 0.970048, 0.981762, 1.23927, 1.19472, 1.951179, 
1.601901, 1.902825, 1.537163, 1.29228, 1.282141, -0.015985, 0.044185, 
-0.045361, -0.04524, 0.042801, -0.015519, 0.025966, -0.029904, 
-0.09063, 0.066278, -0.03439, -0.102582)

expectSE <- c(0.094616, 0.091853, 0.097053, 0.1032, 0.089584, 0.085497, 0.114165, 
0.101867, 0.107596, 0.096892, 0.088672, 0.089026, 0.116041, 0.117149, 
0.123821, 0.134647, 0.129041, 0.121806, 0.205757, 0.166724, 0.192868, 
0.156338, 0.131391, 0.131211, 0.107105, 0.100067, 0.108975, 0.118589, 
0.087647, 0.081824, 0.11288, 0.099279, 0.101171, 0.09197, 0.084342, 
0.085657)

expectMin <- 7823.096

omxCheckCloseEnough(expectVal, threeFactorObliqueRaw3Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(threeFactorObliqueRaw3Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw3Out$output$minimum, 0.001)


