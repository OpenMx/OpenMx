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
# Program: ThreeFactorScale2Test.R
#  Author: Steven M. Boker
#    Date: Sun Mar 14 14:42:16 EDT 2010
#
# This program tests the number of factors in simulated data.
#    using a standard RAM.
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Mar 14 14:42:18 EDT 2010
#      Created ThreeFactorScale2Test.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorScaleExample2)

numberFactors <- 3
indicators <- names(factorScaleExample2)
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
    mxData(observed=factorScaleExample2, type="raw")
    )

threeFactorOrthoRaw1Out <- mxRun(threeFactorOrthoRaw1, suppressWarnings=TRUE)
# NPSOL can nail this one
# omxCheckTrue(any(summary(threeFactorOrthoRaw1Out)[['seSuspect']]))

summary(threeFactorOrthoRaw1Out)

threeFactorObliqueRaw1 <- mxModel(threeFactorOrthoRaw1,
#    mxPath(from=latents,to=latents,all=TRUE,
	 mxPath(from=latents,to=latents,connect="unique.pairs",
           arrows=2, 
           free=TRUE, values=.3),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorOblique"
    )

threeFactorObliqueRaw1Out <- mxRun(threeFactorObliqueRaw1, suppressWarnings=TRUE)
if (mxOption(NULL, 'Default optimizer') != "CSOLNP") {
        omxCheckTrue(any(summary(threeFactorObliqueRaw1Out)[['seSuspect']]))
}

summary(threeFactorObliqueRaw1Out)

threeFactorObliqueRaw2 <- mxModel(threeFactorOrthoRaw1,
    mxPath(from="F2",to="F3",
           arrows=2, 
           free=FALSE, values=1),
    mxPath(from="F1",to="F3",
           arrows=2, 
           free=TRUE, values=.3,labels="C1"),
    mxPath(from="F1",to="F2",
           arrows=2, 
           free=TRUE, values=.3,labels="C1"),
    mxPath(from=latents,
           arrows=2, 
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorOblique2"
    )

threeFactorObliqueRaw2Out <- mxRun(threeFactorObliqueRaw2, suppressWarnings=TRUE)

summary(threeFactorObliqueRaw2Out)

threeFactorObliqueRaw3 <- mxModel(threeFactorOrthoRaw1,
    mxPath(from=latents, to=latents,
#           arrows=2, all=TRUE,
           arrows=2, connect="unique.pairs",
           free=FALSE, values=1, 
           labels=factorVarLabels),
    name="threeFactorAllOne"
    )

threeFactorObliqueRaw3Out <- mxRun(threeFactorObliqueRaw3, suppressWarnings=TRUE)

summary(threeFactorObliqueRaw3Out)


#------------------
# Check values: threeFactorOrthoRaw1Out

expectVal1 <- c(0.703921, 0.642147, 0.95265, 0.985615, 1.218976, 
0.916221, 1.838198, 1.427716, 1.958099, 1.870192, 1.463947, 1.167302, 
1.009401, 1.02327, 0.823179, 0.84571, 1.137994, 1.021853, 0.79516, 
0.96619, 1.066088, 1.172633, 1.09237, 1.294113, -0.06738, 0.019242, 
-0.025738, 0.130504, -0.109013, -0.113767, -0.117153, -0.037544, 
-0.246421, -0.043396, -0.088277, -0.024991)

expectSE1 <- c(0.092156, 0.090827, 0.097117, 0.099476, 0.102465, 0.089168, 
0.118992, 0.104917, 0.127505, 0.12581, 0.108334, 0.103572, 0.12067, 
0.11806, 0.133559, 0.140494, 0.139471, 0.113733, 0.172681, 0.134718, 
0.178728, 0.175763, 0.140762, 0.146253, 0.086742, 0.084724, 0.093026, 
0.095317, 0.114564, 0.096491, 0.144494, 0.122593, 0.156622, 0.152894, 
0.127247, 0.115312)

expectMin <- 8042.699

omxCheckCloseEnough(expectVal1, threeFactorOrthoRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectSE1, 
    as.vector(threeFactorOrthoRaw1Out$output[['standardErrors']]), 0.001)

omxCheckCloseEnough(expectMin, threeFactorOrthoRaw1Out$output$minimum, 0.001)


#------------------
# Check values: threeFactorObliqueRaw1Out

expectVal2 <- c(0.71589, 0.652159, 0.955592, 0.96467, 1.240036, 0.968256, 
1.792033, 1.409752, 1.972189, 1.857734, 1.457215, 1.172104, 0.992406, 
1.010312, 0.817566, 0.886559, 1.086193, 0.923789, 0.962738, 1.017162, 
1.010695, 1.21909, 1.112044, 1.28287, 0.338329, 0.340458, 1.002844, 
-0.06738, 0.019242, -0.025739, 0.130503, -0.109015, -0.113769, 
-0.117156, -0.037547, -0.246422, -0.043399, -0.088279, -0.024991
)

expectSE2 <- c(0.091768, 0.090312, 0.095406, 0.097918, 0.098236, 0.08515, 
0.115792, 0.10241, 0.123501, 0.12343, 0.106157, 0.101291, 0.119642, 
0.116976, 0.129086, 0.135608, 0.118679, 0.098006, 0.128246, 0.117336, 
0.132244, 0.147665, 0.125426, 0.136747, 0.079361, 0.077818, 0.013114, 
0.086761, 0.084739, 0.093047, 0.095348, 0.114681, 0.09655, 0.14468, 
0.122689, 0.156731, 0.15303, 0.127346, 0.115371)

expectMin <- 7699.145
print(max(abs(as.vector(threeFactorObliqueRaw2Out$output[['standardErrors']]))))
omxCheckCloseEnough(expectVal2, threeFactorObliqueRaw1Out$output$estimate, 0.001)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw1Out$output$minimum, 0.001)

omxCheckCloseEnough(c((threeFactorObliqueRaw1Out$output[['standardErrors']] - expectSE2)/expectSE2),
                    rep(0,length(expectSE2)), 0.01)

#------------------
# Check values: threeFactorObliqueRaw2Out

expectVal3 <- c(0.716165, 0.652282, 0.955552, 0.964369, 1.24112, 
0.967984, 1.794554, 1.411694, 1.973517, 1.859021, 1.45814, 1.172644, 
0.992014, 1.01015, 0.817645, 0.887139, 1.083522, 0.92432, 0.953714, 
1.011686, 1.005492, 1.214315, 1.109348, 1.281618, 0.339259, -0.067381, 
0.019241, -0.02574, 0.130502, -0.109017, -0.11377, -0.117159, 
-0.037549, -0.246428, -0.043402, -0.088282, -0.024994)

expectSE3 <- c(0.09167, 0.090289, 0.095411, 0.097799, 0.098073, 0.085161, 
0.114884, 0.101905, 0.123235, 0.123115, 0.106003, 0.101214, 0.119486, 
0.116941, 0.129103, 0.13529, 0.118112, 0.098217, 0.11909, 0.114324, 
0.130268, 0.1455, 0.124857, 0.136625, 0.076153, 0.086673, 0.084666, 
0.092912, 0.095211, 0.114218, 0.096262, 0.143926, 0.122151, 0.15595, 
0.152269, 0.12678, 0.114907)

expectMin <- 7699.195

omxCheckCloseEnough(expectVal3, threeFactorObliqueRaw2Out$output$estimate, 0.001)

omxCheckCloseEnough(c((threeFactorObliqueRaw2Out$output[['standardErrors']] - expectSE3)/expectSE3),
                           rep(0, length(expectSE3)), 0.015)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw2Out$output$minimum, 0.001)

#------------------
# Check values: threeFactorObliqueRaw3Out

expectVal4 <- c(0.323129, 0.29047, 0.34297, 0.282883, 1.243405, 0.96688, 
1.792621, 1.411438, 1.971626, 1.86013, 1.456035, 1.169947, 1.400495, 
1.351251, 1.613098, 1.73713, 1.077842, 0.926457, 0.96065, 1.01241, 
1.012951, 1.210191, 1.11548, 1.287938, -0.067381, 0.019242, -0.025739, 
0.130503, -0.109017, -0.11377, -0.117158, -0.037548, -0.246427, 
-0.043401, -0.088282, -0.024994)

expectSE4 <- c(0.087632, 0.085805, 0.094044, 0.097024, 0.098053, 0.08524, 
0.115056, 0.10198, 0.123396, 0.123167, 0.106145, 0.101346, 0.140644, 
0.13563, 0.162005, 0.174188, 0.117533, 0.098386, 0.119488, 0.114265, 
0.130545, 0.144951, 0.125295, 0.137162, 0.086782, 0.084757, 0.093061, 
0.095339, 0.114828, 0.09669, 0.144951, 0.122938, 0.157055, 0.153328, 
0.127547, 0.115553)

expectMin <- 7837.312

omxCheckCloseEnough(expectVal4, threeFactorObliqueRaw3Out$output$estimate, 0.005)

omxCheckCloseEnough(expectSE4, 
    as.vector(threeFactorObliqueRaw3Out$output[['standardErrors']]), 0.005)

omxCheckCloseEnough(expectMin, threeFactorObliqueRaw3Out$output$minimum, 0.005)




