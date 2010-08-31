# ---------------------------------------------------------------------
# Program: TwoFactorOrthoRaw-OpenMx100221.R
#  Author: Steven M. Boker
#    Date: Sun Feb 21 14:39:50 EST 2010
#
# This program fits a covariance two factor orthogonal FIML model to the 
#     factorExample1.csv simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 21 14:39:47 EST 2010
#      Created TwoFactorOrthoRaw-OpenMx100221.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorExample1)

# ----------------------------------
# Build an OpenMx two factor orthogonal FIML model with fixed variance

indicators <- names(factorExample1)
latents <- c("F1", "F2")
loadingLabels <- c(paste("b_F1", indicators, sep=""), paste("b_F2", indicators, sep=""))
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")
factorCovLabels <- c("Cov_F1F1", "Cov_F2F2", "Cov_F1F2")


twoFactorOrthoRaw1 <- mxModel("Two Factor Orthogonal FIML Model with Fixed Variance",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
           arrows=1, all=TRUE, 
           free=TRUE, values=.2, 
           labels=loadingLabels),
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
    mxData(observed=factorExample1, type="raw")
    )

twoFactorOrthoRaw1Out <- mxRun(twoFactorOrthoRaw1, suppressWarnings=TRUE)

summary(twoFactorOrthoRaw1Out)


# ----------------------------------
# Build an OpenMx two factor orthogonal FIML model with fixed loading


twoFactorOrthoRaw2 <- mxModel("Two Factor Orthogonal FIML Model with Fixed Loading",
    type="RAM",
    manifestVars=indicators,
    latentVars=latents,
    mxPath(from=latents, to=indicators, 
           arrows=1, all=TRUE,
           free=TRUE, values=.2, 
           labels=loadingLabels),
    mxPath(from=indicators, 
           arrows=2, 
           free=TRUE, values=.8, 
           labels=uniqueLabels),
    mxPath(from=latents,
           arrows=2,
           free=TRUE, values=.8, 
           labels=factorVarLabels),
    mxPath(from=latents, to=c("x9", "x7"),
           arrows=1, 
           free=FALSE, values=1),
    mxPath(from="one", to=indicators, 
           arrows=1, free=TRUE, values=.1, 
           labels=meanLabels),
    mxData(observed=factorExample1, type="raw")
    )

twoFactorOrthoRaw2Out <- mxRun(twoFactorOrthoRaw2, suppressWarnings=TRUE)

summary(twoFactorOrthoRaw2Out)

#---------------------
# check values: twoFactorOrthoRaw1Out

expectVal <- c(0.636559, 0.397663, 0.105249, 0.449579, 0.761353, 
1.043744, 1.182638, 0.600451, 0.883266, 0.260461, -0.201878, 
0.008219, 0.088494, -0.438513, 0.446284, 0.526686, 0.231033, 
-0.445119, 0.347544, 0.082811, 0.194243, 0.814987, 0.223275, 
0.334888, 0.250124, 0.239238, 0.392693, -0.011614, -0.006824, 
0.023961, -0.031358, -0.045483, -0.091786, -0.067325, -0.039022, 
-0.059999)

c(NaN, NaN, 0.014064, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 
NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0.024373, 0.007613, 0.012326, 
0.052528, 0.027071, 0.030672, 0.033813, 0.017571, 0.036525, 0.04046, 
0.023705, 0.020264, 0.04525, 0.044552, 0.056887, 0.061969, 0.036093, 
0.05229)

expectMin <-9174.685

omxCheckCloseEnough(expectVal, twoFactorOrthoRaw1Out@output$estimate, 0.001)

#omxCheckCloseEnough(expectSE, as.vector(twoFactorOrthoRaw1Out@output$standardError), 0.001)

omxCheckCloseEnough(expectMin, twoFactorOrthoRaw1Out@output$minimum, 0.001)

#---------------------
# check values: twoFactorOrthoRaw2Out

expectVal <- c(0.439919, 0.45093, 0.089568, 0.35681, 0.888458, 0.712037, 
0.796639, 0.422043, 0.51856, 0.012192, 0.055943, 0.284613, -0.021127, 
0.866696, 0.476616, 0.029882, 0.347544, 0.082811, 0.194243, 0.814986, 
0.223276, 0.33489, 0.250124, 0.239235, 0.392692, 0.977352, 1.055772, 
-0.011614, -0.006824, 0.023961, -0.031358, -0.045484, -0.091786, 
-0.067326, -0.039022, -0.059999)

expectSE <- c(0.846744, 0.020801, 0.091636, 0.460776, 0.089527, 1.415789, 
1.634738, 0.77726, 0.052686, 0.684284, 0.073059, 0.207455, 1.403722, 
0.051649, 0.071861, 1.514011, 0.024374, 0.007613, 0.012326, 0.052533, 
0.027073, 0.030675, 0.033816, 0.017571, 0.036528, 0.13163, 2.607359, 
0.040554, 0.023759, 0.020271, 0.045289, 0.044657, 0.057058, 0.062155, 
0.036188, 0.052422)

expectMin <-  9174.685

omxCheckCloseEnough(expectVal, twoFactorOrthoRaw2Out@output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(twoFactorOrthoRaw2Out@output$standardError), 0.001)

omxCheckCloseEnough(expectMin, twoFactorOrthoRaw2Out@output$minimum, 0.001)




