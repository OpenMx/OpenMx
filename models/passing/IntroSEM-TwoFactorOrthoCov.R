# ---------------------------------------------------------------------
# Program: TwoFactorOrthoCov-OpenMx100221.R
#  Author: Steven M. Boker
#    Date: Sun Feb 21 14:00:09 EST 2010
#
# This program fits a covariance single factor model to the 
#     factorExample1.csv simulated data.
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sun Feb 21 14:00:12 EST 2010
#      Created TwoFactorOrthoCov-OpenMx100221.R.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read the data and print descriptive statistics.

data(factorExample1)

# ----------------------------------
# Build an OpenMx two factor orthogonal covariance model with fixed variance

indicators <- names(factorExample1)
latents <- c("F1", "F2")
loadingLabels <- c(paste("b_F1", indicators, sep=""), paste("b_F2", indicators, sep=""))
uniqueLabels <- paste("U_", indicators, sep="")
meanLabels <- paste("M_", indicators, sep="")
factorVarLabels <- paste("Var_", latents, sep="")

twoFactorOrthoCov1 <- mxModel("Two Factor Orthogonal Covariance Model with Fixed Variance",
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
    mxData(observed=cov(factorExample1), type="cov", numObs=500)
    )

twoFactorOrthoCov1Out <- mxRun(twoFactorOrthoCov1, suppressWarnings=TRUE)

summary(twoFactorOrthoCov1Out)


# ----------------------------------
# Build an OpenMx two factor orthogonal covariance model with fixed loading


twoFactorOrthoCov2 <- mxModel("Two Factor Orthogonal Covariance Model with Fixed Loading",
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
           free=TRUE, values=1, 
           labels=factorVarLabels),
    mxPath(from=latents, to=c("x1", "x1"),
           arrows=1, 
           free=FALSE, values=1),
    mxData(observed=cov(factorExample1), type="cov", numObs=500)
    )

twoFactorOrthoCov2Out <- mxRun(twoFactorOrthoCov2, suppressWarnings=TRUE)

summary(twoFactorOrthoCov2Out)

#---------------------
# check values: twoFactorOrthoCov1Out

expectVal <- c(0.089591, 0.372883, 0.045101, 0.146253, 0.759115, 
0.130195, 0.129255, 0.097246, 0.825376, 0.68262, 0.245451, 0.095567, 
0.434722, 0.444117, 1.128806, 1.289451, 0.636624, 0.54681, 0.34824, 
0.082977, 0.194632, 0.81662, 0.223723, 0.33556, 0.250626, 0.239717, 
0.393481)

expectSE <-c(4.808345, 1.729088, 0.67363, 3.062727, 3.128629, 7.951259, 
9.082674, 4.484453, 3.851983, 0.631535, 2.626168, 0.318282, 1.030758, 
5.346316, 0.917303, 0.910826, 0.685215, 5.812958, 0.024447, 0.007636, 
0.012363, 0.052694, 0.027153, 0.030767, 0.033916, 0.017625, 0.036636
)

expectMin <- 911.4205

omxCheckCloseEnough(expectVal, twoFactorOrthoCov1Out@output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(twoFactorOrthoCov1Out@output$standardError), 0.001)

omxCheckCloseEnough(expectMin, twoFactorOrthoCov1Out@output$minimum, 0.001)

#---------------------
# check values: twoFactorOrthoCov2Out

expectVal <- c(0.490707, 0.28549, 1.485808, 0.880563, 1.640353, 
1.839156, 1.143585, 1.078612, 0.373978, -0.249955, -1.639102, 
0.779338, 1.663724, 1.936191, 0.36078, 0.859068, 0.349497, 0.181883, 
0.169804, -0.08047, 0.653052, 0.344398, 0.282166, 0.193141, 0.880428, 
0.339492, 0.133249)

expectSE <- c(0.28907, 1.309731, 7.643344, 0.260513, 0.103146, 0.251903, 
1.916432, 0.546751, 0.731734, 3.334381, 19.449218, 0.648025, 
0.197105, 0.615349, 4.877527, 1.377921, 0.024518, 0.012159, 0.012296, 
0.19822, 0.043725, 0.030121, 0.032962, 0.017952, 0.058869, 1.661253, 
1.658388)

expectMin <-  1227.254

omxCheckCloseEnough(expectVal, twoFactorOrthoCov2Out@output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
    as.vector(twoFactorOrthoCov2Out@output$standardError), 0.001)

omxCheckCloseEnough(expectMin, twoFactorOrthoCov2Out@output$minimum, 0.001)




