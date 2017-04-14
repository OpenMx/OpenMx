# ===========
# = HISTORY =
# ===========
# 2017-04-14 05:46PM Script no longer runs: 
# Error line 60: could not find function "imxTestPPML"
# Needs updating with imxPPML.Test.Battery???

# One hundred days model
# Each manifest represents a single variable as measured on a different day, days 1-100
# Three predicting latents. Over the trial, one represents a constant part,
# one represents a linear part, and one an exponential part

require(OpenMx)
require(MASS)

# Random covariance matrix
#dataTest <- matrix(rnorm(100*100), 100, 100)
#dataTest <- dataTest %*% t(dataTest)

manifests <- unlist(lapply(1:100, function(n) { paste("Day",n, sep="") } ))
latents   <- c('Const', 'Lin', 'Exp')

amtT          <- 100
amtStep       <- amtT / 99
ConstLoadings <- rep(1,100)
LinLoadings   <- unlist(lapply(1:100, function (n) { n*amtStep/(amtT-1) } ))
EXPFACTOR     <- 0.22
ExpLoadings   <- unlist(lapply(1:100, function (n) { -exp(-EXPFACTOR*n*amtStep) } ))

lambda      <- cbind(ConstLoadings, LinLoadings, ExpLoadings)
dataLatents <- matrix(rnorm(3*3), 3, 3)
dataTest    <- lambda %*% dataLatents %*% t(dataLatents) %*% t(lambda) + diag(rep(1,100))
dataTest    <- (dataTest + t(dataTest)) / 2

dataRaw <- mvrnorm(n=5000, rep(0,100), dataTest)
colnames(dataRaw) <- manifests

colnames(dataTest) <- manifests
rownames(dataTest) <- manifests

factorModel <- mxModel("One Hundred Days Model", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    # A
	  mxPath(from='Const', to = manifests, value = ConstLoadings,	free = FALSE),
	  mxPath(from='Lin', 	 to = manifests, value = LinLoadings,		free = FALSE),
	  mxPath(from='Exp',	 to = manifests, value = ExpLoadings,		free = FALSE),
    # S
	  mxPath(from=manifests, arrows=2,value=rep(1.0, 100), labels = rep("Res", 100)),
	  mxPath(from="Const", to = latents, arrows = 2, values = 1),
	  mxPath(from="Lin"  , to = latents, arrows = 2, values = 1),
	  mxPath(from="Exp"  , to = latents, arrows = 2, values = 1),
	  
	  # Means
	  mxPath(from="one", to=latents, values=0, free=TRUE),
	  
	  # Data
	  mxData(dataRaw, type="raw", numObs=5000)
)

imxTestPPML(factorModel)
