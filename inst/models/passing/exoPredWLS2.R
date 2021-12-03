# ==================================================================================
# = This is a 1-factor model of 2 continuous variables, and 3 factors              =
# = Two definition variables are created (coded with a "c" suffix ?)               =
# = In RAM models, you add a definition variable as a fake latent with no variance =
# = with its mean set to "data.defnVar"                                            =
# ==================================================================================

# TODO: Change the file name: "exoPredWLS.R" is not discoverable

library(OpenMx)

# Make the results repeatable (unclear if this is needed?)
set.seed(1)

data("jointdata", package ="OpenMx", verbose= TRUE)

# ==========================================
# = Here's what our input data looks like: =
# ==========================================
str(jointdata)

# Make z1c = z1 + some noise;
jointdata$z1c <- with(jointdata, z1 * .1 + rnorm(nrow(jointdata)))
# and z2c which is just a column of noise with a mean = the mean factor level of z2 (no clue why as yet...)
jointdata$z2c <- with(jointdata, rnorm(nrow(jointdata), mean=unclass(z2)*.2))

# =============================================
# = Build a function to create WLS RAM models =
# =============================================
# note: This makes the script more complex, but will allow some programatic generation below

buildModel <- function(manifests= paste0('z', 1:5), type = 'WLS', slopes= TRUE) {
  exoFree <- NULL
  if (slopes) {
    exoFree <- matrix(TRUE, length(manifests), 2)
    exoFree[5,] <- FALSE
    exoFree[4,1] <- FALSE
  }
	jr <- mxModel("JointRAM", type= "RAM",
		manifestVars = manifests,
		# TODO: is this legit (to use a var name found in the data as a latent name?)
		latentVars = c('G','z1c','z2c'),
		mxPath('one', c('z1', 'z3'), free= TRUE, labels= c('z1','z3')),
		mxPath(paste0('z', c(1,3)), arrows= 2, free= TRUE, values= .5),
		# ordinal variables have fixed variance and no mean
		# TODO: why var = .5?
		mxPath(paste0('z', c(2,4,5)), arrows= 2, free= FALSE, values= .5),
		# latent scaled to var of 1 (mean is fixed at zero by default)
		mxPath('G', arrows= 2, values= 1, free= FALSE),
		mxPath('G', to = manifests, free= TRUE, values= 1, lbound= 0),
		mxThreshold('z2', 1, free=TRUE, labels="z2_thresh1"),
		mxThreshold('z4', 3, free=TRUE, labels=paste0("z4_thresh",1:3)),
		mxThreshold('z5', 2, free=TRUE, labels=paste0("z5_thresh",1:2)),
		# Note: Data are raw, despite our using WLS
		mxData(jointdata, type = "raw", verbose= 0L, exoFree=exoFree),
		# Note: this call to mxFitFunctionWLS is all that's
		# required to make a model into WLS!
		mxFitFunctionWLS(type)
	)
	# TODO: "slopes" is means?
	if(slopes){
		jr <- mxModel(jr,
			mxPath('one', to = 'z1c', free= FALSE, labels= "data.z1c"),
			mxPath('one', to = 'z2c', free= FALSE, labels= "data.z2c"),
			mxPath('z1c', to = 'z1', labels= "r1"),
			mxPath('z2c', to = 'z2', labels= "r2"))
	}
	# TODO: Shouldn't we have a call to mxExpectationRAM or LISREL??? here???
	# mxExpectationRAM(M = NA, dimnames = NA, thresholds = "T", threshnames = ???)
	return(jr)
}

jointRAM1 <- buildModel()
jointRAM1 <- mxRun(jointRAM1)

omxCheckCloseEnough(jointRAM1$data$observedStats$z1.vcov,
                    vcov(lm(z1 ~ z1c + z2c, jointdata)))
# glm binomial is set up a little different
# omxCheckCloseEnough(jointRAM1$data$observedStats$z2.vcov,
#                     vcov(glm(z2 ~ z1c + z2c, jointdata, family = binomial(link="probit"))),
#                     .01)
omxCheckCloseEnough(jointRAM1$data$observedStats$z3.vcov,
                    vcov(lm(z3 ~ z1c + z2c, jointdata)))
perm <- c(2,3,4,1)
omxCheckCloseEnough(jointRAM1$data$observedStats$z4.vcov,
                    vcov(MASS::polr(z4 ~ z2c, jointdata, method="probit"))[perm,perm], 1e-7)
omxCheckCloseEnough(jointRAM1$data$observedStats$z5.vcov,
                    vcov(MASS::polr(z5 ~ 1, jointdata, method="probit")), 1e-8)

summary(jointRAM1)
# plot(jointRAM1) # (if using umx)

# Where do these come from?
#cat(deparse(round(unname(coef(jointRAM1)),4)))
param <-  c(0.5929, 0.5773, 0.661, 0.6089, 0.1722, 0.0851, 0.0274, 0.5475,
            0.4579, 7.9062, 2.0767, 0.0563, -0.3714, 0.1242, 0.7867, -0.6529,
            -0.2965)
omxCheckCloseEnough(coef(jointRAM1), param, 1e-3)

#cat(deparse(round(unname(c(jointRAM1$output$standardErrors)),4)))
param.se <- c(0.2448, 0.406, 0.2713, 0.374, 0.2614, 0.225, 0.2713, 0.2268,
              0.2583, 0.3, 0.2389, 0.2945, 0.3113, 0.2897, 0.3679, 0.2598,  0.233)
omxCheckCloseEnough(c(jointRAM1$output$standardErrors), param.se, 1e-3)

# ===============================================================
# = Switch jointRAM1 from MxFitFunctionWLS to an ML fitfunction =
# = to compare the estimates of these two approaches            =
# ===============================================================

jointRAM2 <- mxModel(jointRAM1, mxFitFunctionML())
jointRAM2 <- mxRun(jointRAM2)
summary(jointRAM2)

estCmp <- cbind(coef(jointRAM2), coef(jointRAM1))
omxCheckCloseEnough(cor(estCmp)[2,1], 1, 1e-4)

seCmp <- cbind(jointRAM2$output$standardErrors, jointRAM1$output$standardErrors)
omxCheckCloseEnough(cor(seCmp)[2,1], 1, .05)

# ===============================================================
# = Iterate over model types allowed by the buildModel function =
# ===============================================================

mani = paste0('z', 5:1)
for (slopes in c(TRUE,FALSE)) {
	for (type in c('WLS','DWLS','ULS')) {
		jm3 <- buildModel(type=type, slopes=slopes)
		#jm3$data$verbose <- 1L
		jm3 <- mxRun(jm3)
		jm4 <- mxModel(buildModel(manifests = mani, type=type, slopes=slopes), jm3$data)
		jm4 <- mxRun(jm4)
	}
}

