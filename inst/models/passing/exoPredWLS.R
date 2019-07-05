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
		mxData(jointdata, type = "raw", verbose= 0L),
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
summary(jointRAM1)
# plot(jointRAM1) # (if using umx)

# Where do these come from?
param <-  c(0.5808, 0.5903, 0.6526, 0.6066, 0.1745, 0.0871, 0.0504, 0.5476,  0.4639, 7.8323,
            2.0759, 0.0785, -0.3664, 0.1271, 0.7919, -0.6475,  -0.296)
omxCheckCloseEnough(coef(jointRAM1), param, 1e-3)

param.se <- c(0.0613, 0.1056, 0.0684, 0.0942, 0.0665, 0.0541, 0.0648, 0.0559,  0.0644, 0.1054,
              0.0593, 0.0819, 0.0777, 0.0726, 0.0922, 0.0655,  0.0585)
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
omxCheckCloseEnough(cor(seCmp)[2,1], 1, .18)

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

# =========================
# = Make a LISREL version =
# =========================

jointLISREL <- mxModel("JointLISREL", type="LISREL",
    manifestVars = list(endogenous= paste0('z',1:5)),
    latentVars = list(endogenous= c('G','z1c','z2c')),
    mxPath('one', paste0('z', c(1,3)), free= TRUE, labels= c('z1','z3')),
    mxPath(paste0('z', c(1,3)), arrows= 2, free=TRUE, values= .5),
    mxPath(paste0('z', c(2,4,5)), arrows= 2, free=FALSE, values= .5),
    mxPath('G', arrows= 2, values= 1, free= FALSE),
    mxPath('G', paste0('z', 1:5), free= TRUE, values= 1, lbound= 0, ubound= 4),
    mxPath('one', 'z1c', free= FALSE, labels= "data.z1c"),
    mxPath('one', 'z2c', free= FALSE, labels= "data.z2c"),
    mxPath('z1c', 'z1', labels= "r1"),
    mxPath('z2c', 'z2', labels= "r2"),
    mxThreshold('z2', 1, free=TRUE, labels="z2_thresh1"),
    mxThreshold('z4', 3, free=TRUE, labels=paste0("z4_thresh",1:3)),
    mxThreshold('z5', 2, free=TRUE, labels=paste0("z5_thresh",1:2)),
    mxData(jointdata, "raw", verbose=0L),
    mxFitFunctionWLS()
	# Shouldn't we have a call to mxExpectationRAM or LISREL??? here???
	# mxExpectationRAM(M = NA, dimnames = NA, thresholds = "T", threshnames = ???)	
)

# =================================================================
# = TODO: How would a user ever discover this?: What does it do?? =
# =================================================================
# jointLISREL$expectation$verbose <- 1L

jointLISREL <- mxRun(jointLISREL)

# Compare parameter estimates from the RAM and LISREL models
print(max(abs(coef(jointLISREL) - coef(jointRAM1))))
omxCheckCloseEnough(coef(jointLISREL), coef(jointRAM1), 2e-5)

# tmp <- c(jointRAM1$output$standardErrors)
# names(tmp) <- c()
# cat(deparse(round(tmp,4)))

numPeople <- 100
personData <- data.frame(
  snp=rnorm(numPeople),
  isMale=as.numeric(rbinom(numPeople,1,.5)),
  phenotype=rnorm(numPeople),
  snpsex=0)

m1 <- mxModel(
  "gwsem", type="RAM",
  manifestVars = c('phenotype'),
  latentVars = c('snp', 'sex', 'snpsex'),
  
  # residual variances
  mxPath(c('phenotype'), arrows=2, values=1),
  
  mxPath('one', 'sex', free=FALSE, labels="data.isMale"),
  mxPath('one', 'snp', free=FALSE, labels="data.snp"),
  mxAlgebra(data.snp * data.isMale, name="snpsexAlg",
	  dimnames=list(NULL, 'snpsex')),
  mxPath('one', 'snpsex', free=FALSE, labels="data.snpsex"),
  
  mxPath('one', 'phenotype'),
  mxPath(c('snp','sex', 'snpsex'), 'phenotype'),
  
  mxData(personData, type="raw", algebra='snpsexAlg'),
  mxFitFunctionWLS(allContinuousMethod = "marginals"))

m1 <- mxRun(m1)

personData$snpsex <- with(personData, snp * isMale)
c1 <- coef(lm(phenotype ~ isMale + snp + snpsex, personData))

omxCheckCloseEnough(c1[-1], coef(m1)[c('gwsem.A[1,3]', 'gwsem.A[1,2]', 'gwsem.A[1,4]')], 1e-3)
