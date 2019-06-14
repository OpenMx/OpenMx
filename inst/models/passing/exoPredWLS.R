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
# jointdata is a dataframe with 5 variables: "z1", "z2", "z3", "z4", "z5"
# 'data.frame':	250 obs. of  5 variables:
#  $ z1: num  6.83 8.77 8.01 9 8.52 ...
#  $ z2: int  0 0 0 1 0 0 0 0 0 0 ...
#  $ z3: num  -0.0647 2.8983 2.5471 2.9078 3.4518 ...
#  $ z4: int  2 2 1 2 1 1 0 1 0 3 ...
#  $ z5: int  2 1 2 2 0 2 0 2 2 1 ...

round(cov(jointdata),3)
#       z1    z2    z3    z4    z5
# z1 0.927 0.149 0.437 0.368 0.089
# z2 0.149 0.250 0.148 0.202 0.065
# z3 0.437 0.148 0.936 0.454 0.049
# z4 0.368 0.202 0.454 1.279 0.122
# z5 0.089 0.065 0.049 0.122 0.635

# ==============================================================================================================
# = Make data suitable for a joint (continuous and ordinal) model by reformatting some variables as mxFactors  =
# ==============================================================================================================
# z2, z4, and z5 to be mxFactors with 2, 4, and 3 levels respectively
jointdata[,c("z2", "z4", "z5")] <- mxFactor(jointdata[,c("z2", "z4", "z5")], levels= list(c(0, 1), c(0, 1, 2, 3), c(0, 1, 2)))

# Make z1c = z1 + some noise;
jointdata$z1c <- with(jointdata, z1 * .1 + rnorm(nrow(jointdata)))
# and z2c which is just a column of noise with a mean = the mean factor level of z2 (no clue why as yet...)
jointdata$z2c <- with(jointdata, rnorm(nrow(jointdata), mean=unclass(z2)*.2))

# ============================================
# = Create a thresholds matrix for the model =
# ============================================
thresh <- mxMatrix(name="T", "Full", nrow= 3, ncol= 3, free = FALSE, values = 0)

# ... and fill in for columns 1:3 matching vars z2, z4, and z5
thresh$free[,1]   <- c(TRUE  , FALSE, FALSE)
thresh$values[,1] <- c(0     , NA   , NA)
thresh$labels[,1] <- c("z2_thresh1", NA   , NA)

thresh$free[,2]   <- c(TRUE, TRUE, TRUE)
thresh$values[,2] <- c(-1  , 0   , 1)
thresh$labels[,2] <- c("z4_thresh1", "z4_thresh2", "z4_thresh3")

thresh$free[,3]   <- c(TRUE, TRUE, FALSE)
thresh$values[,3] <- c(-1  , 1   , NA)
thresh$labels[,3] <- c("z5_thresh1", "z5_thresh2", NA)
# Add column names
colnames(thresh)  <- paste0('z', c(2, 4, 5))


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
		# include the thresholds matrix made earlier
		thresh,
		# Note: No data are still raw, despite our using WLS
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
	jr$expectation$thresholds <- 'T'
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

# =========================
# = Test permutation code =
# =========================
# TODO: How does this test "permutation" code? Should this be "iterate over model types allowed by the builModel function"

for (slopes in c(TRUE,FALSE)) {
	for (type in c('WLS','DWLS','ULS')) {
		jm3 <- buildModel(type=type, slopes=slopes)
		#jm3$data$verbose <- 1L
		jm3 <- mxRun(jm3)
		jm4 <- mxModel(buildModel(paste0('z', 5:1), type=type, slopes=slopes), jm3$data)
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
    thresh,
	mxData(jointdata, "raw", verbose=0L),
    mxFitFunctionWLS()
	# Shouldn't we have a call to mxExpectationRAM or LISREL??? here???
	# mxExpectationRAM(M = NA, dimnames = NA, thresholds = "T", threshnames = ???)	
)

# =================================================================
# = TODO: How would a user ever discover this?: What does it do?? =
# =================================================================
jointLISREL$expectation$thresholds <- 'T'
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
