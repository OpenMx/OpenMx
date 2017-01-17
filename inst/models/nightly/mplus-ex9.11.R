# MPLUS: Two-level multiple group CFA with continuous factor indicators
# https://www.statmodel.com/usersguide/chapter9.shtml

library(OpenMx)

options(width=120)
ex911 <- suppressWarnings(try(read.table("models/nightly/data/ex9.11.dat")))
if (is(ex911, "try-error")) ex911 <- read.table("data/ex9.11.dat")
colnames(ex911) <- c(paste0('y',1:6), 'g', 'clus')
ex911$clus <- as.integer(ex911$clus)

# For testing
#ex911 <- subset(ex911, clus < 50 | (clus >= 111 & clus < 150))

mkGroup <- function(name, dat) {
    betweenModel <- mxModel(
        paste0('between',name), type='RAM',
        latentVars=c(paste0('y',1:6), 'fb1', 'fb2'),
        mxData(dat[!duplicated(dat$clus),], 'raw', primaryKey='clus'),
        mxPath('fb1', paste0('y',1:3), free=c(FALSE,TRUE,TRUE), values=1, lbound = -1, ubound = 10,
               labels=paste0('bl',1:3)),
        mxPath('fb2', paste0('y',4:6), free=c(FALSE,TRUE,TRUE), values=1, lbound = -1, ubound = 10,
               labels=paste0('bl',4:6)),
        mxPath(c('fb1','fb2'), arrows=2, connect="unique.pairs", values=c(.5,0,.5)),
        mxPath('one', c('fb1','fb2'), free=FALSE),
        mxPath(paste0('y',1:6), arrows=2, values=1))
	
    withinModel <- mxModel(
        paste0('within', name), type='RAM', betweenModel,
        manifestVars=paste0('y',1:6), latentVars=c('fw1','fw2'),
        mxData(dat, 'raw'),
        mxPath('fw1', paste0('y',1:3), free=c(FALSE,TRUE,TRUE), values=1, lbound = -1, ubound = 10),
        mxPath('fw2', paste0('y',4:6), free=c(FALSE,TRUE,TRUE), values=1, lbound = -1, ubound = 10),
        mxPath(paste0('y',1:6), arrows=2, values=1),
        mxPath('one', paste0('y',1:6), labels=paste0('yi',1:6)),
        mxPath(c('fw1','fw2'), arrows=2, connect="unique.pairs", values=c(1,0,1)),
        mxPath(paste0('between',name,'.y',1:6), paste0('y',1:6), free=FALSE, values=1, joinKey='clus'))
}

cfa <- mxModel(
    'cfa',
    mkGroup('Group1', subset(ex911, g==1)),
    mkGroup('Group2', subset(ex911, g==2)),
    mxFitFunctionMultigroup(paste0('withinGroup',1:2)))

cfa$withinGroup2$betweenGroup2$M$free[1,c('fb1','fb2')] <- TRUE

cfa <- mxTryHard(cfa)
omxCheckCloseEnough(cfa$output$fit, 49853.908, 1e-2)
# Mplus -24926.956 * -2 = 49853.91

cfa$withinGroup1$expectation$.rampartCycleLimit <- 0L
cfa$withinGroup2$expectation$.rampartCycleLimit <- 0L
cfa <- mxRun(mxModel(cfa,
		    mxComputeSequence(list(
			mxComputeOnce('fitfunction', 'fit'),
			mxComputeReportExpectation()))))
omxCheckCloseEnough(cfa$output$fit, 49853.908, 1e-2) # same location without Rampart

# Here is the MPlus solution
f1 <- omxSetParameters(cfa, labels=names(coef(cfa)),
                        values=c(1.129, 1.218, 0.969, 1.083,         # group 1 within loadings
                            4.139, 4.089, 3.898, 3.953, 4.715, 3.797,# group 1 y variances
                            1.590, 0.948, 1.795,                     # group 1 fw covariance
                            .046, -.044, .078, -.13, -.223, .076,    # y intercepts (equated)
                            1.382, .587, .799, .942,                 # between loadings (equated)
                            .984, .439, 1.071, .883, .964, 1.329,    # group 1 between y variances
                            .707, .053, .462,                        # group 1 between fb covariance
                            .576, .519, .69, .754,                   # group 2 within loadings
                            3.484, 3.483, 4.072, 3.113, 4.2, 3.808,  # group 2 within y variances
                            2.53, 1.054, 2.663,                      # group 2 fw covariance
                            1.162, .399, 1.244, 1.21, 1.452, 1.052,  # group 2 between y variances
                            .656, .042, .408,                        # group 2 fb covariance
                            -.015, .096                              # group 2 fb means
                                 ))

f1$withinGroup1$expectation$.rampartCycleLimit <- 0L
f1$withinGroup2$expectation$.rampartCycleLimit <- 0L
f1 <- mxRun(mxModel(f1,
		    mxComputeSequence(list(
			mxComputeOnce('fitfunction', 'fit'),
			mxComputeReportExpectation()))))
omxCheckCloseEnough(f1$output$fit, 49853.91, 1e-2)

f1$withinGroup1$expectation$.rampartCycleLimit <- NA_integer_
f1$withinGroup2$expectation$.rampartCycleLimit <- NA_integer_
f1 <- mxRun(mxModel(f1,
		    mxComputeSequence(list(
			mxComputeOnce('fitfunction', 'fit'),
			mxComputeReportExpectation()))))
omxCheckCloseEnough(f1$output$fit, 49853.91, 1e-2)

if (0) { # double check everything
	ed = f1$withinGroup2$expectation$debug
	layout <- ed$layout
	ed$numGroups
	dim(ed$g1$covariance)
	S = f1$withinGroup2$S$values
	A = f1$withinGroup2$A$values
	IA <- solve(diag(8) - A)
	g1cov <- IA %*% S %*% t(IA)
	max(abs(g1cov[1:6,1:6] - ed$g1$covariance))  # OK

	head(layout[layout$group==2,],n=20)
	g2 = ed$g2
	max(abs(f1$withinGroup2$betweenGroup2$S$values - g2$S[1:8,1:8])) #OK
	max(abs(f1$withinGroup2$S$values - g2$S[9:16,9:16])) #OK
	g2$A[1:8,1:8] - f1$withinGroup2$betweenGroup2$A$values #OK
	g2$A[9:16,9:16] - f1$withinGroup2$A$values #OK
	IA <- solve(diag(16) - as.matrix(g2$A))
	g2cov <- IA %*% as.matrix(g2$S) %*% t(IA)
	max(abs(g2cov[9:14, 9:14] - g2$covariance)) #OK

	M = cbind(f1$withinGroup2$betweenGroup2$M$values, f1$withinGroup2$M$values)
	A = as.matrix(g2$A)
	A[9:16, 1:8] <- A[9:16, 1:8] / sqrt(5)
	IA <- solve(diag(16) - A)
	S = g2$S
	IA %*% t(M) - g2$fullMean[1:16]  #OK

	g2$fullMean[9:14] * sqrt(5) - g2$mean[1:6] #OK
}
