demo(GrowthMixtureModel_MatrixRaw)

################################################################################
# Serial Respecification
################################################################################

# how many trials?
trials <- 20

# place all of the parameter names in a vector
parNames <- names(omxGetParameters(gmm))

# make a matrix to hold all of the 
input <- matrix(NA, trials, length(parNames))
dimnames(input) <- list(c(1: trials), c(parNames))

output <- matrix(NA, trials, length(parNames))
dimnames(output) <- list(c(1: trials), c(parNames))

fit <- matrix(NA, trials, 5)
dimnames(fit) <- list(c(1: trials), c("Minus2LL", "Status", "Iterations", "pclass1", "time"))

# populate the class probabilities
input[,"p1"] <- runif(trials, 0.1, 0.9)
input[,"p1"] <- input[,"p1"]/(1-input[,"p1"])

# populate the variances
v <- c("vari1", "vars1", "vari2", "vars2", "residual")
input[,v] <- runif(trials*5, 0, 10)

# populate the means
m <- c("meani1", "means1", "meani2", "means2")
input[,m] <- runif(trials*4, -5, 5)

# populate the covariances
r <- runif(trials*2, -0.9, 0.9)
scale <- c(
    sqrt(input[,"vari1"]*input[,"vars1"]),
    sqrt(input[,"vari2"]*input[,"vars2"]))
input[,c("cov1", "cov2")] <- r * scale


for (i in 1: trials){
	temp1 <- omxSetParameters(gmm,
		labels=parNames,
		values=input[i,]
		)
		
	temp1@name <- paste("Starting Values Set", i)
		
	temp2 <- mxRun(temp1, unsafe=TRUE, suppressWarnings=TRUE, checkpoint=TRUE)
	
	output[i,] <- omxGetParameters(temp2)
	fit[i,] <- c(
		temp2@output$Minus2LogLikelihood,
		temp2@output$status[[1]],
		temp2@output$iterations,
		round(temp2$classProbs@result[1,1], 4),
		temp2@output$wallTime
		)
	}
	
fit
table(round(fit[,1], 3), fit[,2])


################################################################################
# Parallel Respecification
################################################################################
require(snowfall)
sfInit(parallel=TRUE, cpus=4)
sfLibrary(OpenMx)

topModel <- mxModel("Top")	

makeModel <- function(modelNumber){
	temp <- mxModel(gmm, 
		independent=TRUE,
		name=paste("Iteration", modelNumber, sep=""))
	temp <- omxSetParameters(temp,
		labels=parNames,
		values=input[modelNumber,])
	return(temp)
}
	
mySubs <- lapply(1:20, makeModel)
	
topModel@submodels <- mySubs

results <- mxRun(topModel)

fitStats <- function(model){
	retval <- c(
		model@output$Minus2LogLikelihood,
		model@output$status[[1]],
		model@output$iterations,
		round(model$classProbs@result[1,1], 4)
		)	
	return(retval)
}

resultsFit <- t(omxSapply(results@submodels, fitStats))
sfStop()


################################################################################
# Compare Model Estimation Times
################################################################################

# Serial Time, in seconds (ignoring the overhead caused by the for loop)
sum(fit[,5])

# Parallel Time, in seconds
results@output$wallTime

################################################################################
# Reload OpenMx to offload the snowfall OpenMx library
################################################################################
library(OpenMx)