demo(GrowthMixtureModel_MatrixRaw)

trials <- 20

omxGetParameters(gmm)

parNames <- names(omxGetParameters(gmm))
	
input <- matrix(NA, trials, length(parNames))
dimnames(input) <- list(c(1: trials), c(parNames))

output <- matrix(NA, trials, length(parNames))
dimnames(output) <- list(c(1: trials), c(parNames))

fit <- matrix(NA, trials, 4)
dimnames(fit) <- list(c(1: trials), c("Minus2LL", "Status", "Iterations", "pclass1"))
	
for (i in 1: trials){
	cp <- runif(1, 0.1, 0.9) # class probability
	v  <- runif(5, 0.1, 5.0) # variance terms
	cv <- runif(2,-0.9, 0.9) # covariances (as correlations)
	m  <- runif(4,-5.0, 5.0) # means
	cv <- cv*c(sqrt(v[2]*v[3]), sqrt(v[4]*v[5])) #rescale covariances
	
	temp1 <- omxSetParameters(gmm,
		labels=parNames,
		values=c(
			cp, # class probability
			v[1], 
			v[2], cv[1], v[3], m[1], m[2],
			v[4], cv[2], v[5], m[3], m[4]
			)
		)
		
	temp1@name <- paste("Starting Values Set", i)
		
	temp2 <- mxRun(temp1, unsafe=TRUE, suppressWarnings=TRUE)
	
	input[i,] <- omxGetParameters(temp1)
	output[i,] <- omxGetParameters(temp2)
	fit[i,] <- c(
		temp2@output$Minus2LogLikelihood,
		temp2@output$status[[1]],
		temp2@output$iterations,
		temp2@output$estimate[1]
		)
	}
	
fit
