library(OpenMx)

set.seed(1)

start_prob <- c(.2,.4,.4)
transition_prob <- matrix(c(.8, .1, .1,
							.3, .6, .1,
							.1, .3, .6), 3, 3)

# simulate a trajectory
state <- sample.int(3, 1, prob=transition_prob %*% start_prob)
trail <- c(state)
for (rep in 1:1000) {
	state <- sample.int(3, 1, prob=transition_prob[,state])
	trail <- c(trail, state)
}

# add noise
noise <- .2
trailN <- sapply(trail, function(v) rnorm(1, mean=v, sd=sqrt(noise)))

classes <- list()

for (cl in 1:3) {
	classes[[cl]] <- mxModel(paste0("class", cl), type="RAM",
			       manifestVars=c("ob"),
			       mxPath("one", "ob", value=cl, free=FALSE),
			       mxPath("ob", arrows=2, value=noise, free=FALSE),
			       mxFitFunctionML(vector=TRUE))
}

hmm <- mxModel("hmm", classes,
			       mxData(data.frame(ob=trailN), "raw"),
	mxMatrix(values=.1, nrow=1, ncol=2, free=TRUE, name="start"),
	mxMatrix(values=.1, nrow=length(classes)-1, ncol=length(classes),
	         free=TRUE, name="transition"),
	mxAlgebra(cbind(exp(start), 1), "startFull"),
	mxMatrix('Unit', 1, length(classes), name="urow"),
	mxAlgebra(rbind(exp(transition), urow), "transitionFull"),
	mxExpectationHiddenMarkov(paste0("class",1:3), "startFull",
	                          "transitionFull", scale="sum"),
	mxFitFunctionML(),
	mxComputeSequence(list(
		mxComputeGradientDescent(),
		mxComputeReportExpectation())))

hmmFit <- mxRun(hmm)

omxCheckCloseEnough(hmmFit$output$fit, 2263.967, .01)

ex = hmmFit$expectation

omxCheckEquals(order(-ex$output$initial)[1], round(trail[1]))

print(max(abs(ex$output$transition - transition_prob)))
omxCheckCloseEnough(ex$output$transition, transition_prob, .06)
