library(OpenMx)

start_prob <- c(.2,.4,.4)
transition_prob <- matrix(.1,3,3)
diag(transition_prob) <- .8

# simulate a trajectory
state <- sample.int(3, 1, prob=start_prob)
trail <- c(state)
for (rep in 1:5000) {
	state <- sample.int(3, 1, prob=transition_prob[,state])
	trail <- c(trail, state)
}

# add noise
variance <- .2
trailN <- sapply(trail, function(v) rnorm(1, mean=v, sd=sqrt(variance)))

classes <- list()

for (cl in 1:3) {
	classes[[cl]] <- mxModel(paste0("class", cl), type="RAM",
			       manifestVars=c("ob"),
			       mxData(data.frame(ob=trailN), "raw"),
			       mxPath("one", "ob", value=cl, free=FALSE),
			       mxPath("ob", arrows=2, value=variance, free=FALSE),
			       mxFitFunctionML(vector=TRUE))
}

hmm <- mxModel("hmm", classes,
	mxMatrix(values=.1, nrow=1, ncol=2, free=TRUE, name="start"),
	mxMatrix(values=.1, nrow=length(classes)-1, ncol=length(classes),
	         free=TRUE, name="transition"),
	mxAlgebra(cbind(exp(start), 1), "startFull"),
	mxMatrix('Unit', 1, length(classes), name="urow"),
	mxAlgebra(rbind(exp(transition), urow), "transitionFull"),
	mxExpectationHiddenMarkov(paste0("class",1:3), "startFull", "transitionFull"),
	mxFitFunctionML(verbose=0L))
#	mxComputeOnce('fitfunction', 'fit'))

hmmFit <- mxRun(hmm)
summary(hmmFit)

start <- mxEval(startFull, hmmFit, compute=T)
start <- start / sum(start)
print(head(trailN))

transition <- mxEval(transitionFull, hmmFit, compute=T)
transition <- t(t(transition) / colSums(transition))
print(transition)
