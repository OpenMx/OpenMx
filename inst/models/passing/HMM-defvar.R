library(OpenMx)

set.seed(1)

start_prob <- c(.2,.4,.4)
amplitude <- 5

# suppose we can manipulate one of the columns of transition probabilities
transition_prob <- function(tm) {
  tp <- matrix(c(1, 0, 0,
                 0, 1, 0,
                 0, 0, amplitude*sin(tm/10.0)), 3, 3)
  tp <- exp(tp)
  t(t(tp) / colSums(tp))
}

# simulate a trajectory
state <- sample.int(3, 1, prob=transition_prob(0) %*% start_prob)
trail <- c(state)
totalRep <- 5000
for (rep in 1:totalRep) {
	state <- sample.int(3, 1, prob=transition_prob(rep)[,state])
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
			       mxData(data.frame(ob=trailN,
			                         tp1=amplitude*sin(0:totalRep/10.0)), "raw"),
	mxMatrix(nrow=1, ncol=length(classes), name="initial"),
	mxMatrix(nrow=length(classes), ncol=length(classes), name="transition"),
	mxExpectationHiddenMarkov(paste0("class",1:3), "initial",
	                          "transition", scale="softmax"),
	mxFitFunctionML())

hmm$transition$free[c(1,2),] <- TRUE
hmmFit <- mxRun(hmm)
omxCheckCloseEnough(hmmFit$output$fit, 11722.107, .01)

tolerance <- 0.03
ex = hmmFit$expectation
omxCheckTrue(max(abs(ex$output$transition[,1:2] -
                       transition_prob(0)[,1:2])) > tolerance)

hmm$transition$labels[3,3] <- 'data.tp1'
hmm$transition$free[2,3] <- FALSE
hmmFit <- mxRun(hmm)

print(hmmFit$output$fit)
omxCheckCloseEnough(hmmFit$output$fit, 10826.65, .01)

ex = hmmFit$expectation
ex$output$transition[,1:2] - transition_prob(0)[,1:2]
omxCheckCloseEnough(ex$output$transition[,1:2], transition_prob(0)[,1:2], tolerance)
