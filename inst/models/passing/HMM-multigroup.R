library(OpenMx)

set.seed(1)

start_prob <- c(.2,.4,.4)
transition_prob <- matrix(c(.8, .1, .1,
							.3, .6, .1,
							.1, .3, .6), 3, 3)
noise <- .5

people <- list()
for (person in 1:80) {
  # simulate a trajectory
  state <- sample.int(3, 1, prob=transition_prob %*% start_prob)
  trail <- c(state)
  reps <- round(runif(1, 15, 25))
  for (rep in 1:reps) {
    state <- sample.int(3, 1, prob=transition_prob[,state])
    trail <- c(trail, state)
  }
  
  # add noise
  trailN <- sapply(trail, function(v) rnorm(1, mean=v, sd=sqrt(noise)))
  
  name <- paste0("hmm",person)
  
  classes <- list()
  
  for (cl in 1:3) {
    classes[[cl]] <- mxModel(paste0(name, "cl", cl), type="RAM",
                             manifestVars=c("ob"),
                             mxPath("one", "ob", value=cl, free=FALSE),
                             mxPath("ob", arrows=2, value=noise, free=FALSE),
                             mxFitFunctionML(vector=TRUE))
  }

  m1 <-  
    mxModel(name, classes,
            mxData(data.frame(ob=trailN), "raw"),
            mxMatrix(nrow=3, ncol=1,
                     labels=paste0('i',1:3), name="initial"),
            mxMatrix(nrow=length(classes), ncol=length(classes),
                     labels=paste0('t', 1:(length(classes) * length(classes))),
                     name="transition"),
            mxExpectationHiddenMarkov(
              components=sapply(classes, function(m) m$name),
              initial="initial",
              transition="transition", scale="softmax"),
            mxFitFunctionML())
  m1$initial$free[1:(length(classes)-1),1] <- TRUE
  m1$transition$free[1:(length(classes)-1), 1:length(classes)] <- TRUE
  people[[person]] <- m1
}

hmm <- mxModel("hmm", people,
               mxFitFunctionMultigroup(sapply(people, function(m) m$name)))

hmmFit <- mxRun(hmm)
omxCheckTrue(hmmFit$output$status$code <= 1)

ex <- hmmFit$hmm1$expectation

print(max(abs(ex$output$initial - start_prob)))
omxCheckCloseEnough(ex$output$initial,  start_prob, .1)

print(max(abs(ex$output$transition - transition_prob)))
omxCheckCloseEnough(ex$output$transition, transition_prob, .1)

omxCheckCloseEnough(hmmFit$output$fit, 4706.57, .01)
