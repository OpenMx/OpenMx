library(OpenMx)

cubic <- function(marg, state){
    x <- marg@matrices$param@values[1]
    got <- (x-5)*(x-1)*(x+3)
    return(got)
}

model <- mxModel(name="root",
	      mxMatrix(type="Full", ncol=1, nrow=1, name="param", free=TRUE, values=0),
	      mxFitFunctionR(cubic))
model <- mxRun(model, silent=TRUE, suppressWarnings=TRUE)
omxCheckCloseEnough(model@matrices$param@values, 3.309401, 10^-3)

###

infer <- function(marg,state) return(Inf)

model <- mxModel(name="inf", mxFitFunctionR(infer))
omxCheckError(mxRun(model), c("The job for model 'inf' exited abnormally with the error message: Fit function returned inf at iteration 0.1",
			      "The job for model 'inf' exited abnormally with the error message: Fit function returned 1.#INF at iteration 0.1"))

###

NAer <- function(marg,state) return(NA)

model <- mxModel(name="na", mxFitFunctionR(NAer))
omxCheckError(mxRun(model), "The job for model 'na' exited abnormally with the error message: Fit function returned nan at iteration 0.1");

###

counter <<- 1
count <- function(marg,state) {
	counter <<- state[[1]]
	state[[1]] <- state[[1]] + 1
	return(list(1, state))
}

model <- mxModel(name="count", mxFitFunctionR(count, 1))
model <- mxRun(model, silent=TRUE)
omxCheckTrue(counter > 1)

###

toomany <- function(marg,state) {
	return(list(1, state, 5))
}

model <- mxModel(name="toomany", mxFitFunctionR(toomany))
omxCheckError(mxRun(model), "The job for model 'toomany' exited abnormally with the error message: FitFunction returned more than 2 arguments")
