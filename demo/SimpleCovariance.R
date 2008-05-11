
model <- MxModel()
model$A <- FullMatrix(2, 2, free=TRUE)
model$S <- FullMatrix(2, 2, free=TRUE)
model$F <- FullMatrix(2, 2)
model$I <- IdenMatrix(2, 2)
model$cov <- MxAlgebra(model$F %&% (solve(model$I - model$A) %&% model$S))

setValues(model$A,c(0,0,0.2,0))
setValues(model$S,c(0.8,0,0,0.8))
setValues(model$F,c(1,0,0,1))

setParameters(model$A, c(0,0,1,0))
model$S$parameters[1,2] <- model$S$parameters[2,1] <- "0";
model$S$parameters[1,1] <- "banana";
model$S$parameters[2,2] <- "apple";
setParameters(model$F, matrix(0,2,2))

covMatrix <- matrix(c(0.77642931, 0.39590663, 0.39590663, 0.49115615),
	nrow = 2, ncol = 2, byrow = TRUE)

objective <- CovarianceObjective(model$cov, covMatrix)

job <- MxJob(model, objective)

jobClosure <- createMxClosure(job, use_R=TRUE)

jobClosure()
