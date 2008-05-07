


model <- MxModel();
model$A <- FullMatrix(3, 3, free=TRUE);
model$S <- SymmMatrix(3, 3, free=TRUE);
model$F <- IdenMatrix(3, 3);
model$I <- IdenMatrix(3, 3);
model$M <- FullMatrix(3, 1);
model$U <- UnitMatrix(1, 1);

setValues(model$S, c(1.0,0.2,0.0,2.0,0.0,1.0));
setValues(model$A, c(0.0,0.0,0.0,0.0,0.0,0.0,0.3,0.3,0.0));
setValues(model$M, c(0.1,0.1,0.1));

setParameters(model$S, c(1,6,0,2,0,4));
setParameters(model$A, c(0,0,0,0,0,0,8,9,0));
setParameters(model$M, c(13,12,14));

data('SimpleRawData');
observed <- data.matrix(SimpleRawData);
meanAlgebra <- MxAlgebra(model$F %*% solve(model$I - model$A) %*% model$M %*% model$U);
covAlgebra <- MxAlgebra(model$F %&% (solve(model$I - model$A) %&% model$S));

objective <- RawObjective(meanAlgebra, covAlgebra, observed);

job <- MxJob(model, objective);

jobClosure <- createMxClosure(job, use_R=TRUE);
