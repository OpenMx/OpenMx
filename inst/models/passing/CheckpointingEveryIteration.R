#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


require(OpenMx)

#mxOption(NULL, "Default optimizer", "NPSOL")

data(demoOneFactor)
setwd(tempdir())
factorModel <- mxModel("Checkpoint",
                       mxMatrix("Full", 5, 1, values=0.2,
                                free=TRUE, name="A"),
                       mxMatrix("Symm", 1, 1, values=1,
                                free=FALSE, name="L"),
                       mxMatrix("Diag", 5, 5, values=1,
                                free=TRUE, name="U"),
                       mxAlgebra(A %*% L %*% t(A) + U, name="R"),
                       mxExpectationNormal(covariance="R",dimnames=names(demoOneFactor)),
                       mxFitFunctionML(),
                       #mxMLObjective("R", dimnames = names(demoOneFactor)),
                       mxData(cov(demoOneFactor), type="cov", numObs=500))
factorModel <- mxOption(factorModel,"Checkpoint Units",'evaluations')
factorModel <- mxOption(factorModel,"Checkpoint Count",1)
factorFit <- mxRun(factorModel,checkpoint=T)
checkpointdat <- read.table('Checkpoint.omx', as.is=TRUE, sep="\t", header=TRUE, check.names=FALSE)
omxCheckTrue(!is.null(checkpointdat))
omxCheckTrue(nrow(checkpointdat)>1)
omxCheckTrue(all(checkpointdat$OpenMxNumFree == 10))

mask <- checkpointdat$objective < 1000 & checkpointdat$OpenMxContext == mxOption(NULL, "Default optimizer")
traj <- checkpointdat$objective[mask]
omxCheckTrue(length(traj) > 40)
trajDf <- data.frame(y=traj, x=1:length(traj))
m1 <- lm(y ~ I(1/sqrt(x)), trajDf)
omxCheckCloseEnough(summary(m1)$r.squared, .73, .5)

if (0) {
  trajDf$model <- predict(m1)
  require(ggplot2)
  ggplot(trajDf, aes(x, y)) + geom_point(size=2, color="blue") + labs(x="x midpoint") +
    geom_line(aes(x, model), size=.25, color="red")
}
