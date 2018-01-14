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

datawide <- suppressWarnings(try(read.csv("models/passing/data/datawide.csv"),
                              silent=TRUE))
if (is(datawide,"try-error")) datawide <- read.csv("data/datawide.csv")

ylab <- paste("y", 1:5, sep="")
zlab <- paste("z", 1:5, sep="")

colnames(datawide) <- c("id", "x", ylab, zlab)
datawide$x <- as.numeric(datawide$x)
latentLab <- c("intcept", "slope")

Alab <- matrix(NA, 8, 8)
Alab[7, 6] <- "groupdiff"
Alab[1:5, 8] <- paste("data.", zlab, sep="")
Aval <- matrix(0, 8, 8)
Aval[7, 6] <- 0.1 # Need parameter values
Aval[1:5, 7] <- 1
Aval[1:5, 8] <- 1
Afree <- matrix(FALSE, 8, 8)
Afree[7, 6] <- TRUE
colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- c(ylab, "x", latentLab)

Slab <- matrix(NA, 8, 8)
diag(Slab) <- "l1error"
Slab[6, 6] <- "varX"
Slab[7, 7] <- "l2error"
Slab[8, 8] <- NA
Sval <- matrix(0, 8, 8)
diag(Sval) <- 0.8 # Need parameter values
Sval[6, 6] <- 0.25 # Need parameter values
Sval[7, 7] <- 0.2 # Need parameter values
Sval[8, 8] <- 0
Sfree <- matrix(FALSE, 8, 8)
diag(Sfree) <- TRUE
Sfree[8, 8] <- FALSE
colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- c(ylab, "x", latentLab)

Fval <- cbind(diag(6), matrix(0, 6, 2))
Flab <- matrix(NA, 6, 8)
Ffree <- matrix(FALSE, 6, 8)
colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- c(ylab, "x", latentLab)
rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- c(ylab, "x")

Mlab <- matrix(c(rep(NA, 5), "meanX", "mean0", "betaW"), 1, 8)
Mval <- matrix(c(rep(0, 5), 0.5, 0.2, 0.8), 1, 8) # Need parameter values
Mfree <- matrix(c(rep(FALSE, 5), rep(TRUE, 3)), 1, 8)
colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- c(ylab, "x", latentLab)

onecov <- mxModel("With covariate model", type="RAM",
mxData(datawide[,c(ylab, "x", zlab)], type="raw"),
mxMatrix(type="Full", nrow=8, ncol=8, values=Aval, free=Afree, labels=Alab, name="A"),
mxMatrix(type="Symm", nrow=8, ncol=8, values=Sval, free=Sfree, labels=Slab, name="S"),
mxMatrix(type="Full", nrow=6, ncol=8, values=Fval, free=Ffree, labels=Flab, name="F"),
mxMatrix(type="Full", nrow=1, ncol=8, values=Mval, free=Mfree, labels=Mlab, name="M"),
mxFitFunctionML(),mxExpectationRAM("A","S","F","M", dimnames=c(ylab, "x", latentLab)),
mxMatrix(type="Full", nrow=8, ncol=1, values=c(rep(0, 6), 1, 0), free=FALSE, name="s7"),
mxMatrix(type="Full", nrow=8, ncol=1, values=c(rep(0, 5), 1, 0, 0), free=FALSE, name="s6")
) # - end model

onecov <- mxModel(onecov, mxAlgebra(expression = A[7, 6], name = "alpha"), mxCI(c("alpha")))

onecov <- mxRun(onecov, intervals=TRUE)
omxCheckCloseEnough(onecov$output$fit, 775.0977, .01)

d1=onecov$compute$steps[['CI']]$output$detail

#cat(deparse(round(c(onecov$output$confidenceIntervals), 3)))
omxCheckCloseEnough(c(onecov$output$confidenceIntervals),
                    c(-0.253, 0.102, 0.456), .01)
