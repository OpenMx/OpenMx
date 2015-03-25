#
#   Copyright 2007-2010 The OpenMx Project
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

A <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'A')
B <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'B')
C <- mxMatrix(values = runif(30), nrow = 6, ncol = 5, name = 'C')
D <- mxMatrix(values = 1:10, nrow = 2, ncol = 5, name = 'D')
E <- mxMatrix(values = 1:10, nrow = 5, ncol = 2, name = 'E')
# For Mnor and AllInt
F <- mxMatrix("Stand", nrow=2, ncol=2, values=c(.95), free=FALSE, name="F")
G <- mxMatrix("Full", values=rbind(c(0, 0), c(2,3)), free=FALSE, name="G")
J <- mxMatrix("Full", values=rbind(c(1, 0, 0), c(1, 0, 0)), free=FALSE, name="J")
L <- mxMatrix("Full", nrow=1, ncol=3, values=c(0,0,-Inf), free=FALSE, name="L")
M <- mxMatrix("Full", nrow=1, ncol=3, values=c(0,0,0), free=FALSE, name="M")
U <- mxMatrix("Full", nrow=1, ncol=3, values=c(Inf,Inf,Inf), free=FALSE, name="U")
V <- mxMatrix("Stand", nrow=3, ncol=3, values=c(.5, .5, .5), free=FALSE, name="V")

dimnames(A) <- list(letters[1:5], letters[22:26])
dimnames(B) <- dimnames(A)
dimnames(C) <- list(letters[1:6], letters[22:26])

model <- mxModel(A, B, C, D, E, V, M, L, U, name = 'model')
zeta <- 'z'
alpha <- 'a'
two <- 2

# Insert passing tests
model <- mxModel(model, mxAlgebra(omxMnor(J%&%V, M%*%t(J), L%*%t(J), U%*%t(J)), name= 'test35b'))
model <- mxRun(model)

# Check passing tests

omxCheckCloseEnough(model[['test35b']]$result, .33333, 0.001)
