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
loadings <- c(1, -0.625, 0.1953125, 
	1, -0.375, 0.0703125, 
	1, -0.125, 0.0078125,
	1,  0.125, 0.0078125,
	1,  0.375, 0.0703125,
	1,  0.625, 0.1953125)
loadings <- matrix(loadings, 6, 3, byrow = TRUE)
L <- mxMatrix("Full", free=FALSE, values=loadings, name="L", byrow=TRUE)
omxCheckIdentical(loadings, L$values)

mxMatrixFromChar <- function(inputm, name=NA) {
  inputmFixed <- suppressWarnings(matrix(
    as.numeric(inputm),nrow = nrow(inputm), ncol = ncol(inputm)))
  inputmCharacter <- inputm
  inputmCharacter[!is.na(inputmFixed)] <- NA
  mxMatrix(nrow=nrow(inputm), ncol=ncol(inputm),
           free=!is.na(inputmCharacter),
           values=inputmFixed,
           labels=inputmCharacter,
           dimnames=dimnames(inputm), name=name)
}

mat <- mxMatrixFromChar(matrix(c(1,2,"a","b"), nrow=2,ncol=2))
omxCheckEquals(mat$labels[,2], c('a','b'))
omxCheckTrue(all(is.na(mat$labels[,1])))
omxCheckTrue(all(is.na(mat$values[,2])))
omxCheckEquals(mat$values[,1], c(1,2))
omxCheckEquals(c(mat$free), c(FALSE,FALSE,TRUE,TRUE))
