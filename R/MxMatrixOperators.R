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


tr <- function(A) {
	if (is.matrix(A)) {
		return(sum(diag(A)))
	} else {
		return(NA)
	}
}

"%&%" <- function(A, B) {
  return(A %*% B %*% t(A))
}

vech <- function(A) {
	return(A[lower.tri(A, diag=TRUE)])
}

vechs <- function(A) {
	return(A[lower.tri(A, diag=FALSE)])
}