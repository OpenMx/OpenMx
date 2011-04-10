#  Copyright (c) 2001, 2002 Enthought, Inc.
#  All rights reserved.
#
#  Copyright (c) 2003-2009 SciPy Developers.
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#    a. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#    b. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#    c. Neither the name of the Enthought nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
#  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.




# Translated from the function linalg.expm in the SciPy library

omxExponential <- function(x, order = 7) {
	if (nrow(x) != ncol(x)) {
		stop("input to omxExponential() must be a square matrix")
	}
	dimension <- nrow(x)
    normX <- max(apply(abs(x), 1, sum))
	if (normX == 0) {
		return(diag(dimension))
	}
    val <- log2(normX)
    e <- floor(val)
    j <- max(0, e + 1)
    x <- x / (2.0 ^ j)

    A <- x
    const <- 1/2
    N <- diag(dimension) + const * x
    D <- diag(dimension) - const * x
    for(k in 2:order) {
        const <- const * (order - k + 1) / (k * (2 * order - k + 1))
        A <- x %*% A
        cX <- const * A
        N <- N + cX
        if(k %% 2 == 0) {
            D <- D + cX
        } else {
            D <- D - cX
		}
	}
    F <- solve(D, N)
	if (j > 0) {
	    for(k in 1:j) {
    	    F <- F %*% F
		}
	}
    return(F)
}
