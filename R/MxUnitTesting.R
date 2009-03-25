#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


omxCheckEquals <- function(a, b) {
	if (a != b) {
		stop(paste("Error:", a, "and", b, "are not equal"))
	}	
}

omxCheckTrue <- function(a) {	
	if(!a) {
		stop(paste("Error", match.call()$a, "is not true"))
	}
}

omxCheckCloseEnough <- function(a, b, epsilon=10^(-15)) {	
	if(abs(a - b) > epsilon) {
		stop(paste("Error:", a, "and", b, "are equal to within", epsilon))
	}
}

omxCheckWithinPercentError <- function(a, b, epsilon=10^(-15)) {	
	if(abs((a - b)/a*100) > epsilon) {
		stop(paste("Error: ", b, "does not estimate", a, "within", epsilon, "percent error"))
	}
}
