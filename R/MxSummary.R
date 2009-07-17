#
#   Copyright 2007-2009 The OpenMx Project
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

computeOptimizationStatistics <- function(model) {
	retval <- list()
	if(length(model@output) == 0) { return(retval) }
	ptable <- data.frame()
	objective <- model@objective
	parameters <- model@output$estimate
	if (!(is.null(objective) || is(objective, "MxAlgebraObjective"))) {
		retval[['AIC']] <- model@output$minimum + 
			2 * length(parameters)
	}
	if (length(parameters) > 0) {
		for(i in 1:length(parameters)) {
			ptable[i, 'name'] <- names(parameters)[[i]]
			ptable[i, 'parameter estimate'] <- parameters[[i]]
			ptable[i, 'error estimate'] <- model@output$hessian[i, i]			
		}
		retval[['parameters']] <- ptable
	}
	return(retval)
}

setMethod("summary", "MxModel",
	function(object, ...) {	
		retval <- computeOptimizationStatistics(object)
		return(retval)
	}
)

