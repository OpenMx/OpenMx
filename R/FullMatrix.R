#
#   Copyright 2007-2016 The OpenMx Project
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


setClass(Class = "FullMatrix",
	representation = representation(),
	contains = "MxMatrix")

setMethod("imxCreateMatrix", "FullMatrix",
	  function(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
		    condenseSlots, joinKey, joinModel) {
		if (single.na(values)) {
			values <- 0
		}
		if (is.vector(values)) {
			values <- matrix(values, nrow, ncol, byrow = byrow)
		}
    if(condenseSlots && all.false(free) && all.na(labels)){
      labels <- as.character(NA)
      free <- FALSE
    }
    else{
  		if (is.vector(labels)) {labels <- matrix(labels, nrow, ncol, byrow = byrow)}
  		if (is.vector(free)) {free <- matrix(free, nrow, ncol, byrow = byrow)}
    }
    if(condenseSlots && all.na(lbound)){lbound <- as.numeric(NA)}
		else{if (is.vector(lbound)) {
			lbound <- matrix(lbound, nrow, ncol, byrow = byrow)
		}}
		if(condenseSlots && all.na(ubound)){ubound <- as.numeric(NA)}
		else{if (is.vector(ubound)) {
			ubound <- matrix(ubound, nrow, ncol, byrow = byrow)
		}}
		return(callNextMethod(.Object, labels, values, free, lbound, ubound, nrow, ncol, byrow, name,
				      condenseSlots, joinKey, joinModel) )
	}
)

setMethod("imxVerifyMatrix", "FullMatrix",
	function(.Object) {
		callNextMethod(.Object)
	}
)
