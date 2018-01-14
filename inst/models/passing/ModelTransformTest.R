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

fix_matrix <- function(m)
{
 for (i in 1:dim(m$labels)[1])
 {
	for (j in 1:dim(m$labels)[2])
	{
				m$free[i,j] <- FALSE;	
	}
 }
 return(m);
}

# (2)
# generate random normal data  Y~N(0,1)
#
N <- 1000
Y <- rnorm(N,0,1)
data <- data.frame(Y)
names(data) <- c("Y")

# (3)
# create a dummy model which estimates the mean 
# of an observed variable with unit variance
#
model <- mxModel("One Factor", type="RAM",
      manifestVars = c("Y"),
      latentVars = c(),
      mxPath(from=c("Y"), arrows=2,
            free=F, values=1.0, label=c("lat_var")),
mxPath(
        from="one",
        to=c("Y"),
        arrows=1,
        free=TRUE,
        values=c(1),
        labels=c("mean")
    )    ,        
      mxData( data, type="raw", numObs=dim(data)[1]) 
     );

# (4)
# run model (with a single free param)
# and estimate mean
#
run <- mxRun(model);


#################
# rewrite version
#################

m <- run$A$values[1,1]

test <- mxModel("test",
	mxData(data, "raw"),
	mxMatrix("Symm", 1, 1, FALSE, 1, name="S"),
	mxMatrix("Full", 1, 1, name="A"),
	mxMatrix("Iden", 1, name="F"),
	mxMatrix("Full", 1, 1, FALSE, m, name="M"),
	mxFitFunctionML(),mxExpectationRAM("A", "S", "F", "M", dimnames=c("Y"))
	)
	
# runs fine
test2 <- mxRun(test)


# flipping parameters

test3 <- run
test3$M$free[1,1] <- F

test4 <- mxRun(test3)

# test3$output <- list()

test5 <- mxRun(test3)
