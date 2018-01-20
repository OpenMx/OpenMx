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
myDataRaw<-read.table("data/myAutoregressiveData.txt",header=T)

myDataCov<-matrix(
	c(0.672, 0.315, 0.097, -0.037, 0.046,
	  0.315, 1.300, 0.428,  0.227, 0.146,
	  0.097, 0.428, 1.177,  0.568, 0.429,
	 -0.037, 0.227, 0.568,  1.069, 0.468,
	  0.046, 0.146, 0.429,  0.468, 1.031),
	nrow=5,
	dimnames=list(
	c("x1","x2","x3","x4","x5"),
	c("x1","x2","x3","x4","x5"))
	)
	
myDataMeans <- c(3.054, 1.385, 0.680, 0.254, -0.027)
names(myDataMeans) <- c("x1","x2","x3","x4","x5")

model<-mxModel("Autoregressive Model Path", 
      type="RAM",
      mxData(myDataCov,type="cov", means=myDataMeans, numObs=100),
      manifestVars=c("x1","x2","x3","x4","x5"),
      mxPath(from=c("x1","x2","x3","x4"),
            to=c("x2","x3","x4","x5"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1),
            labels=c("beta","beta","beta","beta")
            ),
      mxPath(from=c("x1","x2","x3","x4","x5"),
            arrows=2,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("varx1","e2","e3","e4","e5")
            ),
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("mean1","mean2","mean3","mean4","mean5")
            )
      ) # close model
       
autoregressivePathCov<-mxRun(model)

autoregressivePathCov$output

# Slightly disagrees with old Mx Output, but matches Mplus exactly
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["beta"]], 0.427, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["varx1"]], 0.665, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["e2"]], 1.142, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["e3"]], 1.038, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["e4"]], 0.791, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["e5"]], 0.818, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["mean1"]], 3.054, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["mean2"]], 0.082, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["mean3"]], 0.089, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["mean4"]], -0.036, 0.001)
omxCheckCloseEnough(autoregressivePathCov$output$estimate[["mean5"]], -0.135, 0.001)
