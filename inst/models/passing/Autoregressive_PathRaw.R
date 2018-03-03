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

model<-mxModel("Autoregressive Model Path", 
      type="RAM",
      #mxData(myDataCov,type="cov", means=myDataMeans, numObs=100),
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
       
mxd <- mxData(myDataRaw,type="raw")
dataSize <- object.size(mxd)

model <- mxModel(model, mxd)

autoregressivePathRaw <-mxRun(model)

afterRunSize <- object.size(autoregressivePathRaw)
bareRunSize <- object.size(omxModelDeleteData(autoregressivePathRaw))

omxCheckCloseEnough(afterRunSize - 2*dataSize, bareRunSize, 310)

# Comparing to old Mx Output
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["beta"]], 0.4267, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["varx1"]], 0.6658, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["e2"]], 1.1420, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["e3"]], 1.0379, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["e4"]], 0.7908, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["e5"]], 0.8176, 0.01)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["mean1"]], 3.0537, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["mean2"]], 0.0824, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["mean3"]], 0.0890, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["mean4"]], -0.0362, 0.001)
omxCheckCloseEnough(autoregressivePathRaw$output$estimate[["mean5"]], -0.1356, 0.001)
