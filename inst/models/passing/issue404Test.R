#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# Adapted from script by Timothy C. Bates

library(OpenMx)
# This script doesn't test anything optimization-related, so it doesn't need 
# to be run with anything other than the on-load default:
if(mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")
data(HS.ability.data)  # load the data
Spatial <- c("visual", "cubes", "paper")
m1 <- mxModel(
	"Holzinger_and_Swineford_1939",
	type="RAM",
	manifestVars=Spatial,
	latentVars="vis",
	mxData(observed=HS.ability.data,type="raw"),
	mxPath(from = "vis",  to = Spatial),
	mxPath(from="vis",to="vis",arrows=2,free=F,values=1.0),
	mxPath(from="one",to="vis",arrows=1,free=F,values=0.0),
	mxPath(from=Spatial,to=Spatial,connect="single",arrows=2),
	mxPath(from="one",to=Spatial,connect="single",arrows=1),
	mxAlgebra(name = "myMcols", M, dimnames = list(NULL, c(Spatial, "vis"))),
	mxAlgebra(name = "myMboth", M, dimnames = list("means", c(Spatial, "vis")))
)
m1$expectation$M <- "myMcols"
m2 <- mxRun(m1)
( out2 <- mxStandardizeRAMpaths(m2,SE=T) )

m1$expectation$M <- "myMboth"
m3 <- mxRun(m1)
( out3 <- mxStandardizeRAMpaths(m3,SE=T) )

omxCheckTrue(all(out2$Raw.SE[8:10]!=0))
omxCheckTrue(all(out3$Raw.SE[8:10]!=0))
omxCheckCloseEnough(out2$Raw.SE[8:10],mxSE(M,m2)[1,1:3],1e-8)
omxCheckCloseEnough(out3$Raw.SE[8:10],mxSE(M,m3)[1,1:3],1e-8)
omxCheckCloseEnough(out2$Raw.SE[8:10],out3$Raw.SE[8:10],1e-8)
