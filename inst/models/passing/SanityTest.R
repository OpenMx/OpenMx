#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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

ppml <- getOption('mxOptions')$UsePPML

omxCheckTrue(!is.null(ppml) && ppml == "No")

omxCheckError(mxOption(mxPath('one'), "bar"),
              "The first argument to mxOption must be an MxModel, not 'MxPath'")

mre <- mxOption(key="mvnRelEps")
mxOption(key="mvnRelEps", value= mxOption(key="mvnRelEps") / 5)
omxCheckEquals(mxOption(key="mvnRelEps"), mre / 5)
mxOption(reset=TRUE)
omxCheckEquals(mxOption(key="mvnRelEps"), mre)
