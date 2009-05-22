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

manifest <- c("x1", "x2", "y1", "y2", "y3")
latent <- c("L")

model <- mxModel(type = "RAM", manifestVars = manifest, latentVars = latent)

model <- mxModel(model, 
	mxPath(from = manifest, arrows = 2, free = TRUE))

model <- mxModel(model, 
	mxPath(from = latent, arrows = 2, free = TRUE))
	
model <- mxModel(model,
	mxPath(from = "x1", to = "x2", arrows = 2, free = TRUE))

model <- mxModel(model,
	mxPath(from = c("x1", "x2"), to = "L", arrows = 1, free = TRUE))

model <- mxModel(model,
	mxPath(from = "L", to = c("y1","y2","y3"), 
		arrows = 1, free = TRUE))
