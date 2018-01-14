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


# auto-regression tree
#
# This is not really flexible enough to model a pedigree
# because it is not easy to set the initial conditions.

set.seed(1)

maleEffect <- 13 #cm

mkPerson <- function(df, mother=NA, father=NA) {
  male <- runif(1) > .5
  hm <- 0 #162
  if (!is.na(mother)) {
    hm <- df[mother, 'height']
  }
  hf <- 0 #175
  if (!is.na(father)) {
    hf <- df[father, 'height']
  }
  height <- rnorm(1, mean=(hm + hf) / 2, sd=ifelse(male, 7, 6))
  height <- height + ifelse(male, 1, -1) * maleEffect/2
  nextRow <- 1L+nrow(df)
  if (is.null(df)) nextRow <- 1L
  nperson <- data.frame(personID=nextRow, male=male,
                        motherID=mother, fatherID=father, height=height)
  rbind(df, nperson)
}

pdat <- mkPerson(NULL)
for (px in 1:39) {
  mother <- NA
  father <- NA
  if (any(pdat$male)) {
    father = pdat[sample(which(pdat$male), 1), 'personID']
  }
  if (any(!pdat$male)) {
    mother = pdat[sample(which(!pdat$male), 1), 'personID']
  }
  pdat <- mkPerson(pdat, father, mother)
}

library(OpenMx)

pdat$male <- as.numeric(pdat$male)
pdat[is.na(pdat$motherID) | is.na(pdat$fatherID), 'height'] <- NA

m1 <- mxModel("person", type="RAM",
              manifestVars = c("height"), latentVars = c("maleEffect"),
              mxData(type="raw", observed=pdat, primaryKey = 'personID'),
              mxPath("one", "maleEffect", free=FALSE, labels="data.male"),
              mxPath("maleEffect", "height", values=10),
              mxPath("one", "height"),
              mxPath("height", arrows=2, values=1),
              mxMatrix(name="fromMother", nrow=2, ncol=2,
                       labels=paste0("x", 1:4),
                       dimnames=list(c('height', 'maleEffect'),
                                     c('height', 'maleEffect')),
                       joinKey = "motherID", joinModel = "person"),
              mxMatrix(name="fromFather", nrow=2, ncol=2,
                       labels=paste0("x", 1:4),
                       dimnames=list(c('height', 'maleEffect'),
                                     c('height', 'maleEffect')),
                       joinKey = "fatherID", joinModel = "person"))

m1$fromMother$free['height', 'height'] <- TRUE
m1$fromFather$free['height', 'height'] <- TRUE

m1$expectation$between <- c("fromMother","fromFather")

if (0) {
  m1 <- mxRun(mxModel(m1, mxComputeOnce('fitfunction', 'fit')))
} else {
  m1 <- mxRun(m1)
  summary(m1)  # 102.5036
}

omxCheckCloseEnough(m1$output$fit, 220.699, 1e-2)
omxCheckCloseEnough(coef(m1), c(13.994, 21.607, -6.213, 0.416), 1e-2)
