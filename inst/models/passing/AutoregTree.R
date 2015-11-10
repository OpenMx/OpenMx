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
for (px in 1:19) {
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
              mxPath("one", "height", values=100),
              mxPath("height", arrows=2, values=1),
              mxMatrix(name="Z", nrow=1, ncol=1, dimnames=list('height', 'height'), values=.4),
              mxFitFunctionML(fellner = TRUE))

m1$expectation$join <- list(
  mxJoin(foreignKey = "motherID", expectation = "person", regression = "Z"),
  mxJoin(foreignKey = "fatherID", expectation = "person", regression = "Z")
)

if (0) {
  m1 <- mxRun(mxModel(m1, mxComputeOnce('fitfunction', 'fit')))
} else {
  m1 <- mxRun(m1)
  summary(m1)  # 102.5036
}

omxCheckCloseEnough(m1$output$fit, 102.5036, 1e-2)
