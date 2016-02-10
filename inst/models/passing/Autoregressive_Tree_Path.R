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
  height <- rnorm(1, mean=(hm*2 + hf) / 3, sd=ifelse(male, 7, 6))
  height <- height + ifelse(male, 1, -1) * maleEffect/2
  nextRow <- 1L+nrow(df)
  if (is.null(df)) nextRow <- 1L
  nperson <- data.frame(personID=nextRow, male=male,
                        motherID=mother, fatherID=father, height=height)
  rbind(df, nperson)
}

pdat <- mkPerson(NULL)
for (px in 1:49) {
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
              mxPath("person.height", "height", joinKey='motherID', labels="fromMother"),
              mxPath("person.height", "height", joinKey='fatherID', labels="fromFather"))

#m1 <- mxOption(m1, "RAM Inverse Optimization", "Yes")

if (0) {
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE, iterations=2),
	    mxComputeReportDeriv(),
	    mxComputeReportExpectation()
	))
	m1 <- mxRun(mxModel(m1, plan))

	eo = m1$expectation$output
	ed = m1$expectation$debug
	(solve(diag(nrow(ed$macroA))-ed$macroA) != 0)[ed$layout$numJoins == 0,,drop=FALSE]
	
} else {
  m1 <- mxRun(m1)
  summary(m1)
}

omxCheckCloseEnough(m1$output$fit, 280.489, 1e-2)
# cat(deparse(round(coef(m1),3)))
omxCheckCloseEnough(coef(m1), c(13.864, 21.644, 0.502, 0.474, -3.747), 1e-2)
