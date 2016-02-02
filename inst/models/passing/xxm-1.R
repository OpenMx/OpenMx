# Bivariate random-intercepts model
# Originally from http://xxm.times.uh.edu/learn-xxm/xxm-tutorial/

library(OpenMx)

sData <- suppressWarnings(try(read.csv("models/passing/data/brim.student.csv"), silent=TRUE))
if (is(sData, "try-error")) sData <- read.csv("data/brim.student.csv")

sData$teacher <- as.integer(sData$teacher)

tMod <- mxModel("teacher", type="RAM",
                latentVars = c("eta1", "eta2"),
                mxData(type="raw", observed=data.frame(teacher=unique(sData$teacher)),
                       primaryKey = "teacher"),
                mxPath(c("eta1", "eta2"), connect = "unique.pairs", arrows=2,
                       labels=paste0("psi", 1:3), values=c(.1,0,.1)))

sMod <- mxModel("student", type="RAM", tMod,
                manifestVars = c("y1", "y2"),
                mxData(type="raw", observed=sData, sort=FALSE),
                mxPath("one", paste0("y",1:2)),
                mxPath(c("y1", "y2"), connect = "unique.pairs", arrows=2,
                       labels=paste0("theta", 1:3), values=c(1,0,1)),
                mxPath(paste0("teacher.eta", 1:2), paste0("y",1:2),
                       free=FALSE, values=1, joinKey = "teacher"))

#sMod$expectation$verbose <- 2L

if (0) {
	dist <- mxRun(mxModel(sMod,
			      mxComputeSequence(list(
				  mxComputeOnce('expectation', 'distribution', 'flat'),
				  mxComputeReportExpectation()))))

	str(dist$expectation$output)
	str(dist$expectation$debug)
	eo = dist$expectation$output
	ed = dist$expectation$debug

	ed$layout
	dist$expectation$output$covariance[1:10,1:10]
	dist$expectation$output$covariance[120:130,120:130]
	eigen(dist$expectation$output$covariance)$val

	round(dist$expectation$debug$A[1:20,1:20], 2)
	dist$expectation$debug$S[1:10,1:10]

}

sMod <- mxRun(sMod)

omxCheckCloseEnough(sMod$output$fit, 445.4331, .01)
